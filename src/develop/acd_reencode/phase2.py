#!/usr/bin/env python3
"""Phase 2 of ../plan/acd_cut_reencoding.md: synthesize encoders/decoders for the Phase 1 groups.

For every (circuit, T, side, grouping) that passed the go/no-go gate (K=16), groups are applied
one at a time (largest feasible saving first) to a 2-stage pipeline netlist:
  1. functional dependence: register only S, decode G \\ S from S after the register
  2. general re-encoding:  r code bits = existing G signals where possible + encoder bits
                           computed from the window leaves L; decoder rebuilds G from the code
Encoder/decoder are synthesized by ABC from truth tables (read_truth; strash; resyn2 | dc2 | ...).
ABC's read_truth has no don't-cares, so unused codes are filled with 0 and with 1 and the smaller
result is kept (conservative w.r.t. real don't-care optimization).
A candidate is accepted iff the total area (AND + 10*reg) strictly decreases and both stage depths
stay <= T, recomputed on the whole netlist after the change. Every final netlist is checked with
ABC cec (registers as wires) before any number is reported.

usage: phase2.py [--only circuit[,..]] [--force circuit:tag:side:grouping] [--inject-bad-decoder]
outputs: results/phase2_results.csv, results/phase2_groups.csv, results/phase2/<run>/{pipe,comb}.aig, cec.log
"""
import argparse
import csv
import json
import os
import subprocess
import sys
import tempfile
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
ABC = os.path.join(ROOT, "abc")
sys.path.insert(0, HERE)
import aiger  # noqa: E402
from phase0 import BASELINE  # noqa: E402

REG_W = 10
IGNORE_DELAY = False   # --ignore-delay: accept on area alone (stage depths may exceed T)
RESULTS = os.path.join(HERE, "results")   # --results
RESYN2 = "balance; rewrite; refactor; balance; rewrite; rewrite -z; balance; refactor -z; rewrite -z; balance"
SCRIPTS = ["strash", "strash; " + RESYN2, "strash; dc2", "strash; dc2; " + RESYN2 + "; dc2"]


# ---------------------------------------------------------------------------------------------
# Pipeline netlist: id -> (kind, a, b, stage); kind 'pi' | 'and' | 'reg'. A literal is 2*id + c,
# id 0 is constant false. A reg's input literal is a (a stage-0 signal); its output is in stage 1.

class Pipe:
    def __init__(self):
        self.nodes, self.pos, self.reg_of = {}, [], {}
        self.pis, self.in_names, self.out_names = [], [], []
        self.next_id = 1

    def copy(self):
        p = Pipe()
        p.nodes, p.pos, p.reg_of = dict(self.nodes), list(self.pos), dict(self.reg_of)
        p.pis, p.in_names, p.out_names, p.next_id = self.pis, self.in_names, self.out_names, self.next_id
        return p

    def new(self, kind, a, b, stage):
        i = self.next_id; self.next_id += 1
        self.nodes[i] = (kind, a, b, stage)
        return i

    def fanins(self, i):
        k, a, b, _ = self.nodes[i]
        return (a, b) if k == "and" else ((a,) if k == "reg" else ())

    def topo(self):
        """All nodes reachable from the POs, fanins first (regs are ordinary nodes here)."""
        order, state = [], {}
        for lit in self.pos:
            root = lit >> 1
            if not root or root in state:
                continue
            stack = [(root, 0)]
            while stack:
                v, i = stack.pop()
                if i == 0:
                    if v in state:
                        continue
                    state[v] = 1
                fi = self.fanins(v)
                if i < len(fi):
                    stack.append((v, i + 1))
                    u = fi[i] >> 1
                    if u and u not in state:
                        stack.append((u, 0))
                else:
                    order.append(v)
        return order

    def sweep(self):
        live = set(self.topo())
        for i in [i for i in self.nodes if i not in live and self.nodes[i][0] != "pi"]:
            del self.nodes[i]
        for u in [u for u, r in self.reg_of.items() if r not in self.nodes]:
            del self.reg_of[u]

    def stats(self):
        ands = sum(1 for n in self.nodes.values() if n[0] == "and")
        regs = sum(1 for n in self.nodes.values() if n[0] == "reg")
        return ands, regs, ands + REG_W * regs

    def depths(self):
        """(stage-0 depth, stage-1 depth); raises on an illegal structure."""
        lev, d = {0: 0}, [0, 0]
        for v in self.topo():
            k, a, b, st = self.nodes[v]
            if k == "pi":
                lev[v] = 0
            elif k == "reg":
                src = self.nodes.get(a >> 1)
                if a >> 1 and src[3] != 0:
                    raise AssertionError(f"reg {v} reads a stage-1 signal")
                lev[v] = 0
            else:
                for x in (a >> 1, b >> 1):
                    if x and self.nodes[x][3] != st and not (st == 1 and self.nodes[x][0] == "reg"):
                        raise AssertionError(f"AND {v} (stage {st}) reads node {x} of another stage")
                    if x and st == 0 and self.nodes[x][0] == "reg":
                        raise AssertionError(f"stage-0 AND {v} reads a register")
                lev[v] = 1 + max(lev[a >> 1], lev[b >> 1])
                d[st] = max(d[st], lev[v])
        for lit in self.pos:
            x = lit >> 1
            if x and not (self.nodes[x][0] == "reg" or self.nodes[x][3] == 1):
                raise AssertionError(f"PO driven by unregistered stage-0 node {x}")
        return d

    def instantiate(self, net, inputs, stage):
        """Copy an ABC result (nin, ands, outs) with its inputs bound to literals."""
        nin, ands, outs = net
        m = {0: 0}
        for i in range(nin):
            m[i + 1] = inputs[i] if i < len(inputs) else 0
        lit = lambda l: m[l >> 1] ^ (l & 1)
        for k, (a, b) in enumerate(ands):
            m[nin + 1 + k] = 2 * self.new("and", lit(a), lit(b), stage)
        return [lit(o) for o in outs]

    def substitute_regs(self, subst):
        """subst: reg id -> new literal; rewrites stage-1 fanins and POs, then drops the regs."""
        if not subst:
            return
        f = lambda l: (subst[l >> 1] ^ (l & 1)) if (l >> 1) in subst else l
        for i, (k, a, b, st) in list(self.nodes.items()):
            if k == "and" and st == 1 and ((a >> 1) in subst or (b >> 1) in subst):
                self.nodes[i] = (k, f(a), f(b), st)
        self.pos = [f(l) for l in self.pos]
        for r in subst:
            del self.nodes[r]
        self.sweep()


def build_pipe(g, stage, cut):
    p = Pipe()
    p.in_names, p.out_names = g.in_names, g.out_names
    p.next_id = g.I + len(g.ands) + 1
    for i in range(g.I):
        p.nodes[i + 1] = ("pi", None, None, 0)
        p.pis.append(i + 1)
    for u in cut:
        p.reg_of[u] = p.new("reg", 2 * u, None, 1)
    st = lambda u: stage.get(u, 0)

    def lit_in(l, s):
        u = l >> 1
        if u and s == 1 and st(u) == 0:
            return 2 * p.reg_of[u] + (l & 1)
        return l
    for k, (a, b) in enumerate(g.ands):
        v = g.I + 1 + k
        if v in stage:
            p.nodes[v] = ("and", lit_in(a, stage[v]), lit_in(b, stage[v]), stage[v])
    p.pos = [lit_in(l, 1) for l in g.outputs]
    p.sweep()
    return p


def write_pipe(p, path, comb_path):
    """Pipelined AIG (registers = latches) and the combinational AIG (registers = wires)."""
    order = p.topo()
    regs = [v for v in order if p.nodes[v][0] == "reg"]
    for comb, out in ((False, path), (True, comb_path)):
        g = aiger.Aiger(); g.I = len(p.pis)
        g.in_names, g.out_names = p.in_names, p.out_names
        nl = 0 if comb else len(regs)
        litof = {0: 0}                        # node -> literal of its value in the output AIG
        for i, v in enumerate(p.pis):
            litof[v] = 2 * (i + 1)
        if not comb:
            for j, v in enumerate(regs):
                litof[v] = 2 * (g.I + 1 + j)
        lit = lambda l: litof[l >> 1] ^ (l & 1)
        for v in order:
            k, a, b, _ = p.nodes[v]
            if k == "and":
                litof[v] = 2 * (g.I + nl + 1 + len(g.ands))
                g.ands.append((lit(a), lit(b)))
            elif k == "reg" and comb:
                litof[v] = lit(a)
        if not comb:
            g.latches = [lit(p.nodes[v][1]) for v in regs]
        g.outputs = [lit(l) for l in p.pos]
        aiger.write(g, out)


# ---------------------------------------------------------------------------------------------
# Window simulation and code assignment

def var_pattern(i, P):
    bits = ((np.arange(P) >> i) & 1).astype(np.uint8)
    return int.from_bytes(np.packbits(bits, bitorder="little").tobytes(), "little")


def to_bits(x, P):
    raw = x.to_bytes(max(1, (P + 7) // 8), "little")
    return np.unpackbits(np.frombuffer(raw, np.uint8), bitorder="little")[:P]


def window_values(g, leaves, members):
    """Each member as a P-bit int over all 2^|L| leaf assignments (original AIG logic)."""
    P = 1 << len(leaves)
    mask = (1 << P) - 1
    val = {0: 0}
    for i, l in enumerate(leaves):
        val[l] = var_pattern(i, P)

    def get(v):
        stack = [v]
        while stack:
            u = stack[-1]
            if u in val:
                stack.pop(); continue
            k = u - g.I - 1
            if k < 0:
                raise ValueError(f"PI {u} is not a window leaf")
            a, b = g.ands[k]
            need = [x >> 1 for x in (a, b) if (x >> 1) not in val]
            if need:
                stack += need; continue
            x, y = val[a >> 1], val[b >> 1]
            val[u] = (x ^ mask if a & 1 else x) & (y ^ mask if b & 1 else y)
            stack.pop()
        return val[v]
    return P, [get(m) for m in members]


def image(P, vals):
    M = np.stack([to_bits(v, P) for v in vals], axis=1)            # P x |G|
    packed = np.packbits(M, axis=1)
    _, first, inv = np.unique(packed, axis=0, return_index=True, return_inverse=True)
    return M[first].astype(bool), inv.reshape(-1)                  # image points, pattern -> point


def ceil_log2(n):
    return max(0, (n - 1).bit_length())


def assign_code(U):
    """Greedily reuse G columns as code bits while the total stays r; the rest are new bits."""
    n, m = U.shape
    r = ceil_log2(n)
    cls = np.zeros(n, dtype=np.int64)
    reused = []
    while True:
        best = None
        for j in range(m):
            if j in reused:
                continue
            c2 = cls * 2 + U[:, j]
            _, cnt = np.unique(c2, return_counts=True)
            if len(reused) + 1 + ceil_log2(int(cnt.max())) <= r:
                key = (len(cnt), -int(cnt.max()))
                if best is None or key > best[0]:
                    best = (key, j, c2)
        if best is None:
            break
        reused.append(best[1]); cls = best[2]
    _, cls_id, cnt = np.unique(cls, return_inverse=True, return_counts=True)
    nnew = ceil_log2(int(cnt.max()))
    idx = np.zeros(n, dtype=np.int64)                              # index of a point within its class
    seen = {}
    for p in range(n):
        c = cls_id[p]; idx[p] = seen.get(c, 0); seen[c] = idx[p] + 1
    code = np.zeros(n, dtype=np.int64)
    for k, j in enumerate(reused):
        code |= U[:, j].astype(np.int64) << k
    code |= idx << len(reused)
    return r, reused, nnew, idx, code


# ---------------------------------------------------------------------------------------------
# ABC synthesis from truth tables

def synth(tts, nvars, tmp):
    """tts: list of ints over nvars inputs. Returns candidates [(ands, depth, (nin, ands, outs))]."""
    n = max(nvars, 2)
    P = 1 << n
    lines = []
    for t in tts:
        if nvars < 2:                         # pad: extra inputs are don't-touch
            t = sum(((t >> (p & ((1 << nvars) - 1))) & 1) << p for p in range(P))
        lines.append(format(t, "0%dx" % (P // 4)))
    ttf = os.path.join(tmp, "tt.txt")
    open(ttf, "w").write("\n".join(lines) + "\n")
    cmds, outs = [], []
    for k, s in enumerate(SCRIPTS):
        o = os.path.join(tmp, f"s{k}.aig")
        cmds.append(f"read_truth -f {ttf}; {s}; write_aiger {o}")
        outs.append(o)
    r = subprocess.run([ABC, "-c", "; ".join(cmds)], capture_output=True, text=True)
    if r.returncode:
        raise RuntimeError("ABC failed: " + r.stdout[-500:] + r.stderr[-500:])
    cands = []
    for o in outs:
        g = aiger.read(o)
        assert g.I == n and len(g.outputs) == len(tts)
        lev = [0] * (g.I + 1)
        for a, b in g.ands:
            lev.append(1 + max(lev[a >> 1], lev[b >> 1]))
        depth = max([lev[o >> 1] for o in g.outputs] + [0])
        cands.append((len(g.ands), depth, (g.I, g.ands, g.outputs)))
    return cands


def pick(cands):
    """Min-AND candidate first, then the min-depth one if it differs."""
    a = min(cands, key=lambda c: (c[0], c[1]))
    d = min(cands, key=lambda c: (c[1], c[0]))
    return [a] if d is a or (d[0], d[1]) == (a[0], a[1]) else [a, d]


def decoder_cands(U_code, code, nvars, outputs_cols, tmp):
    """Decoder truth tables over the code bits; unused codes filled with 0 and with 1."""
    N = 1 << nvars
    used = np.zeros(N, dtype=np.uint8); used[code] = 1
    asint = lambda arr: int.from_bytes(np.packbits(arr, bitorder="little").tobytes(), "little")
    dc = asint(1 - used) & ((1 << N) - 1)
    on = []
    for j in outputs_cols:
        t = np.zeros(N, dtype=np.uint8); t[code[U_code[:, j]]] = 1
        on.append(asint(t))
    if not on:
        return [(0, 0, (nvars, [], []))]
    c0 = synth(on, nvars, tmp)
    c1 = synth([t | dc for t in on], nvars, tmp) if dc else []
    return c0 + c1


# ---------------------------------------------------------------------------------------------

def try_candidate(net, T, area0, build):
    """Apply build(copy) and evaluate. Returns (status, new_net, info)."""
    p = net.copy()
    build(p)
    ands, regs, area = p.stats()
    d = p.depths()
    if area >= area0:
        return "area", None, (ands, regs, area, d)
    if max(d) > T and not IGNORE_DELAY:
        return "delay", None, (ands, regs, area, d)
    return "ok", p, (ands, regs, area, d)


notes = []   # per-candidate details of the last process_group call: method:status:area_delta:depths


def process_group(net, g, grp, T, tmp, bad):
    members = [int(x) for x in grp["members"].split()]
    fdset = [int(x) for x in grp["fd_subset"].split()]
    leaves = [int(x) for x in grp["leaves"].split()]
    P, vals = window_values(g, leaves, members)
    U, inv = image(P, vals)
    nimg, m = U.shape
    area0 = net.stats()[2]
    tried = []                                       # statuses of all evaluated candidates
    notes.clear()
    regs = [net.reg_of[u] for u in members]

    # 1. functional dependence: register S only
    sidx = [members.index(s) for s in fdset]
    if len(fdset) < m and len(fdset) <= 16:
        Scode = np.zeros(nimg, dtype=np.int64)
        for k, j in enumerate(sidx):
            Scode |= U[:, j].astype(np.int64) << k
        if len(np.unique(Scode)) == nimg:            # S separates every image point
            rest = [j for j in range(m) if j not in sidx]
            for cand in pick(decoder_cands(U, Scode, len(sidx), rest, tmp)):
                def build(p, cand=cand):
                    outs = p.instantiate(cand[2], [2 * regs[j] for j in sidx], 1)
                    if bad and outs:
                        outs[0] ^= 1
                    p.substitute_regs({regs[j]: o for j, o in zip(rest, outs)})
                st, p2, info = try_candidate(net, T, area0, build)
                tried.append(st); notes.append(f"fd[dec={cand[0]}/d{cand[1]}]:{st}:{info[2] - area0:+d}:{info[3][0]}/{info[3][1]}")
                if st == "ok":
                    return p2, "fd", info, tried, cand[0], 0
    # 2. general re-encoding with r bits
    r, reused, nnew, idx, code = assign_code(U)
    if r < m:
        rest = [j for j in range(m) if j not in reused]
        dec = pick(decoder_cands(U, code, r, rest, tmp))
        if nnew:
            enc_tts = [0] * nnew
            pidx = idx[inv]                          # pattern -> index within class
            for t in range(nnew):
                bits = ((pidx >> t) & 1).astype(np.uint8)
                enc_tts[t] = int.from_bytes(np.packbits(bits, bitorder="little").tobytes(), "little")
            enc = pick(synth(enc_tts, len(leaves), tmp))
        else:
            enc = [(0, 0, (0, [], []))]
        combos = [(enc[0], dec[0])]
        if len(enc) > 1 or len(dec) > 1:
            combos.append((enc[-1], dec[-1]))
        for e, d in combos:
            def build(p, e=e, d=d):
                newbits = []
                if nnew:
                    for o in p.instantiate(e[2], [2 * l for l in leaves], 0):
                        newbits.append(2 * p.new("reg", o, None, 1))
                ins = [2 * regs[j] for j in reused] + newbits
                outs = p.instantiate(d[2], ins, 1)
                if bad and outs:
                    outs[0] ^= 1
                p.substitute_regs({regs[j]: o for j, o in zip(rest, outs)})
            st, p2, info = try_candidate(net, T, area0, build)
            tried.append(st); notes.append(f"enc[r={r},reuse={len(reused)},enc={e[0]}/d{e[1]},dec={d[0]}/d{d[1]}]:{st}:{info[2] - area0:+d}:{info[3][0]}/{info[3][1]}")
            if st == "ok":
                return p2, "encode", info, tried, d[0], e[0]
    return None, None, None, tried, None, None


def run_combo(circuit, tag, side, grouping, groups, args):
    depth, _, trec, _, _ = BASELINE[circuit]
    T = (depth + 1) // 2 if tag == "50%" else trec
    rdir = os.path.join(RESULTS, "runs", f"{circuit}_T{T}_{side}")
    st = {n["id"]: n["s"] for n in json.load(open(os.path.join(rdir, "stage_assign.json")))["nodes"] if n["type"] == "and"}
    cut = [s["node"] for s in json.load(open(os.path.join(rdir, "cuts.json")))["cuts"][0]["signals"]]
    g = aiger.read(os.path.join(args.epfl, circuit + ".aig"))
    net = build_pipe(g, st, cut)
    and0, reg0, area0 = net.stats()
    assert reg0 == len(cut)
    d0 = net.depths()
    assert max(d0) <= T, (d0, T)
    out = os.path.join(RESULTS, "phase2" + args.tag, f"{circuit}_T{T}_{side}_{grouping}")
    os.makedirs(out, exist_ok=True)
    rows, counts = [], {"tried": 0, "accepted": 0, "area": 0, "delay": 0, "skipped": 0}
    order = sorted(groups, key=lambda x: -int(x["saving_feasible"]))
    t0 = time.time()
    with tempfile.TemporaryDirectory() as tmp:
        for grp in order:
            if int(grp["saving"]) <= 0:
                continue
            counts["tried"] += 1
            before = net.stats()
            new, method, info, tried, dec_and, enc_and = process_group(net, g, grp, T, tmp, args.inject_bad_decoder)
            if new is not None:
                net = new
                counts["accepted"] += 1
                status = "accepted"
            elif not tried:
                counts["skipped"] += 1; status = "skipped"
            elif "delay" in tried:
                counts["delay"] += 1; status = "rejected_delay"
            else:
                counts["area"] += 1; status = "rejected_area"
            after = net.stats()
            rows.append({"circuit": circuit, "T": T, "side": side, "grouping": grouping, "group_id": grp["group_id"],
                         "|G|": grp["|G|"], "|L|": grp["|L|"], "r": grp["r"], "fd_subset_size": grp["fd_subset_size"],
                         "upstream_slack": grp["upstream_slack"], "downstream_slack": grp["downstream_slack"],
                         "status": status, "method": method or "", "candidates": " ".join(notes),
                         "enc_and": enc_and if enc_and is not None else "", "dec_and": dec_and if dec_and is not None else "",
                         "reg_delta": after[1] - before[1], "and_delta": after[0] - before[0],
                         "area_delta": after[2] - before[2], "stage_depths": "/".join(map(str, net.depths()))})
    and1, reg1, area1 = net.stats()
    d1 = net.depths()
    pipe, comb = os.path.join(out, "pipe.aig"), os.path.join(out, "comb.aig")
    write_pipe(net, pipe, comb)
    r = subprocess.run([ABC, "-c", f"cec -n {os.path.join(args.epfl, circuit + '.aig')} {comb}"], capture_output=True, text=True)
    open(os.path.join(out, "cec.log"), "w").write(r.stdout + r.stderr)
    eq = "Networks are equivalent" in r.stdout
    res = {"circuit": circuit, "T": T, "T_tag": tag, "side": side, "grouping": grouping, "cec": "equivalent" if eq else "FAIL"}
    if eq:
        res.update({"reg_before": reg0, "reg_after": reg1, "and_before": and0, "and_after": and1,
                    "area_before": area0, "area_after": area1,
                    "area_delta_pct": round(100.0 * (area1 - area0) / area0, 3),
                    "stage_depths_after": f"{d1[0]}/{d1[1]}"})
    else:
        rows = []
        print(f"  cec FAILED for {circuit} {tag} {side} {grouping}: aborting, no numbers reported")
    res.update({"groups_tried": counts["tried"], "groups_accepted": counts["accepted"],
                "groups_rejected_area": counts["area"], "groups_rejected_delay": counts["delay"],
                "groups_skipped": counts["skipped"], "seconds": round(time.time() - t0, 1)})
    return res, rows


RES_COLS = ["circuit", "T", "T_tag", "side", "grouping", "reg_before", "reg_after", "and_before", "and_after",
            "area_before", "area_after", "area_delta_pct", "groups_tried", "groups_accepted",
            "groups_rejected_area", "groups_rejected_delay", "groups_skipped", "stage_depths_after", "cec", "seconds"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--epfl", default=os.path.join(ROOT, "..", "testcases", "aig", "epfl"))
    ap.add_argument("--only", default="")
    ap.add_argument("--force", action="append", default=[], help="circuit:tag:side:grouping, ignores the gate")
    ap.add_argument("--inject-bad-decoder", action="store_true", help="unit test: flip one decoder output")
    ap.add_argument("--ignore-delay", action="store_true", help="accept candidates on area alone")
    ap.add_argument("--all-with-saving", action="store_true",
                    help="run every K=16 (circuit, T, side) of the given grouping with sum_saving > 0, ignoring the gate")
    ap.add_argument("--grouping", default="plan")
    ap.add_argument("--results", default=os.path.join(HERE, "results"),
                    help="directory holding phase0/phase1 outputs; phase2 outputs go there too")
    ap.add_argument("--tag", default="", help="suffix for the output CSV names and directory")
    a = ap.parse_args()
    global RESULTS
    RESULTS = a.results
    summ = list(csv.DictReader(open(os.path.join(RESULTS, "phase1_summary.csv"))))
    allg = list(csv.DictReader(open(os.path.join(RESULTS, "phase1_groups.csv"))))
    global IGNORE_DELAY
    IGNORE_DELAY = a.ignore_delay
    if a.all_with_saving:
        combos = [(s["circuit"], s["T_tag"], s["side"], s["grouping"]) for s in summ
                  if s["K"] == "16" and s["grouping"] == a.grouping and int(s["sum_saving"]) > 0]
    else:
        combos = [(s["circuit"], s["T_tag"], s["side"], s["grouping"]) for s in summ if s["K"] == "16" and s["go"] == "1"]
    combos += [tuple(f.split(":")) for f in a.force]
    if a.only:
        combos = [c for c in combos if c[0] in a.only.split(",")]
    results, grows = [], []
    for c in combos:
        gs = [x for x in allg if (x["circuit"], x["T_tag"], x["side"], x["grouping"]) == c and x["K"] == "16"]
        res, rows = run_combo(*c, gs, a)
        results.append(res); grows += rows
        print(f"{c[0]:10s} {c[1]:3s} {c[2]:6s} {c[3]:5s} cec={res['cec']:10s} reg {res.get('reg_before')}->{res.get('reg_after')} "
              f"and {res.get('and_before')}->{res.get('and_after')} area {res.get('area_delta_pct')}% "
              f"accepted {res['groups_accepted']}/{res['groups_tried']} (area {res['groups_rejected_area']}, "
              f"delay {res['groups_rejected_delay']}) {res['seconds']}s")
    sfx = a.tag
    with open(os.path.join(RESULTS, f"phase2_results{sfx}.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=RES_COLS); w.writeheader(); w.writerows(results)
    if grows:
        with open(os.path.join(RESULTS, f"phase2_groups{sfx}.csv"), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(grows[0])); w.writeheader(); w.writerows(grows)
    return 0 if all(r["cec"] == "equivalent" for r in results) else 1


if __name__ == "__main__":
    sys.exit(main())
