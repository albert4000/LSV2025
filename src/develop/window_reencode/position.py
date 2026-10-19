#!/usr/bin/env python3
"""Where do improving kicks sit? Position / neighbourhood features of the kick root, Phase 2e kicks.

Phase 2e: 200 single kicks per (circuit, sampler) from the raw EPFL AIG x0, R = compress2rs x3, K = 8.
The root ids of results/phase2e/raw/<c>_<m>.csv refer to results/phase2e/kicks/<c>_<m>/x.aig.
Features of the root v in x0:
  rel_pos    lev(v) / (lev(v) + down(v)): 0 = next to the PIs, 1 = next to the POs (longest paths)
  lev, down  longest path from the PIs / to a PO
  fanout     fanout of v
  supp       number of PIs in the transitive fanin of v
  npo        number of POs in the transitive fanout of v
  xor2 / mux2  a structural XOR / MUX root within 2 levels below v (3-AND patterns, structdiff.py)
  maxfo2     largest fanout among the nodes within 2 levels below v (shared select / control signals)
plus the window columns of the CSV (nL, nC, nG, N, kick cost = and_kick - and_x).
Within each circuit, AUC(feature) = P(feature of an improving kick > feature of a non-improving kick)
(ties count 1/2); 0.5 = no relation. Circuits with fewer than 5 improving kicks are skipped. Also the
improving rate per rel_pos bin.
usage: position.py
"""
import collections, csv, os, statistics as st

from structdiff import Aig

HERE = os.path.dirname(os.path.abspath(__file__))
P2E = os.path.join(HERE, "results", "phase2e")
ARITH = ["adder", "bar", "log2", "multiplier", "sin", "sqrt", "square", "max"]          # EPFL arithmetic
CONTROL = ["arbiter", "cavlc", "ctrl", "dec", "i2c", "int2float", "mem_ctrl", "priority", "router", "voter"]
NOISY = {"log2", "mem_ctrl", "i2c", "sin", "cavlc", "max", "voter", "square", "int2float"}   # Phase 2e null kick
FEATS = ["rel_pos", "lev", "down", "fanout", "supp", "npo", "xor2", "mux2", "maxfo2", "nL", "nC", "nG", "N", "cost"]


def node_features(path):
    a = Aig(path)
    n = a.nPi + a.nAnd + 1
    ands = range(a.nPi + 1, n)
    down = [-1] * n
    for o in a.outs:
        down[o >> 1] = max(down[o >> 1], 0)
    for v in reversed(ands):
        if down[v] < 0:
            continue
        for f in (a.f0[v] >> 1, a.f1[v] >> 1):
            down[f] = max(down[f], down[v] + 1)
    supp = [0] * n
    for i in range(1, a.nPi + 1):
        supp[i] = 1 << (i - 1)
    for v in ands:
        supp[v] = supp[a.f0[v] >> 1] | supp[a.f1[v] >> 1]
    po = [0] * n
    for j, o in enumerate(a.outs):
        po[o >> 1] |= 1 << j
    for v in reversed(ands):
        for f in (a.f0[v] >> 1, a.f1[v] >> 1):
            po[f] |= po[v]
    xor, mux = [0] * n, [0] * n
    for v in ands:
        x, y = a.f0[v], a.f1[v]
        if x & 1 and y & 1 and (x >> 1) > a.nPi and (y >> 1) > a.nPi:
            P, Q = {a.f0[x >> 1], a.f1[x >> 1]}, {a.f0[y >> 1], a.f1[y >> 1]}
            k = sum(1 for l in P if l ^ 1 in Q)
            xor[v], mux[v] = k == 2, k == 1
    feats = {}
    for v in ands:
        below = {v}
        fr = [v]
        for _ in range(2):
            fr = [f for u in fr if u > a.nPi for f in (a.f0[u] >> 1, a.f1[u] >> 1)]
            below.update(fr)
        feats[v] = dict(rel_pos=a.lev[v] / max(1, a.lev[v] + max(down[v], 0)), lev=a.lev[v], down=down[v],
                        fanout=a.fo[v], supp=bin(supp[v]).count("1"), npo=bin(po[v]).count("1"),
                        xor2=int(any(xor[u] for u in below)), mux2=int(any(mux[u] for u in below)),
                        maxfo2=max(a.fo[u] for u in below if u != v) if len(below) > 1 else 0)
    return feats


def auc(pos, neg):
    if not pos or not neg:
        return None
    s = 0.0
    for p in pos:
        for q in neg:
            s += 1 if p > q else 0.5 if p == q else 0
    return s / (len(pos) * len(neg))


def main():
    rows = list(csv.DictReader(open(os.path.join(P2E, "phase2e_kicks.csv"))))
    cache = {}
    data = collections.defaultdict(list)   # (circuit, group) -> [(feature dict, better)]
    for r in rows:
        c, m = r["circuit"], r["method"]
        grp = "reenc" if m in ("U", "A", "W4") else m
        if c not in cache:   # x.aig is the same for every sampler of a circuit
            cache[c] = node_features(os.path.join(P2E, "kicks", f"{c}_{m}", "x.aig"))
        f = dict(cache[c].get(int(r["root"]), {}))
        if not f:
            continue
        f.update(nL=int(r["nL"]), nC=int(r["nC"]), nG=int(r["nG"]), N=int(r["N"]),
                 cost=int(r["and_kick"]) - int(r["and_x"]))
        data[(c, grp)].append((f, int(r["d_R"]) < 0))
    for grp in ("reenc", "P0", "N"):
        print(f"\n== {grp}: within-circuit AUC (P(feature of improving kick > others)); circuits with >= 5 improving kicks")
        print(f"{'circuit':11s}{'n+':>5s}" + "".join(f"{k:>8s}" for k in FEATS))
        per = collections.defaultdict(lambda: collections.defaultdict(list))
        for c in ARITH + CONTROL:
            d = data.get((c, grp), [])
            pos = [f for f, b in d if b]
            neg = [f for f, b in d if not b]
            if len(pos) < 5 or len(neg) < 5:
                continue
            vals = {k: auc([f[k] for f in pos], [f[k] for f in neg]) for k in FEATS}
            cls = "arith" if c in ARITH else "control"
            for k, v in vals.items():
                per[cls][k].append(v)
                per["quiet" if c not in NOISY else "noisy"][k].append(v)
            print(f"{c:11s}{len(pos):>5d}" + "".join(f"{vals[k]:>8.2f}" for k in FEATS) + ("  *noisy" if c in NOISY else ""))
        for cls in ("arith", "control", "quiet", "noisy"):
            if per[cls]:
                print(f"{'mean ' + cls:11s}{'':>5s}" + "".join(f"{st.mean(per[cls][k]):>8.2f}" for k in FEATS))
    print("\n== improving rate by rel_pos bin (0-1/3 near PI, 1/3-2/3 middle, 2/3-1 near PO), reenc = U, A, W4 and N")
    for c in ARITH + CONTROL:
        out = []
        for grp in ("reenc", "N"):
            d = data.get((c, grp), [])
            bins = [[b for f, b in d if lo <= f["rel_pos"] < hi] for lo, hi in ((0, 1 / 3), (1 / 3, 2 / 3), (2 / 3, 1.01))]
            out.append("  ".join(f"{st.mean(x):4.0%} (n={len(x):3d})" if x else "   -       " for x in bins))
        print(f"{c:11s} reenc {out[0]}   | N {out[1]}   ({'arith' if c in ARITH else 'control'}{', noisy' if c in NOISY else ''})")

if __name__ == "__main__":
    main()
