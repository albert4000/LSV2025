#!/usr/bin/env python3
"""Residue analysis of one kick: which structure of the kicked AIG y = kick(x) survives R, and which
structure of R(x) is lost.

Files (lsv_reencode -k -W dir, one window per kick): x.aig, Rx.aig, y<k>.aig, y<k>.lab, k<k>.aig (= R(y)).
Every node is keyed by its function (simulation signature up to complement, structdiff.Aig).

Regions of x (from the window ids in y<k>.lab): C = cone minus outputs, G = window outputs,
TFO = transitive fanout of C outside C, OTHER = the rest.
Origins in y (lab): X (copied from x), E (encoder), D (decoder / P0 circuit), T (TFO rebuilt), and the
code bits Y (encoder outputs).

  residue  ANDs of R(y) whose function is not in R(x), by origin of the function:
           code (a code bit that is not a function of x), enc / dec (other new E / D functions),
           old_<region> (a function of x that R(x) had removed but R(y) kept), rnew (made by R from y);
           rnew is split into rnew_down (in the transitive fanout, within R(y), of a node computing a window
           function: C, G, a code bit, E or D) and rnew_far (not downstream of the window at all)
  missing  ANDs of R(x) whose function is not in R(y), by where the function was in x:
           <region>, or rmade (made by R from x, not a function of x), split into rmade_down / rmade_far
           (downstream, within R(x), of a node computing a C or G function of x, or not)
  code bits: support size, algebraic degree (1 = affine, i.e. an XOR of leaves), whether the function
           already exists in x (reused), and whether it survives R (still a function of R(y))
"""
import os

from structdiff import Aig


def load(path):
    a = Aig(path)
    a.simulate()
    return a


def read_lab(path):
    d = {}
    for line in open(path):
        f = line.split()
        d[f[0]] = f[1:]
    lab = d.get("lab", [""])[0]
    ints = lambda k: [int(v) for v in d.get(k, [])]
    return dict(lab=lab, xleaves=ints("xleaves"), xcone=ints("xcone"), xouts=ints("xouts"),
                leaves=ints("leaves"), code=ints("code"), gout=ints("gout"))


def tfo(x, nodes):
    """AND nodes of x in the transitive fanout of nodes (excluding them)."""
    s = set(nodes)
    out = set()
    for v in range(min(s) if s else x.nPi + 1, x.nPi + x.nAnd + 1):
        if v in s:
            continue
        if (x.f0[v] >> 1) in s or (x.f0[v] >> 1) in out or (x.f1[v] >> 1) in s or (x.f1[v] >> 1) in out:
            out.add(v)
    return out


def downstream(a, keys):
    """ANDs of a computing a function in keys, plus their transitive fanout."""
    m = set()
    for v in range(a.nPi + 1, a.nPi + a.nAnd + 1):
        if a.keys[v] in keys or (a.f0[v] >> 1) in m or (a.f1[v] >> 1) in m:
            m.add(v)
    return m


def local_truth(y, leaves, lits, nL):
    """Truth tables (ints over 2^nL leaf patterns) of lits in y, as functions of the leaf literals;
    None for a literal that does not depend on the leaves only."""
    mask = (1 << (1 << nL)) - 1
    val = {0: 0}
    for i, l in enumerate(leaves):
        if l < 0:
            continue
        t = 0
        for p in range(1 << nL):
            if (p >> i) & 1:
                t |= 1 << p
        val[l >> 1] = t ^ (mask if l & 1 else 0)

    def ev(v, depth=0):
        if v in val:
            return val[v]
        if v <= y.nPi or depth > 400:
            return None
        a, b = y.f0[v], y.f1[v]
        x0, x1 = ev(a >> 1, depth + 1), ev(b >> 1, depth + 1)
        if x0 is None or x1 is None:
            val[v] = None
            return None
        val[v] = (x0 ^ (mask if a & 1 else 0)) & (x1 ^ (mask if b & 1 else 0))
        return val[v]

    out = []
    for l in lits:
        t = ev(l >> 1) if l >= 0 else None
        out.append(None if t is None else t ^ (mask if l & 1 else 0))
    return out


def degree(t, n):
    """Algebraic (ANF) degree of a truth table over n variables."""
    a = [(t >> p) & 1 for p in range(1 << n)]
    for i in range(n):
        for p in range(1 << n):
            if (p >> i) & 1:
                a[p] ^= a[p ^ (1 << i)]
    return max((bin(p).count("1") for p in range(1 << n) if a[p]), default=0)


def support(t, n):
    s = 0
    for i in range(n):
        m = 0
        for p in range(1 << n):
            if not (p >> i) & 1 and ((t >> p) & 1) != ((t >> (p | 1 << i)) & 1):
                m = 1
                break
        s += m
    return s


class Base:
    """x and R(x) of one job, shared by its kicks."""

    def __init__(self, d):
        self.x = load(os.path.join(d, "x.aig"))
        self.Rx = load(os.path.join(d, "Rx.aig"))
        self.Fx = set(self.x.keys)
        self.FRx = set(self.Rx.keys)
        self.ands = lambda a: range(a.nPi + 1, a.nPi + a.nAnd + 1)


def analyze(base, d, k):
    x, Rx = base.x, base.Rx
    L = read_lab(os.path.join(d, f"y{k}.lab"))
    y = load(os.path.join(d, f"y{k}.aig"))
    z = load(os.path.join(d, f"k{k}.aig"))
    # regions of x, as function keys
    cone, outs = set(L["xcone"]), set(L["xouts"])
    reg = {}
    for v in base.ands(x):
        r = "G" if v in outs else "C" if v in cone else None
        if r:
            reg.setdefault(x.keys[v], r)
    for v in tfo(x, cone):
        reg.setdefault(x.keys[v], "TFO")
    for v in base.ands(x):
        reg.setdefault(x.keys[v], "OTHER")
    for v in range(0, x.nPi + 1):
        reg[x.keys[v]] = "PI"
    # origins in y
    codeK = {y.keys[l >> 1] for l in L["code"] if l >= 0}
    encK, decK = set(), set()
    lab = L["lab"]
    for v in base.ands(y):
        if lab[v] == "E":
            encK.add(y.keys[v])
        elif lab[v] == "D":
            decK.add(y.keys[v])
    Fz = set(z.keys)
    winK = {f for f, r in reg.items() if r in ("C", "G")}
    downZ = downstream(z, winK | codeK | encK | decK)
    downRx = downstream(Rx, winK)
    row = {}
    res = dict(code=0, enc=0, dec=0, old_C=0, old_G=0, old_TFO=0, old_OTHER=0, rnew_down=0, rnew_far=0)
    for v in base.ands(z):
        f = z.keys[v]
        if f in base.FRx:
            continue
        if f in reg and reg[f] != "PI":
            res["old_" + reg[f]] += 1
        elif f in codeK:
            res["code"] += 1
        elif f in encK:
            res["enc"] += 1
        elif f in decK:
            res["dec"] += 1
        else:
            res["rnew_down" if v in downZ else "rnew_far"] += 1
    miss = dict(C=0, G=0, TFO=0, OTHER=0, rmade_down=0, rmade_far=0)
    for v in base.ands(Rx):
        f = Rx.keys[v]
        if f in Fz:
            continue
        r = reg.get(f)
        if r in miss:
            miss[r] += 1
        else:
            miss["rmade_down" if v in downRx else "rmade_far"] += 1
    row.update({f"res_{a}": b for a, b in res.items()})
    row.update({f"miss_{a}": b for a, b in miss.items()})
    row["res_total"] = sum(res.values())
    row["miss_total"] = sum(miss.values())
    # code bits
    nL = len(L["leaves"])
    tts = local_truth(y, L["leaves"], L["code"], nL) if L["code"] else []
    cb = []
    for l, t in zip(L["code"], tts):
        if l < 0:
            continue
        f = y.keys[l >> 1]
        cb.append(dict(reused=f in base.Fx, survive=f in Fz,
                       deg=degree(t, nL) if t is not None else None,
                       sup=support(t, nL) if t is not None else None))
    row["code_bits"] = len(cb)
    row["code_reused"] = sum(c["reused"] for c in cb)
    row["code_foreign_survive"] = sum(c["survive"] and not c["reused"] for c in cb)
    degs = [c["deg"] for c in cb if c["deg"] is not None]
    row["code_affine"] = sum(dg == 1 for dg in degs)
    row["code_mean_deg"] = round(sum(degs) / len(degs), 2) if degs else ""
    sups = [c["sup"] for c in cb if c["sup"] is not None]
    row["code_mean_sup"] = round(sum(sups) / len(sups), 2) if sups else ""
    # survival of the new E / D logic as a whole
    newE = {f for f in encK | codeK if f not in base.Fx}
    newD = {f for f in decK if f not in base.Fx}
    row["enc_new"] = len(newE)
    row["enc_new_survive"] = len(newE & Fz)
    row["dec_new"] = len(newD)
    row["dec_new_survive"] = len(newD & Fz)
    return row
