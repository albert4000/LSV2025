#!/usr/bin/env python3
"""Phase 3a/3b report: seeds 1-3 of Phase 3a (ILS from the raw EPFL AIG vs &deepsyn) and Phase 3b (ILS
from the &deepsyn result vs &deepsyn continued for the same time). Reads the phase3a_summary.csv files
written by phase3a.py (run phase3a.py --aggregate first) and prints markdown tables.
usage: phase3_report.py
"""
import csv, os, statistics as st

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
CIRC = ["router", "priority", "ctrl", "bar", "voter", "sqrt", "sin", "multiplier", "square", "log2", "max",
        "i2c", "int2float", "cavlc", "arbiter", "mem_ctrl"]
ILS = ["U", "A", "P0", "U-e", "A-e", "P0-e"]


def load(d):
    p = os.path.join(RES, d, "phase3a_summary.csv")
    out = {}
    if os.path.exists(p):
        for r in csv.DictReader(open(p)):
            if r["final"] != "":
                out[(r["circuit"], r["config"])] = r
    return out


def main():
    seeds = {s: load(d) for s, d in ((1, "phase3a"), (2, "phase3a_S2"), (3, "phase3a_S3"))}
    print("## Phase 3a, seeds 1-3 (start = raw EPFL AIG, 20 min each)\n")
    print("| circuit | R(x0) | ILS best S1/S2/S3 | deepsyn S1/S2/S3 | ILS median of per-seed best | deepsyn median | wins |")
    print("|---|---|---|---|---|---|---|")
    cfg_tot = {c: [] for c in ILS}
    for c in CIRC:
        ib, db = [], []
        for s, S in seeds.items():
            v = [int(S[(c, k)]["final"]) for k in ILS if (c, k) in S]
            ib.append(min(v) if len(v) == len(ILS) else None)
            db.append(int(S[(c, "deepsyn")]["final"]) if (c, "deepsyn") in S else None)
        rx = next((S[(c, "U")]["and_Rx"] for S in seeds.values() if (c, "U") in S), "")
        wins = sum(1 for a, b in zip(ib, db) if a is not None and b is not None and a < b)
        n = sum(1 for a, b in zip(ib, db) if a is not None and b is not None)
        f = lambda xs: "/".join(str(x) if x is not None else "-" for x in xs)
        m = lambda xs: st.median([x for x in xs if x is not None]) if any(x is not None for x in xs) else "-"
        print(f"| {c} | {rx} | {f(ib)} | {f(db)} | {m(ib)} | {m(db)} | {wins}/{n} |")
        # per config: relative to deepsyn of the same seed
        for k in ILS:
            for s, S in seeds.items():
                if (c, k) in S and (c, "deepsyn") in S:
                    cfg_tot[k].append(int(S[(c, k)]["final"]) / int(S[(c, "deepsyn")]["final"]))
    print("\nper ILS config (one config fixed for every circuit): geometric-mean AND ratio vs &deepsyn of the same seed")
    for k, v in cfg_tot.items():
        if v:
            g = 1.0
            for x in v:
                g *= x
            print(f"  {k:5s} n={len(v):3d}  ratio {g ** (1 / len(v)):.3f}  wins {sum(x < 1 for x in v)}  ties {sum(x == 1 for x in v)}")

    print("\n## Phase 3b (start = &deepsyn S1 result; ILS seed 2 vs &deepsyn continued with seed 2, 20 min each)\n")
    B = load("phase3b")
    cfgs = ILS + ["N", "deepsyn"]
    print("| circuit | start | " + " | ".join(cfgs) + " | best ILS vs start |")
    print("|---|---|" + "---|" * len(cfgs) + "---|")
    for c in CIRC:
        if not any((c, k) in B for k in cfgs):
            continue
        x0 = next(B[(c, k)]["and_x0"] for k in cfgs if (c, k) in B)
        vals = [B[(c, k)]["final"] if (c, k) in B else "-" for k in cfgs]
        iv = [int(B[(c, k)]["final"]) for k in ILS if (c, k) in B]
        d = f"{min(iv) - int(x0):+d} ({100 * (min(iv) - int(x0)) / int(x0):+.2f}%)" if iv else "-"
        print(f"| {c} | {x0} | " + " | ".join(vals) + f" | {d} |")


if __name__ == "__main__":
    main()
