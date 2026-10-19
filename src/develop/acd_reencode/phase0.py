#!/usr/bin/env python3
"""Phase 0 of ../plan/acd_cut_reencoding.md: reproduce the S=2 min-register baseline.

For every circuit, T in {50%, recommended} and cut side in {source, sink}, run aig_mincut,
compare the register count with the plan's table, and dump every cut signal with its arrival
level and downstream depth.

usage: phase0.py [--epfl DIR] [--out DIR] [-j N]
outputs: <out>/phase0_summary.csv, <out>/phase0_cuts.csv, <out>/runs/<circuit>_T<T>_<side>/{cuts,stage_assign}.json
"""
import argparse
import concurrent.futures as cf
import csv
import json
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
BIN = os.path.join(HERE, "..", "aig_mincut", "aig_mincut")

# circuit: (depth, reg @ T=50%, recommended T, reg @ recommended T, unconstrained)  -- plan section 2
BASELINE = {
    "bar": (12, 259, 7, 131, 128),
    "hyp": (24801, 306, 13641, 257, 128),
    "sqrt": (5058, 166, 2782, 137, 64),
    "sin": (225, 86, 124, 66, 24),
    "mem_ctrl": (114, 1554, 63, 1444, 1193),
    "multiplier": (274, 221, 151, 190, 128),
    "log2": (444, 210, 245, 199, 32),
    "arbiter": (87, 513, 44, 513, 129),
    "max": (287, 258, 144, 258, 130),
    "dec": (3, 16, 2, 16, 8),
    "adder": (255, 193, 192, 161, 129),
    "square": (250, 192, 188, 159, 64),
    "div": (4372, 335, 3279, 299, 128),
    "priority": (250, 72, 188, 40, 8),
    "router": (54, 50, 41, 22, 3),
    "voter": (70, 86, 53, 16, 1),
    "cavlc": (16, 63, 12, 28, 10),
    "int2float": (16, 26, 12, 13, 7),
    "ctrl": (10, 27, 8, 10, 7),
    "i2c": (20, 169, 15, 145, 141),
}
# how the recommended T was chosen (see the S=2 sweep): knee of the curve, flat curve (= 50%),
# or an arbitrary 75% point on a roughly linear trade-off curve
REC_KIND = {c: "knee" for c in ("bar", "hyp", "sqrt", "sin", "mem_ctrl", "multiplier", "log2")}
REC_KIND.update({c: "flat" for c in ("arbiter", "max", "dec")})
REC_KIND.update({c: "linear75" for c in BASELINE if c not in REC_KIND})


def run(job):
    circuit, depth, tag, T, side, expected, epfl, out = job
    rdir = os.path.join(out, "runs", f"{circuit}_T{T}_{side}")
    os.makedirs(rdir, exist_ok=True)
    cmd = [BIN, "--stages", "2", "--period", str(T), "--cut-side", side, "--out-dir", rdir,
           os.path.join(epfl, circuit + ".aig")]
    p = subprocess.run(cmd, capture_output=True, text=True)
    txt = p.stdout
    regs = int(re.search(r"registers\s*: (\d+)", txt).group(1)) if p.returncode == 0 else None
    level = int(re.search(r"level (\d+);", txt).group(1))
    ands = int(re.search(r"\((\d+) live\)", txt).group(1))
    sd = [int(x) for x in re.search(r"stage depth\s*:([\d ]+)", txt).group(1).split()] if regs is not None else [None, None]
    lp_time = float(re.search(r"lp ([\d.]+)s", txt).group(1)) if regs is not None else None
    sigs = []
    if regs is not None:
        cut = json.load(open(os.path.join(rdir, "cuts.json")))["cuts"][0]["signals"]
        for s in cut:
            sigs.append({
                "circuit": circuit, "T_tag": tag, "T": T, "side": side, "node": s["node"],
                "aiger_var": s["aiger_var"], "type": "pi" if s["balancing"] else "and",
                "po_driver": int(s["po_driver"]), "arrival": s["arrival"], "downstream": s["downstream"],
                "upstream_slack": T - s["arrival"], "downstream_slack": T - s["downstream"],
                "global_level": s["level"], "level_to_po": s["level_to_po"],
                "same_stage_uses": s["same_stage_uses"],
            })
    row = {
        "circuit": circuit, "depth": depth, "live_ands": ands, "T_tag": tag, "T": T, "side": side,
        "rec_kind": REC_KIND[circuit], "regs": regs, "expected": expected,
        "match": int(regs == expected), "checks": "FAIL" if ("FAIL" in txt or p.returncode not in (0,)) else "ok",
        "stage0_depth": sd[0], "stage1_depth": sd[1],
        "area": (ands + 10 * regs) if regs is not None else None,
        "pi_signals": sum(1 for s in sigs if s["type"] == "pi"),
        "po_driver_signals": sum(1 for s in sigs if s["po_driver"]),
        "shared_signals": sum(1 for s in sigs if s["same_stage_uses"] > 0),
        "max_arrival": max((s["arrival"] for s in sigs), default=0),
        "max_downstream": max((s["downstream"] for s in sigs), default=0),
        "sigs_upslack_ge2": sum(1 for s in sigs if s["upstream_slack"] >= 2),
        "sigs_downslack_ge2": sum(1 for s in sigs if s["downstream_slack"] >= 2),
        "lp_time_s": lp_time,
    }
    if level != depth:
        row["checks"] = f"depth mismatch {level} vs {depth}"
    return row, sigs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--epfl", default=os.path.join(ROOT, "..", "testcases", "aig", "epfl"))
    ap.add_argument("--out", default=os.path.join(HERE, "results"))
    ap.add_argument("-j", type=int, default=os.cpu_count())
    a = ap.parse_args()
    subprocess.run(["make", "-s", "-C", os.path.dirname(BIN)], check=True)
    jobs = []
    for c, (depth, r50, trec, rrec, _) in BASELINE.items():
        t50 = (depth + 1) // 2
        for tag, T, exp in (("50%", t50, r50), ("rec", trec, rrec)):
            for side in ("source", "sink"):
                jobs.append((c, depth, tag, T, side, exp, a.epfl, a.out))
    rows, sigs = [], []
    with cf.ThreadPoolExecutor(a.j) as ex:
        for row, s in ex.map(run, jobs):
            rows.append(row)
            sigs += s
    with open(os.path.join(a.out, "phase0_summary.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    with open(os.path.join(a.out, "phase0_cuts.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(sigs[0]))
        w.writeheader(); w.writerows(sigs)
    bad = [r for r in rows if not r["match"] or r["checks"] != "ok"]
    print(f"{len(rows)} runs, {len(rows) - len(bad)} match the baseline table and pass all checks")
    for r in bad:
        print("  MISMATCH", r["circuit"], r["T_tag"], r["T"], r["side"], r["regs"], "vs", r["expected"], r["checks"])
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
