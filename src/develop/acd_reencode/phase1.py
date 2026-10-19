#!/usr/bin/env python3
"""Phase 1 of ../plan/acd_cut_reencoding.md: compressibility of the Phase 0 cuts (no circuit change).

usage: phase1.py [--circuits a,b,...] [--K 8,12,16] [--epfl DIR] [--out DIR] [-j N] [--timeout SEC]
reads  <out>/runs/<circuit>_T<T>_<side>/cuts.json   (produced by phase0.py)
writes <out>/phase1_groups.csv, <out>/phase1_summary.csv

Slack-feasible saving of a group (column saving_feasible / sum_saving_with_slack>0):
  upstream_slack >= 1 and downstream_slack >= 1 : |G| - r      (encoder before, decoder after the register)
  only downstream_slack >= 1                     : |G| - |S|    (functional dependence: decoder only)
  otherwise                                      : 0
grouping: "plan" = greedy leaf-overlap grouping of plan 3.2; "slack" = same, but only within the
slack classes A/B (see acd_phase1 --slack). Not part of the plan; reference metric.
Go/no-go (plan 3.4, evaluated at K=16 per (circuit, T, side)):
  sum_saving_with_slack>0 >= max(5% of cut, 4)  and  some group has upstream or downstream slack >= 2
"""
import argparse
import concurrent.futures as cf
import csv
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
BIN = os.path.join(HERE, "acd_phase1")
sys.path.insert(0, HERE)
from phase0 import BASELINE  # noqa: E402

GROUP_COLS = ["circuit", "T_tag", "T", "side", "K", "grouping", "group_id", "|G|", "|L|", "|Img|", "r", "saving",
              "fd_subset_size", "saving_fd", "upstream_slack", "downstream_slack", "r_equals_L",
              "saving_feasible", "members", "fd_subset", "leaves"]
SUM_COLS = ["circuit", "T_tag", "T", "side", "K", "grouping", "cut_size", "groups", "sum_saving", "sum_saving_fd_only",
            "sum_saving_with_slack>0", "groups_slack_ge2", "groups_r_equals_L", "global_log2N", "global_N",
            "saturated", "global_exhaustive", "go", "seconds"]


def run(job):
    circuit, tag, T, side, K, grouping, do_global, epfl, out, timeout = job
    cutfile = os.path.join(out, "runs", f"{circuit}_T{T}_{side}", "cut.txt")
    cmd = [BIN, "--K", str(K)] + (["--global"] if do_global else []) + (["--slack", str(T)] if grouping == "slack" else []) + [os.path.join(epfl, circuit + ".aig"), cutfile]
    t0 = time.time()
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return job, None, None, time.time() - t0, "timeout"
    if p.returncode:
        return job, None, None, time.time() - t0, p.stderr.strip()
    groups, glob = [], None
    for line in p.stdout.splitlines():
        f = line.split(",")
        if f[0] == "global":
            glob = {"global_N": int(f[1]), "global_log2N": int(f[2]), "saturated": int(f[3]), "global_exhaustive": int(f[5])}
        elif f[0] == "group":
            G, L, img, r, sav, fd, marr, mdown, req = map(int, f[2:11])
            up, down = T - marr, T - mdown
            feas = sav if (up >= 1 and down >= 1) else (G - fd if down >= 1 else 0)
            groups.append(dict(zip(GROUP_COLS, [circuit, tag, T, side, K, grouping, int(f[1]), G, L, img, r, sav, fd, G - fd,
                                                up, down, req, feas, f[11], f[12], f[13]])))
    return job, groups, glob, time.time() - t0, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--circuits", default=",".join(BASELINE))
    ap.add_argument("--K", default="8,12,16")
    ap.add_argument("--epfl", default=os.path.join(ROOT, "..", "testcases", "aig", "epfl"))
    ap.add_argument("--out", default=os.path.join(HERE, "results"))
    ap.add_argument("-j", type=int, default=os.cpu_count())
    ap.add_argument("--timeout", type=int, default=1800, help="per run; plan 9: stop and report beyond 30 min")
    a = ap.parse_args()
    subprocess.run(["make", "-s", "-C", HERE], check=True)
    Ks = [int(k) for k in a.K.split(",")]
    jobs = []
    for c in a.circuits.split(","):
        depth, _, trec, _, _ = BASELINE[c]
        for tag, T in (("50%", (depth + 1) // 2), ("rec", trec)):
            for side in ("source", "sink"):
                for K in Ks:
                    for grouping in ("plan", "slack"):
                        jobs.append((c, tag, T, side, K, grouping, K == max(Ks) and grouping == "plan",
                                     a.epfl, a.out, a.timeout))
    # skip (circuit, T, side) where phase0 found T infeasible (no cuts.json)
    jobs = [j for j in jobs if os.path.exists(os.path.join(a.out, "runs", f"{j[0]}_T{j[2]}_{j[3]}", "cuts.json"))]
    for c, tag, T, side, *_ in jobs:   # write each cut file once, before the parallel runs
        rdir = os.path.join(a.out, "runs", f"{c}_T{T}_{side}")
        sigs = json.load(open(os.path.join(rdir, "cuts.json")))["cuts"][0]["signals"]
        with open(os.path.join(rdir, "cut.txt"), "w") as f:
            for s in sigs:
                f.write(f"{s['node']} {s['arrival']} {s['downstream']}\n")
    groups, summary, failures, globs = [], [], [], {}
    with cf.ThreadPoolExecutor(a.j) as ex:
        for job, gs, glob, sec, err in ex.map(run, jobs):
            c, tag, T, side, K, grouping = job[:6]
            if err:
                failures.append((c, tag, side, K, round(sec), err))
                continue
            groups += gs
            if glob:
                globs[(c, tag, side)] = glob
            cut = sum(g["|G|"] for g in gs)
            feas = sum(g["saving_feasible"] for g in gs)
            ge2 = sum(1 for g in gs if g["upstream_slack"] >= 2 or g["downstream_slack"] >= 2)
            summary.append({"circuit": c, "T_tag": tag, "T": T, "side": side, "K": K, "grouping": grouping, "cut_size": cut,
                            "groups": len(gs), "sum_saving": sum(g["saving"] for g in gs),
                            "sum_saving_fd_only": sum(g["saving_fd"] for g in gs),
                            "sum_saving_with_slack>0": feas, "groups_slack_ge2": ge2,
                            "groups_r_equals_L": sum(g["r_equals_L"] for g in gs),
                            "go": int(feas >= max(0.05 * cut, 4) and ge2 > 0), "seconds": round(sec, 2)})
    for s in summary:
        s.update(globs.get((s["circuit"], s["T_tag"], s["side"]),
                           {"global_N": "", "global_log2N": "", "saturated": "", "global_exhaustive": ""}))
    order = {c: i for i, c in enumerate(BASELINE)}
    key = lambda r: (order[r["circuit"]], r["T_tag"], r["side"], r["K"], r["grouping"])
    summary.sort(key=key)
    groups.sort(key=lambda r: key(r) + (r["group_id"],))
    with open(os.path.join(a.out, "phase1_groups.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=GROUP_COLS); w.writeheader(); w.writerows(groups)
    with open(os.path.join(a.out, "phase1_summary.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=SUM_COLS); w.writeheader(); w.writerows(summary)
    for s in summary:
        if s["K"] == max(Ks):
            print(f"{s['circuit']:10s} {s['T_tag']:3s} T={s['T']:<5} {s['side']:6s} {s['grouping']:5s} cut={s['cut_size']:<5} "
                  f"saving={s['sum_saving']:<4} fd={s['sum_saving_fd_only']:<4} feasible={s['sum_saving_with_slack>0']:<4} "
                  f"global_log2N={s['global_log2N']}{' (sat)' if s['saturated'] == 1 else ''} go={s['go']} {s['seconds']}s")
    for fl in failures:
        print("FAILED", *fl)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
