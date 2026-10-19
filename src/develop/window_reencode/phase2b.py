#!/usr/bin/env python3
"""Phase 2b (diagnostic, outside the plan's go/no-go): Phase 2 with a stronger re-optimization script R.

Phase 2 found that a single dc2 neither undoes nor exploits most kicks. This reruns the single-kick
experiment of phase2.py on the circuits where Phase 2 escaped (d_R < 0), with R replaced by
  dc2x5   "&put; dc2 x5; &get -n"               (the same operator, closer to convergence)
  rs3     "&put; compress2rs x3; &get -n"       (adds resub -K 6..12 -N 2, i.e. multi-node resubstitution;
                                                 compress2rs as in &deepsyn, giaDeep.c)
One &deepsyn iteration (&dch -f; &if -a -K 6; &mfs; compress2rs) was dropped: on the start AIGs it is
4-6% worse than x and takes up to 65 s per call, i.e. it is a restart, not a local optimizer.
R(x) is recomputed with the same R, so d_R compares like with like; thr is the Phase 2 threshold.

outputs (results/phase2b/<R>/): raw/, best/, phase2_kicks.csv, phase2_summary.csv, phase2_go.csv
usage: phase2b.py [-j jobs] [-n kicks] [-K K ...] [-R name ...] [--aggregate] [circuit ...]
"""
import argparse, concurrent.futures as cf, csv, os, subprocess, sys

import phase2
from phase1 import start_aig

C2RS = ("balance -l; resub -K 6 -l; rewrite -l; resub -K 6 -N 2 -l; refactor -l; resub -K 8 -l; balance -l; "
        "resub -K 8 -N 2 -l; rewrite -l; resub -K 10 -l; rewrite -z -l; resub -K 10 -N 2 -l; balance -l; "
        "resub -K 12 -l; refactor -z -l; resub -K 12 -N 2 -l; rewrite -z -l; balance -l")
SCRIPTS = {"dc2x5": "&put; " + "dc2; " * 5 + "&get -n",
           "rs3": "&put; " + (C2RS + "; ") * 3 + "&get -n"}
CIRCUITS = ["log2", "sin", "sqrt", "priority"]


def run(job):
    c, m, K, R, path, n = job
    tag = f"{c}_{m}_K{K}"
    out = os.path.join(phase2.HERE, "results", "phase2b", R)
    dst, best = os.path.join(out, "raw", tag + ".csv"), os.path.join(out, "best", tag + ".aig")
    p = subprocess.run([phase2.ABC, "-c", f'&r {path}; lsv_reencode -k -M {m} -K {K} -n {n} -S 1 '
                                          f'-R "{SCRIPTS[R]}" -o {dst} -w {best}'],
                       cwd=phase2.ROOT, capture_output=True, text=True)
    log = [l for l in p.stdout.splitlines() if l.startswith("lsv_reencode")]
    return f"{R}/{tag}", "\n".join(log) if log else p.stdout[-400:] + p.stderr[-400:]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-j", type=int, default=40)
    ap.add_argument("-n", type=int, default=200)
    ap.add_argument("-K", type=int, nargs="*", default=[8])
    ap.add_argument("-R", nargs="*", default=list(SCRIPTS))
    ap.add_argument("--aggregate", action="store_true")
    ap.add_argument("circuits", nargs="*", default=CIRCUITS)
    a = ap.parse_args()
    phase2.KS = a.K
    for R in a.R:
        for d in ("raw", "best"):
            os.makedirs(os.path.join(phase2.HERE, "results", "phase2b", R, d), exist_ok=True)
    starts = {c: start_aig(c) for c in a.circuits}
    jobs = [(c, m, K, R, starts[c][3], a.n) for R in a.R for c in a.circuits for m in phase2.SAMPLERS for K in a.K]
    jobs.sort(key=lambda j: -starts[j[0]][1] * (3 if j[3] == "rs3" else 1))   # slowest first
    if a.aggregate:
        jobs = []
    with cf.ThreadPoolExecutor(a.j) as ex:
        for tag, log in ex.map(run, jobs):
            print(f"{tag:<28} {log}", flush=True)
    sig = phase2.sigma()
    for R in a.R:
        phase2.OUT = os.path.join(phase2.HERE, "results", "phase2b", R)
        rows = []
        for c in a.circuits:
            for m in phase2.SAMPLERS:
                for K in a.K:
                    for r in csv.DictReader(open(os.path.join(phase2.OUT, "raw", f"{c}_{m}_K{K}.csv"))):
                        r["circuit"] = c
                        rows.append(r)
        phase2.write("phase2_kicks.csv", rows)
        summ = phase2.summarize(rows, a.circuits, sig)
        phase2.write("phase2_summary.csv", summ)
        phase2.write("phase2_go.csv", phase2.go(summ))
        print(R, len(rows), "kicks")


if __name__ == "__main__":
    sys.exit(main())
