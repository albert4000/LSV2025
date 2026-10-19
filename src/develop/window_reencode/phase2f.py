#!/usr/bin/env python3
"""Phase 2f: how many code bits survive R, and what happens if R must keep them (lsv_reencode -b).

Setting of Phase 2e / 2c mode m1: raw EPFL start, R = compress2rs x3, K = 8, 200 single kicks, seed 1, samplers
U, A, W4, base. Each job runs twice, without and with -b; -b does not use the random generator, so both
runs make the same kicks (same root, window and code) and differ only in R:
  normal  R on the kicked AIG (Phase 2e)
  keep    R with every code bit cut (exported as a CO, read back as a CI), reconnected afterwards: E and D are
          optimized separately and the r code lines stay
Columns from lsv_reencode: code_bits (distinct non-constant code functions), code_fresh (not a function of
x0), code_surv / code_fresh_surv (still computed by some node of R(kick); 2048-bit simulation signatures),
ncut (cut code nodes, -b).
Uses ./abc_dev (built with `make PROG=abc_dev`), so that ./abc can stay in use by other experiments.

outputs (results/phase2f/): raw/<circuit>_<sampler>_<mode>.csv, phase2f_kicks.csv, phase2f_summary.csv,
  phase2f_paired.csv (per circuit and sampler: d_R keep - d_R normal over the same kicks)
usage: phase2f.py [-j jobs] [--aggregate] [circuit ...]
"""
import argparse, concurrent.futures as cf, csv, glob, os, statistics as st, subprocess, sys

from phase2c import CIRCUITS, EPFL, R, ROOT

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results", "phase2f")
ABC = os.path.join(ROOT, "abc_dev")
SAMPLERS = ["U", "A", "W4", "base"]
MODES = {"normal": "", "keep": "-b"}


def run(job):
    c, m, mode = job
    dst = os.path.join(OUT, "raw", f"{c}_{m}_{mode}.csv")
    p = subprocess.run([ABC, "-c", f'&r {EPFL}/{c}.aig; lsv_reencode -k {MODES[mode]} -M {m} -K 8 -n 200 -S 1 '
                                   f'-R "{R}" -o {dst}'], cwd=ROOT, capture_output=True, text=True)
    log = [l for l in p.stdout.splitlines() if l.startswith("lsv_reencode")]
    return f"{c}_{m}_{mode}", "\n".join(log) if log else p.stdout[-400:] + p.stderr[-400:]


def write(name, rows):
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(os.path.join(OUT, name), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-j", type=int, default=40)
    ap.add_argument("--aggregate", action="store_true")
    ap.add_argument("circuits", nargs="*", default=CIRCUITS)
    a = ap.parse_args()
    os.makedirs(os.path.join(OUT, "raw"), exist_ok=True)
    if not a.aggregate:
        jobs = [(c, m, mode) for c in a.circuits for m in SAMPLERS for mode in MODES]
        with cf.ThreadPoolExecutor(a.j) as ex:
            for tag, log in ex.map(run, jobs):
                print(f"{tag:<24} {log}", flush=True)
    rows, summ, paired = [], [], []
    for c in a.circuits:
        for m in SAMPLERS:
            runs = {}
            for mode in MODES:
                p = os.path.join(OUT, "raw", f"{c}_{m}_{mode}.csv")
                if os.path.exists(p):
                    runs[mode] = list(csv.DictReader(open(p)))
                    for r in runs[mode]:
                        r.update(circuit=c, mode=mode)
                    rows += runs[mode]
            for mode, R_ in runs.items():
                d = [int(r["d_R"]) for r in R_]
                tot = lambda k: sum(int(r[k]) for r in R_)
                summ.append(dict(circuit=c, method=m, mode=mode, kicks=len(R_), and_Rx=R_[0]["and_Rx"],
                                 code_bits=tot("code_bits"), code_fresh=tot("code_fresh"),
                                 surv_frac=round(tot("code_surv") / max(1, tot("code_bits")), 3),
                                 fresh_surv_frac=round(tot("code_fresh_surv") / max(1, tot("code_fresh")), 3),
                                 kicks_all_fresh_gone=round(sum(int(r["code_fresh_surv"]) == 0 for r in R_) / len(R_), 3),
                                 frac_better=round(sum(x < 0 for x in d) / len(d), 3),
                                 frac_worse=round(sum(x > 0 for x in d) / len(d), 3),
                                 best_dR=min(d), med_dR=st.median(d)))
            if len(runs) == 2:
                k = {r["kick"]: r for r in runs["keep"]}
                diff = [int(k[r["kick"]]["d_R"]) - int(r["d_R"]) for r in runs["normal"] if r["kick"] in k]
                same = all(k[r["kick"]]["root"] == r["root"] for r in runs["normal"] if r["kick"] in k)
                paired.append(dict(circuit=c, method=m, kicks=len(diff), same_kicks=same,
                                   med_keep_minus_normal=st.median(diff) if diff else "",
                                   frac_keep_better=round(sum(x < 0 for x in diff) / len(diff), 3) if diff else "",
                                   frac_keep_worse=round(sum(x > 0 for x in diff) / len(diff), 3) if diff else ""))
    write("phase2f_kicks.csv", rows)
    write("phase2f_summary.csv", summ)
    write("phase2f_paired.csv", paired)
    print(len(rows), "kicks")


if __name__ == "__main__":
    sys.exit(main())
