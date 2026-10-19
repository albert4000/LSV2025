#!/usr/bin/env python3
"""Phase 2 of src/develop/plan/window_reencode_pareto.md: single kicks followed by R (go/no-go #1).

Start point per circuit: the Phase 1 start AIG x (smallest cec-verified AIG of Phase 0). For every
(circuit, sampler, K), `lsv_reencode -k` makes 200 independent kicks from x: pick a random window with
N >= 5 classes, re-encode it with the sampler (U, A, W1, W4; base and P0 are controls), cec the kicked
AIG against x, run R = "&put; dc2; &get -n", and compare with R(x).

  undo    R(kick(x)) has the same (AND, level) as R(x)
  noop    the kick itself left (AND, level) unchanged (structural hashing rebuilt the same logic)
  d_R     AND(R(kick(x))) - AND(R(x)); < 0 is an escape to a smaller AIG
  thr     max(sigma, 0.5% of AND(R(x))), sigma = &deepsyn seed spread at T=600 (Phase 0)

Go (per sampler, plan section 5): undo rate < 80%, and d_R < 0 on >= 3 circuits with best -d_R > thr.

outputs (results/phase2/):
  raw/<circuit>_<sampler>_K<K>.csv, best/<circuit>_<sampler>_K<K>.aig (best R(kick), cec-verified)
  phase2_kicks.csv     all kicks
  phase2_summary.csv   per (circuit, sampler, K)
  phase2_go.csv        go/no-go per (sampler, K)
usage: phase2.py [-j jobs] [-n kicks] [--only TAG ...] [--aggregate] [circuit ...]
  --only <circuit>_<sampler>_K<K>  rerun just these jobs; --aggregate  only rebuild the CSV summaries
"""
import argparse, concurrent.futures as cf, csv, os, statistics as st, subprocess, sys

from phase1 import CIRCUITS, start_aig

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
ABC = os.path.join(ROOT, "abc")
P0 = os.path.join(HERE, "results", "phase0")
OUT = os.path.join(HERE, "results", "phase2")
SAMPLERS = ["U", "A", "W1", "W4", "base", "P0"]
CONTROLS = {"base", "P0"}
KS = [6, 8]


def run(job):
    c, m, K, path, n = job
    tag = f"{c}_{m}_K{K}"
    dst, best = os.path.join(OUT, "raw", tag + ".csv"), os.path.join(OUT, "best", tag + ".aig")
    p = subprocess.run([ABC, "-c", f"&r {path}; lsv_reencode -k -M {m} -K {K} -n {n} -S 1 -o {dst} -w {best}"],
                       cwd=ROOT, capture_output=True, text=True)
    log = [l for l in p.stdout.splitlines() if l.startswith("lsv_reencode")]
    return tag, "\n".join(log) if log else p.stdout[-400:] + p.stderr[-400:]


def sigma():
    s = {}
    for r in csv.DictReader(open(os.path.join(P0, "phase0_noise.csv"))):
        if r["T"] == "600":
            s[r["circuit"]] = float(r["sigma"])
    return s


def summarize(rows, circuits, sig):
    out = []
    for c in circuits:
        for m in SAMPLERS:
            for K in KS:
                R = [r for r in rows if r["circuit"] == c and r["method"] == m and int(r["K"]) == K]
                if not R:
                    continue
                d = [int(r["d_R"]) for r in R]
                real = [r for r in R if (r["and_kick"], r["lev_kick"]) != (r["and_x"], r["lev_x"])]
                aRx = int(R[0]["and_Rx"])
                thr = max(sig.get(c, 0.0), 0.005 * aRx)
                out.append(dict(
                    circuit=c, method=m, K=K, kicks=len(R), and_x=R[0]["and_x"], and_Rx=aRx, lev_Rx=R[0]["lev_Rx"],
                    thr=round(thr, 1),
                    noop=round(1 - len(real) / len(R), 3),
                    undo=round(sum(r["undo"] == "1" for r in R) / len(R), 3),
                    undo_real=round(sum(r["undo"] == "1" for r in real) / len(real), 3) if real else "",
                    kick_cost_med=st.median(int(r["and_kick"]) - int(r["and_x"]) for r in real) if real else "",
                    kick_dlev_med=st.median(int(r["lev_kick"]) - int(r["lev_x"]) for r in real) if real else "",
                    frac_better=round(sum(x < 0 for x in d) / len(R), 3),
                    frac_better_lev_ok=round(sum(int(r["d_R"]) < 0 and int(r["lev_R"]) <= int(r["lev_Rx"]) for r in R) / len(R), 3),
                    frac_worse=round(sum(x > 0 for x in d) / len(R), 3),
                    best_dR=min(d), median_dR=st.median(d),
                    best_dR_pct=round(100.0 * min(d) / aRx, 3),
                    beats_thr=int(-min(d) > thr),
                    secs=round(sum(float(r["secs"]) for r in R), 1)))
    return out


def go(summ):
    out = []
    for m in SAMPLERS:
        for K in KS:
            S = [s for s in summ if s["method"] == m and s["K"] == K]
            if not S:
                continue
            kicks = sum(s["kicks"] for s in S)
            undo = sum(s["undo"] * s["kicks"] for s in S) / kicks
            esc = [s["circuit"] for s in S if s["best_dR"] < 0]
            sig = [s["circuit"] for s in S if s["beats_thr"]]
            out.append(dict(method=m, K=K, control=int(m in CONTROLS), kicks=kicks, undo=round(undo, 3),
                            circuits_escape=len(esc), circuits_beat_thr=len(sig),
                            go=int(undo < 0.8 and len(sig) >= 3), go_without_thr=int(undo < 0.8 and len(esc) >= 3),
                            escape_list=" ".join(esc), beat_thr_list=" ".join(sig)))
    return out


def write(name, rows):
    with open(os.path.join(OUT, name), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-j", type=int, default=40)
    ap.add_argument("-n", type=int, default=200)
    ap.add_argument("--only", nargs="*", default=[])
    ap.add_argument("--aggregate", action="store_true")
    ap.add_argument("circuits", nargs="*", default=CIRCUITS)
    a = ap.parse_args()
    for d in ("raw", "best"):
        os.makedirs(os.path.join(OUT, d), exist_ok=True)
    starts = {c: start_aig(c) for c in a.circuits}
    jobs = [(c, m, K, starts[c][3], a.n) for c in a.circuits for m in SAMPLERS for K in KS]
    jobs.sort(key=lambda j: -starts[j[0]][1])   # big circuits first
    if a.only:
        jobs = [j for j in jobs if f"{j[0]}_{j[1]}_K{j[2]}" in a.only]
    if a.aggregate:
        jobs = []
    with cf.ThreadPoolExecutor(a.j) as ex:
        for tag, log in ex.map(run, jobs):
            print(f"{tag:<24} {log}", flush=True)
    rows = []
    for c in a.circuits:
        for m in SAMPLERS:
            for K in KS:
                for r in csv.DictReader(open(os.path.join(OUT, "raw", f"{c}_{m}_K{K}.csv"))):
                    r["circuit"] = c
                    rows.append(r)
    write("phase2_kicks.csv", rows)
    summ = summarize(rows, a.circuits, sigma())
    write("phase2_summary.csv", summ)
    write("phase2_go.csv", go(summ))
    print(len(rows), "kicks")


if __name__ == "__main__":
    sys.exit(main())
