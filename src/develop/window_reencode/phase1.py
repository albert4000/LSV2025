#!/usr/bin/env python3
"""Phase 1 of src/develop/plan/window_reencode_pareto.md: window statistics and direct gain of re-encoding.

Start point per circuit: the smallest cec-verified AIG of Phase 0 (ties -> fewer levels), i.e. a circuit
that is already at ABC's local minimum. For every K, `lsv_reencode -a` samples up to 1000 random window
roots and costs, per window, P0 (direct resynthesis), base (greedy assignCode), U, A, W1, W4 (16 codes
each) and X (all canonical codes, N <= 8, first 30 such windows).

outputs (results/phase1/):
  raw/<circuit>_K<K>.csv   lsv_reencode output
  phase1_windows.csv       all rows
  phase1_summary.csv       per (circuit, K): window statistics + per-sampler gain statistics
  phase1_methods.csv       per (K, sampler) over all circuits
usage: phase1.py [-j jobs] [-n roots] [circuit ...]
"""
import argparse, concurrent.futures as cf, csv, os, statistics as st, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
ABC = os.path.join(ROOT, "abc")
P0 = os.path.join(HERE, "results", "phase0")
OUT = os.path.join(HERE, "results", "phase1")
CIRCUITS = ["adder", "bar", "log2", "multiplier", "sqrt", "square", "voter", "sin", "max", "priority",   # group B
            "int2float", "cavlc"]                                                                         # group A controls
KS = [6, 8, 10, 12]
METHODS = ["P0", "base", "U", "A", "W1", "W4", "X"]


def start_aig(c):
    rows = [r for r in csv.DictReader(open(os.path.join(P0, "phase0_baseline.csv")))
            if r["circuit"] == c and r["cec"] == "ok" and r["relax"] == "0" and r["flow"] == "nf"]
    best = min(rows, key=lambda r: (int(r["ands"]), int(r["levels"])))
    return best["variant"], int(best["ands"]), int(best["levels"]), os.path.join(P0, "aig", c, best["variant"] + ".aig")


def run(job):
    c, K, path, n = job
    dst = os.path.join(OUT, "raw", f"{c}_K{K}.csv")
    p = subprocess.run([ABC, "-c", f"&r {path}; lsv_reencode -a -K {K} -n {n} -t 16 -S 1 -X 30 -o {dst}"],
                       cwd=ROOT, capture_output=True, text=True)
    log = [l for l in p.stdout.splitlines() if l.startswith("lsv_reencode")]
    return c, K, dst, log[-1] if log else p.stdout[-300:]


def num(x):
    return float(x) if x not in ("", None) else None


def summarize(rows, circuits):
    summ = []
    for c in circuits:
        for K in KS:
            R = [r for r in rows if r["circuit"] == c and int(r["K"]) == K]
            if not R:
                continue
            wins = {}
            for r in R:
                wins.setdefault(r["win"], r)
            W = list(wins.values())
            N = [int(w["N"]) for w in W]
            s = dict(circuit=c, K=K, windows=len(W),
                     mean_nL=round(st.mean(int(w["nL"]) for w in W), 1),
                     mean_nC=round(st.mean(int(w["nC"]) for w in W), 1),
                     mean_nG=round(st.mean(int(w["nG"]) for w in W), 1),
                     median_N=st.median(N),
                     frac_N_le2=round(sum(x <= 2 for x in N) / len(W), 3),
                     frac_N_ge5=round(sum(x >= 5 for x in N) / len(W), 3),
                     frac_G_gt_r=round(sum(int(w["nG"]) > int(w["r"]) for w in W) / len(W), 3),
                     frac_injective=round(sum(int(w["N"]) == 2 ** int(w["nL"]) for w in W) / len(W), 3),
                     frac_trunc=round(sum(w["trunc"] == "1" for w in W) / len(W), 3),
                     frac_shrunk=round(sum(w.get("shrunk") == "1" for w in W) / len(W), 3))
            for m in METHODS:
                M = [r for r in R if r["method"] == m]
                ok = [r for r in M if r["best_gain"] != ""]
                s[f"{m}_wins"] = len(M)
                s[f"{m}_pos"] = round(sum(num(r["best_gain"]) > 0 for r in ok) / len(M), 3) if M else ""
                s[f"{m}_pos_safe"] = round(sum(num(r["best_gain"]) > 0 and num(r["best_ddep"]) == 0 for r in ok) / len(M), 3) if M else ""
                s[f"{m}_med_best"] = st.median(num(r["best_gain"]) for r in ok) if ok else ""
                s[f"{m}_fail"] = round(sum(int(r["fail"]) for r in M) / max(1, sum(int(r["samples"]) for r in M)), 3) if M else ""
            summ.append(s)
    return summ


def by_method(rows):
    out = []
    for K in KS:
        for m in METHODS:
            M = [r for r in rows if int(r["K"]) == K and r["method"] == m]
            if not M:
                continue
            ok = [r for r in M if r["best_gain"] != ""]
            p0 = {(r["circuit"], r["win"]): num(r["best_gain"]) for r in rows
                  if int(r["K"]) == K and r["method"] == "P0" and r["best_gain"] != ""}
            beats = [r for r in ok if (r["circuit"], r["win"]) in p0 and num(r["best_gain"]) > p0[(r["circuit"], r["win"])]]
            out.append(dict(K=K, method=m, windows=len(M),
                            frac_best_gain_pos=round(sum(num(r["best_gain"]) > 0 for r in ok) / len(M), 4),
                            frac_best_gain_pos_safe=round(sum(num(r["best_gain"]) > 0 and num(r["best_ddep"]) == 0 for r in ok) / len(M), 4),
                            frac_best_gain_ge0=round(sum(num(r["best_gain"]) >= 0 for r in ok) / len(M), 4),
                            frac_beats_P0=round(len(beats) / len(M), 4),
                            median_best_gain=st.median(num(r["best_gain"]) for r in ok) if ok else "",
                            median_mean_gain=round(st.median(num(r["mean_gain"]) for r in ok), 2) if ok else "",
                            mean_frac_pos_samples=round(st.mean(num(r["frac_pos"]) for r in ok), 4) if ok else "",
                            median_min_dlev=st.median(num(r["min_dlev"]) for r in ok) if ok else "",
                            fail_rate=round(sum(int(r["fail"]) for r in M) / max(1, sum(int(r["samples"]) for r in M)), 4),
                            total_gain_pos=sum(num(r["best_gain"]) for r in ok if num(r["best_gain"]) > 0)))
    return out


def write(name, rows):
    with open(os.path.join(OUT, name), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-j", type=int, default=36)
    ap.add_argument("-n", type=int, default=1000)
    ap.add_argument("circuits", nargs="*", default=CIRCUITS)
    a = ap.parse_args()
    os.makedirs(os.path.join(OUT, "raw"), exist_ok=True)
    starts = {c: start_aig(c) for c in a.circuits}
    for c, (v, ands, lev, _) in starts.items():
        print(f"start {c:<10} {v:<18} and={ands} lev={lev}")
    jobs = [(c, K, starts[c][3], a.n) for c in a.circuits for K in KS]
    with cf.ThreadPoolExecutor(a.j) as ex:
        for c, K, dst, log in ex.map(run, jobs):   # log reports duplicate / no-window / cyclic roots
            print(f"{c:<10} K={K:<2} {log}", flush=True)
    rows = []
    for c in a.circuits:
        for K in KS:
            for r in csv.DictReader(open(os.path.join(OUT, "raw", f"{c}_K{K}.csv"))):
                r["circuit"] = c
                r["start"] = starts[c][0]
                rows.append(r)
    write("phase1_windows.csv", rows)
    write("phase1_summary.csv", summarize(rows, a.circuits))
    write("phase1_methods.csv", by_method(rows))
    print(len(rows), "rows")


if __name__ == "__main__":
    sys.exit(main())
