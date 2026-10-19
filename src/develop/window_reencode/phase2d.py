#!/usr/bin/env python3
"""Phase 2d (diagnostic): larger windows with random class numbering.

Same as phase2c.py mode m1 (one window per kick, 200 kicks, R = compress2rs x3, raw EPFL start), with
K in {6, 8, 10, 12} for U (uniform random class -> code) and the control P0, on the circuits where
Phase 2c found escapes. Outputs: results/phase2d/K<K>/raw/*_struct.csv, results/phase2d/phase2d_summary.csv.
--budget B relaxes the AND budget of each synthesized circuit to B * (|C| + 10) (default 3; 0 = none) and
writes to results/phase2d_B<B>/. --nocec skips the per-kick cec against x (every window is still checked
exhaustively, E o D = G, and the best R(kick) of each job is cec'd); used for B = 0, where kicked AIGs
grow to ~6x and cec dominates the run time.
usage: phase2d.py [-j jobs] [--budget B] [--samplers U P0] [--nocec] [--aggregate] [circuit ...]
"""
import argparse, concurrent.futures as cf, csv, glob, os, statistics as st, sys

import phase2c

OUT = os.path.join(phase2c.HERE, "results", "phase2d")
CIRCUITS = ["bar", "router", "priority", "ctrl", "int2float", "cavlc", "i2c", "sin"]
KS = [6, 8, 10, 12]
SAMPLERS = ["U", "P0"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-j", type=int, default=40)
    ap.add_argument("--aggregate", action="store_true")
    ap.add_argument("--budget", type=int, default=3)
    ap.add_argument("--samplers", nargs="*", default=SAMPLERS)
    ap.add_argument("--nocec", action="store_true")
    ap.add_argument("circuits", nargs="*", default=CIRCUITS)
    a = ap.parse_args()
    out_dir = OUT if a.budget == 3 else OUT + f"_B{a.budget}"
    extra = f"-B {a.budget}" + (" -c" if a.nocec else "")
    jobs = [((c, m, "m1"), K) for c in a.circuits for m in a.samplers for K in KS]
    jobs.sort(key=lambda j: -j[1])
    if not a.aggregate:
        with cf.ThreadPoolExecutor(a.j) as ex:
            futs = {ex.submit(phase2c.run, j, K, os.path.join(out_dir, f"K{K}"), extra): (j, K) for j, K in jobs}
            for f in cf.as_completed(futs):
                tag, log = f.result()
                print(f"K{futs[f][1]:<3} {tag:<22} {log}", flush=True)
    out = []
    for c in a.circuits:
        for m in a.samplers:
            for K in KS:
                p = os.path.join(out_dir, f"K{K}", "raw", f"{c}_{m}_m1_struct.csv")
                if not os.path.exists(p):
                    continue
                R = list(csv.DictReader(open(p)))
                if not R:
                    continue
                d = [int(r["d_R"]) for r in R]
                aRx = int(R[0]["and_Rx"])
                out.append(dict(circuit=c, method=m, K=K, budget=a.budget, kicks=len(R), and_Rx=aRx,
                                mean_nL=round(st.mean(int(r["nL"]) for r in R), 1),
                                mean_nC=round(st.mean(int(r["nC"]) for r in R), 1),
                                frac_nL_eq_K=round(sum(int(r["nL"]) == K for r in R) / len(R), 3),
                                mean_tries=round(st.mean(int(r["tries"]) for r in R), 1),
                                med_N=st.median(int(r["N"]) for r in R),
                                med_kick_dAND=st.median(int(r["and_kick"]) - int(r["and_x"]) for r in R),
                                frac_better=round(sum(x < 0 for x in d) / len(R), 3),
                                frac_worse=round(sum(x > 0 for x in d) / len(R), 3),
                                best_dR=min(d), best_dR_pct=round(100.0 * min(d) / aRx, 2),
                                med_dR=st.median(d),
                                med_new_frac=round(st.median(float(r["new_frac"]) for r in R), 4),
                                secs=round(sum(float(r["secs"]) for r in R), 1)))
    with open(os.path.join(out_dir, "phase2d_summary.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)


if __name__ == "__main__":
    sys.exit(main())
