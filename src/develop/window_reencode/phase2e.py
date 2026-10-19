#!/usr/bin/env python3
"""Phase 2e (diagnostic): null-kick baseline and residue analysis.

The null kick (sampler N) puts the window back unchanged: the kicked AIG is the same AIG, only its nodes
are renumbered the way every real kick renumbers them (outside the window's fanout first, then the
window, then its fanout). R is order dependent, so N measures how much of any kick's effect is just
node order. Setting = Phase 2c mode m1: raw EPFL start, R = compress2rs x3, K = 8, 200 kicks, seed 1
(the U ... P0 kicks are the Phase 2c kicks again).

Every kick is traced (lsv_reencode -W) and analyzed with residue.py: which functions of R(kick) are not
in R(x) and where they come from, which functions of R(x) are lost and where they were.

outputs (results/phase2e/): raw/<circuit>_<sampler>.csv (kicks + residue columns),
  phase2e_kicks.csv, phase2e_summary.csv (per circuit and sampler)
usage: phase2e.py [-j jobs] [--samplers ...] [--aggregate] [circuit ...]
"""
import argparse, concurrent.futures as cf, csv, glob, os, statistics as st, subprocess, sys

import phase2
import residue
from phase2c import CIRCUITS, EPFL, R, ROOT

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results", "phase2e")
SAMPLERS = ["N", "P0", "base", "U", "A", "W4"]


def run(job):
    c, m = job
    tag = f"{c}_{m}"
    kd = os.path.join(OUT, "kicks", tag)
    os.makedirs(kd, exist_ok=True)
    dst = os.path.join(OUT, "raw", tag + "_lsv.csv")
    p = subprocess.run([phase2.ABC, "-c", f'&r {EPFL}/{c}.aig; lsv_reencode -k -M {m} -K 8 -n 200 -S 1 '
                                          f'-R "{R}" -o {dst} -W {kd}'],
                       cwd=ROOT, capture_output=True, text=True)
    log = [l for l in p.stdout.splitlines() if l.startswith("lsv_reencode")]
    rows = list(csv.DictReader(open(dst))) if os.path.exists(dst) else []
    if rows:
        base = residue.Base(kd)
        for r in rows:
            k = r["kick"]
            r.update(residue.analyze(base, kd, int(k)))
            for f in (f"y{k}.aig", f"y{k}.lab", f"k{k}.aig"):
                if int(r["d_R"]) >= 0 or f.startswith("y"):
                    os.remove(os.path.join(kd, f))
        with open(os.path.join(OUT, "raw", tag + ".csv"), "w", newline="") as fo:
            w = csv.DictWriter(fo, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
    return tag, "\n".join(log) if log else p.stdout[-400:] + p.stderr[-400:]


def summarize(rows):
    out = []
    for c in CIRCUITS:
        for m in SAMPLERS:
            R_ = [r for r in rows if r["circuit"] == c and r["method"] == m]
            if not R_:
                continue
            d = [int(r["d_R"]) for r in R_]
            aRx = int(R_[0]["and_Rx"])
            med = lambda k: round(st.median(float(r[k]) for r in R_), 1)
            out.append(dict(circuit=c, method=m, kicks=len(R_), and_Rx=aRx,
                            frac_better=round(sum(x < 0 for x in d) / len(d), 3),
                            frac_same=round(sum(x == 0 for x in d) / len(d), 3),
                            frac_worse=round(sum(x > 0 for x in d) / len(d), 3),
                            best_dR=min(d), best_pct=round(100 * min(d) / aRx, 2),
                            worst_dR=max(d), med_dR=st.median(d),
                            sd_dR=round(st.pstdev(d), 2),
                            **{f"med_{k}": med(k) for k in
                               ("res_total", "res_code", "res_enc", "res_dec", "res_old_C", "res_old_G",
                                "res_old_TFO", "res_old_OTHER", "res_rnew_down", "res_rnew_far",
                                "miss_total", "miss_C", "miss_G", "miss_TFO", "miss_OTHER",
                                "miss_rmade_down", "miss_rmade_far")}))
    return out


def write(name, rows):
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(os.path.join(OUT, name), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-j", type=int, default=40)
    ap.add_argument("--samplers", nargs="*", default=SAMPLERS)
    ap.add_argument("--aggregate", action="store_true")
    ap.add_argument("circuits", nargs="*", default=CIRCUITS)
    a = ap.parse_args()
    os.makedirs(os.path.join(OUT, "raw"), exist_ok=True)
    size = {c: residue.Aig(f"{EPFL}/{c}.aig").nAnd for c in a.circuits}
    jobs = sorted([(c, m) for c in a.circuits for m in a.samplers], key=lambda j: -size[j[0]])
    if not a.aggregate:
        with cf.ThreadPoolExecutor(a.j) as ex:
            for tag, log in ex.map(run, jobs):
                print(f"{tag:<18} {log}", flush=True)
    rows = []
    for p in sorted(glob.glob(os.path.join(OUT, "raw", "*.csv"))):
        if p.endswith("_lsv.csv"):
            continue
        c, m = os.path.basename(p)[:-4].rsplit("_", 1)
        for r in csv.DictReader(open(p)):
            r["circuit"] = c
            rows.append(r)
    write("phase2e_kicks.csv", rows)
    write("phase2e_summary.csv", summarize(rows))
    print(len(rows), "kicks")


if __name__ == "__main__":
    sys.exit(main())
