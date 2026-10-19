#!/usr/bin/env python3
"""Phase 2c (diagnostic): re-encoding kicks on the unoptimized EPFL AIGs, and how much they change structure.

Start x0 = the original EPFL AIG (no &deepsyn). R = compress2rs x3 (phase2b.SCRIPTS["rs3"]). A kick
re-encodes, one after another, either one window (mode m1, 200 kicks) or windows until their cones add up
to 10% / 50% of AND(x0) (modes P10 / P50, 20 kicks each), then runs R; compared with R(x0).

Structure (structdiff.py), every R(kick) against R(x0):
  new_frac   fraction of ANDs of R(kick) whose function (up to complement) does not occur in R(x0)
  kept_frac  fraction of ANDs of R(x0) whose function still occurs in R(kick)
  xor / mux / frac_shared / depth, and new_lev (mean level / depth of the new nodes)
Reference rows (kind = ref) compare x0, R(x0), R*(x0) (R iterated 4x) and the Phase 0 best AIG (mostly
&deepsyn) with R(x0): the structural change that ordinary optimization makes.

outputs (results/phase2c/):
  raw/<tag>.csv, kicks/<tag>/{Rx,k<i>}.aig (only R(kick) with d_R < 0 are kept)
  phase2c_kicks.csv     every kick + structure columns
  phase2c_ref.csv       reference rows
  phase2c_summary.csv   per (circuit, sampler, mode)
usage: phase2c.py [-j jobs] [--modes m1 P10 P50] [--samplers ...] [--aggregate] [circuit ...]
"""
import argparse, concurrent.futures as cf, csv, glob, os, statistics as st, subprocess, sys

import phase2
from phase1 import start_aig
from phase2b import C2RS, SCRIPTS
from structdiff import Aig, compare

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
EPFL = os.path.abspath(os.path.join(ROOT, "..", "testcases", "aig", "epfl"))
OUT = os.path.join(HERE, "results", "phase2c")
CIRCUITS = ["adder", "bar", "log2", "multiplier", "sqrt", "square", "voter", "sin", "max", "priority",
            "int2float", "cavlc", "arbiter", "mem_ctrl", "i2c", "ctrl", "router", "dec"]
SAMPLERS = ["U", "A", "W4", "base", "P0"]
MODES = {"m1": ("-m 1", 200), "P10": ("-P 10", 20), "P50": ("-P 50", 20)}
R = SCRIPTS["rs3"]


def struct_row(A, path):
    B = Aig(path)
    B.simulate()
    s = B.stats()
    s.update(compare(A, B))
    return s


def run(job, K=8, out=None, extra=""):
    c, m, mode = job
    out = out or OUT
    tag = f"{c}_{m}_{mode}"
    kd = os.path.join(out, "kicks", tag)
    os.makedirs(kd, exist_ok=True)
    os.makedirs(os.path.join(out, "raw"), exist_ok=True)
    dst = os.path.join(out, "raw", tag + ".csv")
    flag, n = MODES[mode]
    p = subprocess.run([phase2.ABC, "-c", f'&r {EPFL}/{c}.aig; lsv_reencode -k -M {m} -K {K} -n {n} -S 1 {flag} {extra} '
                                          f'-R "{R}" -o {dst} -W {kd}'],
                       cwd=ROOT, capture_output=True, text=True)
    log = [l for l in p.stdout.splitlines() if l.startswith("lsv_reencode")]
    # structure of every R(kick) against R(x0); keep only the improving AIGs on disk
    rows = list(csv.DictReader(open(dst))) if os.path.exists(dst) else []
    if rows:
        A = Aig(os.path.join(kd, "Rx.aig"))
        A.simulate()
        for r in rows:
            f = os.path.join(kd, f"k{r['kick']}.aig")
            r.update(struct_row(A, f))
            if int(r["d_R"]) >= 0:
                os.remove(f)
        with open(os.path.join(out, "raw", tag + "_struct.csv"), "w", newline="") as fo:
            w = csv.DictWriter(fo, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
    return tag, "\n".join(log) if log else p.stdout[-400:] + p.stderr[-400:]


def refs(c):
    """x0, R(x0), R*(x0), Phase 0 best, each against R(x0)."""
    d = os.path.join(OUT, "ref", c)
    os.makedirs(d, exist_ok=True)
    cmd = f"&r {EPFL}/{c}.aig; &put; {C2RS}; {C2RS}; {C2RS}; write_aiger -s {d}/Rx.aig"
    cmd += "".join(f"; {C2RS}; {C2RS}; {C2RS}" for _ in range(3)) + f"; write_aiger -s {d}/Rstar.aig"
    subprocess.run([phase2.ABC, "-c", cmd], cwd=ROOT, capture_output=True, text=True)
    A = Aig(os.path.join(d, "Rx.aig"))
    A.simulate()
    out = []
    p0 = start_aig(c)
    for name, path in [("x0", f"{EPFL}/{c}.aig"), ("Rx", f"{d}/Rx.aig"), ("Rstar", f"{d}/Rstar.aig"),
                       ("phase0_best:" + p0[0], p0[3])]:
        r = dict(circuit=c, ref=name)
        r.update(struct_row(A, path))
        out.append(r)
    return out


def num(x):
    return float(x) if x not in ("", None) else None


def summarize(rows):
    out = []
    keys = sorted({(r["circuit"], r["method"], r["mode"]) for r in rows},
                  key=lambda k: (CIRCUITS.index(k[0]), SAMPLERS.index(k[1]), list(MODES).index(k[2])))
    for c, m, mode in keys:
        R_ = [r for r in rows if (r["circuit"], r["method"], r["mode"]) == (c, m, mode)]
        d = [int(r["d_R"]) for r in R_]
        better = [r for r in R_ if int(r["d_R"]) < 0]
        med = lambda xs: round(st.median(xs), 4) if xs else ""
        out.append(dict(
            circuit=c, method=m, mode=mode, kicks=len(R_), and_x0=R_[0]["and_x"], and_Rx=R_[0]["and_Rx"],
            lev_Rx=R_[0]["lev_Rx"], mean_nwin=round(st.mean(int(r["nwin"]) for r in R_), 1),
            med_kick_dAND=med([int(r["and_kick"]) - int(r["and_x"]) for r in R_]),
            frac_better=round(len(better) / len(R_), 3), frac_worse=round(sum(x > 0 for x in d) / len(R_), 3),
            undo=round(sum(r["undo"] == "1" for r in R_) / len(R_), 3),
            best_dR=min(d), best_dR_pct=round(100.0 * min(d) / int(R_[0]["and_Rx"]), 3), med_dR=med(d),
            med_dlev=med([int(r["lev_R"]) - int(r["lev_Rx"]) for r in R_]),
            med_new_frac=med([num(r["new_frac"]) for r in R_]),
            max_new_frac=max(num(r["new_frac"]) for r in R_),
            med_new_frac_better=med([num(r["new_frac"]) for r in better]),
            best_new_frac=[r["new_frac"] for r in R_ if int(r["d_R"]) == min(d)][0],
            med_dxor=med([int(r["xor"]) - int(R_[0]["xor_Rx"]) for r in R_]) if "xor_Rx" in R_[0] else "",
            secs=round(sum(float(r["secs"]) for r in R_), 1)))
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
    ap.add_argument("--modes", nargs="*", default=list(MODES))
    ap.add_argument("--samplers", nargs="*", default=SAMPLERS)
    ap.add_argument("--aggregate", action="store_true")
    ap.add_argument("circuits", nargs="*", default=CIRCUITS)
    a = ap.parse_args()
    os.makedirs(os.path.join(OUT, "raw"), exist_ok=True)
    size = {c: Aig(f"{EPFL}/{c}.aig").nAnd for c in a.circuits}
    cost = {"m1": 1, "P10": 1, "P50": 3}
    jobs = [(c, m, mode) for c in a.circuits for m in a.samplers for mode in a.modes]
    jobs.sort(key=lambda j: -size[j[0]] * cost[j[2]])
    if a.aggregate:
        jobs = []
    with cf.ThreadPoolExecutor(a.j) as ex:
        futs = [ex.submit(refs, c) for c in a.circuits]
        for tag, log in ex.map(run, jobs):
            print(f"{tag:<26} {log}", flush=True)
        ref = [r for f in futs for r in f.result()]
    write("phase2c_ref.csv", ref)
    rx = {r["circuit"]: r for r in ref if r["ref"] == "Rx"}
    rows = []
    for p in sorted(glob.glob(os.path.join(OUT, "raw", "*_struct.csv"))):
        c, m, mode = os.path.basename(p)[:-len("_struct.csv")].rsplit("_", 2)
        for r in csv.DictReader(open(p)):
            r.update(circuit=c, mode=mode, xor_Rx=rx[c]["xor"] if c in rx else "")
            rows.append(r)
    write("phase2c_kicks.csv", rows)
    write("phase2c_summary.csv", summarize(rows))
    print(len(rows), "kicks")


if __name__ == "__main__":
    sys.exit(main())
