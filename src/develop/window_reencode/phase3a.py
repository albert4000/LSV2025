#!/usr/bin/env python3
"""Phase 3a: iterated local search (lsv_reencode -i) against &deepsyn, area (AND count) only.

Start x0 = the original EPFL AIG, R = compress2rs x3 (as Phase 2c/2e). Each run: cur = R(x0); repeat
{ kick cur with one window (K = 8), run R, keep the result if it has fewer ANDs (or as many and fewer
levels) } for T seconds (wall clock). Samplers: U, A, P0 (window resynthesis without the code
bottleneck), N (null kick: only the node order changes, i.e. R rerun on a reshuffled AIG); "-e" variants
also accept results with as many ANDs and no more levels (plateau moves). Reference: &deepsyn -T T from
x0 (its improvements over time from -v).

outputs (results/phase3a/):
  raw/<circuit>_<config>.csv    one row per ILS iteration (lsv_reencode -i -o)
  raw/<circuit>_deepsyn.txt     &deepsyn -v log
  best/<circuit>_<config>.aig   final AIG of each run (cec-checked by lsv_reencode)
  phase3a_curve.csv             best AND at checkpoints (seconds) per circuit and config
  phase3a_summary.csv           per circuit and config: final AND, vs R(x0), iterations, accepted
--start: start AIG per circuit as a template with {c} (default: the EPFL AIG); --out: output directory
(default results/phase3a). Phase 3b = --start results/phase3a/best/{c}_deepsyn.aig (ILS on top of
&deepsyn; the deepsyn config then continues &deepsyn from that AIG with another seed).
usage: phase3a.py [-j jobs] [-T sec] [-S seed] [--start tmpl] [--out dir] [--configs ...] [--aggregate] [circuit ...]
"""
import argparse, concurrent.futures as cf, csv, glob, os, re, subprocess, sys

import phase2
from phase2c import CIRCUITS, EPFL, R, ROOT
from structdiff import Aig

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results", "phase3a")
CONFIGS = {"U": "-M U", "A": "-M A", "P0": "-M P0", "N": "-M N",
           "U-e": "-M U -e", "A-e": "-M A -e", "P0-e": "-M P0 -e", "N-e": "-M N -e",
           # Phase 3d: no per-iteration cec (-c; accepted results are still cec'd), R on the whole AIG or
           # (-L h) only on the region within h fanin/fanout steps of the replacement
           "U-e-c": "-M U -e -c", "U-e-L2": "-M U -e -c -L 2", "U-e-L4": "-M U -e -c -L 4", "U-e-L8": "-M U -e -c -L 8",
           "A-e-L4": "-M A -e -c -L 4"}
CHECKPOINTS = [10, 30, 60, 120, 300, 600, 1200, 1800, 3600]


START = EPFL + "/{c}.aig"


def run(job):
    c, cfg, T, seed = job
    tag = f"{c}_{cfg}"
    src = START.format(c=c)
    if cfg == "deepsyn":
        log = os.path.join(OUT, "raw", tag + ".txt")
        cmd = f"&r {src}; &deepsyn -T {T} -S {seed} -v; &ps; &w {OUT}/best/{tag}.aig"
        p = subprocess.run([phase2.ABC, "-c", cmd], cwd=ROOT, capture_output=True, text=True)
        open(log, "w").write(p.stdout + p.stderr)
        return tag, (p.stdout.strip().splitlines() or [""])[-1][:160]
    dst = os.path.join(OUT, "raw", tag + ".csv")
    cmd = (f'&r {src}; lsv_reencode -i {CONFIGS[cfg]} -K 8 -n 1000000 -T {T} -S {seed} '
           f'-R "{R}" -o {dst} -w {OUT}/best/{tag}.aig')
    p = subprocess.run([phase2.ABC, "-c", cmd], cwd=ROOT, capture_output=True, text=True)
    log = [l for l in p.stdout.splitlines() if l.startswith("lsv_reencode")]
    return tag, "\n".join(log) if log else p.stdout[-400:] + p.stderr[-400:]


def ils_curve(path):
    """[(seconds, best AND)] of one ILS run, and (iterations, accepted, R(x0) AND, x0 AND)."""
    rows = list(csv.DictReader(open(path)))
    pts, it, acc = [], 0, 0
    aRx = int(rows[0]["and_R"]) if rows else None
    for r in rows:
        if r["iter"] == "-1":
            pts.append((float(r["elapsed"]), aRx))
            continue
        it += 1
        if r["acc"] == "1":
            acc += 1
            pts.append((float(r["elapsed"]), int(r["and_R"])))
    return pts, dict(iters=it, accepted=acc, and_Rx=aRx, and_x0=int(rows[0]["and_cur"]) if rows else None)


def deepsyn_curve(path):
    pts = []
    for line in open(path):
        m = re.search(r"Time\s+([\d.]+) sec : And =\s*(\d+)", line)
        if m:
            pts.append((float(m.group(1)), int(m.group(2))))
    m = re.findall(r"and =\s*(\d+)", open(path).read())
    final = int(m[-1]) if m else None
    return pts, final


def at(pts, t):
    v = [a for s, a in pts if s <= t]
    return min(v) if v else ""


def aggregate(circuits, configs, T):
    curve, summ = [], []
    for c in circuits:
        x0 = Aig(START.format(c=c)).nAnd
        aRx = None
        rows = {}
        for cfg in configs:
            p = os.path.join(OUT, "raw", f"{c}_{cfg}.csv")
            if os.path.exists(p):
                pts, info = ils_curve(p)
                aRx = aRx or info["and_Rx"]
                rows[cfg] = (pts, info)
        p = os.path.join(OUT, "raw", f"{c}_deepsyn.txt")
        if os.path.exists(p):
            pts, final = deepsyn_curve(p)
            pts = [(0.0, x0)] + pts
            rows["deepsyn"] = (pts, dict(iters="", accepted=len(pts) - 1, final=final))
        for cfg, (pts, info) in rows.items():
            final = info.get("final") or (min(a for _, a in pts) if pts else "")
            curve.append(dict(circuit=c, config=cfg, **{f"t{t}": at(pts, t) for t in CHECKPOINTS if t <= T}))
            summ.append(dict(circuit=c, config=cfg, and_x0=x0, and_Rx=aRx, final=final,
                             d_vs_Rx=final - aRx if aRx and final != "" else "",
                             pct_vs_Rx=round(100 * (final - aRx) / aRx, 2) if aRx and final != "" else "",
                             iters=info["iters"], accepted=info["accepted"],
                             last_impr_s=round(pts[-1][0], 1) if pts else ""))
    for name, rows in (("phase3a_curve.csv", curve), ("phase3a_summary.csv", summ)):
        keys = list(dict.fromkeys(k for r in rows for k in r))
        with open(os.path.join(OUT, name), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            w.writerows(rows)


def main():
    global START, OUT
    ap = argparse.ArgumentParser()
    ap.add_argument("-j", type=int, default=40)
    ap.add_argument("-T", type=int, default=1200)
    ap.add_argument("-S", type=int, default=1)
    ap.add_argument("--configs", nargs="*", default=list(CONFIGS) + ["deepsyn"])
    ap.add_argument("--start", default=START)
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--aggregate", action="store_true")
    ap.add_argument("circuits", nargs="*", default=CIRCUITS)
    a = ap.parse_args()
    START, OUT = os.path.abspath(a.start), os.path.abspath(a.out)
    for d in ("raw", "best"):
        os.makedirs(os.path.join(OUT, d), exist_ok=True)
    if not a.aggregate:
        jobs = [(c, cfg, a.T, a.S) for c in a.circuits for cfg in a.configs]
        # resume: a job is complete once its final AIG is written (ILS: only if it passed cec)
        jobs = [j for j in jobs if not os.path.exists(os.path.join(OUT, "best", f"{j[0]}_{j[1]}.aig"))]
        print(len(jobs), "jobs to run", flush=True)
        with cf.ThreadPoolExecutor(a.j) as ex:
            for tag, log in ex.map(run, jobs):
                print(f"{tag:<20} {log}", flush=True)
    aggregate(a.circuits, a.configs, a.T)


if __name__ == "__main__":
    sys.exit(main())
