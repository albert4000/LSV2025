#!/usr/bin/env python3
"""Phase 3c: alternating &deepsyn and ILS (lsv_reencode -i) under the same total wall-clock budget.

Start x0 = the original EPFL AIG. Every job has a budget of 3600 s, split into steps run one after another,
each step starting from the AIG left by the previous one:
  D    &deepsyn -T 3600 (one run)
  DR   &deepsyn -T 600, six times (restarted from its own result with a new seed): the effect of restarts
  DI   (&deepsyn -T 600; ILS U-e -T 600) x 3: the alternating flow
  DN   (&deepsyn -T 600; ILS N-e -T 600) x 3: the same with null kicks (R + node order only)
  DLU  (&deepsyn -T 600; ILS U-e with local R -L 4 -T 600) x 3: DI with R only around the kick (Phase 3d)
  DLA  the same with sampler A (the best local-R config of Phase 3d)
ILS = lsv_reencode -i -K 8 -R "compress2rs x3" (Phase 3a); &deepsyn steps use seed 10*S + step (D: 10*S, so
step 0 of every config is the same &deepsyn call). Steps overrun their -T by at most one iteration; the
times recorded are wall clock. The final AIG of every job is checked with &cec against x0.

outputs (results/phase3c/):
  work/<job>/s<k>.aig, s<k>.log / s<k>.csv   AIG after each step, &deepsyn -v log or ILS CSV
  work/<job>/steps.csv                       the steps of the job (written when it finishes)
  steps.csv                                  one row per step: and / lev in and out, seconds, cumulative
  phase3c_summary.csv                        per job: final and / lev, cec, AND at the step boundaries
usage: phase3c.py [-j jobs] [-S seeds] [--configs ...] [--scale f] [--out dir] [--aggregate] [circuit ...]
(--scale multiplies every step time, for smoke tests)
"""
import argparse, concurrent.futures as cf, csv, glob, os, re, subprocess, sys, time

import phase2
from phase2c import EPFL, R, ROOT

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results", "phase3c")
CIRCUITS = ["bar", "log2", "multiplier", "sqrt", "square", "voter", "sin", "max", "priority", "int2float",
            "cavlc", "arbiter", "mem_ctrl", "i2c", "ctrl", "router"]
TOTAL = 3600
SCALE = 1.0
CONFIGS = {"D": [("deepsyn", 3600)],
           "DR": [("deepsyn", 600)] * 6,
           "DI": [("deepsyn", 600), ("ils:-M U -e", 600)] * 3,
           "DN": [("deepsyn", 600), ("ils:-M N -e", 600)] * 3,
           "DLU": [("deepsyn", 600), ("ils:-M U -e -c -L 4", 600)] * 3,
           "DLA": [("deepsyn", 600), ("ils:-M A -e -c -L 4", 600)] * 3}


def stats(path):
    p = subprocess.run([phase2.ABC, "-c", f"&r {path}; &ps"], cwd=ROOT, capture_output=True, text=True)
    s = re.sub(r"\x1b\[[0-9;]*m", "", p.stdout)
    a, l = re.search(r"and =\s*(\d+)", s), re.search(r"lev =\s*(\d+)", s)
    return (int(a.group(1)), int(l.group(1))) if a and l else (None, None)


def run(job):
    try:
        return run1(job)
    except Exception as e:   # one broken job (e.g. abc not executable) must not stop the others
        return f"{job[0]}_{job[1]}_S{job[2]}", [], f"error: {e!r}", None


def run1(job):
    c, cfg, S = job
    tag = f"{c}_{cfg}_S{S}"
    wd = os.path.join(OUT, "work", tag)
    os.makedirs(wd, exist_ok=True)
    x0 = f"{EPFL}/{c}.aig"
    cur = x0
    rows, cum = [], 0.0
    a0, l0 = stats(x0)
    for k, (kind, T) in enumerate(CONFIGS[cfg]):
        dst = os.path.join(wd, f"s{k}.aig")
        T = max(1, round(T * SCALE))
        t0 = time.time()
        if kind == "deepsyn":
            cmd = f"&r {cur}; &deepsyn -T {T} -S {(10 * S + k) % 101} -v; &w {dst}"
            p = subprocess.run([phase2.ABC, "-c", cmd], cwd=ROOT, capture_output=True, text=True)
            open(os.path.join(wd, f"s{k}.log"), "w").write(p.stdout + p.stderr)
        else:
            cmd = (f'&r {cur}; lsv_reencode -i {kind[4:]} -K 8 -n 1000000 -T {T} -S {100 * S + k} '
                   f'-R "{R}" -o {wd}/s{k}.csv -w {dst}')
            p = subprocess.run([phase2.ABC, "-c", cmd], cwd=ROOT, capture_output=True, text=True)
            open(os.path.join(wd, f"s{k}.log"), "w").write(p.stdout + p.stderr)
        secs = time.time() - t0
        cum += secs
        if not os.path.exists(dst):
            return tag, rows, f"step {k} ({kind}) wrote no AIG", None
        ai, li = stats(cur)
        ao, lo = stats(dst)
        rows.append(dict(job=tag, circuit=c, config=cfg, seed=S, step=k, kind=kind.split(":")[0],
                         T=T, and_in=ai, lev_in=li, and_out=ao, lev_out=lo, secs=round(secs, 1), cum=round(cum, 1)))
        cur = dst
    p = subprocess.run([phase2.ABC, "-c", f"&r {x0}; &cec {cur}"], cwd=ROOT, capture_output=True, text=True)
    cec = "ok" if "Networks are equivalent" in p.stdout else "FAIL"
    for r in rows:
        r["cec"] = cec
    write(os.path.join(wd, "steps.csv"), rows)
    return tag, rows, f"{a0} -> {rows[-1]['and_out']} (lev {rows[-1]['lev_out']}), cec {cec}, {cum:.0f}s", cec


def write(path, rows):
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def from_disk(tag):
    """Steps of a finished job: its steps.csv, or (jobs of the first run, which wrote no steps.csv) rebuilt
    from the step AIGs; secs then = difference of file times (step 0: its planned time)."""
    wd = os.path.join(OUT, "work", tag)
    p = os.path.join(wd, "steps.csv")
    if os.path.exists(p):
        return list(csv.DictReader(open(p)))
    c, cfg, S = tag.rsplit("_", 2)
    plan = CONFIGS[cfg]
    files = [os.path.join(wd, f"s{k}.aig") for k in range(len(plan))]
    if not all(os.path.exists(f) for f in files):
        return []
    x0 = f"{EPFL}/{c}.aig"
    rows, prev, cum = [], x0, 0.0
    for k, ((kind, T), f) in enumerate(zip(plan, files)):
        secs = T if k == 0 else os.path.getmtime(f) - os.path.getmtime(files[k - 1])
        cum += secs
        ai, li = stats(prev)
        ao, lo = stats(f)
        rows.append(dict(job=tag, circuit=c, config=cfg, seed=int(S[1:]), step=k, kind=kind.split(":")[0], T=T,
                         and_in=ai, lev_in=li, and_out=ao, lev_out=lo, secs=round(secs, 1), cum=round(cum, 1)))
        prev = f
    p = subprocess.run([phase2.ABC, "-c", f"&r {x0}; &cec {prev}"], cwd=ROOT, capture_output=True, text=True)
    for r in rows:
        r["cec"] = "ok" if "Networks are equivalent" in p.stdout else "FAIL"
    write(os.path.join(wd, "steps.csv"), rows)
    return rows


def main():
    global SCALE, OUT
    ap = argparse.ArgumentParser()
    ap.add_argument("-j", type=int, default=43)
    ap.add_argument("-S", type=int, nargs="*", default=[1, 2])
    ap.add_argument("--configs", nargs="*", default=list(CONFIGS))
    ap.add_argument("--scale", type=float, default=1.0)
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--aggregate", action="store_true", help="only rebuild the summaries from the work directories")
    ap.add_argument("circuits", nargs="*", default=CIRCUITS)
    a = ap.parse_args()
    SCALE, OUT = a.scale, os.path.abspath(a.out)
    os.makedirs(OUT, exist_ok=True)
    jobs = [(c, cfg, S) for S in a.S for c in a.circuits for cfg in a.configs]
    tags = [f"{c}_{cfg}_S{S}" for c, cfg, S in jobs]
    if not a.aggregate:
        # resume: run only the jobs whose work directory is missing or empty (a non-empty one finished or is
        # still being run by another driver)
        wdir = lambda t: os.path.join(OUT, "work", t)
        todo = [j for j, t in zip(jobs, tags) if not os.path.isdir(wdir(t)) or not os.listdir(wdir(t))]
        print(len(todo), "jobs to run", flush=True)
        with cf.ThreadPoolExecutor(a.j) as ex:
            for tag, rows, msg, cec in ex.map(run, todo):
                print(f"{tag:<24} {msg}", flush=True)
    steps, summ = [], []
    for t in tags:
        rows = from_disk(t) if os.path.isdir(os.path.join(OUT, "work", t)) else []
        if not rows:
            print(f"{t:<24} incomplete", flush=True)
            continue
        steps += rows
        r0, rl = rows[0], rows[-1]
        cfg = r0["config"]
        b = {}
        for r in rows:   # AND at each step boundary, by planned cumulative time
            b[sum(T for _, T in CONFIGS[cfg][:int(r["step"]) + 1])] = r["and_out"]
        summ.append(dict(job=t, circuit=r0["circuit"], config=cfg, seed=r0["seed"], and_x0=r0["and_in"],
                         final=rl["and_out"], lev=rl["lev_out"], cec=rl["cec"], secs=rl["cum"],
                         **{f"t{x}": b.get(x, "") for x in range(600, TOTAL + 1, 600)}))
    write(os.path.join(OUT, "steps.csv"), steps)
    write(os.path.join(OUT, "phase3c_summary.csv"), summ)

if __name__ == "__main__":
    sys.exit(main())
