#!/usr/bin/env python3
"""Phase 0 of src/develop/plan/window_reencode_pareto.md: baseline area/delay Pareto fronts and noise.

For every EPFL circuit build AIG variants with ABC's own scripts (iterated to convergence where it
makes sense, &deepsyn with several time limits and seeds, delay-oriented scripts), cec each one against
the original, and map each one with asap7 (&nf -p, plain and after &dch -f, sweeping the delay
relaxation -R; &nf's -D is not wired to the mapper, so -R is the only delay knob).

M1 = (AND count, AIG level); M2 = (stime area, stime delay in ps).

outputs (results/phase0/):
  aig/<circuit>/<variant>.aig
  phase0_baseline.csv   one row per (variant, mapping flow, R)
  phase0_front_M1.csv   non-dominated (ands, levels) per circuit
  phase0_front_M2.csv   non-dominated (map_area, map_delay) per circuit
  phase0_noise.csv      &deepsyn seed spread (sigma) per circuit and time limit
usage: phase0.py [-j jobs] [-T 60,600] [-o dir] [circuit ...]
"""
import argparse, concurrent.futures as cf, csv, os, re, shutil, statistics, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
ABC = os.path.join(ROOT, "abc")
EPFL = os.path.abspath(os.path.join(ROOT, "..", "testcases", "aig", "epfl"))
LIB = os.path.abspath(os.path.join(ROOT, "..", "Downloads", "gradmap_libs", "asap7.lib"))
OUT = os.path.join(HERE, "results", "phase0")

CIRCUITS = ["adder", "arbiter", "bar", "cavlc", "ctrl", "dec", "div", "hyp", "i2c", "int2float", "log2", "max",
            "mem_ctrl", "multiplier", "priority", "router", "sin", "sqrt", "square", "voter"]

# name -> (script, iterate until the AND count stops decreasing)
SCRIPTS = {
    "orig": ("strash", False),
    "dc2": ("strash; dc2", True),
    "resyn2": ("strash; resyn2", True),
    "compress2rs": ("strash; compress2rs", True),
    "balance": ("strash; balance", False),
    "syn2": ("&get -n; &syn2; &put", False),
    "bd": ("&get -n; &st; &b -d; &put", False),
}
DEEPSYN_T = [60, 600]
DEEPSYN_S = [1, 2, 3]
MAX_ITERS, ITER_BUDGET = 30, 1800      # iterated scripts: max rounds / seconds
RELAX = [0, 10, 25, 50, 100]
MAP_FLOWS = {"nf": "", "dch": "&dch -f; "}
CEC_TIMEOUT = 3600

STATS = re.compile(r"and\s*=\s*(\d+)\s+lev\s*=\s*(\d+)")
STIME = re.compile(r"Area\s*=\s*([\d.]+).*?Delay\s*=\s*([\d.]+)\s*ps")


def abc(cmd, timeout=None):
    p = subprocess.run([ABC, "-c", cmd], cwd=ROOT, capture_output=True, text=True, timeout=timeout)
    return re.sub(r"\x1b\[[0-9;]*m", "", p.stdout + p.stderr)


def stats(path):
    m = STATS.findall(abc(f"read {path}; print_stats"))
    return (int(m[-1][0]), int(m[-1][1])) if m else (None, None)


def build(job):
    """Build one AIG variant; returns (circuit, variant, seed, path, seconds, rounds, note)."""
    c, v, seed = job
    src = os.path.join(EPFL, c + ".aig")
    d = os.path.join(OUT, "aig", c)
    os.makedirs(d, exist_ok=True)
    dst = os.path.join(d, v + ".aig")
    t0 = time.time()
    rounds, note = 1, ""
    try:
        if v.startswith("deepsyn"):
            T = int(v.split("_")[1][1:])
            abc(f"&r {src}; &deepsyn -T {T} -S {seed}; &w {dst}", timeout=4 * T + 1800)
        else:
            script, iterate = SCRIPTS[v]
            abc(f"read {src}; {script}; write_aiger -s {dst}", timeout=ITER_BUDGET)
            if iterate:
                best = stats(dst)[0]
                tmp = dst + ".tmp.aig"
                while rounds < MAX_ITERS and time.time() - t0 < ITER_BUDGET:
                    abc(f"read {dst}; {script}; write_aiger -s {tmp}", timeout=ITER_BUDGET)
                    a = stats(tmp)[0]
                    if a is None or a >= best:
                        break
                    shutil.move(tmp, dst)
                    best, rounds = a, rounds + 1
                if os.path.exists(tmp):
                    os.remove(tmp)
                if rounds >= MAX_ITERS or time.time() - t0 >= ITER_BUDGET:
                    note = "not converged"
    except subprocess.TimeoutExpired:
        note = "timeout"
    if not os.path.exists(dst):
        note = note or "no output"
    return c, v, seed, dst, time.time() - t0, rounds, note


def check(job):
    """cec against the original, AIG stats, and all mapping points of one variant."""
    c, v, seed, path, secs, rounds, note = job
    if note in ("timeout", "no output") and not os.path.exists(path):
        return [dict(circuit=c, variant=v, seed=seed, time=round(secs, 1), rounds=rounds, ands="", levels="",
                     cec="", flow="", relax="", map_area="", map_delay="", note=note)]
    src = os.path.join(EPFL, c + ".aig")
    try:
        log = abc(f"read {path}; cec -n {src}", timeout=CEC_TIMEOUT)
        cec = "ok" if "are equivalent" in log else ("FAIL" if "are NOT EQUIVALENT" in log else "undecided")
    except subprocess.TimeoutExpired:
        cec = "timeout"
    ands, lev = stats(path)
    rows = []
    for flow, pre in MAP_FLOWS.items():
        cmd = f"read_lib {LIB}; " + "; ".join(f"&r {path}; {pre}&nf -p -R {r}; &put; topo; stime -p" for r in RELAX)
        try:
            pts = STIME.findall(abc(cmd, timeout=7200))
        except subprocess.TimeoutExpired:
            pts = []
        for i, r in enumerate(RELAX):
            a, dl = pts[i] if i < len(pts) else ("", "")
            rows.append(dict(circuit=c, variant=v, seed=seed, time=round(secs, 1), rounds=rounds, ands=ands,
                             levels=lev, cec=cec, flow=flow, relax=r, map_area=a, map_delay=dl, note=note))
    return rows


def front(points):
    """Non-dominated subset (minimize both coordinates) of [(x, y, row)]."""
    out, best_y = [], float("inf")
    for x, y, row in sorted(points, key=lambda p: (p[0], p[1])):
        if y < best_y:
            out.append(row)
            best_y = y
    return out


def write(name, rows, fields=None):
    with open(os.path.join(OUT, name), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields or list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main():
    global OUT, DEEPSYN_T
    ap = argparse.ArgumentParser()
    ap.add_argument("-j", type=int, default=36)
    ap.add_argument("-T", default=",".join(map(str, DEEPSYN_T)), help="&deepsyn time limits (comma list)")
    ap.add_argument("-o", default=OUT, help="output directory")
    ap.add_argument("circuits", nargs="*", default=CIRCUITS)
    a = ap.parse_args()
    OUT, DEEPSYN_T = os.path.abspath(a.o), [int(t) for t in a.T.split(",")]
    os.makedirs(OUT, exist_ok=True)
    jobs = [(c, v, 0) for c in a.circuits for v in SCRIPTS]
    jobs += [(c, f"deepsyn_T{T}_S{s}", s) for c in a.circuits for T in DEEPSYN_T for s in DEEPSYN_S]
    jobs.sort(key=lambda j: -os.path.getsize(os.path.join(EPFL, j[0] + ".aig")))   # big ones first
    rows = []
    with cf.ThreadPoolExecutor(a.j) as ex:
        built = [ex.submit(build, j) for j in jobs]
        checked = []
        for f in cf.as_completed(built):
            b = f.result()
            print(f"built {b[0]:<10} {b[1]:<18} {b[4]:8.1f}s rounds={b[5]} {b[6]}", flush=True)
            checked.append(ex.submit(check, b))
        for f in cf.as_completed(checked):
            r = f.result()
            rows += r
            print(f"checked {r[0]['circuit']:<10} {r[0]['variant']:<18} cec={r[0]['cec']}", flush=True)
    rows.sort(key=lambda r: (r["circuit"], r["variant"], r["flow"], r["relax"]))
    write("phase0_baseline.csv", rows)

    good = [r for r in rows if r["cec"] == "ok"]
    f1, f2 = [], []
    for c in a.circuits:
        aigs = {r["variant"]: r for r in good if r["circuit"] == c}
        f1 += front([(r["ands"], r["levels"], r) for r in aigs.values()])
        f2 += front([(float(r["map_area"]), float(r["map_delay"]), r)
                     for r in good if r["circuit"] == c and r["map_area"] != ""])
    m1 = ["circuit", "variant", "ands", "levels", "time"]
    write("phase0_front_M1.csv", [{k: r[k] for k in m1} for r in f1], m1)
    m2 = ["circuit", "variant", "flow", "relax", "map_area", "map_delay", "ands", "levels"]
    write("phase0_front_M2.csv", [{k: r[k] for k in m2} for r in f2], m2)

    noise = []
    for c in a.circuits:
        for T in DEEPSYN_T:
            v = {r["variant"]: r for r in good if r["circuit"] == c and r["variant"].startswith(f"deepsyn_T{T}_")}
            ands = [r["ands"] for r in v.values()]
            if len(ands) >= 2:
                noise.append(dict(circuit=c, T=T, n=len(ands), mean=round(statistics.mean(ands), 1),
                                  sigma=round(statistics.stdev(ands), 1), min=min(ands), max=max(ands)))
    if noise:
        write("phase0_noise.csv", noise)
    bad = [r for r in rows if r["cec"] not in ("ok",)]
    print(len(rows), "rows;", len({(r['circuit'], r['variant']) for r in bad}), "variants without cec ok")


if __name__ == "__main__":
    sys.exit(main())
