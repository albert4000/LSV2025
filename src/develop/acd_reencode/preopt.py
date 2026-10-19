#!/usr/bin/env python3
"""Does logic optimization before the min-cut change the min-register result?

For every EPFL circuit build optimized variants with ABC, check them with cec, then run
aig_mincut (S=2, source side) at
  - the original absolute T (50% and recommended), i.e. the same delay budget, and
  - 50% of the variant's own depth.
area = live ANDs + 10 * registers.
outputs: results/preopt/<variant>/<circuit>.aig, results/preopt_results.csv
"""
import concurrent.futures as cf, csv, os, re, subprocess, sys
from phase0 import BASELINE, BIN, ROOT, HERE

ABC = os.path.join(ROOT, "abc")
EPFL = os.path.join(ROOT, "..", "testcases", "aig", "epfl")
OUT = os.path.join(HERE, "results", "preopt")
VARIANTS = {
    "orig": "strash",
    # level-preserving resubstitution (resub default keeps levels)
    "rs": "strash; rs -K 6; rs -K 8 -N 2; rs -K 10; rs -K 12 -N 2; rs -K 12 -N 2",
    # area-oriented, levels not preserved (compress2rs = rs/rw/rf with -l)
    "area": "strash; compress2rs; dc2; compress2rs; compress2rs",
}


def opt(job):
    c, v = job
    os.makedirs(os.path.join(OUT, v), exist_ok=True)
    dst = os.path.join(OUT, v, c + ".aig")
    src = os.path.join(EPFL, c + ".aig")
    p = subprocess.run([ABC, "-c", f"read {src}; {VARIANTS[v]}; write_aiger {dst}; cec {src}; print_stats"],
                       cwd=ROOT, capture_output=True, text=True)
    ok = "are equivalent" in p.stdout
    m = re.search(r"and\s*=\s*(\d+)\s+lev\s*=\s*(\d+)", p.stdout)
    return c, v, dst, ok, int(m.group(1)), int(m.group(2))


def mincut(path, T):
    p = subprocess.run([BIN, "--stages", "2", "--period", str(T), "--cut-side", "source", "--no-json",
                        "--sim", "0", path], capture_output=True, text=True)
    if p.returncode == 3:
        return None
    m = re.search(r"registers\s*: (\d+)", p.stdout)
    return int(m.group(1)) if m else None


def main():
    jobs = [(c, v) for c in BASELINE for v in VARIANTS]
    with cf.ThreadPoolExecutor(os.cpu_count()) as ex:
        opts = list(ex.map(opt, jobs))
    runs = []
    for c, v, path, ok, ands, lev in opts:
        depth, _, trec, _, _ = BASELINE[c]
        for tag, T in (("50%", (depth + 1) // 2), ("rec", trec), ("own50%", (lev + 1) // 2)):
            runs.append((c, v, path, ok, ands, lev, tag, T))
    with cf.ThreadPoolExecutor(os.cpu_count()) as ex:
        regs = list(ex.map(lambda r: mincut(r[2], r[7]), runs))
    rows = []
    for (c, v, _, ok, ands, lev, tag, T), r in zip(runs, regs):
        rows.append({"circuit": c, "variant": v, "cec": "ok" if ok else "FAIL", "ands": ands, "depth": lev,
                     "T_tag": tag, "T": T, "regs": r if r is not None else "infeasible",
                     "area": ands + 10 * r if r is not None else ""})
    with open(os.path.join(HERE, "results", "preopt_results.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print(len(rows), "rows;", sum(r["cec"] != "ok" for r in rows), "cec failures")


if __name__ == "__main__":
    sys.exit(main())
