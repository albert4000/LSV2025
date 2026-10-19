# aig_mincut

Fixed-stage pipeline register minimization on a combinational AIG — implementation of
`../plan/aig_mincut_spec.md` (stages A–D, milestones 1–6). Standalone C++17, no dependencies.

```sh
make
./aig_mincut --stages 4 --period 70 --out-dir out sin.aig     # writes out/stage_assign.json, out/cuts.json
./aig_mincut --stages 3 --write-blif p.blif --write-ref-blif r.blif dec.aig
../../../abc -c "dsec p.blif r.blif"                            # sequential equivalence
./run_tests.sh                                                  # full regression (section 9)
```

Options: `--stages S`, `--mode count|area`, `--reg-weight W`, `--cut-side source|sink`,
`--period T`, `--method auto|flow|lp`, `--brute N`, `--spec-graph`, `--sim N`, `-q` (see `--help`).
Exit code 3 = infeasible for the given T.

## Files

| file | content |
|---|---|
| `aig.{h,cpp}` | AIGER reader (`.aag`/`.aig`), topological renumbering, CSR fanout |
| `flow.h` | residual graph, Dinic (iterative), primal-dual min-cost flow, Dijkstra with potentials |
| `mincut.cpp` | stages A/B/C, brute force, verification, JSON / BLIF output, report |
| `gen_aig.py` | adder / parity / comparator / random AIG generators |
| `tests/fanout_share.aag` | minimal case where the original spec 4.1 graph (`v_in -> u_out`) is suboptimal |

## Method

- **A (S=2)**: node split `v_in -> v_out` (cap W), `u_out -> v_in` ∞, closure arc **`v_in -> u_in`** ∞,
  `src -> PI_in`, `PO_out -> snk`. Source side = stage 0; v is registered iff `v_in ∈ X, v_out ∉ X`.
- **B (any S)**: the LP of spec 5.1 (plus the redundant bound `m(u) <= S-1`). Each constraint
  `x_b - x_a <= c` is an uncapacitated arc `a -> b` of cost c; supply W at `m(u)`, demand W at `s(u)`.
  Solved by primal-dual min-cost flow (Dijkstra + Dinic on zero-reduced-cost arcs; the number of
  phases is O(S), 7 for S=4 on every EPFL circuit). The initial feasible potential is the ASAP schedule.
  Optimal `s, m` = shortest distances from the reference node in the final residual graph
  (`source` side = pointwise-max s, cut nearest the PIs) or distances to it (`sink` side).
- **C (depth T)**: level windows `ceil(L/T)-1 <= s(v) <= S-ceil((R+1)/T)`, then lazy constraints
  `s(v) >= s(u)+1` for every same-stage path of T+1 ANDs ending at a violating node. Every lazy
  constraint is *necessary*, so the result is exactly optimal and "infeasible" is exact
  (not heuristic as spec 7/11 assumes); the brute-force tests confirm both.
- Registers reset to the node's value under all-zero inputs, which makes the pipeline exactly
  sequentially equivalent to the original with S-1 zero-reset delays on every PI (used for `dsec`).

## Deviations from the spec

1. **Spec 4.1 reverse arc was wrong (fixed in the spec since).** `v_in -> u_out` forbids a registered node u from also feeding a
   stage-0 node directly (register sharing), so the min-cut over-counts. Correct arc: `v_in -> u_in`.
   With the original arc (`--spec-graph`) the result is suboptimal on 68/300 random AIGs and on EPFL
   `max` (512 vs 130), `priority` (128 vs 8), `i2c` (147 vs 141). See `tests/fanout_share.aag` (3 vs 2).
2. No scipy prototype (scipy is not installed); the LP is validated by brute force and the S=2
   max-flow cross-check instead.
3. `--cut-side balanced` is not implemented (spec allows deferring it).
4. "balancing" registers in the report / `cuts.json` = registers on PI chains (spec 5.3).

## Verification (`run_tests.sh`)

1. hand-written adder/parity/comparator, S=2..4, both cut sides: legal, A == B, == brute force
2. 2700 random instances (≤14 ANDs, S=2..4, T∈{none,2,3}): LP == brute force, including infeasibility
3. spec-graph comparison (informational)
4. EPFL: S=2 A==B with identical extreme cuts, monotone in S, T-constrained runs
5. ABC `dsec` on pipeline vs PI-delayed original (negative tests: a flipped gate and a bypassed
   register are both reported NOT EQUIVALENT)

Every run also checks legality/balance and a 64-way cycle-accurate latency simulation.

## Scale (S=4)

| circuit | ANDs | depth | no T | T = ⌈depth/S⌉ | T ×1.2 | time (T ×1.2) |
|---|---|---|---|---|---|---|
| adder | 1020 | 255 | 387 | 579 | 543 | <0.01s |
| sin | 5416 | 225 | 72 | 234 | 184 | 0.2s |
| multiplier | 27062 | 274 | 384 | 687 | 616 | 7.0s |
| log2 | 32060 | 444 | 96 | 611 | 485 | 17.3s |
| mem_ctrl | 46836 | 114 | 3579 | 4530 | 4321 | 17.8s |
| hyp | 214335 | 24801 | 384 | 858 | 771 | 10.3s |

Without T, every EPFL circuit runs in ≤ 3 s (hyp) and the optimum is (S−1) × the S=2 min-cut,
as expected: all cuts stack on the narrowest boundary.
