#!/usr/bin/env bash
# Regression tests for the ABC commands lsv_pipe / lsv_pipe_acd.
#   usage: src/ext-lsv/tests/lsv_pipe_test.sh [epfl_dir]
# Needs ./abc built, python3, and src/develop/aig_mincut (for the PI-delayed reference used by dsec).
set -u
ROOT=$(cd "$(dirname "$0")/../../.." && pwd)
EPFL=${1:-$ROOT/../testcases/aig/epfl}
ABC=$ROOT/abc
REF=$ROOT/src/develop/aig_mincut/aig_mincut
GEN="python3 $ROOT/src/develop/aig_mincut/gen_aig.py"
# ABC cannot read ASCII AIGER (.aag); convert with the helper from src/develop/acd_reencode
to_aig() { python3 -c "import sys; sys.path.insert(0, '$ROOT/src/develop/acd_reencode'); import aiger; aiger.write(aiger.read('$1'), '$2')"; }
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
make -s -C "$(dirname "$REF")" >/dev/null || exit 1
fail=0
ok()  { printf '  %-60s ok\n' "$1"; }
bad() { printf '  %-60s FAIL %s\n' "$1" "${2:-}"; fail=1; }
abc() { "$ABC" -c "$1" 2>&1; }

echo "== lsv_pipe: register counts (S, T) vs. known optima, all checks, dsec of -p output"
for spec in dec:2:0:8 dec:3:1:48 ctrl:4:3:67 int2float:3:6:48 router:2:41:22 cavlc:4:5:173 i2c:3:8:340 adder:4:64:579; do
  IFS=: read -r c S T want <<<"$spec"
  out=$(abc "&r $EPFL/$c.aig; lsv_pipe -S $S -T $T -p; &w $TMP/p.aig")
  got=$(sed -n 's/.*registers *: \([0-9]*\).*/\1/p' <<<"$out")
  grep -q "checks *: all passed" <<<"$out" || { bad "$c S=$S T=$T checks" "$(grep FAIL <<<"$out")"; continue; }
  [ "$got" = "$want" ] || { bad "$c S=$S T=$T registers" "$got, want $want"; continue; }
  "$REF" -q --no-json --sim 0 --stages "$S" --period "$T" --write-ref-blif "$TMP/r.blif" "$EPFL/$c.aig" >/dev/null
  abc "dsec $TMP/p.aig $TMP/r.blif" | grep -q "Networks are equivalent" && ok "$c S=$S T=$T: $got registers, dsec equivalent" || bad "$c S=$S T=$T dsec"
done

echo "== lsv_pipe: brute-force optimality and both cut sides on random AIGs"
n=0
for seed in $(seq 1 60); do
  $GEN random $((2 + seed % 4)) $((6 + seed % 9)) $((1 + seed % 4)) "$seed" > "$TMP/r.aag"
  to_aig "$TMP/r.aag" "$TMP/r.aig"
  for S in 2 3; do for T in 0 3; do for k in "" "-k"; do
    out=$(abc "&r $TMP/r.aig; lsv_pipe -S $S -T $T -b 20 $k")
    n=$((n + 1))
    grep -q "INFEASIBLE" <<<"$out" && continue
    grep -q "checks *: all passed" <<<"$out" || { bad "random seed=$seed S=$S T=$T $k" "$(grep FAIL <<<"$out" | head -1)"; }
  done; done; done
done
ok "$n random instances checked"

echo "== lsv_pipe_acd: results, cec, dsec of -p output"
for spec in "ctrl:5:-g:23" "dec:2:-g:16" "log2:222:-g -k:204" "div:3279:-g -k:273" "mem_ctrl:63:-g:1406"; do
  IFS=: read -r c T f want <<<"$spec"
  out=$(abc "&r $EPFL/$c.aig; lsv_pipe_acd -T $T $f -p; &w $TMP/a.aig")
  got=$(sed -n 's/.*registers [0-9]* -> \([0-9]*\),.*/\1/p' <<<"$out")
  grep -q "cec (registers as wires) vs. original: equivalent" <<<"$out" || { bad "$c T=$T $f cec" "$(tail -2 <<<"$out")"; continue; }
  [ "$got" = "$want" ] || { bad "$c T=$T $f registers" "$got, want $want"; continue; }
  "$REF" -q --no-json --sim 0 --stages 2 --period "$T" --write-ref-blif "$TMP/r.blif" "$EPFL/$c.aig" >/dev/null
  abc "dsec $TMP/a.aig $TMP/r.blif" | grep -q "Networks are equivalent" && ok "$c T=$T $f: $got registers, cec + dsec equivalent" || bad "$c T=$T dsec"
done

[ $fail -eq 0 ] && echo "ALL PASS" || echo "SOME TESTS FAILED"
exit $fail
