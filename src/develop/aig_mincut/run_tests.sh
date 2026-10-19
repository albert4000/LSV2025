#!/usr/bin/env bash
# Regression tests for aig_mincut (spec section 9).
#   usage: run_tests.sh [epfl_dir]     (default ../../../../testcases/aig/epfl)
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../../.." && pwd)
EPFL=${1:-$ROOT/../testcases/aig/epfl}
BIN=$HERE/aig_mincut
ABC=$ROOT/abc
GEN="python3 $HERE/gen_aig.py"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
make -s -C "$HERE" >/dev/null || exit 1
fail=0
ok()  { printf '  %-58s ok\n' "$1"; }
bad() { printf '  %-58s FAIL %s\n' "$1" "${2:-}"; fail=1; }
regs() { "$BIN" -q --no-json "$@" | sed -n 's/.* regs=\([0-9]*\).*/\1/p'; }

echo "== 1. hand-written circuits: legality, A==B, brute force, sink side"
$GEN adder 3   > "$TMP/adder3.aag"
$GEN parity 8  > "$TMP/parity8.aag"
$GEN compare 3 > "$TMP/cmp3.aag"
for f in adder3 parity8 cmp3; do
  for S in 2 3 4; do
    for side in source sink; do
      out=$("$BIN" --no-json --brute 24 --stages $S --cut-side $side "$TMP/$f.aag")
      if grep -q FAIL <<<"$out"; then bad "$f S=$S $side" "$(grep FAIL <<<"$out")"; else ok "$f S=$S $side regs=$(sed -n 's/.*registers *: \([0-9]*\).*/\1/p' <<<"$out")"; fi
    done
  done
done

echo "== 2. random AIGs (<= 14 ANDs): LP == brute force, flow == LP, with and without T"
nrand=0; nfail=0
for seed in $(seq 1 300); do
  $GEN random $((2 + seed % 4)) $((6 + seed % 9)) $((1 + seed % 4)) $seed > "$TMP/r.aag"
  for S in 2 3 4; do
    for T in 0 2 3; do
      out=$("$BIN" -q --no-json --brute 20 --stages $S --period $T "$TMP/r.aag")
      rc=$?
      nrand=$((nrand + 1))
      [ $rc -eq 3 ] && continue   # infeasible, already confirmed by brute force (rc 1 otherwise)
      if [ $rc -ne 0 ]; then nfail=$((nfail + 1)); [ $nfail -le 3 ] && bad "random seed=$seed S=$S T=$T" "$out"; cp "$TMP/r.aag" "$TMP/fail_$seed.aag"; fi
    done
  done
done
[ $nfail -eq 0 ] && ok "$nrand random instances" || bad "$nfail of $nrand random instances"

echo "== 3. original spec 4.1 graph (v_in->u_out arc) vs exact optimum"
nbad=0
for seed in $(seq 1 300); do
  $GEN random $((2 + seed % 4)) $((6 + seed % 9)) $((1 + seed % 4)) $seed > "$TMP/r.aag"
  out=$("$BIN" --no-json --spec-graph --stages 2 "$TMP/r.aag")
  if grep -q FAIL <<<"$out"; then nbad=$((nbad + 1)); [ -f "$TMP/spec_cex.aag" ] || cp "$TMP/r.aag" "$TMP/spec_cex.aag"; fi
done
echo "  spec graph is suboptimal on $nbad of 300 random AIGs (informational)"
r_ok=$(regs --brute 20 "$HERE/tests/fanout_share.aag"); r_spec=$(regs --method flow --spec-graph "$HERE/tests/fanout_share.aag")
[ "$r_ok" = 2 ] && [ "$r_spec" = 3 ] && ok "tests/fanout_share.aag: fixed graph 2, spec graph 3" || bad "tests/fanout_share.aag" "$r_ok / $r_spec"

echo "== 4. EPFL: S=2 cross-check, monotonicity in S, T constraint"
for f in dec ctrl int2float router cavlc priority i2c adder bar max; do
  [ -f "$EPFL/$f.aig" ] || continue
  prev=0; line=""
  for S in 2 3 4 5; do
    out=$("$BIN" -q --no-json --stages $S "$EPFL/$f.aig")
    r=$(sed -n 's/.* regs=\([0-9]*\).*/\1/p' <<<"$out")
    grep -q ' ok$' <<<"$out" || bad "$f S=$S" "$out"
    [ "$r" -ge "$prev" ] || bad "$f monotone S=$S" "$r < $prev"
    prev=$r; line="$line S$S=$r"
  done
  ok "$f:$line"
done
for f in dec ctrl int2float router cavlc; do
  [ -f "$EPFL/$f.aig" ] || continue
  out=$("$BIN" -q --no-json --stages 3 --period 6 "$EPFL/$f.aig"); rc=$?
  [ $rc -eq 0 ] || [ $rc -eq 3 ] && ok "$f S=3 T=6: $(sed 's/.*S=3 //' <<<"$out")" || bad "$f S=3 T=6" "$out"
done

echo "== 5. sequential equivalence with ABC dsec (pipeline vs PI-delayed original)"
if [ -x "$ABC" ]; then
  for spec in "$TMP/adder3.aag:3:0" "$TMP/cmp3.aag:4:2" "$EPFL/dec.aig:2:0" "$EPFL/dec.aig:3:1" \
              "$EPFL/ctrl.aig:3:0" "$EPFL/int2float.aig:4:4" "$EPFL/router.aig:3:0" "$EPFL/cavlc.aig:3:6"; do
    IFS=: read -r f S T <<<"$spec"
    [ -f "$f" ] || continue
    "$BIN" -q --no-json --stages $S --period $T --write-blif "$TMP/p.blif" --write-ref-blif "$TMP/r.blif" "$f" >/dev/null || { bad "$(basename $f) S=$S T=$T" "solver"; continue; }
    res=$("$ABC" -c "dsec $TMP/p.blif $TMP/r.blif" 2>&1)
    if grep -q "Networks are equivalent" <<<"$res"; then ok "dsec $(basename $f) S=$S T=$T"; else bad "dsec $(basename $f) S=$S T=$T" "$(tail -2 <<<"$res")"; fi
  done
else
  echo "  (skipped: $ABC not built)"
fi

[ $fail -eq 0 ] && echo "ALL PASS" || echo "SOME TESTS FAILED"
exit $fail
