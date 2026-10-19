#!/usr/bin/env bash
# Unit tests of plan section 5 that apply to Phase 1.
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
TOOL=$HERE/../acd_phase1
MINCUT=$HERE/../../aig_mincut/aig_mincut
EPFL=${1:-$HERE/../../../../../testcases/aig/epfl}
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
make -s -C "$HERE/.." && make -s -C "$(dirname "$MINCUT")" || exit 1
fail=0
check() { if [ "$2" = "$3" ]; then printf '  %-55s ok\n' "$1"; else printf '  %-55s FAIL (got %s, want %s)\n' "$1" "$2" "$3"; fail=1; fi; }
# group line fields: 3=|G| 4=|L| 5=|Img| 6=r 7=saving 8=fd_size 12=members 13=fd_subset
field() { cut -d, -f"$1"; }

out=$("$TOOL" "$HERE/dec3to8.aag" "$HERE/dec3to8.cut" | grep ^group)
check "3-to-8 one-hot decoder: image"  "$(field 5 <<<"$out")" 8
check "3-to-8 one-hot decoder: r"      "$(field 6 <<<"$out")" 3

out=$("$TOOL" "$HERE/fd.aag" "$HERE/fd_ab.cut" | grep ^group)
check "cut {a, b, a&b}: functional-dependence subset" "$(field 13 <<<"$out")" "1 2"

out=$("$TOOL" "$HERE/fd.aag" "$HERE/fd_pis.cut" | grep ^group)
check "independent PIs: total saving" "$(field 7 <<<"$out" | paste -sd+ | bc)" 0

"$MINCUT" -q --stages 2 --period 2 --out-dir "$TMP" "$EPFL/dec.aig" >/dev/null
python3 -c "
import json,sys
for s in json.load(open('$TMP/cuts.json'))['cuts'][0]['signals']: print(s['node'], s['arrival'], s['downstream'])" > "$TMP/cut.txt"
out=$("$TOOL" --global "$EPFL/dec.aig" "$TMP/cut.txt")
g=$(grep ^group <<<"$out")
check "dec @ T=2: cut width"                 "$(field 3 <<<"$g" | paste -sd+ | bc)" 16
check "dec @ T=2: compressible to (sum r)"   "$(field 6 <<<"$g" | paste -sd+ | bc)" 8
check "dec @ T=2: global ceil(log2 N)"       "$(grep ^global <<<"$out" | field 3)" 8
check "dec @ T=2: max downstream depth = T (no slack)" "$(cut -d' ' -f3 "$TMP/cut.txt" | sort -n | tail -1)" 2

# Phase 2 tests (need results/phase1_*.csv and results/runs from phase0.py / phase1.py)
P2="python3 $HERE/../phase2.py"
BADC=div; BADT=rec; BADS=sink; BADG=slack   # a combo with accepted groups
if [ -f "$HERE/../results/phase1_groups.csv" ]; then
  export TMPDIR=${TMPDIR:-/tmp}
  $P2 --only dec --force dec:50%:source:plan --tag _unit_dec >/dev/null 2>&1
  r=$(python3 -c "
import csv; r=list(csv.DictReader(open('$HERE/../results/phase2_results_unit_dec.csv')))[0]
print(r['groups_tried'], r['groups_accepted'], r['groups_rejected_delay'], r['cec'])")
  check "dec @ T=2 Phase 2: tried accepted rejected_delay cec" "$r" "4 0 4 equivalent"
  $P2 --only "$BADC" --force "$BADC:$BADT:$BADS:$BADG" --inject-bad-decoder --tag _unit_bad >/dev/null 2>&1; rc=$?
  r=$(python3 -c "
import csv; r=[x for x in csv.DictReader(open('$HERE/../results/phase2_results_unit_bad.csv')) if x['T_tag']=='$BADT' and x['side']=='$BADS' and x['grouping']=='$BADG'][0]
print(r['cec'], 'reported' if r['area_after'] else 'not-reported')")
  check "injected bad decoder: cec fails, no area reported" "$r $rc" "FAIL not-reported 1"
fi

[ $fail -eq 0 ] && echo "ALL PASS" || echo "SOME TESTS FAILED"
exit $fail
