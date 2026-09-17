#!/bin/sh
# Validation harness: runs the reference C `euclid` tool and the new Chapel
# pipeline on the same battery of fields and reports both results side by
# side. The C tool's answer is treated as ground truth (exact rational
# minimum); the Chapel pipeline reports a rigorous numeric bracket that
# should contain it.
#
# Usage: tests/validate.sh   (run from the euclid-1.0 directory)
set -eu

cd "$(dirname "$0")/.."

if [ ! -x ./euclid ]; then
  echo "Reference C tool ./euclid not found; run 'make' first." >&2
  exit 1
fi
if [ ! -x ./euclid_chpl ]; then
  echo "Chapel binary ./euclid_chpl not found; run:" >&2
  echo "  chpl -M . main.chpl -o euclid_chpl" >&2
  exit 1
fi

# name : polynomial : chapel extra args
#
# q19 and q57 regression-guard two specific bugs found via a wider real
# quadratic battery (see CHAPEL.md "Bugs found and fixed"): q19 caught the
# bracket-finding phase locking the lower bound above the true minimum
# (converged to exactly 1.0 instead of 170/171), and q57 caught the
# small-elements candidate filter excluding needed large-coefficient/
# small-norm absorbers (converged to exactly 2x the true minimum, 28/19
# instead of 14/19).
FIELDS='
x2-2     x^2-2      --initialK=1.0 --tolerance=0.001
x2-61    x^2-61     --initialK=2.0 --tolerance=0.001
q19      x^2-19     --initialK=1.0 --tolerance=0.001
q57      x^2-57     --initialK=1.0 --tolerance=0.001
x3+x2-1  x^3+x^2-1  --initialK=0.5 --tolerance=0.001
x3-3x-1  x^3-3*x-1  --initialK=0.6 --tolerance=0.001
x5-x-1   x^5-x-1    --initialK=0.6 --tolerance=0.01 --refineDepth=16 --maxProblems=20000
'

echo "$FIELDS" | while read -r name poly extra; do
  [ -z "$name" ] && continue
  echo "=========================================================="
  echo "Field: $poly"
  echo "--- reference C tool (exact) ---"
  ./euclid "$poly" || echo "(C tool failed/timed out)"
  echo "--- Chapel pipeline (rigorous numeric bracket) ---"
  fixture="tests/fixtures/$name.txt"
  if [ ! -f "$fixture" ]; then
    echo "(missing fixture $fixture, generating)"
    POLY="$poly" gp -q generate_field.gp >/dev/null
    cp field_data.txt "$fixture"
  fi
  ./euclid_chpl --fieldFile="$fixture" $extra
  echo
done
