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
FIELDS='
x2-2     x^2-2      --initialK=1.0 --tolerance=0.001
x2-61    x^2-61     --initialK=2.0 --tolerance=0.001
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
