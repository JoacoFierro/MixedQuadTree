#!/usr/bin/env bash
# Regression test for Issue A (first-subquad Unknown) + Issue (WindingState
# not derived from winding numbers).
#
# Verifies that after every bridge split the new sub-quads (replacement +
# appended) have a proper WindingState, not Unknown.

set -e

BIN="${BIN:-./build/mesher_roi}"
DATA="${DATA:-./data}"

if [ ! -x "$BIN" ]; then
    echo "FAIL: binary $BIN not found" >&2
    exit 2
fi

failed=0
for f in a.poly tusqh_bridge_small_gap.poly tusqh_boundary_filter.poly tusqh_small_feature.poly; do
    name=$(basename "$f" .poly)
    out="reg_$name"
    rm -f "${out}"_*.vtk "${out}"_*.oct 2>/dev/null
    log=$("$BIN" -p "$DATA/$f" -u "$out" -T -J -a 3 2>&1)
    unknown=$(echo "$log" | grep -o "Unknown=[0-9]*" | head -1 | cut -d= -f2)
    if [ "$unknown" != "0" ]; then
        echo "FAIL: $f has Unknown=$unknown (expected 0)"
        echo "$log" | grep -E "postbridge|bridges added" | head -5
        failed=$((failed + 1))
    else
        kept=$(echo "$log" | grep -o "[0-9]* will keep" | head -1 | awk '{print $1}')
        echo "OK:   $f Unknown=0, will keep=$kept"
    fi
done

if [ "$failed" -gt 0 ]; then
    echo
    echo "$failed case(s) failed"
    exit 1
fi
echo
echo "All cases passed."
