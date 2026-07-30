#!/bin/sh
# Regression check for multiple pair_style runner instances.
#
# Runs the single-instance reference input and the hybrid/overlay input
# that stacks the same 2G water potential twice, then checks that the
# overlaid pair energy is exactly twice the single-instance one (within
# a small relative tolerance for the log formatting). With shared global
# state, the second instance would have clobbered the first and this
# invariant would not hold.
#
# Usage: ./check_multi_instance.sh /path/to/lmp
set -e

LMP=${1:?"usage: $0 /path/to/lmp"}

epair () {
    # Pull epair of step 0 from a run-0 log: the thermo data line follows
    # the header line that starts with 'Step'.
    awk '/^ *Step /{getline; print $3; exit}' "$1"
}

$LMP -in in.ml-runner.H2O.single -log log.single.tmp > /dev/null
$LMP -in in.ml-runner.H2O.multi -log log.multi.tmp > /dev/null

E1=$(epair log.single.tmp)
E2=$(epair log.multi.tmp)

echo "single instance : epair = $E1"
echo "two instances   : epair = $E2"

ok=$(awk -v e1="$E1" -v e2="$E2" 'BEGIN {
    diff = e2 - 2.0 * e1
    if (diff < 0.0) diff = -diff
    ref = 2.0 * e1
    if (ref < 0.0) ref = -ref
    print (diff <= 1.0e-8 * ref) ? "yes" : "no"
}')

rm -f log.single.tmp log.multi.tmp

if [ "$ok" = "yes" ]; then
    echo "OK: two pair_style runner instances are independent (E2 = 2*E1)."
else
    echo "FAILED: expected epair(two instances) = 2 * epair(single instance)."
    exit 1
fi
