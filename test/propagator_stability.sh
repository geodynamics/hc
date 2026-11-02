#!/bin/bash
#
# diagnostic for the precision loss of the propagator matrix solution
# with increasing maximum degree. inserts a viscosity layer boundary
# without a viscosity change, which should not change the solution;
# the relative difference of the two geoid solutions measures the
# round-off amplification. this is not a pass/fail test.
#
# usage: ./propagator_stability.sh [bindir] [viscosity file] [lmax list] [extra options]
#
# e.g. ./propagator_stability.sh ../bin ../example_data/visc.sh08 "63 127 255" "-solver 1"
#
bindir=${1-../bin}
vfile=${2-../example_data/visc.sh08}
lmaxs=${3-"15 31 63 127 255"}
extra=${4-""}

bindir=`cd $bindir; pwd`
here=`pwd`
vfile=`cd \`dirname $vfile\`; pwd`/`basename $vfile`
prem=$here/../prem/prem.dat
tmp=tmp.$$
if command -v gawk > /dev/null 2>&1;then
    awk=gawk
else
    awk=awk
fi
mkdir -p $tmp
cd $tmp
# insert a layer boundary 0.01 above the third viscosity radius
$awk '{print($0);if(NR==3)print($1+0.01,$2)}' $vfile > visc.split
for L in $lmaxs;do
    $awk -v L=$L -f $here/gen_synth_dens.awk > dens.ab
    $bindir/hc -dens dens.ab -dshs -vf $vfile -prem $prem -fs $extra > /dev/null 2>&1
    mv geoid.ab geoid.1.ab
    $bindir/hc -dens dens.ab -dshs -vf visc.split -prem $prem -fs $extra > /dev/null 2>&1
    printf "lmax=%3i  " $L
    $awk -f $here/cmp_ab.awk -v tol=1 geoid.1.ab geoid.ab
done
cd ..
rm -rf $tmp
