#!/bin/bash
#
# regression tests for hc, using only the data in example_data/
#
# usage: ./run_test.sh [bindir] [mode] [tol]
#
#   bindir: directory with the hc binaries (../bin)
#   mode:   compare (default): compare with stored reference output
#           create: (re)write the reference output
#   tol:    relative tolerance for the comparison (1e-6)
#   extra:  additional options passed to hc and hc_visc_scan for all
#           runs, e.g. "-solver 1" to test the multiple shooting
#           solver against the propagator matrix references ("")
#
bindir=${1-../bin}
mode=${2-compare}
tol=${3-1e-6}
extra=${4-""}

# absolute paths, since the runs take place in a temporary directory
bindir=`cd $bindir; pwd`
here=`pwd`
hc=$bindir/hc
scan=$bindir/hc_visc_scan
data=$here/../example_data
prem=$here/../prem/prem.dat
ref=reference
tmp=tmp.$$
if command -v gawk > /dev/null 2>&1;then
    awk=gawk
else
    awk=awk
fi
awkcmp="$awk -f $here/cmp_ab.awk -v tol=$tol"
extract=$bindir/hc_extract_sh_layer

common="-dens $data/smean.31.m.ab -dshs -vf $data/visc.sh08 -prem $prem $extra"

mkdir -p $tmp $ref
nfail=0
ntest=0

# compare or store one output file
#
# $1: case name, $2: file produced in $tmp, $3: file name in $ref
check(){
    name=$1
    out=$tmp/$2
    rfile=$ref/$3
    ((ntest=ntest+1))
    # remove trailing white space (keeps the reference files clean)
    sed 's/[[:space:]]*$//' $out > $out.tmp
    mv $out.tmp $out
    if [ $mode = create ];then
	cp $out $rfile
	echo "$name: stored $rfile"
    else
	if [ ! -s $rfile ];then
	    echo "$name: reference $rfile missing"
	    ((nfail=nfail+1))
	elif [ ! -s $out ];then
	    echo "$name: output $out missing"
	    ((nfail=nfail+1))
	else
	    printf "%-28s %-14s " $name $3
	    if ! $awkcmp $rfile $out;then
		((nfail=nfail+1))
	    fi
	fi
    fi
}
# extract layers 2 (above CMB), 12 and 22 (surface) in r and pol/tor
# mode from a binary solution file $1 into $1.l.dat
extract_layers(){
    rm -f $1.l.dat
    for layer in 2 12 22;do
	$extract $1 $layer 3 >> $1.l.dat 2> /dev/null
    done
}

cd $tmp
#
# 1) free slip, velocities and geoid at all layers
#
$hc $common -fs -ag > /dev/null 2>&1
extract_layers vel.sol.bin
cd ..
check free_slip        vel.sol.bin.l.dat  fs.vel.sol.dat
check free_slip        geoid.ab     fs.geoid.ab
#
# 2) free slip, radial tractions
#
cd $tmp; rm -f *
$hc $common -fs -rtrac > /dev/null 2>&1
extract_layers rtrac.sol.bin
cd ..
check free_slip_rtrac  rtrac.sol.bin.l.dat fs.rtrac.sol.dat
#
# 3) no slip
#
cd $tmp; rm -f *
$hc $common -ns > /dev/null 2>&1
extract_layers vel.sol.bin
cd ..
check no_slip          vel.sol.bin.l.dat  ns.vel.sol.dat
check no_slip          geoid.ab     ns.geoid.ab
#
# 4) prescribed plate velocities, with and without density anomalies
#
cd $tmp; rm -f *
$hc $common -pvel $data/pvel.sh.dat > /dev/null 2>&1
extract_layers vel.sol.bin
cd ..
check plate_vel        vel.sol.bin.l.dat  pv.vel.sol.dat
check plate_vel        geoid.ab     pv.geoid.ab
cd $tmp; rm -f *
$hc $common -pvel $data/pvel.sh.dat -ds 0 -rtrac > /dev/null 2>&1
extract_layers rtrac.sol.bin
cd ..
check plate_vel_only_rtrac rtrac.sol.bin.l.dat pvo.rtrac.sol.dat
#
# 5) depth dependent density scaling, mean density instead of PREM
#
cd $tmp; rm -f *
$hc $common -fs -dsf $here/scale.2.dat -dnp > /dev/null 2>&1
extract_layers vel.sol.bin
cd ..
check dens_scale_file  vel.sol.bin.l.dat  dsf.vel.sol.dat
check dens_scale_file  geoid.ab     dsf.geoid.ab
#
# 6) geoid kernels
#
cd $tmp; rm -f *
$hc $common -fs -pk > /dev/null 2>&1
cd ..
check kernels          geoidk.dat   kernel.geoidk.dat
#
# 7) CMB solver kludge (only applies to the propagator matrix solver,
#    the multiple shooting solver ignores -cbckl)
#
if [ "`echo $extra | grep -c 'solver 1'`" -eq 0 ];then
    cd $tmp; rm -f *
    $hc $common -fs -cbckl 15 -ascii > /dev/null 2>&1
    head -1500 vel.sol.dat > vel.sol.head.dat	# first few layers only
    cd ..
    check cmb_kludge       geoid.ab     kludge.geoid.ab
    check cmb_kludge_ascii vel.sol.head.dat  kludge.vel.sol.dat
fi
#
# 8) viscosity scan, two layers, coarse
#
cd $tmp; rm -f *
$scan -gref $here/../geoid_ref.ab -dens $data/smean.31.m.ab -dshs -prem $prem \
      -vs_n 2 -vs_dv 1 -vs_em 1 $extra > scan.dat 2> /dev/null
cd ..
check visc_scan        scan.dat     scan.dat
#
# 9) consistency: inserting a viscosity layer boundary without a
#    viscosity change must not change the solution (no reference
#    needed)
#
cd $tmp; rm -f *
$awk '{print($0);if(NR==3)print($1+0.01,$2)}' $data/visc.sh08 > visc.split
$hc $common -fs > /dev/null 2>&1
mv geoid.ab geoid.orig.ab
$hc -dens $data/smean.31.m.ab -dshs -vf visc.split -prem $prem -fs $extra > /dev/null 2>&1
printf "%-28s %-14s " visc_layer_split geoid.ab
if ! $awkcmp geoid.orig.ab geoid.ab;then
    ((nfail=nfail+1))
fi
((ntest=ntest+1))
cd ..
#
# 10) same for a synthetic lmax = 127 density model with the
#     multiple shooting solver, which has to be stable at high degrees
#     (the propagator solver gives order unity differences here)
#
cd $tmp; rm -f *
$awk '{print($0);if(NR==3)print($1+0.01,$2)}' $data/visc.sh08 > visc.split
$awk -v L=127 -f $here/gen_synth_dens.awk > dens.ab
$hc -dens dens.ab -dshs -vf $data/visc.sh08 -prem $prem -fs -solver 1 > /dev/null 2>&1
mv geoid.ab geoid.orig.ab
$hc -dens dens.ab -dshs -vf visc.split -prem $prem -fs -solver 1 > /dev/null 2>&1
printf "%-28s %-14s " ms_layer_split_l127 geoid.ab
if ! $awk -f $here/cmp_ab.awk -v tol=1e-8 geoid.orig.ab geoid.ab;then
    ((nfail=nfail+1))
fi
((ntest=ntest+1))
cd ..

rm -rf $tmp
if [ $mode = compare ];then
    echo "$ntest tests, $nfail failed"
    exit $nfail
fi
