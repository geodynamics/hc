#!/bin/bash
#
# build hc with address and undefined behavior sanitizers in a
# separate directory and run the regression tests with those
# binaries. requires gcc or clang with libasan.
#
# usage: ./run_asan.sh [build_dir] [extra options for the tests]
#
bdir=${1-/tmp/hc_asan_build}
extra=${2-""}
mkdir -p $bdir
make -C .. all_no_gmt HC_HOME=$bdir ADD_FLAGS=-DNO_GMT \
     CFLAGS="-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer" \
     LDFLAGS="-fsanitize=address,undefined" > $bdir/make.log 2>&1 || \
    { echo "$0: build failed, see $bdir/make.log"; exit 1; }
# leak reports at exit are not useful here
ASAN_OPTIONS=detect_leaks=0:abort_on_error=1 UBSAN_OPTIONS=print_stacktrace=1 \
	    ./run_test.sh $bdir/bin compare 1e-5 "-solver 0" && \
ASAN_OPTIONS=detect_leaks=0:abort_on_error=1 UBSAN_OPTIONS=print_stacktrace=1 \
	    ./run_test.sh $bdir/bin compare 1e-5 "-solver 1 $extra"
