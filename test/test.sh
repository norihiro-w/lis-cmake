#!/bin/sh

if test -z "$srcdir"; then
  srcdir=`echo "$0" | sed 's,[^/]*$,,'`
  test "$srcdir" = "$0" && srcdir=.
  test -z "$srcdir" && srcdir=.
  test "${VERBOSE+set}" != set && VERBOSE=1
fi
  VERBOSE=1
. $srcdir/defs

OMPNUM=1
if test "$enable_mpi" = "yes"; then
  if test -z "$MPIRUN"; then
    MPIRUN="mpirun -np 2"
  else
    MPIRUN="$MPIRUN $MPINP 2"
  fi
else
  MPIRUN=""
fi
if test "$enable_omp" = "yes"; then
  OMPNUM=2
fi

echo ' '
echo 'checking linear solvers...'
$MPIRUN $srcdir/test1 $srcdir/testmat.mtx 0 $srcdir/sol.txt $srcdir/res.txt -omp_num_threads $OMPNUM

echo 'checking eigensolvers...'
$MPIRUN $srcdir/etest1 $srcdir/testmat.mtx 0 $srcdir/sol.txt $srcdir/res.txt -omp_num_threads $OMPNUM

if test "$enable_fortran" = "yes"; then
  echo 'checking Fortran interface...'
#  $MPIRUN $srcdir/test1f $srcdir/testmat.mtx 0 $srcdir/sol.txt $srcdir/res.txt -omp_num_threads $OMPNUM
  cd $srcdir; 
  $MPIRUN ./test4f
fi

if test "$enable_quad" = "yes"; then
  echo ' '
  echo 'cheking quad precision...'
  $MPIRUN $srcdir/test5 200 2.0 -f double -omp_num_threads $OMPNUM
  $MPIRUN $srcdir/test5 200 2.0 -f quad -omp_num_threads $OMPNUM
fi

OMPNUM=1
if test "$enable_saamg" = "yes"; then
  echo ' '
  echo 'checking SAAMG preconditioner...'
  $MPIRUN $srcdir/test2 10 10 1 $srcdir/sol.txt $srcdir/res.txt -i cg -p saamg -omp_num_threads $OMPNUM
#  $MPIRUN $srcdir/test1 $srcdir/add32.mtx $srcdir/sol.txt $srcdir/res.txt -i bicgstab -p saamg -saamg_unsym true -omp_num_threads $OMPNUM
fi

