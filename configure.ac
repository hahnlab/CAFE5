AC_INIT(cafexp, 1.0)
 
AC_CONFIG_HEADER([config.h])

AC_LANG(C++)
AC_PROG_CXX([icc g++])
AC_OPENMP

AC_CHECK_LIB([m], [floor])

ax_blas_ok=no

AC_SEARCH_LIBS(sgemm, mkl_intel_lp64,
	[ax_blas_ok=yes;BLAS_LIBS="-lmkl_sequential -lmkl_core"],,
	[-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread])

AC_SEARCH_LIBS(cblas_dgemm, satlas,
	[ax_blas_ok=yes;BLAS_LIBS="-lsatlas";AC_DEFINE([HAVE_ATLAS],[1], [Atlas for matrix multiplication])],,
	[])

if test $ax_blas_ok = no; then
	unset ac_cv_search_cblas_dgemm

	AC_SEARCH_LIBS(cblas_dgemm, openblas,
		[ax_blas_ok=yes;BLAS_LIBS="-lopenblas";AC_DEFINE([HAVE_OPENBLAS],[1], [OpenBLAS for matrix multiplication])],,
		[])
fi 

if test $ax_blas_ok = yes; then
  save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS"
  AC_DEFINE([HAVE_BLAS],[1], [If a matrix multiplication library exists])
fi

AC_SEARCH_LIBS(vdExp, mkl_intel_lp64,
		[AC_DEFINE([HAVE_VECTOR_EXP],[1], [Have a vector exponential function])],,
		[])

AC_DEFINE(NUM_OPTIMIZER_INITIALIZATION_ATTEMPTS, 100, [Number of times optimizer will restart if it fails to find legal values])
AC_DEFINE(LAMBDA_PERTURBATION_STEP_SIZE, 2, [Number of values simulator will generate before modifying lambda value])
AC_DEFINE(OPTIMIZER_STRATEGY_SIMILARITY_CUTOFF, 15, [Optimizer will stop after 15 iterations with no significant change in -ln likelihood]
)
AC_DEFINE(PHASED_OPTIMIZER_PHASE1_ATTEMPTS, 4, [Number of attempts optimizer will make to initialize to a good value])
AC_DEFINE(OPTIMIZER_LOW_PRECISION, 1e-3, Precision of values optimizer will use before abandoning a set of values)
AC_DEFINE(OPTIMIZER_HIGH_PRECISION, 1e-6, Precision of values optimizer will use before abandoning a set of values)

dnl Process Makefile.in to create Makefile
AC_OUTPUT(Makefile)

