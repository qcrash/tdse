AC_PREREQ([2.68])
AC_INIT([tdse], [1.0])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_PROG_FC
AC_CONFIG_FILES([Makefile])
AC_LANG_(Fortran)
AC_FC_SRCEXT(f90)
## The AC_ FC_FREEFORM tries to ensure that the Fortran compiler ($FC) allows free-format source code (as opposed to the older fixed-format style from Fortran 77). If necessary, it may add some additional flags to FCFLAGS.
AC_FC_FREEFORM
## The AC_FC_LINE_LENGTH macro tries to ensure that the Fortran compiler ($FC) accepts long source code lines
AC_FC_LINE_LENGTH
#AC_ FC_LIBRARY_LDFLAGS
# Add -m64 -fdefault-real-8 flags
AX_CHECK_COMPILE_FLAG([-m64], [FCFLAGS="$FCFLAGS -m64"])
AX_CHECK_COMPILE_FLAG([-fdefault-integer-8], [FCFLAGS="$FCFLAGS -fdefault-integer-8"])
## Detect system LAPACK library, or use internal functions
AC_ARG_WITH([lapack],
            [AS_HELP_STRING([--without-lapack], [ignore LAPACK/BLAS libs and use internal])])
AS_IF([test "x$with_lapack" != "xno"], [], [
  with_blas="no"
  with_lapack="no"
  ])
## Looks like openblas has pkg-config file, so use it
## and intel mkl has them too, but hidden, so don't expect them.
AS_CASE([$with_blas],
    [no],[AM_CONDITIONAL([HAVE_openblas],[false])],
    [PKG_HAVE_WITH_MODULES([openblas],[openblas >= 0.2])])
##PKG_HAVE_WITH_MODULES([mkl-dynamic-lp64-seq],[mkl-dynamic-lp64-seq >= 2018], [])
AS_IF([test "x$with_openblas" = "xyes"],[
       ax_package_requires="openblas >= 0.2" # have pkg-config support
       ax_blas_ok=yes # disables testing later
       ax_lapack_ok=yes # disables testing later
      ])
## If pkg-config openblas failed, discover BLAS/LAPACK
AS_IF([test "x$with_openblas" = "xno"],[AX_BLAS()])
AS_IF([test "x$with_openblas" = "xno"],[AX_LAPACK()])
## Add whatever AX_BLAS discovered to pc file
## and what the Intel MKL Link Line Advisor says (for a few cases)
AS_CASE([$BLAS_LIBS],
  [-lmkl_gf_*lp64*],[
	package_extra_libs="-Wl,--no-as-needed $LAPACK_LIBS $BLAS_LIBS $LIBS"
	package_extra_flags="-m64"
	],
  [-lmkl*lp64*],[
	package_extra_libs="$LAPACK_LIBS $BLAS_LIBS $LIBS"
	package_extra_flags=""
	],
  [-l*],  [package_extra_libs="$LAPACK_LIBS $BLAS_LIBS $LIBS"],
  [""],   [],
  )
AM_CONDITIONAL([HAVE_BLAS], [test x"$ax_blas_ok" = xyes])
AM_CONDITIONAL([HAVE_LAPACK], [test x"$ax_lapack_ok" = xyes])
