dnl ######################################################################
dnl Basic autotools settings.
dnl ######################################################################
AC_PREREQ([2.64])
dnl Update the corresponding entry in doc/Doxyfile when changing this
AC_INIT([yaxt],[0.4.4],[hanke@dkrz.de],,
  [http://www.dkrz.de/redmine/projects/yaxt])
AC_CONFIG_SRCDIR([src/xt_idxlist.c])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([m4])
dnl forbid acx macro names from remaining unexpanded
m4_pattern_forbid([^ACX_])
AC_CANONICAL_HOST
AM_INIT_AUTOMAKE([1.10 foreign subdir-objects parallel-tests])
AM_MAINTAINER_MODE
dnl ######################################################################
dnl Diagnose incompatible tool versions
dnl ######################################################################
m4_if(m4_version_compare(AC_AUTOCONF_VERSION,[2.67]),1,
  [m4_if(m4_version_compare(LT_PACKAGE_VERSION,[2.4.2]),[-1],[m4_fatal([autoconf versions 2.68 and newer require using libtool 2.4.2 or newer])])])
dnl ######################################################################
dnl Checks for programs.
dnl ######################################################################
AC_PROG_CC
AC_PROG_CC_C99
AX_PTHREAD([AC_DEFINE([HAVE_PTHREAD],1,[Define if you have POSIX threads libraries and header files.])
   LIBS="$PTHREAD_LIBS $LIBS"
   CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
   CC="$PTHREAD_CC"])
AX_TLS(,[ac_cv_tls=`$CC -qversion 2>&1 | sed -n '/^IBM XL C/{
n
s/^Version: \(@<:@0-9@:>@*\).*/\1/
t print
b
: print
p
}'`
   AS_IF([test x"$ac_cv_tls" == x],
       [ac_cv_tls=none],
     [test "$ac_cv_tls" -gt 7],
       [saved_CFLAGS=$CFLAGS
        CFLAGS=`echo "$CFLAGS" | sed -n '/.*-qtls\(=@<:@^ @:>@*\)\{0,1\}/{
p
q
}
s/$/ -qtls=initial-exec/
p
q
'`
dnl unless the user already set the -qtls option, add it and retry test
       AS_IF([test x"$CFLAGS" = x"$saved_CFLAGS"],
             [ac_cv_tls=none],
             [AC_MSG_NOTICE([retrying with -qtls=initial-exec added to CFLAGS])
              AS_UNSET([ac_cv_tls])
              AX_TLS(,[CFLAGS=$saved_CFLAGS])])
],[ac_cv_tls=none])])
AC_PROG_FC
AC_FC_SRCEXT([f90])
AC_FC_FREEFORM
AC_LANG_PUSH([Fortran])
ACX_PROG_FC_CHECK_FPP
FCFLAGS="$FCFLAGS $FC_FPP_FLAG"
AC_PROG_FPP
AC_LANG_POP([Fortran])
AC_PROG_MKDIR_P
AX_WITH_PERL([/bin/false], $PATH:$prefix/bin:/usr/bin:/usr/local/bin)
AX_PROG_PERL_VERSION([5.0005],[],
        AC_MSG_ERROR([Perl (version>=5.0005) is required]))
# Set up libtool.
dnl before switching on libtool, identify compilers that prevent us from
dnl certain build configurations
ACX_LT_PROBLEMS
LT_INIT
dnl _KPSE_USE_LIBTOOL ensures libtool is also used for configure-time tests,
dnl which deduces dependent libraries automatically
_KPSE_USE_LIBTOOL
AC_SUBST([lt_cv_nm_interface])
dnl ######################################################################
dnl Properties of the compiler
dnl ######################################################################
AM_PROG_CC_C_O
AC_PROG_FC_C_O
ACX_SL_FC_MOD_SUFFIX(,[AC_MSG_ERROR(dnl
[Fortran module file extension must be known to install F90 modules.])])
AM_CONDITIONAL([FC_MOD_UPPERCASE],
  [test x$ac_fc_mod_uppercase = xyes])
AC_LANG_PUSH([Fortran])
ACX_SL_FC_CHECK_MOD_PATH_FLAG
AC_LANG_POP([Fortran])
ACX_FC_VOLATILE([FCFLAGS="${FCFLAGS} ${FPP_DEFOPT}HAVE_VOLATILE"])
AC_FC_DUMMY_MAIN(,[:])
AC_FC_WRAPPERS
dnl
dnl Check wether char is unsigned by default
AC_C_CHAR_UNSIGNED
AS_IF([test x$ac_cv_c_char_unsigned = xyes],
  [C_CHAR_IS_UNSIGNED=1],
  [C_CHAR_IS_UNSIGNED=0])
AC_SUBST([C_CHAR_IS_UNSIGNED])
ACX_CHECK_SIZEOF_RELATION([long], [int])
ACX_CHECK_SIZEOF_RELATION([long long], [long])
AC_CHECK_SIZEOF([unsigned long])
dnl
AS_FOR([builtin_macro],[builtin],[__builtin_ctzl __builtin_clzl],
  [AS_VAR_PUSHDEF([builtin_cache],[acx_cv_have_decl_]builtin_macro)
   AC_CACHE_CHECK([whether ]builtin_macro[ is declared],
     [builtin_cache],
     [AC_LINK_IFELSE([AC_LANG_PROGRAM(,[  long lbz = builtin_macro][[(56L)]])],
        [AS_VAR_SET([builtin_cache],[yes])],
        [AS_VAR_SET([builtin_cache],[no])])])
   AS_VAR_IF([builtin_cache],[yes],
     [AC_DEFINE_UNQUOTED([HAVE_DECL_]AS_TR_CPP([builtin_macro]),[1])
      break],
     [AC_DEFINE_UNQUOTED([HAVE_DECL_]AS_TR_CPP([builtin_macro]),[0])])])
AH_TEMPLATE([HAVE_DECL___BUILTIN_CTZL],
  [Define to 1 if __builtin_ctzl is available, 0 if not])
AH_TEMPLATE([HAVE_DECL___BUILTIN_CLZL],
  [Define to 1 if __builtin_clzl is available, 0 if not])
dnl ######################################################################
dnl Test whether cfortran.h works correctly
dnl ######################################################################
ACX_FIND_CFORTRAN_DEF
ACX_XLF_QEXTNAME_ADD_APPENDUS
ACX_CHECK_CFORTRAN

dnl ######################################################################
dnl Checks for header files.
dnl ######################################################################
dnl
AC_HEADER_STDC

#  ----------------------------------------------------------------------
# Compile with MPI support
dnl
dnl ######################################################################
dnl Checks for MPI.
dnl ######################################################################
ACX_C_PACKAGE([MPI],[mpi.h],,,
  [AC_MSG_ERROR([Required header mpi.h not found or not compilable.])],
  [MPI_Waitall],[mpi mpich],,,
  [AC_MSG_ERROR([Cannot link C MPI programs.])])
ACX_FORTRAN_PACKAGE([MPI],[mpif.h],,,
  [AC_MSG_ERROR([Required include mpif.h not found or not compilable.])],
  [mpi_waitall],[mpi mpi_f90 mpi_f77 mpich],[[-lmpi_f77 -lmpi],[-lmpi]],,
  [AC_MSG_ERROR([Cannot link Fortran MPI programs.])],[])
ACX_F90_PACKAGE([MPI],[mpi],,,
  [AC_MSG_ERROR([Required Fortran 90 module not found or incompatible.])],
  [mpi_waitall],[mpi mpi_f90 mpi_f77 mpich],[[-lmpi_f77 -lmpi],[-lmpi]],,
  [AC_MSG_ERROR([Cannot link Fortran 90 MPI programs.])],,
[      INTEGER :: req(1), stat(mpi_status_size, 1), ierror
],[(1, req, stat, ierror)])

save_CFLAGS=$CFLAGS
CFLAGS="$CFLAGS $MPI_C_INCLUDE"
TJ_FIND_INTEGRAL_TYPE([MPI_Fint],[MPI_FINT],[@%:@include <mpi.h>])
CFLAGS=$save_CFLAGS
AS_CASE(["$MPI_FINT"],
  [int],[XT_MPI_FINT_MAX=INT_MAX XT_MPI_FINT_FC_KIND=c_int],
  [long],[XT_MPI_FINT_MAX=LONG_MAX XT_MPI_FINT_FC_KIND=c_long],
  ["long long"],[XT_MPI_FINT_MAX=LLONG_MAX XT_MPI_FINT_FC_KIND=c_long_long],
  [short],[XT_MPI_FINT_MAX=SHRT_MAX XT_MPI_FINT_FC_KIND=c_short])
AC_DEFINE_UNQUOTED([XT_MPI_FINT_MAX],[$XT_MPI_FINT_MAX],
  [Maximum value of MPI_Fint.])
# Match MPI_Fint in reverse to make sure Fortran/C-APIs match
FCFLAGS="$FCFLAGS ${FPP_DEFOPT}XT_MPI_FINT_FC_KIND=$XT_MPI_FINT_FC_KIND"
# the following is needed for source parts from libraries with
# optional MPI dependencies
AC_DEFINE([USE_MPI],[1],[use MPI parallel process environment])
saved_CFLAGS=$CFLAGS
saved_FCFLAGS=$FCFLAGS
saved_LIBS=$LIBS
CFLAGS="$CFLAGS $MPI_C_INCLUDE"
FCFLAGS="$FCFLAGS $MPI_FC_INCLUDE"
LIBS="$LIBS $MPI_C_LIB"
ACX_MPIRUN
   # check MPI defects
   AC_ARG_WITH([regard-for-quality],
     [AS_HELP_STRING([--without-regard-for-quality],
        [ignore results of checks for known defects @<:@default: abort on error@:>@])],
     [],[with_regard_for_quality=yes])
ACX_MPI_DEFECTS(,,
  [ACX_MPI_DEFECTS_DOCUMENT
   AS_IF([test "$with_regard_for_quality" != no],
     [AC_MSG_FAILURE([test for known defect $subtestname failed, re-configure with --without-regard-for-quality in case you must use the given MPI implementation, but some tests might fail])])])
AS_IF([test "x$MPI_LAUNCH" = xtrue],
  [AC_MSG_WARN([MPI launch command unavailable])])
FCFLAGS=$saved_FCFLAGS
CFLAGS=$saved_CFLAGS
LIBS=$saved_LIBS

# enable user to select an alternative enumeration type in case
# 2^31 elements isn't enough or too much
AC_ARG_WITH([idxtype],
  AS_HELP_STRING([--with-idxtype],[select type to use for global index values @<:@default=int@:>@]),
  [XT_INT=$withval
   AC_COMPILE_IFELSE([AC_LANG_PROGRAM(,[[char a[((($XT_INT)-1)<0)];
char b[($XT_INT)((-($XT_INT)-0.5)+1.5)]]])],,
     [AC_MSG_ERROR([Selected index type not compilable or unsigned])])],
  [XT_INT=int])
m4_foreach([substtype],[[XT_INT]],
  [AS_CASE(["$substtype"],
     [int],
     substtype[_MPIDT=MPI_INT
     AC_CHECK_SIZEOF([int])
     ]substtype[_FC_KIND=c_int
     ACX_FORTRAN_TYPE_KIND([integer],[c_int],[int],
       [      use iso_c_binding],,substtype[_FC_KIND=$acx_fortran_kind_subst])
     ]substtype[_MAX=INT_MAX
     ]substtype[_MIN="(-INT_MAX)"
     ]substtype[_FMT="\"d\""
     ]substtype[_abs='(abs((v)))'],
     [long],
     substtype[_MPIDT=MPI_LONG
     AC_CHECK_SIZEOF([long])
     ]substtype[_FC_KIND=c_long
     ACX_FORTRAN_TYPE_KIND([integer],[c_long],[long],
       [      use iso_c_binding],,substtype[_FC_KIND=$acx_fortran_kind_subst])
     ]substtype[_MAX=LONG_MAX
     ]substtype[_MIN="(-LONG_MAX)"
     ]substtype[_FMT="\"ld\""
     ]substtype[_abs='(labs((v)))'],
     [short],
     substtype[_MPIDT=MPI_SHORT
     AC_CHECK_SIZEOF([short])
     ]substtype[_FC_KIND=c_short
     ACX_FORTRAN_TYPE_KIND([integer],[c_short],[short],
       [      use iso_c_binding],,substtype[_FC_KIND=$acx_fortran_kind_subst])
     ]substtype[_MAX=SHRT_MAX
     ]substtype[_MIN="(-SHRT_MAX)"
     ]substtype[_FMT="\"hd\""
     ]substtype[_abs='((short)abs((v)))'],
     ["long long"],
     substtype[_MPIDT=MPI_LONG_LONG
     AC_CHECK_SIZEOF([long long])
     ]substtype[_FC_KIND=c_long_long
     ACX_FORTRAN_TYPE_KIND([integer],[c_long_long],[long long],
       [      use iso_c_binding],,substtype[_FC_KIND=$acx_fortran_kind_subst])
     ]substtype[_MAX=LLONG_MAX
     ]substtype[_MIN="(-LLONG_MAX)"
     ]substtype[_FMT="\"lld\""
     ]substtype[_abs='(llabs((v)))'],
     [AC_MSG_ERROR([Invalid type selected for substtype])])
   save_FCFLAGS=$FCFLAGS
   AC_LANG_PUSH([Fortran])
   FCFLAGS="$MPI_FC_INCLUDE $FCFLAGS"
   for mpi_int2test in MPI_INTEGER2 MPI_INTEGER4 MPI_INTEGER8; do
     AC_MSG_CHECKING([if Fortran MPI bindings define $mpi_int2test])
     AC_COMPILE_IFELSE(AC_LANG_PROGRAM([
        IMPLICIT NONE
        INCLUDE 'mpif.h'
        INTEGER :: i
        i = $mpi_int2test
]),[AC_MSG_RESULT([yes])
    eval HAVE_FC_$mpi_int2test=1],
   [AC_MSG_RESULT([no])
    ASX_VAR_UNSET([HAVE_FC_$mpi_int2test])])
   done
   FCFLAGS=$save_FCFLAGS
   AC_LANG_POP([Fortran])
   AS_IF([eval eval test x\\\"\\\${HAVE_FC_MPI_INTEGER[]\$ac_cv_sizeof_]AS_TR_SH([$substtype])[+set}\\\" = xset],
     [eval substtype[]_FC_MPIDT=\"MPI_INTEGER[]\$ac_cv_sizeof_]AS_TR_SH([$substtype])[\"],
     substtype[_FC_MPIDT=MPI_DATATYPE_NULL])
   FCFLAGS="$FCFLAGS ${FPP_DEFOPT}substtype[]_FC_KIND=$substtype[]_FC_KIND ${FPP_DEFOPT}substtype[]_FC_MPIDT=$substtype[]_FC_MPIDT"
   AC_DEFINE_UNQUOTED([CONF_]substtype[_MAX],[$]substtype[_MAX],
     [Maximal value for ]substtype)
   AC_DEFINE_UNQUOTED([CONF_]substtype[_MIN],[$]substtype[_MIN],
     [Minimal value for ]substtype)
   AC_DEFINE_UNQUOTED([CONF_]substtype[_FMT],[$]substtype[_FMT],
     [printf(3) format string for ]subsstype)
   AC_DEFINE_UNQUOTED(substtype[_ABS(v)],[$]substtype[_abs],
     [return absolute value of v])
  ])

AS_CASE([$XT_INT],
  [int],[XT_UINT_MPIDT=MPI_UNSIGNED],
  [short],[XT_UINT_MPIDT=MPI_UNSIGNED_SHORT],
  [long],[XT_UINT_MPIDT=MPI_UNSIGNED_LONG],
  ["long long"],[XT_UINT_MPIDT=MPI_UNSIGNED_LONG_LONG])

AC_DEFINE_UNQUOTED([XT_INT],[$XT_INT],[Type to use for enumeration of distributed elements.])
AC_DEFINE_UNQUOTED([XT_INT_MPIDT],[$XT_INT_MPIDT],[MPI Datatype for type for enumeration of distributed elements.])
AC_DEFINE_UNQUOTED([XT_UINT_MPIDT],[$XT_UINT_MPIDT],[MPI Datatype for type for counts of distributed elements.])


AC_ARG_VAR([BUILD_CFLAGS],
  [append to CFLAGS during build but not in configure phase])
AC_CONFIG_COMMANDS_PRE([CFLAGS="$CFLAGS${BUILD_CFLAGS+ $BUILD_CFLAGS}"])
AC_ARG_VAR([BUILD_FCFLAGS],
  [append to FCFLAGS during build but not in configure phase])
AC_CONFIG_COMMANDS_PRE([FCFLAGS="$FCFLAGS${BUILD_FCFLAGS+ $BUILD_FCFLAGS}"])
AC_ARG_VAR([BUILD_LDFLAGS],
  [append to LDFLAGS during build but not in configure phase])
AC_CONFIG_COMMANDS_PRE([LDFLAGS="$LDFLAGS${BUILD_LDFLAGS+ $BUILD_LDFLAGS}"])
AC_ARG_VAR([BUILD_LIBS],
  [prepend to LIBS during build but not in configure phase])
AC_CONFIG_COMMANDS_PRE([LIBS="${BUILD_LIBS+$BUILD_LIBS }$LIBS"])
AC_ARG_VAR([BUILD_CC],
  [replace CC with expansion of $BUILD_CC during build but not in configure phase])
AC_CONFIG_COMMANDS_PRE([CC="${BUILD_CC-$CC}"])
AC_ARG_VAR([BUILD_FC],
  [replace FC with expansion of $BUILD_FC during build but not in configure phase])
AC_CONFIG_COMMANDS_PRE([FC="${BUILD_FC-$FC}"])

AC_CONFIG_HEADER([include/config.h])

AC_CONFIG_FILES([						\
	perf/perf_idxsection_get_positions_of_indices_run	\
	perf/perf_sparse_array_gather_run			\
	tests/test_exchanger_parallel_run			\
	tests/test_redist_p2p_parallel_run			\
	tests/test_redist_collection_parallel_run		\
	tests/test_redist_collection_static_parallel_run	\
	tests/test_redist_repeat_parallel_run			\
	tests/test_xmap_all2all_parallel_run			\
	tests/test_xmap_dist_dir_parallel_run			\
	tests/test_xmap_intersection_parallel_run		\
        tests/test_initialized_finalized_run			\
	tests/test_idxempty_run					\
	tests/test_idxvec_run					\
	tests/test_idxlist_collection_run			\
	tests/test_idxmod_run					\
	tests/test_idxsection_run				\
	tests/test_idxstripes_run				\
	tests/test_redist_collection_run			\
	tests/test_redist_collection_static_run			\
	tests/test_redist_repeat_run				\
	tests/test_redist_p2p_run				\
	tests/test_exported_symbols				\
	tests/test_xmap_dist_dir_run				\
	tests/test_xmap_all2all_run				\
	tests/test_xmap_all2all_fail_run			\
	tests/test_handles_run					\
	tests/test_ut_run					\
        tests/test_perf_run					\
        tests/test_perf_stripes_run				\
        tests/test_sort_run					\
        tests/test_uid_run					\
        tests/test_yaxt_run					\
	tests/test_mpi_generate_datatype_run			\
	tests/test_mpi_smartdedup_run				\
	examples/rrobin_run					\
	examples/row2col_run					\
	examples/tr_gp2fs_run					\
	examples/unstruct_halo_exchange_run			\
	util/serialrun],[chmod a+x "$ac_file"])
AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile examples/Makefile \
  inst_headers/Makefile inst_headers/f90/Makefile \
  src/pkgconfig/yaxt.pc perf/Makefile])
# force make to regenerate FC.deps file (because it depends on Makefile)
AC_CONFIG_COMMANDS([force_dep_regeneration],[sleep 1
touch tests/Makefile src/Makefile])

AC_OUTPUT
