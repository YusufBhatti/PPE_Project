dnl Process this file with autoconf to produce a configure script.
dnl
dnl Luis Kornblueh, MPI-M
dnl Uwe Schulzweida, MPI-M
dnl Thomas Jahns, DKRZ
dnl _______________________________________________________________________________________________

AC_INIT([echam],[6.3.00],[luis.kornblueh@mpimet.mpg.de],[echam-6.3.00])
AC_PREREQ(2.69)

dnl _______________________________________________________________________________________________

AC_PREFIX_DEFAULT([$(pwd)])

AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([config/config.h])
AC_CONFIG_HEADERS([config/config.inc:config/config.h.in],
                  [sed 's%^   %! %
                        s%^/\* %! %
                        s% \*/$%%' config/config.inc > $tmp/config.inc
                   rm -f config/config.inc
                   AS_IF([diff config/config.inc.check $tmp/config.inc > /dev/null 2>&1],
                         [AC_MSG_NOTICE(config/config.inc is unchanged)
                          mv config/config.inc.check config/config.inc],
                         [mv $tmp/config.inc config/config.inc])
                   rm -f config/config.inc.check],
		  [AS_IF([test -f config/config.inc],
		         [mv config/config.inc config/config.inc.check])])
AC_CANONICAL_SYSTEM

AC_PROG_AWK
AC_PROG_SED
AC_PROG_GREP

dnl _______________________________________________________________________________________________
dnl get host- and sitename
AS_CASE([$ac_hostname],
        [plogin*], [ac_sitename=dkrz.de],
        [btc*], [ac_sitename=dkrz.de],
        [
        AS_IF([test "x$(hostname)" = "x$ac_hostname"],
              [fqdn=$(host -n $ac_hostname | $AWK '/has address/ {print $1}')],
              [fqdn=$(hostname)])
        AS_IF([test -n "$fqdn"],
          [ac_sitename=$(echo $fqdn | $AWK -F. '{i=NF-1; printf("%s.%s",$i,$NF);}')],
          [ac_sitename=])
])

AC_MSG_CHECKING([for hostname])
AC_MSG_RESULT([$[ac_hostname]])
AC_MSG_CHECKING([for site name])
AC_MSG_RESULT([$[ac_sitename]])

dnl _______________________________________________________________________________________________
dnl Checks for package options.

AC_ARG_WITH([fortran], 
  [AS_HELP_STRING([--with-fortran=COMPILER],
    [select_available compiler (see README.compiler)])],
  [],
  [])
AS_IF([test x$with_fortran != x],
      [fortran_compiler=$with_fortran],
      [fortran_compiler="default"])

dnl _______________________________________________________________________________________________
dnl read OS, compiler and site dependent information 

AX_MACHINE_CONFIG

dnl _______________________________________________________________________________________________
dnl init automake and libtool

AM_INIT_AUTOMAKE([1.11]) 
AM_MAINTAINER_MODE([disable])

dnl the following needs to be tested for NEC SX:
ax_is_nec_sx=no
AS_CASE([$host],[sx*-*-*],[ax_is_nec_sx=yes])
AM_CONDITIONAL([NEC_SX_BUILD],[test x$ax_is_nec_sx = xyes])
dnl _______________________________________________________________________________________________
dnl required tools for further search of system type/name

dnl for later installation

AC_PROG_MKDIR_P
AS_IF([test "x$MKDIR_P" = "xconfig/install-sh -c -d"],
      [MKDIR_P="\$(abs_top_srcdir)/$MKDIR_P"],
      [])
AC_PROG_INSTALL
AS_IF([test "x$INSTALL" = "xconfig/install-sh -c"],
      [INSTALL="\$(abs_top_srcdir)/$INSTALL"],
      [])

dnl TODO: OpenMP support temporarily disabled due to possible race condition
dnl AC_ARG_WITH([openmp],
dnl   [AS_HELP_STRING([--with-openmp],
dnl     [enable OpenMP, if available for selected compiler])],[],
dnl     [withval=no])
dnl AS_IF([test "x$withval" == "xno"],
dnl       [use_openmp="not_enabled"],
dnl       [use_openmp="enabled"])

dnl _______________________________________________________________________________________________
dnl
dnl AS_IF([test -z "$CROSS"],
dnl       [AC_CONFIG_HEADER(config/config.h)],
dnl       [])    
dnl
dnl reduced_CFLAGS=$(echo $CFLAGS | sed 's/-g/ /')
dnl if test "x$CFLAGS" != "x$reduced_CFLAGS" ; then
dnl    ac_cflags_debug_option="-g"
dnl    CFLAGS=$reduced_CFLAGS
dnl else
dnl    ac_cflags_debug_option=""
dnl fi
dnl
dnl _______________________________________________________________________________________________
dnl for handling old configurations
dnl
dnl FC=${F90:+$F90}
dnl
dnl _______________________________________________________________________________________________
dnl save the variables disturbing tests for later use 

ac_config_set_libs=$LIBS
LIBS=""

dnl _______________________________________________________________________________________________
dnl look for compiler and their tool sets and charateristics

AC_PROG_CC
AC_PROG_FC
AC_PROG_F77

AC_FC_SRCEXT([f90])
AC_FC_PP_SRCEXT([f90])
AC_FC_PP_DEFINE
AC_FC_MODULE_FLAG
AC_FC_MODULE_EXTENSION

dnl TODO: OpenMP support temporarily disabled due to possible race condition
dnl AC_LANG_PUSH([Fortran])
dnl AC_OPENMP
dnl dnl this is a fix required by xlf so that directives get evaluated
dnl AS_IF([test "x$ac_cv_prog_fc_openmp" = "x-qsmp=omp"],
dnl       [ac_cv_prog_fc_openmp="-qsmp=omp -qdirective"],
dnl       [])	    
dnl AC_LANG_POP([Fortran])

AX_CC_VERSION
AX_FC_VERSION

AS_IF([test "x$ac_cv_fc_pp_srcext_f90" = "xnone"],
      [FPPFLAG=],	    
      [FPPFLAG=$ac_cv_fc_pp_srcext_f90])
FPPFLAG="${FPPFLAG} -DINPUT_IN_ECHAM"
AC_SUBST([FPPFLAG])

dnl TODO: OpenMP support temporarily disabled due to possible race condition
dnl AS_IF([test "x$use_openmp" = "xenabled"],
dnl       [AS_IF([test "x$ac_cv_prog_fc_openmp" = "xnone needed"],
dnl              [OMPFLAG=],	    
dnl              [OMPFLAG=$ac_cv_prog_fc_openmp])],
dnl       [OMPFLAG=])
AC_SUBST([OMPFLAG])

dnl case $host in	
dnl powerpc-ibm-aix*)     
dnl     check_version=$(eval $FC -qversion 2>&1 | awk '/^Version:/ {print $2}' )
dnl     AS_VERSION_COMPARE([$check_version], [13.01.0000.0002],[
dnl         AS_BOX([ERROR: only IBM XLF 13.1.0.2 (EFIXed) or above is supported.])
dnl         rm -f Makefile 
dnl         AS_EXIT(1)])
dnl     ;;  
dnl esac

dnl _______________________________________________________________________________________________
dnl initialize libtool
dnl before switching on libtool, identify compilers that prevent us from
dnl certain build configurations
ACX_LT_PROBLEMS
dnl do libtool tests
LT_INIT
dnl _KPSE_USE_LIBTOOL ensures libtool is also used for configure-time tests,
dnl which deduces dependent libraries automatically
_KPSE_USE_LIBTOOL

dnl _______________________________________________________________________________________________
dnl check for libraries

AX_LIB_ZLIB
AX_LIB_SZIP
AX_LIB_HDF5
AX_LIB_NETCDF4
AX_LIB_POSTGRESQL

dnl _______________________________________________________________________________________________
dnl check for YAXT

AX_LIB_YAXT

dnl _______________________________________________________________________________________________
dnl check for CDI

AX_LIB_CDI

dnl _______________________________________________________________________________________________
dnl check for SCT

SCT_LIB=""
SCT_INCLUDE=""

AC_ARG_WITH([sct],
            [AC_HELP_STRING([--with-sct=DIR],
                            [Enable external timer library])],,
                            withval=no)

if test x$withval != xno; then
  if test x$with_sct != x; then
    SCTROOT=$with_sct
  fi

  AC_CHECK_FILE([$SCTROOT/include/sct.$FC_MODEXT],,
    [AC_MSG_ERROR([SCTROOT not properly defined])])
  SCT_INCLUDE="$FC_MODINC$SCTROOT/include"
  if test -n "$SCTROOT" ; then
    if test -r "$SCTROOT/lib/libsct.a" ; then
       SCT_LIB='-L$(SCTROOT)/lib -lsct'
    else
       echo "cannot read $SCTROOT/lib/libsct.a"
       exit 1
    fi
  else
    echo "error: the variable SCTROOT is not set in $host_frag"
    echo "       set SCTROOT to the SCT root directory"
    exit 1
  fi
fi

dnl _______________________________________________________________________________________________
dnl check for MPI

AC_ARG_WITH([mpi],
  [AS_HELP_STRING([--with-mpi=DIR],
	          [MPI library root directory])],,
                  [with_mpi=yes])
AS_IF([test "x$with_mpi" = "xno"],
      [use_mpi=no],
      [use_mpi=yes
      AS_IF([test "x$with_mpi" != "xyes"],
            [use_MPIROOT=$with_mpi],
            [use_MPIROOT=""])])

AS_IF([test "x$use_mpi" = "xno"],
      [MPI="${FC_DEFINE}NOMPI"
       MPI_LIB=""
       MPI_VERSION=none
       MPIROOT=""
       use_mpi_compiler=no], 
      [AS_CASE([$host],
               [powerpc-ibm-aix*],
               [MPIROOT=/usr/lpp/ppe.poe    
                q64_result="no"
                mt_result="no"
                q64_result=$(echo $FCFLAGS | awk -v b="-q64" '$0 ~ b { print "yes" }')
                mt_result=$(echo $FC | awk -v b="_r" '$0 ~ b { print "yes" } ')
                if [[ x$mt_result != xyes ]]; then
                    echo "error: IBM MPI requires thread safe compiler versions."  
                    echo "       add _r to compiler selected by FC."
                    exit 1
                fi  
                if test x$q64_result = xyes; then
                    MPI_INCLUDE_THREAD="thread64/"      
                else
                    MPI_INCLUDE_THREAD="thread/"    
                fi],
               [*-*-darwin*],
#lk++
           [if test -n "$MPIROOT" ; then
                mpifh=$(find $MPIROOT/include -name mpif.h)
                mpifh=${mpifh%\/mpif.h}
                mpifh=${mpifh##*/}
                if test -n $mpifh; then
                    MPI_INCLUDE_THREAD="$mpifh/"
                fi
            fi
            case $MPI_INCLUDE_THREAD in #(
            include*)
                MPI_INCLUDE_THREAD=""
                ;;
            esac],
#lk++
#lk--
#                [MPI_INCLUDE_THREAD=""],		  
#lk--
                [MPI_INCLUDE_THREAD=""])])

  MPI_INCLUDE=""

  AS_IF([test "x$with_mpi" != x -a "x$with_mpi" != xno -a "x$with_mpi" != xyes],
        [MPIROOT=$with_mpi],
        [])

  AS_IF([test -n "$MPIROOT"],
        [AC_CHECK_FILE([$MPIROOT/include/${MPI_INCLUDE_THREAD}mpif.h],
                       [MPI_C_INCLUDE=-I$MPIROOT/include/${MPI_INCLUDE_THREAD}
                        MPI_FC_INCLUDE=-I$MPIROOT/include/${MPI_INCLUDE_THREAD}],
                       [AC_MSG_ERROR([MPIROOT not properly defined])])
         AS_CASE([$host],		      
                  [*-*-darwin*],		      
		  [MPI_LIB="-L$MPIROOT/lib/mpich-gcc49 "$MPI_LIB],
		  [MPI_LIB="-L$MPIROOT/lib "$MPI_LIB])
         AS_CASE([$host],
                 [powerpc-ibm-aix*],
                 [use_mpi_compiler=yes
                  mod_file=$(find $MPIROOT -name mpi.mod -print 2>/dev/null | grep mpi.mod | grep ${MPI_INCLUDE_THREAD})
                 ],
                 [sx?-nec-superux],
                 [real_path=$(cd $MPIROOT; pwd -P)
                  mod_file=$(find $real_path -name mpi.mod -print | grep -v modulew)
                  mod_file="${MPIROOT}/lib${mod_file##*lib0}"
                 ],
                 [AS_IF([ test "x${FC}" = "xifort"],
                        [real_path=$(cd $MPIROOT; pwd -P)
                         mod_file=$(find $real_path -not \( -path */gfortran/* -prune \) -name mpi.mod -print 2> /dev/null)],
                        [real_path=$(cd $MPIROOT; pwd -P)
                         mod_file=$(find $real_path -name mpi.mod -print 2> /dev/null)])])
         AS_IF([test -r $mod_file],
               [mod_file=${mod_file#$real_path}
                mod_file=${mod_file#$MPIROOT}       
                mod_file="${MPIROOT}${mod_file}"
                MPI_INCLUDE=${FC_MODINC}${mod_file%/*}
                MPI_FC_MOD=${FC_MODINC}${mod_file%/*}],
               [echo "cannot read $mod_file/mpi.mod"  
                exit 1])],    
        [MPI="${FC_DEFINE}NOMPI"
         MPI_LIB=""])

dnl _______________________________________________________________________________________________
dnl Switch off unwanted input routines in JSBACH

FCDEFS="${FCDEFS} ${FC_DEFINE}INPUT_IN_ECHAM"
AC_SUBST([FCDEFS])

dnl _______________________________________________________________________________________________
dnl  check wether we are building with coupling code and oasis3 support

AC_ARG_WITH([coupler],
  [AS_HELP_STRING([--with-coupler],
   [compile for coupled model run @<:@default: off@:>@])],[],[with_coupler=off])

AS_IF([test "x$with_coupler" = "xoasis3mct"],
      [AC_MSG_WARN(["the value for --with-coupler is oasis3-mct - correct this next time!])
       with_coupler="oasis3-mct"],
      [])

COUPLERROOT=
COUPLER_INCLUDE=	
COUPLER_LIB=		

AC_ARG_VAR([OASIS3ROOT],	
  [directory oasis3 is installed in])
AC_ARG_VAR([OASIS3_FC_LIB],
  [libraries needed for linking	with oasis3])

AC_ARG_VAR([OASIS3MCTROOT],
  [directory oasis3-mct is installed in])
AC_ARG_VAR([OASIS3MCT_FC_LIB],
  [libraries needed for linking with oasis3-mct])

AC_ARG_WITH([oasis3mode],
  [AS_HELP_STRING([--with-oasis3mode],
     [MPI communicator mode of oasis3 @<:@default: MPI1@:>@])],
  [],
  [with_oasis3mode=MPI1])

AS_IF([test "x$with_coupler" = "xoasis3" -o "x$with_coupler" = "xoasis3-mct"],
  [AS_CASE([$with_oasis3mode],
           [MPI1|MPI2],[],
           [AC_MSG_FAILURE([Unsupported oasis3 communicator mode $with_oasis3mode!])])])

AS_CASE([$with_coupler],
        [oasis3],[AS_IF([test "x$OASIS3ROOT" = x],
                        [AC_MSG_FAILURE([Option --with-coupler=oasis3 requires OASIS3ROOT to be known!])])
                  ac_oasis3_comm_mode_fcflag="${FC_DEFINE}use_comm_${with_oasis3mode}"		  
   	          ac_oasis3_coupling_fcflag="${FC_DEFINE}__prism ${FC_DEFINE}__cpl_mpiom"
                  ac_oasis3_mct_fcflag=""	
                  ac_oasis3_modinc_fcflag="${FC_MODINC}\$(COUPLERROOT)/build/lib/psmile.${with_oasis3mode}"
		  COUPLERROOT=$OASIS3ROOT
                  COUPLER_INCLUDE="$ac_oasis3_comm_mode_fcflag $ac_oasis3_coupling_fcflag $ac_oasis3_modinc_fcflag"
                  COUPLER_LIB="-L\$(COUPLERROOT)/lib -lpsmile.${with_oasis3mode} -lmpp_io ${OASIS3_FC_LIB}"],
        [oasis3-mct],[AS_IF([test "x$OASIS3MCTROOT" = x],
                            [AC_MSG_FAILURE([Option --with-coupler=oasis3-mct requires OASIS3MCTROOT to be known!])])
 	          AS_IF([test "x$with_oasis3mode" = "xMPI2"],
  	                [AC_MSG_WARN([Option --with-coupler=oasis3-mct requires oasis3 mode MPI1 (default)!])])
                  with_oasis3mode=MPI1
                  ac_oasis3_comm_mode_fcflag="${FC_DEFINE}use_comm_${with_oasis3mode}"		  
   	          ac_oasis3_coupling_fcflag="${FC_DEFINE}__prism ${FC_DEFINE}__cpl_mpiom"
   	          ac_oasis3_mct_fcflag="${FC_DEFINE}__oa3mct"
                  ac_oasis3_modinc_fcflag="${FC_MODINC}\$(COUPLERROOT)/include"
		  COUPLERROOT=$OASIS3MCTROOT
                  COUPLER_INCLUDE="$ac_oasis3_comm_mode_fcflag $ac_oasis3_coupling_fcflag $ac_oasis3_mct_fcflag $ac_oasis3_modinc_fcflag"
                  COUPLER_LIB="-L\$(COUPLERROOT)/lib -lpsmile -lscrip -lmct -lmpeu $OASIS3_FC_LIB"],
        [])	

AC_SUBST([COUPLERROOT])
AC_SUBST([COUPLER_LIB])
AC_SUBST([COUPLER_INCLUDE])

dnl _______________________________________________________________________________________________
dnl check for support 

dnl _______________________________________________________________________________________________
dnl check for BLAS/LAPACK
AC_SUBST([LAPACKROOT])
AC_ARG_VAR([LAPACK_LIB],[Set LAPACK_LIB to the link time arguments
needed for an alternative instead of the internally distributed LAPACK
and BLAS libraries.  Setting this variable causes the internal LAPACK
not to be built.])

AS_IF([test x${LAPACK_LIB+set} != xset],
      [BUILD_INTERNAL_LAPACK=yes
       LAPACK_LIB='$(top_builddir)/lapack/liblapack.la'],
      [test -z "$LAPACK_LIB"],
      [BUILD_INTERNAL_LAPACK=yes
       LAPACK_LIB='$(top_builddir)/lapack/liblapack.la'],
      [BUILD_INTERNAL_LAPACK=no])
AM_CONDITIONAL([BUILD_INTERNAL_LAPACK],[test $BUILD_INTERNAL_LAPACK = yes])

dnl _______________________________________________________________________________________________
dnl settings for quickplots' table printer

FCDEFS="${FCDEFS} ${FC_DEFINE}HAVE_CONFIG_INC"
AC_SUBST([FCDEFS])
AC_ARG_VAR([QPLOT_TABLE_LINES],
           [Number of lines to get page-aligned quickplot tables @<:@53@:>@])
AC_DEFINE_UNQUOTED([QPLOT_TABLE_LINES], [${QPLOT_TABLE_LINES:-53}],
                   [Number of lines to get page-aligned quickplot tables])

dnl _______________________________________________________________________________________________
dnl Checks for header files.

AC_HEADER_STDC   
AC_CHECK_HEADERS_ONCE(fcntl.h fortran.h limits.h malloc.h netdb.h pwd.h unistd.h execinfo.h ucontext.h sys/param.h sys/time.h sys/unistd.h sys/utsname.h)

dnl Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_HEADER_TIME

dnl Checks for library functions.
AC_CHECK_FUNCS(getrusage gettimeofday sysconf uname valloc)

dnl Check for byte ordering

AC_C_BIGENDIAN

dnl check for C datatypes

AC_CHECK_SIZEOF(char)
AC_CHECK_SIZEOF(short) 
AC_CHECK_SIZEOF(int) 
AC_CHECK_SIZEOF(long) 
AC_CHECK_SIZEOF(long long) 
AC_CHECK_SIZEOF(float) 
AC_CHECK_SIZEOF(double) 
AC_CHECK_SIZEOF(long double)
AC_CHECK_SIZEOF([int *])

dnl Check for Fortran naming convention

AC_FC_WRAPPERS

LIBS=$ac_config_set_libs

AS_IF([test -z "$AR"],
      [AR="ar"],[])
AS_IF([test -z "$AS"],
      [AS="as"],[])
AS_IF([test -z "$ARFLAGS"],
      [ARFLAGS="crv"],[])

AC_ARG_ENABLE([profile],
  [AS_HELP_STRING([--enable-profile],
    [add compile options and libraries for profiling (IBM pSeries and NEC-SX, IBM requires additionally --with-hpm=DIR)])],
    [],
    [withval=no])

if test x$withval != xno; then
  case $host in
  powerpc-ibm-aix*)
    AC_ARG_WITH([hpm],
      [AS_HELP_STRING([--with-hpm=DIR],
        [enable experimental support for profiling (IBM hpm])],
      [],
      [withval=no])

    if test x$withval != xno; then
      HPMROOT=$with_hpm
      AC_CHECK_FILE($HPMROOT/include/libhpm.h,[],AC_MSG_ERROR([--with-hpm=DIR not properly defined]))
      HPM_INCLUDE="-I$HPMROOT/include"
      if test x$q64_result = xyes; then
        HPM_LIB="-L$HPMROOT/lib64 -lhpm -lpmapi -lm" 
      else
        HPM_LIB="-L$HPMROOT/lib -lhpm -lpmapi -lm" 
      fi
    else
      AC_MSG_ERROR([--with-hpm=DIR needs to be added])
    fi
    CFLAGS="${CFLAGS} ${HPM_INCLUDE} -WF,-D__hpm__ -WF,-D_PROFILE" 
    FFLAGS="${FFLAGS} ${HPM_INCLUDE} -WF,-D__hpm__ -WF,-D_PROFILE" 
    FCFLAGS="${FCFLAGS} ${HPM_INCLUDE} -WF,-D__hpm__ -WF,-D_PROFILE" 
    LIBS="$LIBS $HPM_LIB"
    ;;      
  sx*-*-*)
    CFLAGS="${CFLAGS} -ftrace -D_PROFILE" 
    FFLAGS="${FFLAGS} -ftrace -D_PROFILE" 
    FCFLAGS="${FCFLAGS} -ftrace -D_PROFILE" 
    ;;
  esac
fi

AC_SUBST(NOMPI_FC,[${NOMPI_FC:-$FC}])

AS_CASE(["$host"],
  [powerpc-ibm-aix*],
  [AS_IF([test "$use_mpi_compiler" = yes],
     [AS_IF([test -x /sw/aix61/smartmpxlf/bin/smartmpcc_r],
        [CC=/sw/aix61/smartmpxlf/bin/smartmpcc_r],
        [CC=mpcc_r])
      AS_IF([test -x /sw/aix61/smartmpxlf/bin/smartmpxlf95_r],
        [FC=/sw/aix61/smartmpxlf/bin/smartmpxlf95_r],
        [FC=mpxlf95_r])
      AS_IF([test -x /sw/aix61/smartmpxlf/bin/smartmpxlf_r],
        [F77=/sw/aix61/smartmpxlf/bin/smartmpxlf_r],
        [F77=mpxlf_r])
      MPIROOT=""
      MPI_INCLUDE=""
      MPI_LIB=""])])

AC_SUBST(AR,$AR)
AC_SUBST(AS,$AS)
AC_SUBST(ARFLAGS,$ARFLAGS)

AC_SUBST(SCTROOT,$SCTROOT)
AC_SUBST(SCT_LIB,$SCT_LIB)
AC_SUBST(SCT_INCLUDE,$SCT_INCLUDE)

AC_SUBST(MPIROOT,$MPIROOT)
AC_SUBST(MPI_LIB,$MPI_LIB)
AC_SUBST(MPI_INCLUDE,$MPI_INCLUDE)
AC_SUBST(MPI,$MPI)

AC_SUBST(SUPPORTROOT,$SUPPORTROOT)
AC_SUBST(SUPPORT_LIB,$SUPPORT_LIB)

AC_SUBST(SRCDIRS,$SRCDIRS)

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([blas/Makefile])
AC_CONFIG_FILES([lapack/Makefile])
AC_CONFIG_FILES([support/Makefile])
AC_CONFIG_FILES([src/Makefile])
AC_CONFIG_FILES([util/Makefile])
AC_CONFIG_FILES([util/quickplots/Makefile])

AC_CONFIG_COMMANDS_POST([ac_configure_args=$ac_sub_configure_args])dnl
dnl
dnl
ac_sub_configure_args=`echo " $ac_configure_args " | sed 's/ '\''FC=[^'\'']*//
s/ '\''FCFLAGS=[^'\'']*//
s/ '\''CC=[^'\'']*//
s/ '\''CXX=[^'\'']*//
s/ '\''CFLAGS=[^'\'']*//
s/ '\''LDFLAGS=[^'\'']*//
s/ '\''LIBS=[^'\'']*//
s/ '\''F77=[^'\'']*//
s/ '\''FFLAGS=[^'\'']*//
s/^ //
s/ $//
'`
ac_sub_configure_args="$ac_sub_configure_args 'CC=$CC' 'CFLAGS=$CFLAGS' 'LDFLAGS=$LDFLAGS' 'LIBS=$LIBS' 'FC=$FC' 'F77=$F77' 'FFLAGS=$FFLAGS' 'CXX=$CXX'"
dnl
dnl when building with internal yaxt, sub-configure FCFLAGS need to account
dnl for the different place of yaxt headers vs. Fortran module files
AS_IF([test x"${enable_cdi_pio}" = xno],
  [ac_sub_configure_args="$ac_sub_configure_args 'FCFLAGS=$FCFLAGS'"],
  [AM_COND_IF([BUILD_INTERNAL_YAXT],
     [AM_COND_IF([BUILD_INTERNAL_CDI],
        [ac_sub_configure_args="$ac_sub_configure_args FCFLAGS='$FCFLAGS ${FC_MODINC}${YAXTROOT}/inst_headers/f90'"])])])
m4_foreach([lib_build_type],[shared,static],
  [AS_CASE([x"$enable_]lib_build_type["],
    [xyes],[ac_sub_configure_args="$ac_sub_configure_args --enable-]lib_build_type",
    [xno],[ac_sub_configure_args="$ac_sub_configure_args --disable-]lib_build_type")
])
dnl

AS_IF([test "x$use_mpi" != "xno" -a x"$with_yaxt" != xno],
  [ac_sub_configure_args="$ac_sub_configure_args 'MPI_FC_LIB=${MPI_FC_LIB-$MPI_LIB}' 'MPI_C_LIB=${MPI_C_LIB-$MPI_LIB}' 'MPI_C_INCLUDE=${MPI_C_INCLUDE-${MPIROOT+-I$MPIROOT/include}}' 'MPI_FC_INCLUDE=${MPI_FC_INCLUDE-${MPIROOT+-I$MPIROOT/include}}' 'MPI_FC_MOD=${MPI_FC_MOD-${MPI_FC_INCLUDE-${MPIROOT+$FC_MODINC$MPIROOT/include $FC_MODINC$MPIROOT/lib}}}' '${MPI_LAUNCH+MPI_LAUNCH=$MPI_LAUNCH}'"])
dnl
AS_IF([test x"${enable_cdi_pio}" != xno],
  [ac_sub_configure_args="$ac_sub_configure_args --enable-mpi"])
dnl
dnl
AS_IF([test x${YAXT_CFLAGS+set} = xset -a x${YAXT_LIBS+set} = xset],
  [ac_sub_configure_args="$ac_sub_configure_args 'YAXT_CFLAGS=${YAXT_CFLAGS}' 'YAXT_LIBS=${YAXT_LIBS}'"])

AM_COND_IF([BUILD_INTERNAL_YAXT],
           [AC_CONFIG_SUBDIRS([yaxt])],[])

AS_CASE([x"$NETCDFROOT"],
  [x],[],
  [x/usr],[ac_sub_configure_args="$ac_sub_configure_args --with-netcdf"],
  [ac_sub_configure_args="$ac_sub_configure_args '--with-netcdf=$NETCDFROOT'"])
AS_CASE([x"$HDF5ROOT"],
  [x],[],
  [x/usr],[ac_sub_configure_args="$ac_sub_configure_args --with-hdf5"],
  [ac_sub_configure_args="$ac_sub_configure_args '--with-hdf5=$HDF5ROOT'"])
AS_CASE([x"$SZIPROOT"],
  [x],[],
  [x/usr],[ac_sub_configure_args="$ac_sub_configure_args --with-szlib"],
  [ac_sub_configure_args="$ac_sub_configure_args '--with-szlib=$SZIPROOT'"])
AS_CASE([x"$ZLIBROOT"],
  [x],[],
  [x/usr],[ac_sub_configure_args="$ac_sub_configure_args --with-zlib"],
  [ac_sub_configure_args="$ac_sub_configure_args '--with-zlib=$ZLIBROOT'"])

ac_sub_configure_args="$ac_sub_configure_args $extra_cdi_lib_configure_args"

AM_COND_IF([BUILD_INTERNAL_CDI],
           [AC_CONFIG_SUBDIRS([cdi])])

AC_OUTPUT

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
AS_IF([test x$MPIROOT = x],
      [AS_CASE([$host],
               [powerpc-ibm-aix*],[
                   if test "$use_mpi_compiler" = yes ; then
                       mpi_version_number=$(lslpp -l ppe.poe | grep ppe.poe | uniq | awk '{print $2}')
                       MPI_VERSION="IBM MPI ${mpi_version_number}"
                   fi],
               [sx?-nec-superux],[
                   mpi_version_number=$(strings $MPIROOT/lib/libmpi.a | grep "MPI/SX: library Version" | awk '{print $4}')
                   MPI_VERSION="SX MPI ${mpi_version_number}"],
               [MPI_VERSION=none])],
      [AS_CASE([$host],
               [*-*-darwin*],[
                   mpi_version_number=$(mpiexec --version | awk 'NR==2{print $2}')
                   MPI_VERSION="mpich ${mpi_version_number} (MacPorts)"],
               [
                   MPI_VERSION="${MPIROOT##*/}"
               ])
])
dnl _______________________________________________________________________________________________
dnl
AS_ECHO([])
AS_ECHO(["Selected compiler:"])
AS_ECHO([])

AS_ECHO(["    Fortran   : $[]ax_fc_version"])
AS_ECHO(["    C         : $[]ax_cc_version"])
dnl TODO: OpenMP support temporarily disabled due to possible race condition
dnl AS_ECHO([])
dnl AS_ECHO(["    OpenMP    : $[]use_openmp"])

AS_ECHO([])
AS_ECHO(["Selected libraries:"])
AS_ECHO([])

AS_ECHO(["    MPI       : $[]MPI_VERSION"])
AS_ECHO(["                installation  $[]MPIROOT"])

AS_IF([test x$NETCDFROOT = x],
      [AS_ECHO(["    NETCDF    : none (error)"])],
      [AS_ECHO(["    NETCDF    : netcdf-$[]NETCDF_VERSION"])
       AS_ECHO(["                installation  $[]NETCDFROOT"])])

AS_IF([test x$NETCDFFROOT = x],
      [AS_ECHO(["    NETCDFF   : none (error if NETCDF version >= 4.2)"])],
      [AS_ECHO(["    NETCDFF   : netcdff-$[]NETCDF_FVERSION"])
       AS_ECHO(["                installation  $[]NETCDFFROOT"])])

AS_IF([test x$HDF5ROOT = x],
      [AS_ECHO(["    HDF5      : none (error)"])],
      [AS_ECHO(["    HDF5      : hdf5-$[]HDF5_VERSION"])
       AS_ECHO(["                installation  $[]HDF5ROOT"])])

AS_IF([test x$SZIPROOT = x],
      [AS_ECHO(["    SZIP      : none (error)"])],
      [AS_ECHO(["    SZIP      : szip-$[]SZIP_VERSION"])
       AS_ECHO(["                installation  $[]SZIPROOT"])])

AS_IF([test x$ZLIBROOT = x],
      [AS_ECHO(["    ZLIB      : none (error)"])],
      [AS_ECHO(["    ZLIB      : zlib-$[]ZLIB_VERSION"])
       AS_ECHO(["                installation  $[]ZLIBROOT"])])

AS_ECHO([])

AS_IF([test x$CDIROOT = x],
      [AS_ECHO(["    CDI       : none"])],
      [AS_ECHO(["    CDI       : cdi-$[]CDI_VERSION"])
       AS_ECHO(["                installation  $[]CDIROOT"])])

AS_IF([test x$YAXTROOT = x],
      [AS_ECHO(["    YAXT      : none"])],
      [AS_ECHO(["    YAXT      : yaxt-$[]YAXT_VERSION"])
       AS_ECHO(["                installation  $[]YAXTROOT"])])

AS_IF([test x$SCTROOT = x],
      [AS_ECHO(["    SCT       : none"])],
      [AS_ECHO(["    SCT       : $[]SCTROOT"])])

AS_IF([test $BUILD_INTERNAL_LAPACK = yes],
      [AS_ECHO(["    LAPACK    : internal"])],
      [AS_ECHO(["    LAPACK    : $[]LAPACKROOT"])])

AS_ECHO([])

AS_IF([test x$POSTGRESQLROOT = x],
      [AS_ECHO(["    POSTGRESQL: none"])],
      [AS_ECHO(["    POSTGRESQL: postgresql-$[]POSTGRESQL_VERSION"])
       AS_ECHO(["                installation $[]POSTGRESQLROOT"])])

AS_ECHO([])

AS_CASE([$with_coupler],
        [oasis3],
        [AS_ECHO(["    OASIS3    : $[]OASIS3ROOT"])],
      	[oasis3-mct],
        [AS_ECHO(["    OASIS3-MCT: $[]OASIS3MCTROOT"])],
        [AS_ECHO(["    Coupling  : not enabled"])])

AS_IF([test $FC = "pgfortran"],
      [AS_ECHO([])
       AS_ECHO(['To prevent problems with pgfortran/libtool make sure'])
       AS_ECHO(['   @S|@HOME/.mypgfortranrc'])
       AS_ECHO(['is available and is containing'])
       AS_ECHO(['   switch -pthread is replace@{:@-lpthread@:}@ positional@{:@linker@:}@;'])     
      ],[])

AS_ECHO([])
AS_ECHO(['Configuration done. Run make (GNU make) to build ECHAM.'])
AS_ECHO([])
dnl _______________________________________________________________________________________________
