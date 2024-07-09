AU_ALIAS([CHECK_NETCDF4_LIB], [AX_LIB_NETCDF4])
AC_DEFUN([AX_LIB_NETCDF4], [

AC_REQUIRE([AC_PROG_SED])
AC_REQUIRE([AC_PROG_AWK])

CFLAGSsave=$CFLAGS
CPPFLAGSsave=$CPPFLAGS
LIBSsave=$LIBS

ax_netcdf_environment_defined=false
AS_IF([test -n "$NETCDF_LIB"],
      [ax_netcdf_environment_defined=true],
      [])

AS_IF([test -z "$NETCDF_INCLUDE" -a $ax_netcdf_environment_defined],
      [ax_netcdf_environment_defined=false],
      [])

AS_IF([test -z "$NETCDF_FLIB" -a $ax_netcdf_environment_defined],
      [ax_netcdf_environment_defined=false],
      [])

AS_IF([test -z "$NETCDF_FINCLUDE" -a $ax_netcdf_environment_defined],
      [ax_netcdf_environment_defined=false],
      [])

AS_IF([test "$ax_netcdf_environment_defined" == "false"],
      [

NETCDF_VERSION=""
NETCDF_INCLUDE=""
NETCDF_LIB=""

NETCDF_FVERSION=""
NETCDF_FINCLUDE=""
NETCDF_FLIB=""

NETCDFROOT=${NETCDFROOT:-""}
NETCDFFROOT=${NETCDFFROOT:-""}

AC_ARG_WITH([netcdf],
            AS_HELP_STRING([--with-netcdf=[PATH]], [location of netCDF installation directory]))
AS_IF([test -d "$with_netcdf"],
      [NETCDFROOT="$with_netcdf"],
      [])    	    
PATH=$NETCDFROOT/bin:$PATH 
AC_PATH_PROGS([NC_CONFIG],
              [nc-config], 
	      [])

AS_IF([test -z "$NC_CONFIG"],
      [ax_netcdf_version=3],
      [ax_netcdf_version=4])

AC_MSG_CHECKING([for netCDF library])

AS_IF([test "$ax_netcdf_version" -eq 3],
       [
        NETCDF_VERSION=$(${NETCDFROOT}/bin/ncdump 2>&1 | $AWK -F"\"" 'END{print $[]2}')
        AC_MSG_RESULT([yes (version $[NETCDF_VERSION])])

	# C interface
        AC_CHECK_FILE([$NETCDFROOT/include/netcdf.h],
                      [],
                      [AC_MSG_ERROR([NETCDFROOT not properly defined])])
        NETCDF_INCLUDE="-I$NETCDFROOT/include"
        NETCDF_LIB="-L$NETCDFROOT/lib -lnetcdf"
        CFLAGS=$NETCDF_INCLUDE
        LIBS=$NETCDF_LIB
        CPPFLAGS=$CFLAGS
        AC_CHECK_HEADER([netcdf.h],[ax_netcdf_header_ok=yes],[ax_netcdf_header_ok=no],[])
        AC_CHECK_LIB([netcdf],[nc_inq_libvers],[ax_netcdf_lib_ok=yes],[ax_netcdf_lib_ok=no],[])

        # Fortran interface
        NETCDF_FVERSION=$NETCDF_VERSION
        AC_CHECK_FILE([$NETCDFROOT/include/netcdf.inc],
                      [],
                      [AC_MSG_ERROR([NETCDFROOT not properly defined, Fortran not available])])
        AC_LANG_PUSH([Fortran])
        AC_CHECK_LIB([netcdf],[nf_inq_libvers],[ax_netcdff_lib_ok=yes],[ax_netcdff_lib_ok=no],[])
        AC_LANG_POP([Fortran])
        AS_IF([test -r $NETCDFROOT/lib/libnetcdff.a],
              [NETCDF_FLIB="-L/$NETCDFROOT/lib -lnetcdff -lnetcdf"],
              [NETCDF_FLIB="-L/$NETCDFROOT/lib -lnetcdf"]) 	
       ],[
        NETCDF_VERSION=$(eval $NC_CONFIG --version | $AWK '{print $[]2}')
	AC_MSG_RESULT([yes (version $[NETCDF_VERSION])])

	# C interface
        ax_nc_cflags=$(eval $NC_CONFIG --cflags | $AWK '{for(i=1;i<=NF;i++) a[[$i]]++} END{for(i in a) printf i" "}')
        NETCDF_INCLUDE=$(echo $ax_nc_cflags | $AWK '{for(i=0;++i<=NF;){if ($[]i ~ /^-I/){printf "%s ",$[]i}}}')
        NETCDF_LIB=$(eval $NC_CONFIG --libs)
        CFLAGS=$NETCDF_INCLUDE
        LIBS=$NETCDF_LIB
        CPPFLAGS=$CFLAGS
        AC_CHECK_HEADER([netcdf.h],[ax_netcdf_header_ok=yes],[ax_netcdf_header_ok=no],[])
        AC_CHECK_LIB([netcdf],[nc_inq_libvers],[ax_netcdf_lib_ok=yes],[ax_netcdf_lib_ok=no],[])
      
        # Fortran interface 
        ax_netcdf_version=$(eval $NC_CONFIG --version | $AWK '{print $[]2}')
        ax_netcdf_major=$(echo $ax_netcdf_version | $AWK -F. '{printf "%d",$[]1}')
        ax_netcdf_minor=$(echo $ax_netcdf_version | $AWK -F. '{printf "%d",$[]2}')

        AS_IF([test $ax_netcdf_minor -lt 2],
              [
    	       # pre netcdf 4.2
               NETCDF_FVERSION=$NETCDF_VERSION
               ax_nc_fcflags=$(eval $NC_CONFIG --fflags | $AWK '{for(i=1;i<=NF;i++) a[[$i]]++} END{for(i in a) printf i" "}')
               NETCDF_FINCLUDE=$(echo $ax_nc_fcflags | $AWK '{for(i=0;++i<=NF;){if ($[]i ~ /^-I|^-L/){printf "%s ",$[]i}}}')
               NETCDF_FLIB=$(eval $NC_CONFIG --flibs)
               LIBS=$NETCDF_FLIB
               AC_LANG_PUSH([Fortran])
               AC_CHECK_LIB([netcdff],[nf_inq_libvers],[ax_netcdff_lib_ok=yes],[ax_netcdff_lib_ok=no],[])
               AC_LANG_POP([Fortran])
              ],[
               # netcdf 4.2
               AC_ARG_WITH([netcdff],
                           AS_HELP_STRING([--with-netcdff=[PATH]], [location of netCDF Fortran installation directory]))
               AS_IF([test -d "$with_netcdff"],
                     [NETCDFFROOT="$with_netcdff"],
                     [AS_IF([test -z $NETCDFFROOT],
                            [NETCDFFROOT="$NETCDFROOT"],
                            [])])
               PATH=$NETCDFFROOT/bin:$PATH 
               AC_PATH_PROGS([NF_CONFIG],
                             [nf-config], 
	                     [])
               AS_IF([test -z $NF_CONFIG],
                     [AC_MSG_ERROR([nf-config not available])],
                     [])
               NETCDF_FVERSION=$(eval $NF_CONFIG --version | $AWK '{print $[]2}')
               ax_nc_fcflags=$(eval $NF_CONFIG --fflags | $AWK '{for(i=1;i<=NF;i++) a[[$i]]++} END{for(i in a) printf i" "}')
               NETCDF_FINCLUDE=$(echo $ax_nc_fcflags | $AWK '{for(i=0;++i<=NF;){if ($[]i ~ /^-I|^-L/){printf "%s ",$[]i}}}')

               NETCDF_FLIB=$(eval $NF_CONFIG --flibs | $AWK '$[]2~/lnetcdff/{print $[]1,$[]2}')
               LIBS=$NETCDF_FLIB
               AC_LANG_PUSH([Fortran])
               AC_CHECK_LIB([netcdff],[nf_inq_libvers],[ax_netcdff_lib_ok=yes],[ax_netcdff_lib_ok=no],[])
               AC_LANG_POP([Fortran])
              ])
       ])
       ],
       [
           NETCDF_VERSION=${NETCDF_VERSION:-"unknown"}
	   NETCDFROOT=${NETCDFROOT:-""}
           NETCDF_FVERSION=${NETCDF_FVERSION:-"unknown"}
	   NETCDFFROOT=${NETCDFFROOT:-""}
       ])

AS_CASE([$FC],
        [nagfor|gfortran|ifort|pgfortran],[
            NETCDF_LIB=$(echo $NETCDF_LIB | $SED -e "s;\(-L\)\(/sw/\S*\);\1\2 $WLFLAG,-rpath,\2;g")
            NETCDF_FLIB=$(echo $NETCDF_FLIB | $SED -e "s;\(-L\)\(/sw/\S*\);\1\2 $WLFLAG,-rpath,\2;g")
        ],[])  

AC_SUBST([NETCDFROOT])
AC_SUBST([NETCDF_VERSION])
AC_SUBST([NETCDF_INCLUDE])
AC_SUBST([NETCDF_LIB])
AC_SUBST([NETCDFFROOT])
AC_SUBST([NETCDF_FVERSION])
AC_SUBST([NETCDF_FINCLUDE])
AC_SUBST([NETCDF_FLIB])
AC_DEFINE([HAVE_NETCDF], [1], [Defined if you have netCDF support available])

CFLAGS=$CFLAGSsave
CPPFLAGS=$CPPFLAGSsave
LIBS=$LIBSsave
])
