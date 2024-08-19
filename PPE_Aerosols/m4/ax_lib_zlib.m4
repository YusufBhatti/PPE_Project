AU_ALIAS([CHECK_ZLIB], [AX_LIB_ZLIB])
AC_DEFUN([AX_LIB_ZLIB],[

ax_zlib_search_paths="${ZLIBROOT:-""} /usr /usr/local /opt/local /sw /client"

AC_ARG_WITH([zlib],
            [AS_HELP_STRING([--with-zlib=[PATH]],[location of zlib installation directory])]
            [],
            [withval=no])

AS_IF([test "x$withval" != "xno"],
      [AS_IF([test -d "$withval"],
             [ax_zlib_search_paths="$withval $ax_zlib_search_paths"],
             [])],
      [])

AS_IF([test -n "$ax_zlib_search_paths"],
      [for ax_zlibroot in $ax_zlib_search_paths
       do
          if test -f "$ax_zlibroot/include/zlib.h"
          then 
              break
          fi
          ax_zlibroot=""
       done

       LIBSsave=$LIBS
       CPPFLAGSsave=$CPPFLAGS

       AC_CHECK_HEADER([zlib.h],[ax_zlib_header_ok=yes],[ax_zlib_header_ok=no],[])
       AC_CHECK_LIB([z],[inflateEnd],[ax_zlib_lib_ok=yes],[ax_zlib_lib_ok=no],[])

       AS_IF([test "$ax_zlib_header_ok" = "yes" && test "$ax_zlib_lib_ok" = "yes"],
             [ZLIBROOT="$ax_zlibroot"
              ZLIB_INCLUDE="-I$ax_zlibroot/include"
              ZLIB_LIB="-L$ax_zlibroot/lib -lz"
	      ZLIB_VERSION=$(grep "define ZLIB_VERSION" ${ZLIBROOT}/include/zlib.h | $AWK '{gsub(/@<:@:"@:>@*/,""); print $[]3}')
              AC_DEFINE([HAVE_LIBZ],[1],[Defined if you have zlib support])],
             [ZLIBROOT=""
              ZLIB_INCLUDE=""
              ZLIB_LIB=""
	      ZLIB_VERSION=""])

AS_CASE([$FC],
        [nagfor|gfortran|ifort|pgfortran],[
            ZLIB_LIB=$(echo $ZLIB_LIB | $SED -e "s;\(-L\)\(/sw/\S*\);\1\2 $WLFLAG,-rpath,\2;g")
        ],[])  

       AC_SUBST([ZLIBROOT])
       AC_SUBST([ZLIB_INCLUDE])
       AC_SUBST([ZLIB_LIB])
       AC_SUBST([ZLIB_VERSION])
      ],[
       AC_MSG_ERROR([specify a valid zlib installation with --with-zlib=PATH])
      ])

])
