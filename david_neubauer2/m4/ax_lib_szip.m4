AU_ALIAS([CHECK_SZIP], [AX_LIB_SZIP])
AC_DEFUN([AX_LIB_SZIP],[

AC_REQUIRE([AC_PROG_GREP])
AC_REQUIRE([AC_PROG_AWK])

ax_szip_search_paths="${SZIPROOT:-""} /usr /usr/local /opt/local /sw /client"

AC_ARG_WITH([szip],
            [AS_HELP_STRING([--with-szip=[PATH]],[location of szip installation directory])]
            [],
            [withval=no])

AS_IF([test "x$withval" != "xno"],
      [AS_IF([test -d "$withval"],
             [ax_szip_search_paths="$withval $ax_szip_search_paths"],
             [])],
      [])

AS_IF([test -n "$ax_szip_search_paths"],
      [for ax_sziproot in $ax_szip_search_paths
       do
          if test -f "$ax_sziproot/include/szlib.h"
          then 
              break
          fi
          ax_sziproot=""
       done
       
       CPPFLAGSsave=$CPPFLAGS
       CPPFLAGS="$CPPFLAGS -I$ax_sziproot/include"

       LIBSsave=$LIBS
       LIBS="$LIBS -L$ax_sziproot/lib -lsz"

       AC_CHECK_HEADER([szlib.h],[ax_szip_header_ok=yes],[ax_szip_header_ok=no],[])
       AC_CHECK_LIB([sz],[SZ_encoder_enabled],[ax_szip_lib_ok=yes],[ax_szip_lib_ok=no],[])

       CPPFLAGS=$CPPFLAGSsave

       LIBS=$LIBSsave

       AS_IF([test "$ax_szip_header_ok" = "yes" && test "$ax_szip_lib_ok" = "yes"],
             [SZIPROOT="$ax_sziproot"
              SZIP_CPPFLAGS="-I$ax_sziproot/include"
              SZIP_INCLUDE="-I$ax_sziproot/include"
              SZIP_LIB="-L$ax_sziproot/lib -lsz"
	      SZIP_VERSION=$(grep "define SZLIB_VERSION" ${SZIPROOT}/include/szlib.h | $AWK '{gsub(/@<:@:"@:>@*/,""); print $[]3}')
              AC_DEFINE([HAVE_LIBZ],[1],[Defined if you have szip support])],
             [SZIPROOT=""
              SZIP_CPPFLAGS=""
              SZIP_INCLUDE=""
              SZIP_LIB=""
              SZIP_VERSION=""])

AS_CASE([$FC],
        [nagfor|gfortran|ifort|pgfortran],[
            SZIP_LIB=$(echo $SZIP_LIB | $SED -e "s;\(-L\)\(/sw/\S*\);\1\2 $WLFLAG,-rpath,\2;g")
        ],[])  

       AC_SUBST([SZIPROOT])
       AC_SUBST([SZIP_INCLUDE])
       AC_SUBST([SZIP_LIB])
       AC_SUBST([SZIP_VERSION])
      ],[
       AC_MSG_ERROR([specify a valid szip installation with --with-szip=PATH])
      ])

])
