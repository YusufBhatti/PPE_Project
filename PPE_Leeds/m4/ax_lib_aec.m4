AU_ALIAS([CHECK_AEC], [AX_LIB_AEC])
AC_DEFUN([AX_LIB_AEC],[

ax_aec_search_paths="${AECROOT:-""} /usr /usr/local /opt/local /sw /client"

AC_ARG_WITH([aec],
            [AS_HELP_STRING([--with-aec=[PATH]],[location of aec installation directory])]
            [],
            [withval=no])

AS_IF([test "x$withval" != "xno"],
      [AS_IF([test -d "$withval"],
             [ax_aec_search_paths="$withval $ax_aec_search_paths"],
             [])],
      [])

AS_IF([test -n "$ax_aec_search_paths"],
      [for ax_aecroot in $ax_aec_search_paths
       do
          if test -f "$ax_aecroot/include/libaec.h"
          then 
              break
          fi
          ax_aecroot=""
       done

       LIBSsave=$LIBS
       CPPFLAGSsave=$CPPFLAGS

       AS_IF([test -n "$ax_aecroot"],
             [LIBS="-L$ax_aecroot/lib -laec"
              CPPFLAGS="-I$ax_aecroot/include"],
             [])

       AC_CHECK_HEADER([libaec.h],[ax_aec_header_ok=yes],[ax_aec_header_ok=no],[])
       AC_CHECK_LIB([aec],[ae_decode],[ax_aec_lib_ok=yes],[ax_aec_lib_ok=no],[])

       AS_IF([test "$ax_aec_header_ok" = "yes" && test "$ax_aec_lib_ok" = "yes"],
             [AECROOT="$ax_aecroot"
              AEC_CPPFLAGS="-I$ax_aecroot/include"
              AEC_LIB="-L$ax_aecroot/lib -laec"
              AC_DEFINE([HAVE_LIBZ],[1],[Defined if you have aec support])],
             [AECROOT=""
              AEC_CPPFLAGS=""
              AEC_LIB=""])
       AC_SUBST([AECROOT])
       AC_SUBST([AEC_INCLUDE])
       AC_SUBST([AEC_LIBS])
      ],[
       AC_MSG_ERROR([specify a valid aec installation with --with-aec=PATH])
      ])

])
