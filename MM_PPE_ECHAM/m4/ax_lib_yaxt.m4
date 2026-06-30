AC_DEFUN([AX_LIB_YAXT], [

AC_ARG_WITH([yaxt],
            [AS_HELP_STRING([--with-yaxt=DIR],
            [Required communication library, use internal, if not set])],[],[])

AC_ARG_VAR([YAXTROOT],[directory YAXT is installed in])

AS_IF([test x"$with_yaxt" != x],
  [YAXTROOT=$with_yaxt],
  [test x${YAXTROOT+set} != xset],
  [with_yaxt=internal
   YAXTSRC="$(cd ${srcdir} ; pwd -P)/yaxt"
   YAXTROOT="$(pwd -P)/yaxt"])

#AS_ECHO([])
#AS_ECHO([with_yaxt: $with_yaxt])
#AS_ECHO([])
#AS_ECHO([YAXTROOT: $YAXTROOT])
#AS_ECHO([])

AS_IF([test "x$with_yaxt" = "xinternal"],
      [AC_CHECK_FILE([$YAXTSRC/src/yaxt.f90],
                     [YAXT_LIB="-L${YAXTROOT}/src -lyaxt"
                      YAXT_INCLUDE="${FC_MODINC}${YAXTROOT}/inst_headers/f90"
                      YAXT_CFLAGS="-I${YAXTROOT}/inst_headers"
                      YAXT_LIBS="${YAXT_LIB}"
                      YAXT_VERSION="internal"
		      FCDEFS="${FCDEFS+$FCDEFS }${FC_DEFINE}HAVE_YAXT"],
                     [AC_MSG_ERROR([internal YAXT sources not found])])],
      [AC_CHECK_FILE([$YAXTROOT/include/yaxt.$FC_MODEXT],
                     [YAXT_LIB="-L${YAXTROOT}/lib -lyaxt"
                      YAXT_INCLUDE="${FC_MODINC}${YAXTROOT}/include"
                      YAXT_CFLAGS="-I${YAXTROOT}/include"
                      YAXT_LIBS="-L${YAXTROOT}/lib -lyaxt"
                      FCDEFS="${FCDEFS+$FCDEFS }${FC_DEFINE}HAVE_YAXT"],
                     [AC_MSG_ERROR([YAXTROOT or --with-yaxt value not properly defined])])])

#AS_ECHO([])
#AS_ECHO([YAXT_LIB: $YAXT_LIB])
#AS_ECHO([YAXT_CFLAGS: $YAXT_CFLAGS])
#AS_ECHO([YAXT_INCLUDE: $YAXT_INCLUDE])
#AS_ECHO([FCDEFS: $FCDEFS])
#AS_ECHO([])

# AC_ARG_WITH([yaxt],
#             [AS_HELP_STRING([--with-yaxt=DIR|yes|no],
#             [Set directory to search for YAXT headers and library (DIR) or use version
#              of YAXT distributed with ECHAM (yes). Make use of YAXT library for
#              transpositions or I/O when set to yes or DIR. @<:@default=no@:>@])],
#             [],
#             [AS_IF([test x${YAXTROOT+set} = xset],
#                    [with_yaxt=yes],
#                    [with_yaxt=no])])

# AC_ARG_VAR([YAXTROOT],[directory YAXT is installed in])

# AS_IF([test x"$with_yaxt" = xno],
#       [YAXTROOT=""
#        YAXT_VERSION=""
#        YAXT_LIB=""
#        YAXT_INCLUDE=""],
#       [AS_IF([test x"$with_yaxt" = xyes -o x"$with_yaxt" = x],
#              [AS_IF([test x${YAXTROOT+set} != xset],
#                     [YAXTROOT=`pwd -P`/yaxt
#                      SRCDIRS="yaxt ${SRCDIRS}"
#                      YAXT_LIB='-L$(YAXTROOT)/src -lyaxt'
#                      YAXT_INCLUDE='$(FC_MODINC)$(YAXTROOT)/inst_headers/f90'
#                      YAXT_CFLAGS="-I$YAXTROOT/inst_headers"
#                      YAXT_LIBS="-L${YAXTROOT}/src -lyaxt"
#                      YAXT_VERSION=`sed -n '/^A@<:@C@:>@_INIT(/{
# s/A@<:@C@:>@_INIT(\@<:@@<:@^@:>@@:>@*\@:>@,//
# s/A@<:@C@:>@_INIT(@<:@^,@:>@*,//
# s/\@<:@\(@<:@^@:>@@:>@*\)\@:>@.*/\1/
# t found
# s/\(@<:@^,@:>@\).*/\1/
# : found
# p
# }' "$srcdir/yaxt/configure.ac"`
#                      AS_IF([test -d "$YAXTROOT" -o -d "$srcdir/yaxt" ],
#                            [AC_MSG_NOTICE([Using internal YAXT library])],
#                            [AC_MSG_ERROR([YAXT directory not found in echam6 distribution.
# Specify a valid YAXT installation using YAXTROOT or option --with-yaxt.])])])],
#              [YAXTROOT="$with_yaxt"
#               AC_CHECK_FILE([$YAXTROOT/include/yaxt.$FC_MODEXT],,
#                             [AC_MSG_ERROR([YAXTROOT or --with-yaxt value not properly defined])])])

#        AS_IF([test x"$YAXTROOT" != x],
#              [FCDEFS="${FCDEFS+$FCDEFS }${FC_DEFINE}HAVE_YAXT"
#               AS_IF([test x${YAXT_LIB+set} != xset],
#                     [YAXT_LIB='-L$(YAXTROOT)/lib -lyaxt'
#                      AS_IF([test "x$ax_rpath" != "xno"],
#                            [AS_IF([test "x$ax_rpath_split" = "xyes"],
#                                   [YAXT_LIB="${YAXT_LIB} $ax_rpath_prefix$ax_rpath$ax_rpath_prefix"'$(YAXTROOT)/lib'],
#                                   [YAXT_LIB="${YAXT_LIB} $ax_rpath_prefix$ax_rpath"'$(YAXTROOT)/lib'])])])
#               AS_IF([test x${YAXT_INCLUDE+set} != xset],
#                     [YAXT_INCLUDE='$(FC_MODINC)$(YAXTROOT)/include'],
#                     [])])

#        PKG_CONFIG_PATH="$YAXTROOT/lib/pkgconfig${PKG_CONFIG_PATH+:$PKG_CONFIG_PATH}"
#        export PKG_CONFIG_PATH

#        AS_IF([pkg-config --exists yaxt],
#              [YAXT_VERSION=`pkg-config --modversion yaxt`],
#              [])
# ])

AM_CONDITIONAL([BUILD_INTERNAL_YAXT],
               [test "$with_yaxt" = internal])

#AS_ECHO([BUILD_INTERNAL_YAXT: $BUILD_INTERNAL_YAXT]) 

# AM_CONDITIONAL([BUILD_INTERNAL_YAXT],[test x"${YAXTROOT}" = x"`pwd -P`/yaxt"])

AC_SUBST([YAXTROOT])
AC_SUBST([YAXT_LIB])
AC_SUBST([YAXT_INCLUDE])
AC_SUBST([YAXT_VERSION])

])
