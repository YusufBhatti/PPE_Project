AC_DEFUN([AX_LIB_CDI],[
AC_ARG_WITH([cdi],
            [AS_HELP_STRING([--with-cdi=DIR],
                            [Required I/O package, use internal, if not set])],[],[])
AC_ARG_VAR([CDIROOT],[directory CDI is installed in])

AS_IF([test x"$with_cdi" != x],
  [CDIROOT=$with_cdi],
  [test x${CDIROOT+set} != xset],
  [with_cdi=internal
   CDISRC="$(cd ${srcdir} ; pwd -P)/cdi"
   CDIROOT="$(pwd -P)/cdi"])

AC_CHECK_FILE([$CDIROOT/include/cdi.inc],
  [AS_IF([test x${CDI_INCLUDE+set} = x],
     [CDI_INCLUDE='-I$(CDIROOT)/include'])
      AS_IF([test x${CDI_LIB+set} = x],
        [CDI_LIB='-L$(CDIROOT)/lib -lcdi'])],
  [AC_CHECK_FILE([$CDISRC/src/cdi.inc],
    [AS_IF([test x${CDI_INCLUDE+set} = x],
       [CDI_INCLUDE="-I${CDISRC}/src"])
     AS_IF([test x${CDI_LIB+set} = x],
       [CDI_LIB="-L$ac_pwd/cdi/src -lcdi"])],
    [AC_MSG_ERROR([CDIROOT not properly defined])])])

AC_ARG_ENABLE([cdi-pio],
  [AS_HELP_STRING([--enable-cdi-pio],
    [Make use of CDI-PIO library for I/O @<:@default=no@:>@])],
    [],
    [enable_cdi_pio=no])

AS_IF([test x"$enable_cdi_pio" = x"yes"],
  [AS_IF([test x"${YAXTROOT}" = x],
     [AC_MSG_ERROR([YAXT is required by CDI-PIO but unavailable.
Specify a valid YAXT installation using YAXTROOT or option --with-yaxt.])])
   AC_CHECK_FILE([$CDIROOT/include/cdipio.inc],
     [PKG_CONFIG_PATH="$CDIROOT/lib/pkgconfig${PKG_CONFIG_PATH+:$PKG_CONFIG_PATH}"],
     [AC_CHECK_FILE([$CDISRC/src/cdipio.inc],,
        [AC_MSG_ERROR([CDI-PIO include file 'cdipio.inc' is unavailable.
Specify a valid CDI installation using CDIROOT or option --with-cdi.])])])

   CDI_LIB=`echo $CDI_LIB | sed -e 's%-lcdi%-lcdipio -lcdi%'`
   FCDEFS="${FCDEFS+$FCDEFS }${FC_DEFINE}HAVE_CDIPIO"])

PKG_CONFIG_PATH="$CDIROOT/lib/pkgconfig${PKG_CONFIG_PATH+:$PKG_CONFIG_PATH}"
export PKG_CONFIG_PATH
AS_IF([pkg-config --exists cdi],
      [CDI_VERSION=`pkg-config --modversion cdi`],
      [CDI_VERSION='internal'])

AM_CONDITIONAL([BUILD_INTERNAL_CDI],
               [test "$with_cdi" = internal])

AS_CASE([$FC],
        [nagfor|gfortran|ifort|pgfortran],[
            CDI_LIB=$(echo $CDI_LIB | $SED -e "s;\(-L\)\(/sw/\S*\);\1\2 $WLFLAG,-rpath,\2;g")
        ],[])

AC_SUBST([CDIROOT])
AC_SUBST([CDI_LIB])
AC_SUBST([CDI_INCLUDE])

])

dnl
dnl Local Variables:
dnl mode: autoconf
dnl End:
dnl
