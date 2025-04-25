! config/config.inc.  Generated from config.h.in by configure. 
! config/config.h.in.  Generated from configure.ac by autoheader. 

! Define if building universal (internal helper macro)
! #undef AC_APPLE_UNIVERSAL_BUILD

! Define to dummy `main' function (if any) required to link to the Fortran
! libraries.
! #undef FC_DUMMY_MAIN

! Define if F77 and FC dummy `main' functions are identical.
! #undef FC_DUMMY_MAIN_EQ_F77

! Define to a macro mangling the given C identifier (in lower and upper
! case), which must not contain underscores, for linking with Fortran.
#define FC_FUNC(name,NAME) name ## _

! As FC_FUNC, but for C identifiers containing underscores.
#define FC_FUNC_(name,NAME) name ## _

! Define to 1 if you have the <dlfcn.h> header file.
#define HAVE_DLFCN_H 1

! Define to 1 if you have the <execinfo.h> header file.
#define HAVE_EXECINFO_H 1

! Define to 1 if you have the <fcntl.h> header file.
#define HAVE_FCNTL_H 1

! Define to 1 if you have the <fortran.h> header file.
! #undef HAVE_FORTRAN_H

! Define to 1 if you have the `getrusage' function.
#define HAVE_GETRUSAGE 1

! Define to 1 if you have the `gettimeofday' function.
#define HAVE_GETTIMEOFDAY 1

! Defined if you have HDF5 support available
#define HAVE_HDF5 1

! Define to 1 if you have the <inttypes.h> header file.
#define HAVE_INTTYPES_H 1

! Defined if you have szip support
#define HAVE_LIBZ 1

! Define to 1 if you have the <limits.h> header file.
#define HAVE_LIMITS_H 1

! Define to 1 if you have the <malloc.h> header file.
#define HAVE_MALLOC_H 1

! Defined if you have netCDF support available
#define HAVE_NETCDF 1

! Define to 1 if you have the <netdb.h> header file.
#define HAVE_NETDB_H 1

! Define to 1 if PostgreSQL libraries are available
! #undef HAVE_POSTGRESQL

! Define to 1 if you have the <pwd.h> header file.
#define HAVE_PWD_H 1

! Define to 1 if you have the <stdint.h> header file.
#define HAVE_STDINT_H 1

! Define to 1 if you have the <stdio.h> header file.
#define HAVE_STDIO_H 1

! Define to 1 if you have the <stdlib.h> header file.
#define HAVE_STDLIB_H 1

! Define to 1 if you have the <strings.h> header file.
#define HAVE_STRINGS_H 1

! Define to 1 if you have the <string.h> header file.
#define HAVE_STRING_H 1

! Define to 1 if you have the `sysconf' function.
#define HAVE_SYSCONF 1

! Define to 1 if you have the <sys/param.h> header file.
#define HAVE_SYS_PARAM_H 1

! Define to 1 if you have the <sys/stat.h> header file.
#define HAVE_SYS_STAT_H 1

! Define to 1 if you have the <sys/time.h> header file.
#define HAVE_SYS_TIME_H 1

! Define to 1 if you have the <sys/types.h> header file.
#define HAVE_SYS_TYPES_H 1

! Define to 1 if you have the <sys/unistd.h> header file.
#define HAVE_SYS_UNISTD_H 1

! Define to 1 if you have the <sys/utsname.h> header file.
#define HAVE_SYS_UTSNAME_H 1

! Define to 1 if you have the <ucontext.h> header file.
#define HAVE_UCONTEXT_H 1

! Define to 1 if you have the `uname' function.
#define HAVE_UNAME 1

! Define to 1 if you have the <unistd.h> header file.
#define HAVE_UNISTD_H 1

! Define to 1 if you have the `valloc' function.
#define HAVE_VALLOC 1

! Define to the sub-directory where libtool stores uninstalled libraries.
#define LT_OBJDIR ".libs/"

! Name of package
#define PACKAGE "echam-6.3.00"

! Define to the address where bug reports for this package should be sent.
#define PACKAGE_BUGREPORT "luis.kornblueh@mpimet.mpg.de"

! Define to the full name of this package.
#define PACKAGE_NAME "echam"

! Define to the full name and version of this package.
#define PACKAGE_STRING "echam 6.3.00"

! Define to the one symbol short name of this package.
#define PACKAGE_TARNAME "echam-6.3.00"

! Define to the home page for this package.
#define PACKAGE_URL ""

! Define to the version of this package.
#define PACKAGE_VERSION "6.3.00"

! Number of lines to get page-aligned quickplot tables
#define QPLOT_TABLE_LINES 53

! The size of `char', as computed by sizeof.
#define SIZEOF_CHAR 1

! The size of `double', as computed by sizeof.
#define SIZEOF_DOUBLE 8

! The size of `float', as computed by sizeof.
#define SIZEOF_FLOAT 4

! The size of `int', as computed by sizeof.
#define SIZEOF_INT 4

! The size of `int *', as computed by sizeof.
#define SIZEOF_INT_P 8

! The size of `long', as computed by sizeof.
#define SIZEOF_LONG 8

! The size of `long double', as computed by sizeof.
#define SIZEOF_LONG_DOUBLE 16

! The size of `long long', as computed by sizeof.
#define SIZEOF_LONG_LONG 8

! The size of `short', as computed by sizeof.
#define SIZEOF_SHORT 2

! Define to 1 if all of the C90 standard headers exist (not just the ones
! required in a freestanding environment). This macro is provided for
! backward compatibility; new code need not use it.
#define STDC_HEADERS 1

! Define to 1 if you can safely include both <sys/time.h> and <time.h>. This
! macro is obsolete.
#define TIME_WITH_SYS_TIME 1

! Version number of package
#define VERSION "6.3.00"

! Define WORDS_BIGENDIAN to 1 if your processor stores words with the most
! significant byte first (like Motorola and SPARC, unlike Intel).
#if defined AC_APPLE_UNIVERSAL_BUILD
# if defined __BIG_ENDIAN__
#  define WORDS_BIGENDIAN 1
# endif
#else
# ifndef WORDS_BIGENDIAN
! #  undef WORDS_BIGENDIAN
# endif
#endif

! Define to empty if `const' does not conform to ANSI C.
! #undef const
