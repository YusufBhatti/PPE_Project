AC_DEFUN([AX_LIB_POSTGRESQL],
[
    AC_REQUIRE([AC_PROG_AWK])	
    AC_REQUIRE([AX_RPATH_HANDLING])

    AC_ARG_WITH([postgresql],
        AS_HELP_STRING([--with-postgresql=@<:@ARG@:>@],
            [use PostgreSQL library @<:@default=no@:>@, optionally specify path to pg_config]
        ),
        [AS_IF([test "$withval" = "no"],[want_postgresql="no"],
               [test "$withval" = "yes"], [want_postgresql="yes"],
               [want_postgresql="yes"
                PG_CONFIG="$withval"])],
        [want_postgresql="no"]
    )

    POSTGRESQL_INCLUDE=""
    POSTGRESQL_LIB=""
    POSTGRESQL_VERSION=""

    AS_IF([test "$want_postgresql" = "yes"],[

        AS_IF([test -z "$PG_CONFIG" -o test],
              [AC_PATH_PROG([PG_CONFIG], [pg_config], [])],
              [])

        AS_IF([test ! -x "$PG_CONFIG"],
              [AC_MSG_ERROR([$PG_CONFIG does not exist or it is not an exectuable file])
               PG_CONFIG="no"
               found_postgresql="no"],
              [])

        AS_IF([test "$PG_CONFIG" != "no"],
	      [AC_MSG_CHECKING([for PostgreSQL libraries])

               POSTGRESQL_INCLUDE="-I$($PG_CONFIG --includedir)"
	       AS_IF([test "x$ax_rpath" != "xno"],
	             [ax_postgresql_lib_path=$($PG_CONFIG --libdir)
                      AS_IF([test "x$ax_rpath_split" = "xyes"],
                            [ax_postgresql_lib="$ax_rpath_prefix$ax_rpath$ax_rpath_prefix$ax_postgresql_lib_path -lpq"],
                            [ax_postgresql_lib="$ax_rpath_prefix$ax_rpath$ax_postgresql_lib_path -lpq"])],  
                     [ax_postgresql_lib=""])	       
	       POSTGRESQL_LIB="-L$($PG_CONFIG --libdir) -lpq $ax_postgresql_lib"
               POSTGRESQL_VERSION=$($PG_CONFIG --version | $AWK '{print $2}')

               AC_DEFINE([HAVE_POSTGRESQL], 
	                 [1], 
                         [Define to 1 if PostgreSQL libraries are available])

               found_postgresql="yes"
               AC_MSG_RESULT([yes])],
              [found_postgresql="no"
               AC_MSG_RESULT([no])])],
    [])

    AC_SUBST([POSTGRESQL_VERSION])
    AC_SUBST([POSTGRESQL_INCLUDE])
    AC_SUBST([POSTGRESQL_LIB])
])
