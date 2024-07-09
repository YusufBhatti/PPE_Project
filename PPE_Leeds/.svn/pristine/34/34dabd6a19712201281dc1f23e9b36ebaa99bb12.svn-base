# ===========================================================================
#
# SYNOPSIS
#
# AX_MACHINE_CONFIG([file-to-read])
#
# DESCRIPTION
#
# Load machine specific compiler options
#
# LICENSE
#
#   Copyright (c) 2012 Luis Kornblueh <luis.kornblueh@zmaw.de>
#
#   This program is free software; you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation; either version 2 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

AU_ALIAS([MACHINE_CONFIG], [AX_MACHINE_CONFIG])
AC_DEFUN([AX_MACHINE_CONFIG], [

AC_MSG_CHECKING([for machine dependent configuration])

AC_REQUIRE([AC_PROG_SED])

AC_ARG_ENABLE([mh-file],      
              [AC_HELP_STRING([--enable-mh-file=<name>],
                              [select specific build and host configuration @<:@default: auto@:>@])],
              [],
              [enable_mh_file=auto])

AS_IF([test x"$enable_mh_file" = xyes],
      [enable_mh_file=auto],
      [])      	    

AS_IF([test x"$enable_mh_file" = xauto],
      [host_frag=""
       AS_CASE([$host],
               [*-ibm-aix*],      [host_frag=$srcdir/config/mh-aix],
               [*86*-*-linux-*],  [host_frag=$srcdir/config/mh-linux],
               [*-apple-darwin*], [host_frag=$srcdir/config/mh-darwin],
               [AC_MSG_WARN([This configuration is not supported. Please create a valid config/mh-* file.])])],
      [host_frag="$srcdir/config/$enable_mh_file"])

AS_IF([test x"$enable_mh_file" != xno],
      [AS_IF([test x"$host_frag" != x -a ! -f $host_frag],
             [AC_MSG_FAILURE([machine dependent configuration file $host_frag does not exist!], 1)],
             [])

ax_mc_cc=${CC:-""}		       
ax_mc_cflags=${CFLAGS:-""}	       
ax_mc_fc=${FC:-""}		       
ax_mc_fcflags=${FCFLAGS:-""}	       
ax_mc_f77=${F77:-""}		       
ax_mc_fflags=${FFLAGS:-""}	       
ax_mc_mpiroot=${MPIROOT:-""}	       
ax_mc_netcdfroot=${NETCDFROOT:-""}     
ax_mc_hdf5root=${HDF5ROOT:-""}	       
ax_mc_sziproot=${SZIPROOT:-""}	       
ax_mc_zlibroot=${ZLIBROOT:-""}         

changequote(,)
cat > confsed <<EOF
/^[ ]*\#/d
/^[ ]*\$/d
s/^[[:blank:]]*/ /
s/[ ]*\$//
s/[ ]*=[ ]*/=/
/^[ ]*[A-Za-z0-9_]*=/ {
    /\$(/! {	   
	s/=/="/ 
	s/\$/"/    
    }
}
p
EOF
changequote([,])
$SED -n -f confsed $host_frag > conftest
. ./conftest
/bin/rm -f confsed conftest
if test "$host_frag" != 0 ; then
    AC_MSG_RESULT([$host_frag])
else
    AC_MSG_RESULT([unavailable])
fi

AS_IF([test -n "$ax_mc_cc"],
      [CC=$ax_mc_cc], 
      [])
AS_IF([test -n "$ax_mc_cflags"],
      [CFLAGS=$ax_mc_cflags],
      [])
AS_IF([test -n "$ax_mc_fc"], 
      [FC=$ax_mc_fc],
      [])
AS_IF([test -n "$ax_mc_fcflags"], 
      [FCFLAGS=$ax_mc_fcflags],
      [])
AS_IF([test -n "$ax_mc_f77"], 
      [F77=$ax_mc_f77],
      [])
AS_IF([test -n "$ax_mc_f77flags"], 
      [FFLAGS=$ax_mc_fflags],
      [])
AS_IF([test -n "$ax_mc_mpiroot"], 
      [MPIROOT=$ax_mc_mpiroot],
      [])
AS_IF([test -n "$ax_mc_netcdfroot"], 
      [NETCDFROOT=$ax_mc_netcdfroot],
      [])
AS_IF([test -n "$ax_mc_hdf5root"], 
      [HDF5ROOT=$ax_mc_hdf5root],
      [])
AS_IF([test -n "$ax_mc_sziproot"], 
      [SZIPROOT=$ax_mc_sziproot],
      [])
AS_IF([test -n "$ax_mc_zlibroot"], 
      [ZLIBROOT=$ax_mc_zlibroot],
      [])
],
[AS_ECHO([disabled])])

])
