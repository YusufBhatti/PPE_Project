dnl acx_mpi_defects.m4 --- check whether MPI has one or more of
dnl                        several known defects
dnl
dnl Copyright  (C)  2014  Thomas Jahns <jahns@dkrz.de>
dnl
dnl Keywords: configure configure.ac autoconf MPI mpirun mpiexec
dnl Author: Thomas Jahns <jahns@dkrz.de>
dnl Maintainer: Thomas Jahns <jahns@dkrz.de>
dnl URL: https://www.dkrz.de/redmine/projects/show/scales-ppm
dnl
dnl Redistribution and use in source and binary forms, with or without
dnl modification, are  permitted provided that the following conditions are
dnl met:
dnl
dnl Redistributions of source code must retain the above copyright notice,
dnl this list of conditions and the following disclaimer.
dnl
dnl Redistributions in binary form must reproduce the above copyright
dnl notice, this list of conditions and the following disclaimer in the
dnl documentation and/or other materials provided with the distribution.
dnl
dnl Neither the name of the DKRZ GmbH nor the names of its contributors
dnl may be used to endorse or promote products derived from this software
dnl without specific prior written permission.
dnl
dnl THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
dnl IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
dnl TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
dnl PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
dnl OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
dnl EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
dnl PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
dnl PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
dnl LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
dnl NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
dnl SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
dnl
dnl
dnl ACX_MPI_DEFECTS([TEST-SOURCE-DIR=config/checksrc],
dnl                 [ACTION-IF-CHECK-SUCCEEDS],
dnl                 [ACTION-IF-CHECK-FAILED=AC_MSG_FAILURE])
dnl
dnl Requires MPI_LAUNCH program. Also CC/FC must be setup to
dnl build MPI programs.
dnl Builds and runs simple programs from TEST-SOURCE-DIR, each of which
dnl should represent a test for a known defect that affects the library
dnl code.
dnl Each test is built according to it's file suffix as either Fortran
dnl or C MPI program.
dnl Within ACTION-IF-CHECK-SUCCEEDS and ACTION-IF-CHECK-FAILED,
dnl the following variables are set to test-specific values:
dnl acx_subtestname = base file name of the test
dnl acx_mpi_check_src = path to check source file
dnl acx_suffix = file suffix of source file
dnl
dnl Each test source may contain zero or more of the following stanzas
dnl acx_mpirun_num_tasks = N
dnl   specify number of tasks N (positive integer) to run this test with
dnl TODO: extend for F77 and C++
AC_DEFUN([ACX_MPI_DEFECTS],
  [AS_IF([test x"$MPI_LAUNCH" = xtrue],
     [AC_MSG_NOTICE([Skipping tests for known MPI defects: MPI launcher unavailable])],
     [AC_MSG_CHECKING([MPI for known defects])
      AC_MSG_RESULT([])
      for acx_mpi_check_src in "$srcdir/m4_ifval([$1],[$1],[config/checksrc])/"* ; do
        acx_suffix=`echo "$acx_mpi_check_src" | sed 's/^.*\.\(@<:@^.@:>@*\)$/\1/'`
        acx_subtestname=`echo "$acx_mpi_check_src" | sed 's/^.*\/\(@<:@^\/@:>@*\)\.@<:@^.@:>@*/\1/'`
        AS_CASE([$acx_suffix],
          [c],
          [cat confdefs.h "$acx_mpi_check_src" >conftest."$acx_suffix"
           AC_LANG_PUSH([C])],
          [f90|F90],[cat "$acx_mpi_check_src" >conftest."$acx_suffix"
           AC_LANG_PUSH([Fortran])],
          [AC_MSG_FAILURE([Unexpected language in MPI check: ${acx_subtestname}.${acx_suffix}])])
        AC_MSG_CHECKING([$acx_subtestname])
        acx_mpirun_num_tasks=`sed -n '/acx_mpirun_num_tasks *= *\(@<:@0-9@:>@*\)/{
s/.*acx_mpirun_num_tasks *= *\(@<:@0-9@:>@*\).*/\1/
p
q
}
' "$acx_mpi_check_src"`
        AS_IF([test `expr "$acx_mpirun_num_tasks" : "@<:@0-9@:>@@<:@0-9@:>@*$"` -gt 0 \
               && test "$acx_mpirun_num_tasks" -gt 0],,
          [acx_mpirun_num_tasks=1])
        AC_LINK_IFELSE(,
          [acx_mpirun_num_tasks="$MPI_LAUNCH -n $acx_mpirun_num_tasks ./conftest$EXEEXT"
           AC_TRY_EVAL([acx_mpirun_num_tasks])
           AS_IF([test $ac_status -eq 0],
             [AC_MSG_RESULT([okay])m4_ifval([$2],[
              $2])],
             [AC_MSG_RESULT([error])
              m4_ifval([$3],[$3],
                [AC_MSG_FAILURE([chosen MPI has known error $acx_subtestname])])])],
          [AC_MSG_RESULT([error])
           m4_ifval([$3],[$3],
             [AC_MSG_RESULT([chosen MPI has known error $acx_subtestname])])])
        AS_CASE([$acx_suffix],
          [f90|F90],[AC_LANG_POP([Fortran])],
          [c],[AC_LANG_POP([C])])
      done
      ASX_VAR_UNSET([acx_mpirun_num_tasks])
      ASX_VAR_UNSET([acx_mpi_check_src])
      ASX_VAR_UNSET([acx_suffix])
      ASX_VAR_UNSET([acx_subtestname])
     ])])
dnl dump text documentation of defect test to stderr
dnl ACX_MPI_DEFECTS_DOCUMENT([TEST-DOC-DIR=config/checkdoc])
AC_DEFUN([ACX_MPI_DEFECTS_DOCUMENT],
  [AS_IF([test -r "$srcdir/m4_ifval([$1],[$1],[config/checkdoc])/${acx_subtestname}.txt"],
             [cat "$srcdir/m4_ifval([$1],[$1],[config/checkdoc])/${acx_subtestname}.txt" >&2])])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://www.dkrz.de/redmine/projects/show/scales-ppm"
dnl license-default: "bsd"
dnl End:
