dnl
dnl Copyright  (C)  2010  Thomas Jahns <jahns@dkrz.de>
dnl
dnl Version: 1.0
dnl Keywords: configure configure.ac autotools
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
dnl Commentary:
dnl
dnl
dnl
dnl Code:
dnl
dnl ACX_CHECK_SIZEOF_RELATION(TYPE_A, TYPE_B)
dnl sets C_SIZEOF_TYPE_A_IS_GREATER_THAN_SIZEOF_TYPE_B to 1 or 0
dnl depending on actual relation
AC_DEFUN([ACX_CHECK_SIZEOF_RELATION],
  [AS_VAR_PUSHDEF([sizeof_type_a], [ac_cv_sizeof_$1])
   AS_VAR_PUSHDEF([sizeof_type_b], [ac_cv_sizeof_$2])
   AC_CHECK_SIZEOF([$1])
   AC_CHECK_SIZEOF([$2])
   AC_MSG_CHECKING([if sizeof($1) is greater than sizeof($2)])
   AS_IF([test AS_VAR_GET([sizeof_type_a]) -gt AS_VAR_GET([sizeof_type_b])],
     [AS_TR_CPP([C_$1_IS_LARGER_THAN_$2])=1
      AC_MSG_RESULT([yes])],
     [AS_TR_CPP([C_$1_IS_LARGER_THAN_$2])=0
      AC_MSG_RESULT([no])])
   AC_SUBST(AS_TR_CPP([C_$1_IS_LARGER_THAN_$2]))])
dnl
dnl Local Variables:
dnl license-project-url: "https://www.dkrz.de/redmine/projects/show/scales-ppm"
dnl license-default: "bsd"
dnl mode: autoconf
dnl End:
