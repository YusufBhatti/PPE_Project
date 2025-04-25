AC_DEFUN([AX_FC_VERSION],[

AC_REQUIRE([AC_PROG_FC])
AC_REQUIRE([AC_PROG_AWK])

AS_CASE([$FC],
	[gfortran],[ax_fc_version=$(eval $FC --version | $AWK 'NR==1')],
	[pgfortran],[ax_fc_version=$(eval $FC -V | $AWK 'NR==2')],
	[ifort],[ax_fc_version=$(eval $FC --version 2>&1 | $AWK 'NR==1')],
	[nagfor],[ax_fc_version=$(eval $FC -V 2>&1 | $AWK 'NR==1')],
        [crayftn],[ax_fc_version=$(eval $FC -V 2>&1 | $AWK 'NR==1{printf "%s %s %s\n", $[]1, $[]2, $[]5}')],
	[xlf*],[ax_fc_version=$(eval $FC -qversion | $AWK 'NR%2{printf "%s, ",$[]0;next}{printf "(build %s)\n",$[]2;}')],
	[sxf*],[ax_fc_version=$(eval $FC -V 2>&1 | awk 'NR==1;NR==5' | tr -d "\n")],
	[ax_fc_version="unknown"])
])

AC_DEFUN([AX_CC_VERSION],[

AC_REQUIRE([AC_PROG_CC])
AC_REQUIRE([AC_PROG_AWK])

AS_CASE([$CC],
	[gcc],[ax_cc_version=$(eval $CC --version | $AWK 'NR==1')],
	[pgcc],[ax_cc_version=$(eval $CC --version | $AWK 'NR==2')],
	[icc],[ax_cc_version=$(eval $CC --version 2>&1 | $AWK 'NR==1')],
	[craycc],[ax_cc_version=$(eval $CC -V 2>&1 | $AWK 'NR==1{printf "%s %s %s\n", $[]1, $[]2, $[]5}')],
	[xlc*],[ax_cc_version=$(eval $CC -qversion | $AWK 'NR%2{printf "%s, ",$[]0;next}{printf "(build %s)\n",$[]2;}')],
	[sxc*],[ax_cc_version=$(eval sxc++ -V 2>&1 | sed -e 's/.c.*$//' | sed -e 's/C\/C++ Compiler //' | awk 'NR==1;NR==3' | tr -d "\n")],
	[ax_cc_version="unknown"])
])
