AC_DEFUN([AX_RPATH_HANDLING],[

AC_REQUIRE([AC_CANONICAL_SYSTEM])

AS_CASE([$target_os],
	[aix*],
	[ax_rpath=no
         ax_rpath_split=no
         ax_rpath_prefix=
	],
        [darwin*],
	[ax_rpath="-rpath,"
         ax_rpath_split=no
         ax_rpath_prefix="-Wl,"
	],
	[linux*],
	[ax_rpath="-rpath,"
         ax_rpath_split=no
         ax_rpath_prefix="-Wl,"
	],
	[ax_rpath=no
         ax_rpath_split=no
         ax_rpath_prefix=
	])

AS_IF([test -n $FC],
      [AS_CASE([$FC],
 	       [nagfor], [ax_rpath=",-rpath "
	                  ax_rpath_split=yes
                          ax_rpath_prefix="-Wl,-Wl,"
                         ],
               [])],
      []) 
])
