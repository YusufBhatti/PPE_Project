require 'mkmf'
load "extconf.rb"
$libs = append_library($libs, "stdc++")
$srcs = %w[cdilib_wrap.c]
$objs = %w[cdilib_wrap.o]
create_makefile('CdiLib')
