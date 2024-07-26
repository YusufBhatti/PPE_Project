require 'mkmf'
load "extconf.rb"
$libs = append_library($libs, "stdc++")
$LDFLAGS  += " ../.libs/libcdipp.a ../../src/.libs/libcdi.a"
$srcs = %w[cdiobj_wrap.cpp]
$objs = %w[cdiobj_wrap.o]
create_makefile('CdiObj')
