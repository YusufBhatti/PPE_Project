lib_LTLIBRARIES = libsupport.la
libsupport_la_SOURCES = \
	util_reshape.c \
	util_convert.c \
	util_pbio.c \
	util_sysinfo.c \
	util_system.c \
	util_timer.c \
	util_backtrace.c \
	util_logging.c \
	util_file.c \
	util_buffer_pool.c \
	util_mem.c \
	util_db_timings.c

libsupport_la_LIBADD =
# if anyone has an SX to test this, feedback would be appreciated
if NEC_SX_BUILD
libsupport_la_LIBADD += rtc_sx.o
endif

rtc_sx.o: rtc_sx.s
	$(AS) -c rtc_sx.s
