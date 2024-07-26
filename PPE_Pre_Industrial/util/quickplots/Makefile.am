FC = @NOMPI_FC@

noinst_PROGRAMS = table_diff.x table.x

table_diff_x_SOURCES = mo_tables_diff.f90 mo_util_string.f90 table_diff.f90
table_x_SOURCES = mo_tables.f90 mo_util_string.f90 table.f90

table.$(OBJEXT): $(CONFIG_HEADER) mo_tables.$(OBJEXT) \
	mo_util_string.$(OBJEXT) $(top_builddir)/config/config.inc
table_diff.$(OBJEXT) table.$(OBJEXT): $(CONFIG_HEADER) \
	mo_tables_diff.$(OBJEXT) mo_util_string.$(OBJEXT) \
	$(top_builddir)/config/config.inc

$(top_builddir)/config/config.inc:
	$(am__cd) $(top_builddir) && $(MAKE) $(AM_MAKEFLAGS) config/config.inc

CLEANFILES = mo_tables.mod mo_tables_diff.mod mo_util_string.mod

AM_FCFLAGS = $(DEFAULT_INCLUDES) $(FCDEFS)
