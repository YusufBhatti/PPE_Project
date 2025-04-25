#!/usr/bin/env perl
#
# Usage: createMakefiles.pl
#
# Generate a Makefile in the modules and src directory from ECHAM.
# To run createMakefiles, it will be necessary to modify the first line of
# this script to point to the actual location of Perl on your system.
#
# Written by Uwe Schulzweida <schulzweida@dkrz.de> June, 1999
#
use strict;
use warnings;

my $EXECUTABLE="echam6";
#
my $PROG=`basename "$0"`;
#
# Create Makefile in src
#
my $dir = 'src';
#
chdir $dir;
#
rename "Makefile.am", "Makefile.am.old";
open(MAKEFILE, "> Makefile.am");
print "create new Makefile.am in: $dir\n";
print MAKEFILE <<"EOF"
# Generated automatically by $PROG

AM_FCFLAGS = \$(FCDEFS) \\
	-I\$(abs_top_srcdir)/include \$(FC_MODINC)./ -I\$(abs_top_builddir)/config \\
	\$(CDI_INCLUDE) \\
	\$(COUPLER_INCLUDE) \\
	\$(NETCDF_FINCLUDE) \$(NETCDF_INCLUDE) \$(HDF5_INCLUDE) \$(SZIP_INCLUDE) \$(ZLIB_INCLUDE) \\
	\$(YAXT_INCLUDE) \$(SCT_INCLUDE) \$(MPI_INCLUDE)


CLEANFILES = version.c

.PHONY: create_version_c

create_version_c:
	\$(top_srcdir)/config/pvcs.pl --srcdir \$(top_srcdir)

version.c: | create_version_c

BUILT_SOURCES = version.c

bin_PROGRAMS = $EXECUTABLE

# FIXME: the direct linking to liblapack.la/libblas.la needs to be
# fixed in favour of optimized versions
echam6_LDADD = \\
        \$(top_builddir)/support/libsupport.la \\
	\$(LAPACK_LIB) \\
	\$(CDI_LIB) \\
	\$(COUPLER_LIB) \\
	\$(NETCDF_FLIB) \\
	\$(NETCDF_LIB) \$(HDF5_LIB) \$(SZIP_LIB) \$(ZLIB_LIB) \\
	\$(YAXT_LIB) \$(SCT_LIB) \$(MPI_LIB) \$(POSTGRESQL_LIB)

INCLUDE = \$(top_srcdir)/include

EOF
    ;

#
# Source listing
#
my @srcs = &uniq(<*.f90 *.f *.F *.c>, 'version.c');
print MAKEFILE "${EXECUTABLE}_SOURCES = ";
&PrintWords(8, 0, @srcs);
print MAKEFILE "\n\n";
#
# Add version file
#
#print MAKEFILE "all: create_version_c \$(PROG)\n\n";
#print MAKEFILE "\$(PROG): \$(OBJS) version.o ../lib/libsupport.a\n";
#print MAKEFILE "\t\$(", &LanguageCompiler($ARGV[1], @srcs);
#print MAKEFILE ") \$(LDFLAGS) -o \$@ \$(OBJS) version.o \$(LIBS)\n\n";
#print MAKEFILE ".PHONY: create_version_c\n\n";
#print MAKEFILE "create_version_c:\n";
#print MAKEFILE "\t../config/pvcs.pl --srcdir ..\n\n";
#print MAKEFILE "version.c: | create_version_c\n\n";
#print MAKEFILE "version.o: version.c\n\n";
#
# make install
#
#print MAKEFILE "install: all\n";
#print MAKEFILE "\t\$(MKDIR_P) \$(bindir)\n";
#print MAKEFILE "\t\$(INSTALL) -m 755 ../bin/$EXECUTABLE \$(bindir)\n\n";
#
# make clean
#
#print MAKEFILE "distclean:\n";
#print MAKEFILE "\trm -f \$(PROG) \$(OBJS) *.mod version.c version.o\n\n";
#print MAKEFILE "mostlyclean:\n";
#print MAKEFILE "\trm -f \$(PROG) \$(OBJS) *.mod version.c version.o\n\n";
#print MAKEFILE "clean:\n";
#print MAKEFILE "\trm -f \$(PROG) \$(OBJS) *.mod version.c version.o\n\n";
#
# Make .f90 a valid suffix
#
#print MAKEFILE ".SUFFIXES: \$(SUFFIXES) .f90\n\n";
#
# .f90 -> .o
#
#print MAKEFILE "%.o: %.f90\n";
#print MAKEFILE "\t\$(FC) \$(FCFLAGS) -c \$<\n\n";
#
# Dependency listings
#
&MakeDependsf90($ARGV[1]);
#
# mo_transpose:
#
# Note: PGI compiler seems to interrogate $FCFLAGS from the environment
# and a -fastsse flag specified there will override the command line option
#print MAKEFILE "\nifeq (\$(strip \$(ARCH)), CRAY_XT3)\n";
#print MAKEFILE "mo_transpose.o: mo_transpose.f90\n";
#print MAKEFILE "\t( FCFLAGS=\" \" ; \$(FC) -O2  -c mo_transpose.f90 )\n";
#print MAKEFILE "endif\n\n";
#
close MAKEFILE;
#
exit;
#
#
# &PrintWords(current output column, extra tab?, word list); --- print words
#    nicely
#
sub PrintWords {
   my($columns) = 78 - shift(@_);
   my($extratab) = shift(@_);
   my($wordlength);
   #
   print MAKEFILE $_[0];
   $columns -= length(shift(@_));
   foreach my $word (@_) {
      $wordlength = length($word);
      if ($wordlength + 1 < $columns) {
         print MAKEFILE " $word";
         $columns -= $wordlength + 1;
         }
      else {
         #
         # Continue onto a new line
         #
         if ($extratab) {
            print MAKEFILE " \\\n\t\t$word";
            $columns = 62 - $wordlength;
            }
         else {
            print MAKEFILE " \\\n\t$word";
            $columns = 70 - $wordlength;
            }
         }
      }
   }

#
# &LanguageCompiler(compiler, sources); --- determine the correct language
#    compiler
#
sub LanguageCompiler {
   my ($compiler) = lc(shift(@_));
   my (@srcs) = @_;
   #
   if (length($compiler) > 0) {
      CASE: {
         grep(/^$compiler$/, ("fc", "f77")) &&
            do { $compiler = "F77"; last CASE; };
         grep(/^$compiler$/, ("cc", "c"))   &&
            do { $compiler = "CC"; last CASE; };
         $compiler = "FC";
         }
      }
   else {
      CASE: {
         grep(/\.f90$/, @srcs)   && do { $compiler = "FC"; last CASE; };
         grep(/\.(f|F)$/, @srcs) && do { $compiler = "F77";  last CASE; };
         grep(/\.c$/, @srcs)     && do { $compiler = "CC";  last CASE; };
         $compiler = "???";
         }
      }
   $compiler;
   }

#
# &uniq(sorted word list); --- remove adjacent duplicate words
#
sub uniq {
  my (@words);
  if (@_) {
    @words = $_[0];
    foreach my $word (@_) {
      if ($word ne $words[$#words]) {
        push(@words, $word);
      }
    }
  }
  @words;
}
#
# &MakeDependsf90(f90 compiler); --- FORTRAN 90 dependency maker
#

# Map include file names to path needed for Makefile
sub map_incs {
   my $inc = shift;
   if (-e "../include/$inc") {
      return "\$(INCLUDE)/$inc";
   }
   else {
      print "  $inc: include file $inc not in \$(INCLUDE)\n";
      return ();
   }
}

# Map module names to their corresponding files. Skip self references.
sub map_modules {
   my $obj = shift;
   my $file = shift;
   if($file && $obj ne $file) {
      return $file;
   }
   else {
      return ();
   }
}


sub MakeDependsf90 {
   my(%filename);
   my(@incs);
   my(@modules);
   my($objfile);
   #
   # Associate each module with the name of the file that contains it
   #
   foreach my $file (<*.f90>) {
      open(FILE, $file) || warn "Cannot open $file: $!\n";
      while (<FILE>) {
         /^\s*module\s+([^\s!]+)/i &&
            ($filename{lc($1)} = $file) =~ s/\.f90$/.o/;
      }
   }
   #
   # Print the dependencies of each file that has one or more include's or
   # references one or more modules
   #
   foreach my $file (<*.f90>) {
      my @incs = ();
      my @modules = ();
      ($objfile = $file) =~ s/\.f90$/.o/;
      open(FILE, $file);
      while (<FILE>) {
         /^\s*#*include\s+["\']([^"\']+)["\']/i && push(@incs, $1);
         /^\s*use\s+([^\s,!]+)/i && push(@modules, lc($1));
      }

      @incs = uniq(sort(map(map_incs($_), @incs)));
      @modules = uniq(sort(map(map_modules($objfile, $filename{$_}), @modules)));

      if (@incs || @modules) {
#          print "'$objfile':\n".
#                "    incs: '".join("','", @incs)."'\n".
#                "    mods: '".join("','", @modules)."'\n".
         print MAKEFILE "$objfile: ";
         &PrintWords(length($objfile) + 2, 0, @modules, @incs);
         print MAKEFILE "\n";
      }
   } # for file in *.f90
} # MakeDependsf90
