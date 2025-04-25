#!/usr/bin/perl

use strict;
use warnings;

sub match_version_from_ac_init {
  my $file = shift || die "sub match_version_from_ac_init: missing argument";
  open(F, "<$file") || die "Cannot open $file.";
  while(my $line = <F>) {
    if ($line =~ /^\s*AC_INIT/) {
      $line =~ s/\s+//g;
      if ($line =~ /^AC_INIT\(\[yaxt\],\[([\d\.\w]+)\]/) {
        close(F);
        return $1;
      }
    }
  }
  close(F);
  return undef;
}

sub create_doxyfile {
  my ($version, $in, $out) = @_;
  die "sub create_doxyfile: missing argument"
    unless defined($version) and defined($in) and defined($out);
  open(F, "<$in") || die "Cannot open $in.";
  open(G, ">$out") || die "Cannot open $out.";
  my $match = 0;
  while(my $line = <F>) {
    if ($line =~ /^(\s*PROJECT_NUMBER\s*=\s*)/) {
      print G "PROJECT_NUMBER = $version\n";
      $match++;
    } else {
      print G $line;
    }
  }
  close(G);
  close(F);
  return $match;
}

my ($conf, $doxy_in, $doxy_out) = @ARGV;

die "$0 usage: your_configure.ac Doxygen_template outputfile\n"
  unless defined($conf) and defined($doxy_in) and defined($doxy_out);


my $version = match_version_from_ac_init($conf);
die "Cannot find yaxt version in $conf." unless defined $version;

create_doxyfile($version, $doxy_in, $doxy_out) || die "Found no match in $doxy_in\n";


