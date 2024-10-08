#!/usr/bin/env perl
#
# f2html, html source converter for fortran 77/90 code
# Copyright (C) 1997,98  Beroud Jean-Marc

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.    

# $Id: f2html.pl.in,v 1.1 1999/07/02 08:30:32 m214003 Exp $

# need version >=5.000
require 5.000;

# modules
# use diagnostics;
use strict;
use Getopt::Long;
use Cwd;

# strip path from the calling program name
$0 =~ s/.*\///;

# force buffer flush
$| = 1;

#==============================================================================

# global vars
use vars qw(@dirs $fgenrc $htmldir $f77header $help);
use vars qw($fgen_datadir $copyright $calling_dir $date $incs);

use vars qw(%rcvars %db %vars %programs %modules %includes %entries
	    %subroutines %functions %interfaces %targets);

use vars qw($keywords1 $keywords2 $keywords3 $keywords4 $keywords5
	    $keywords6);

#==============================================================================


$fgenrc       = "_fgenrc_";
$fgen_datadir = "/afs/dkrz.de/pf/m/m214030/share/fgen";
$copyright    = "$0 v0.3 (C) 1997,98 Beroud Jean-Marc";
$calling_dir  = cwd();
$date         = `date`; chomp($date);

&parse_cmdline();

$rcvars{fc}{suffixes} |= ".h";

$incs = join("|", split(/\s+/, ".h $rcvars{fc}{suffixes}"));
$incs =~ s/\.//og;

if ($htmldir) {
    if ($htmldir !~ /^s*\//) {
	$htmldir = "$calling_dir/$htmldir";
    }
    else {
	$htmldir =~ s/\/$//; # remove trailing / if necessary
    }
}
else {
    if ($rcvars{html}{dir} !~ /^s*\//) {
	$htmldir = "$calling_dir/$rcvars{html}{dir}";
    }
    else {
	$htmldir = "$rcvars{html}{dir}";
    }

    $htmldir =~ s/\/$//;
}

 VARS:{
     my($type, $dir, $var);

     foreach $type (@dirs) {
	 # reset counter
	 my $num = 0;

	 foreach $dir (@dirs) {
	     # build varname
	     $var = "$type";
		 
	     # change dir and store absolute pathname
	     chdir("$dir") || die "$0: can't chdir to $dir: $!\n";

	     # store absolute pathdir
	     $dir = cwd();

	     # store varname (srcdir[.], libdir[.])
	     $vars{$dir} = "$var";

	     # back to calling dir
	     chdir("$calling_dir") ||
		 die "$0: can't chdir to $calling_dir: $!\n";
	 }
     }
}

 SCAN:{
     my($type, $dir);

     &init_html();

     foreach $dir (@dirs) {
	 print "scanning ", abs2rel($dir, $calling_dir), "\n";

	 # change dir and store absolute pathname
	 chdir("$dir") || die "$0: can't chdir to $dir: $!\n";
	 
	 &html_scan($dir);

	 # back to calling dir
	 chdir("$calling_dir") || die "$0: can't chdir to $calling_dir: $!\n";
     }

     print "\n";
 }

 CREATE:{
     my($type, $dir, $file);

     if (! -e $htmldir) {
	 mkdir ("$htmldir", 0755) || die "can't mkdir $htmldir: $!\n";
     }

     # change dir and store absolute pathname
     chdir("$htmldir") || die "$0: can't chdir to $htmldir: $!\n";

     $htmldir = cwd();

     # back to calling dir
     chdir("$calling_dir") || die "$0: can't chdir to $calling_dir: $!\n";

     &create_index($htmldir);
     &create_source_index($htmldir);

     foreach $dir (@dirs) {
	 print "htmlzing ", abs2rel($dir, $calling_dir), " ...\n";

	 foreach $file (sort keys %{ $db{$dir} }) {
	     &create_html_pages($dir, $file);
	 }
     }

     print "\nhtml conversion completed\n";
     print "URL: file:$htmldir/index.html\n";
 }

# shell return value
exit 0;

#==============================================================================


sub parse_cmdline {

    # parse command line switches

    &GetOptions("d=s"       => \$htmldir  , # htmldir
		"h"         => \$help     , # show help
                "77"        => \$f77header, # *.h in f77 dialect
                "f=s"       => \$fgenrc   ) # alternate config file
	|| (&show_usage() && exit 1);
    
    # help
    if ($help) { &show_usage(); exit 0 }

    # make or depend mode is mandatory
    unless ($htmldir) {
	print "$0: htmldir not specified. Try `$0 -h' for more information.\n";
	exit 1;
    }

    # test the avilability of the config file
    if ("$fgenrc" eq "_fgenrc_") {
	$fgenrc = "fgenrc";

	if (-e "$ENV{HOME}/.$fgenrc") {
	    $fgenrc = "$ENV{HOME}/.$fgenrc";
	}
	elsif (-e "$fgen_datadir/$fgenrc") {
	    $fgenrc = "$fgen_datadir/$fgenrc";
	}
	else {
	    print "$0: no config file found\n";
	    print "please take one from the directory $fgen_datadir\n";
	    exit 2
	}
    }
    elsif (! -e "$fgenrc") {
	print "$0: no such config file: ($fgenrc)\n";
	exit 2;
    }

    # parse config file
    &parse_rcfile($fgenrc);

    # remaining arguments ar search directories
    @dirs = @ARGV if (@ARGV);

    # if none dir specified
    push(@dirs, $calling_dir) unless (@dirs);
}

#==============================================================================

sub show_usage {

    # print usage on stdout if called without arguments

    print <<EOF;
$copyright
$0 -d dir [-h] [-77] [-f fgenrc] [dir [dir] ...]

-d dir   : html output directory
-h       : show this page
-77      : header files (*.h ...) are in f77 (default is f90)
-f fgenrc: specify an alternate config file for $0

dir      : search dirs (default is the current working directory)

$0 scans for *.F90 *.f90 *.F *.f *.h files  + user defined file extensions
EOF
}

#==============================================================================

sub parse_rcfile() {

    # open and parse user configuration file.

    my($file) = @_;

    # slurp file into memory
    open (RCFILE, "$file") || die "cannot open $file:$!\n";

    my @file_content = <RCFILE>; close RCFILE;

    # parse file
    my($part, $key, $value);

    foreach (@file_content) {
	# skip those
	next if /^\s*(\#|$)/;

	chomp;

	# trim white spaces
	s/^\s+|\s+$//og;

	if (/^\[(.*)\]$/) {
	    $part = $1;
	}
	else {
	    ($key, $value) = split(/\s*=\s*/, $_, 2);

	    # store peers
	    $rcvars{$part}{$key} = $value;
	}
    }
}

#==============================================================================

sub which {

    # search execs in $PATH (different names for a same program e.g:
    # gcc, cc, lcc) passed as arguments. Return the first occurence found.
    # If not successfull, return empty

    # arguments
    my @bins = @_;

    # return if arguments are empty
    return undef if ($#bins == 0 && $bins[0] eq "");

    # local variables
    my @paths = split (":", $ENV{PATH});

    my($path, $bin);

    foreach $bin (@bins) {
	# absolute path
	if ($bin =~ /^\s*\//) {
	    if (-e "$bin" && -x "$bin" && ! -d "$bin") {
		return "$bin";
	    }
	}
	# relative/without path 
	else {
	    foreach $path (@paths) {
		if (-e "$path/$bin"  && -x "$path/$bin" &&
		    ! -d "$path/$bin") {
		    return "$path/$bin";
		}
	    }
	}
    }

    # none found
    return undef;
}

#==============================================================================

sub html_scan {

    # argument
    my($dir) = @_;

    my $file;

    my $exts = "\*\." . join(" \*\.", split(/\|/, $incs));

    # scan directory
    foreach $file (glob "*.[fF] *.f90 *.F90 $exts") {
	my(@includes, @uses, @modules, @subroutines, @programs, @entries,
	   @labels, @functions, @calls, @externals, @interfaces,
	   $interface_body);
 
	my $f77 = 0;

	# open file for reading
	open (FILE, "$file") || die "$0: can't open $file: $!\n";

	# slurp file in memory
	my @file_content = <FILE>; close FILE;

	if    ($file =~ /\.(f|F)$/o) {
	    $f77 = 1;
	}
	elsif ($file =~ /\.($incs)$/o) {
	    $f77 = 1 if $f77header;
	}

	my @state;

	foreach (@file_content) {
	    # skip those
	    next if /^\s*(!|$)/o;

	    # comment lines (fixed format)
	    next if ($f77 && /^[^\d^\s]/o);

	    # labels
	    push(@labels, /^\s*(\w+)\s*:[^:]/o);

	    # assume: one statement, at the beginning of a line
	    # include, program, entry, subroutine and function keywords
	    if (/^\s*include\s*(?:\'(.*?)\'|\"(.*?)\")/oi) {
		push(@includes, $1 || $2); next;
	    }

	    if (/^\s*module\s+(\w+)/oi && $1 !~ /procedure/oi) {
		push(@state, "m");
		push(@modules, $1); next;
	    }

	    if (/^\s*program\s+(\w+)/oi) {
		push(@state, "p");
		push(@programs, $1); next;
	    }

	    if (/^\s*entry\s+(\w+)/oi) {
		push(@entries, $1);
	    }

	    # FUNCTION/SUBROUTINE CODE (IN MODULE/CONTAINS ...)
	    if (/^\s*contains\b/oi) { push(@state, "c"); next }

	    if(/^\s*end\s*(do|if|interface|select|type|where)\b/oi) {
		next;
	    }

	    if (/^\s*end\b/oi) { pop(@state) }
		    
	    if (/^\s*interface\b/oi) { push(@state, "i") }

	    if (/^\s*interface\s+(\w+)/oi &&
		$1 !~ /operator|assignment/oi) {
		push(@interfaces, $1); push(@state, "i"); next;
	    }

	    if (/^\s*end\s+interface\b/oi) {
		print "error\n" if (pop(@state) ne "i"); next;
	    }
		 
	    # WRONG (interface)
	    if (! @state[-1] ne "i") {
		if (/^\s*subroutine\s+(\w+)/oi) {
		    push(@subroutines, $1); push(@state, "s"); next;
		}

		# remove inline comments to be sure
		if (/!/o) { $_ = skip_inline_comments($_) }

		if (/\bfunction\s+(\w+)/oi) {
		    push(@functions, $1); push(@state, "f"); next;
		}
	    }

	    if (/\bcall\s/oi) {
		# remove inline comments to be sure
		if (/!/o) { $_ = skip_inline_comments($_) }
		    
		push(@calls, /\bcall\s+(\w+)/oig); next;
	    }

	    if (/\buse\s/oi) {
		# remove inline comments to be sure
		if (/!/o) { $_ = skip_inline_comments($_) }

		push(@uses, /\buse\s+(\w+)/oig);
	    }
	}

	$db{$dir}{$file}{includes   } = [ &uniq  (@includes) ];
	$db{$dir}{$file}{uses       } = [ &uniqtr(@uses    ) ];
	$db{$dir}{$file}{calls      } = [ &uniqtr(@calls   ) ];
	$db{$dir}{$file}{labels     } = [ &uniqtr(@labels  ) ];
	$db{$dir}{$file}{programs   } = [ @programs  = &uniqtr(@programs ) ];
 	$db{$dir}{$file}{entries    } = [ @entries   = &uniqtr(@entries  ) ];
	$db{$dir}{$file}{functions  } = [ @functions = &uniqtr(@functions) ];
	$db{$dir}{$file}{modules    } = [ @modules   = &uniqtr(@modules  ) ];

	$db{$dir}{$file}{subroutines} =
	    [ @subroutines = &uniqtr(@subroutines) ];
	$db{$dir}{$file}{interfaces } =
	    [ @interfaces  = &uniqtr(@interfaces ) ];

	foreach (@programs   ) { $programs{$_}    = $file }
	foreach (@entries    ) { $entries{$_}     = $file }
	foreach (@functions  ) { $functions{$_}   = $file }
	foreach (@modules    ) { $modules{$_}     = $file }
	foreach (@subroutines) { $subroutines{$_} = $file }
	foreach (@interfaces ) { $interfaces{$_}  = $file }

	if ($file =~ /\.($incs)$/) { $includes{$file} = $file }
    }
}

#==============================================================================

sub skip_inline_comments {

    # argument
    my($code) = @_;

    if ($code =~ /\'|\"/o) {
	my @chunks = split(/(\".*?\"|\'.*?\'|\s*!.*$)/o, $code);

	pop @chunks if ($chunks[-1] =~ /^!/o);

	$code = join("", @chunks);
    }
    else {
	$code =~ s/!.*$//o;
    }

    return $code;
}

#==============================================================================

sub uniqtr {

    # return only uniq keywords

    # argument
    my @list = @_;

    my %uniqs;

    foreach (@list) {
	tr [A-Z] [a-z];
	$uniqs{$_} = undef;
    }

    return sort keys %uniqs;
}

#==============================================================================

sub uniq {

    # return only uniq keywords

    # argument
    my @list = @_;

    my %uniqs;

    foreach (@list) {
	$uniqs{$_} = undef;
    }

    return sort keys %uniqs;
}

#==============================================================================

sub init_html {

    my @types       = ('integer(?:\s*\*\d+)?', 'character(?:\s*\*\d+)?',
		       'complex(?:\s*\*\d+)?', 'double\s+precision',
		       'logical(?:\s*\*\d+)?', 'real(?:\s*\*\d+)?');

    # not on word boundary: ")" at the end
    my @types2      = ('character\s*(?:\*?\([\s\w=*]+\))',
		       'type\s*\([\s\w]+\)');

    my @attribs_f77 = qw(common data dimension equivalence external
			 implicit(?:\s*none)? namelist parameter save);
    my @attribs_f90  = qw(allocatable intrinsic optional pointer private
			  public sequence target);

    # remark: call cpu_time appears in f95
    my @f90_call     = qw(cpu_time date_and_time mvbits
			  random_(?:number|seed) system_clock);
    my @attributes  = (@attribs_f77, @attribs_f90);

    # not on word boundary: ")" at the end
    my @attributes2 = ('dimension\s*\([\s\w,:]+\)',
		       'intent\s*\([\s\w]+\)');

    my @anchors     = ("use", "call", "include", "entry");

    my @intrinsics  = ("[cdi]?abs", "(?:i?a|i?)char", "(?:da|[acd]?)cos",
		       "d?cosh", "adjust[l|r]", "aimag", "[acd]?log",
		       "(?:idn?|[ad]?n|[ad]?)int", "all(?:ocated)",
		       "[ad]?log10", "amax[0|1]", "max[01]?", "amin[0|1]",
		       "min[01]?", "[ad]?mod", "(?:da|[acd]?)sin", "d?sinh",
		       "[ul]bound", "[di]?dim", "[di]?sign", "d?atan2?",
		       "d?tanh", "tanh?", "[cd]?exp", "[cd]?sqrt", "any",
		       "associated", "bit_size", "btest", "ceiling", "cmplx",
		       "conjg", "count", "cshift", "dmax1", "dmin1", "dprod",
		       "date_and_time", "dble", "dot_product", "dprod",
		       "eoshift", "epsilon", "exponent", "float", "floor",
		       "fraction", "huge", "iand", "ibclr", "ibits", "ibset",
		       "ifix", "ieor", "index", "ior", "ishftc?", "kind",
		       "len", "len_trim", "lge", "lgt", "lle", "llt",
		       "logical", "matmul", "maxeponent", "maxloc", "maxval",
		       "merge", "minexponent", "minloc", "minval", "modulo",
		       "mvbits", "nearest", "not", "pack", "precision",
		       "present", "product", "radix", "random_number",
		       "random_seed", "range", "real", "repeat", "reshape",
		       "rrspacing", "scale", "scan", "selected_int_kind",
		       "selected_real_kind", "set_exponent", "shape", "sngl",
		       "size", "spacing", "spread", "sum", "system_clock",
		       "tiny", "transfer", "transpose", "trim", "unpack",
		       "verify");

    my @intrinsics2 = ("close", "(?:de)?allocate", "format", "inquire",
		       "nullify", "open", "print", "read", "result",
		       "write", "assign");

    # not on word boundary
    my @separators  = ('\.eq\.', '==', '\.ne\.', '/=', '\.lt\.', '&lt=',
		       '\.le\.', '&lt', '\.gt\.', '&gt=', '\.ge\.', '&gt',
		       '\.and\.', '\.or\.', '\.eqv\.', '\.neqv\.', '\.not\.');

    my @blocks      = ("case", "contains", "continue", "cycle", "else",
		       'end\s*blockdata(?:\s+\w+)?', "blockdata",
		       'end\s*do(?:\s+\w+)?', "do", 
		       'end\s*if(?:\s+\w+)?', "if", 
		       '(?:end\s*)?select(?:\s+\w+)?', 'end\s*$', "exit", #'
		       'go\s*to', "return", "then", "stop", "then", "while",
		       '(?:\w+\s*:\s+)', '(?:end\s*)?(?:else)?where');

    $keywords1 = join("|", sort(@types2, @attributes2));
    $keywords2 = join("|", sort(@blocks, @types, @attributes));
    $keywords3 = join("|", sort @anchors);
    $keywords4 = join("|", sort @intrinsics);
    $keywords5 = join("|", sort @intrinsics2);
    $keywords6 = join("|", @separators);


}

#==============================================================================

sub create_index {

    my $indexfile = "$htmldir/index.html";
    my $frame2file;

    foreach (sort keys %programs) { # take the first program
	$frame2file = "$programs{$_}.html"; last;
    }
    
    if (! $frame2file) {
	foreach (sort keys %modules) { # take the first program
	    $frame2file = "$modules{$_}.html"; last;
	}

	if (! $frame2file) {
	    foreach (sort keys %subroutines) { # take the first subroutine
		$frame2file = "$subroutines{$_}.html"; last;
	    }
	}

	if (! $frame2file) {
	    foreach (sort keys %functions) { # take the first function
		$frame2file = "$functions{$_}.html"; last;
	    }
	}
    }

    # create htmlfile
    open(HTMLFILE, ">$indexfile") || die "$0: can't open $indexfile: $!\n";

    print HTMLFILE <<HTML;
<!DOCTYPE HTML PUBLIC-//IETF//DTD HTML 3.2 Strict Level 1//EN>
<HTML>
<!-- $copyright -->
<HEAD>
<TITLE>Index</TITLE>
<BASE TARGET=source_code>
</HEAD>
<FRAMESET COLS="20%,*">
  <FRAME SRC=source_index.html NAME=source_index>
  <FRAME SRC=$frame2file NAME=source_code>
</FRAMESET>
</HTML>
HTML
    close HTMLFILE;
}

#==============================================================================

sub create_source_index {

    my $indexfile = "$htmldir/source_index.html";

    # create htmlfile
    open(HTMLFILE, ">$indexfile") || die "$0: can't open $indexfile: $!\n";

    print HTMLFILE <<HTML;
<!DOCTYPE HTML PUBLIC-//IETF//DTD HTML 3.2 Strict Level 1//EN>
<HTML>
<!-- $copyright -->
<HEAD>
<TITLE>Source Index</TITLE>
<BASE TARGET=source_code>
</HEAD>
<BODY>
<H3>Index</H3>
HTML

    my @items = ([\%programs   , 'programs'   , 'p'],
		 [\%modules    , 'modules'    , 'm'],
		 [\%subroutines, 'subroutines', 's'],
		 [\%functions  , 'functions'  , 'f'],
		 [\%entries    , 'entries'    , 'e'],
		 [\%interfaces , 'interfaces' , 'w'],
		 [\%includes   , 'includes'   , 'i']
		 );

    foreach (@items) {
	my ($ref, $title, $key) = @{ $_ };

	# list
	print HTMLFILE "<UL>\n<LI TYPE=square><B>$title</B>\n";
	print HTMLFILE "  <UL>\n";

	foreach (sort keys %$ref) {
	    print HTMLFILE "  <LI TYPE=disc>",
	          &hrefs("$key", $_, "source_index.html", 1), "\n";
	}

	print HTMLFILE "  </UL>\n</UL>\n";
    }

    print HTMLFILE "</BODY>\n</HTML>";

    close HTMLFILE;
}

#==============================================================================

sub create_html_pages {

    # arguments
    my ($dir, $file) = @_;

    # open source file
    open (FILE, "$dir/$file") || die "$0: can't open $dir/$file: $!\n";

    # slurp file in memory
    my @file_content = <FILE>; close FILE;

    # colors
    my $c_anchors    = "$rcvars{html}{c_anchors}";
    my $c_strings    = "$rcvars{html}{c_strings}";
    my $c_blocks     = "$rcvars{html}{c_blocks}";
    my $c_intrinsics = "$rcvars{html}{c_intrinsics}";
    my $c_comments   = "$rcvars{html}{c_comments}";

    my $htmlfile = "$htmldir/$file.html";

    # create htmlfile
    open(HTMLFILE, ">$htmlfile") || die "$0: can't open $htmlfile: $!\n";

    # hreferences
    #------------
    my $includes    = join("|", sort @{ $db{$dir}{$file}->{includes   } } );
    #my $externals   = join("|", sort @{ $db{$dir}{$file}->{externals  } } );
    # FGEN03

    # anchors
    #--------
    my $programs    = join("|", sort @{ $db{$dir}{$file}->{programs   } });
    my $modules     = join("|", sort @{ $db{$dir}{$file}->{modules    } });
    my $subroutines = join("|", sort @{ $db{$dir}{$file}->{subroutines} });
    my $functions   = join("|", sort @{ $db{$dir}{$file}->{functions  } });
    my $entries     = join("|", sort @{ $db{$dir}{$file}->{entries    } });
    my $interfaces  = join("|", sort @{ $db{$dir}{$file}->{interfaces } });
    my $labels      = join("|", sort @{ $db{$dir}{$file}->{labels     } });


    # create html header
    print HTMLFILE &header($file);

    if ($file =~ /\.($incs)$/) {
	print HTMLFILE "<A NAME=$file></A>";
    }

    # utility
    my($line, $interface_body, $f77) = (0, 0, 0);

    $f77 = 1 if ($file =~ /\.(f|F|h)$/);
    $f77 = 0 if ($1 eq "h" && ! $f77header);

    # main loop
    foreach (@file_content) {
	
	# line counter (each 5 lines)
	if ((++$line)%5 == 0) {
	    printf HTMLFILE ("%4d: ", $line);
	}
	else {
	    print HTMLFILE " "x6;
	}

	# to avoid confusing netscape
	s/</\&lt/g;
	s/>/\&gt/g;
	
	# skip empty line
	if (/^\s*$/) { print HTMLFILE; next }

	# remove line feed
	chomp;

	# comment lines (fixed format)
	if ($f77 && /^[^\d^\s]/) {
	    print HTMLFILE &color($c_comments, $_), "\n"; next;
	}

	# comment line
	if (/^\s*!/) {
	    print HTMLFILE &color($c_comments, $_), "\n"; next;
	}

	# includes
	s/\b(include\s*)([\"\'].*[\"\'])\s*$/&hrefs($1, $2, $file)/ie;

	# split line
	my @chunks  = split(/(\".*?\"|\'.*?\'|\s*!.*$)/);
	my $comment = pop @chunks if ($chunks[-1] =~ /^\s*!/);

	foreach (@chunks) {
	    next if ! defined;
		
	    if (/^(\"|\')/) {
		print HTMLFILE &color($c_strings, $_) ; next;
	    }

	    if (/\bend\b\s*interface\s*$/i) {
		$interface_body = 0;
	    }
	    elsif (/\binterface\b\s*$/i) {
		$interface_body = 1;
	    }

	    # hreferences
	    #------------

	    # FGEN03
	    # external
	    #if (! $interface_body) {
		#if ($externals && ! /\bfunction\b/i) {
		 #   s/\b($externals)\b/&hrefexts($1)/ieg;
		#}
	    #}

	    # labels
	    if ($labels && /\b(end\s*\w+\s*($labels)|goto|cycle|exit)\b/i) {
		s/\b($labels)\b/&hrefexts($1, $file)/ieg;
	    }

	    # call, use
	    s/(\b(?:call|use|module\s+procedure)\b\s+)\b(\w+)\b/
		&hrefs($1, $2, $file)/gei;

	    # includes
	    s/\b(include\s*)([\"\'].*[\"\'])\s*$/&hrefs($1, $2, $file)/ie;

	    # colorize and put anchors (assume only one per code line)
	    #--------
	    # subroutines

	    if (s/\b(end\b\s*)?(subroutine\s+)(\w+)/
		&color($c_anchors, "$1$2$3")/ie) {
		    

		unless ($1 =~ /end/i) {
		    s/\b(subroutine\s+)($3)\b/&anchor($1, $2)/ie;
		}
	    }
	    # functions
	    elsif (s/\b(end\s+)?(function\s+)(\w+)\b/
		   &color( $c_anchors, "$1$2$3" )/ie) {

		unless ($1 =~ /end/i) {
		    s/\b(function\s+)($3)\b/&anchor($1, $2)/ie;
		}
	    }
	    # modules
	    elsif (s/\b(end\s+)?(module\s+)(\w+)\b/
		   &color($c_anchors, "$1$2$3")/ie) {
		
		unless ($1 =~ /end/i) {
		    s/\b(module\s+)($3)\b/&anchor($1, $2)/ie;
		}
	    }
	    # type
	    elsif (s/\b(end\s+)?(type\s+\w+)\b/
		   &color($c_anchors, "$1$2")/ie) {
	    }
	    # interface 1
	    elsif (s/\b(end\s+interface)\b/&color($c_anchors, $1)/ie) {
	    }
	    # interface 2
	    elsif (s/\b(interface)(\s+\w+)?\b/
		   &color($c_anchors, "$1$2")/ie) {
		
		unless ($2 =~ /\s+operator/i || ! $2) {
		    s/\b(interface)($2)\b/&anchor($1, $2)/ie;
		}
	    }
	    # programs
	    elsif (s/\b(end\s+)?(program\s+)(\w+)\b/
	       &color($c_anchors, "$1$2$3")/ie) {
		
		unless ($1 =~ /end/i) {
		    s/\b(program\s+)($3)\b/&anchor($1, $2)/ie;
		}
	    }
	    # entry
	    elsif (s/\b(end\s+)?(entry\s+)(\w+)\b/
		   &color($c_anchors, "$1$2$3")/ie) {
		
		unless ($1 =~ /end/i) {
		    s/\b(entry\s+)($3)\b/&anchor($1, $2)/ie;
		}
	    }
	    # labels
	    elsif (s/(^\s*\w+\s*:)[^:]/
		   &color($c_anchors, $1)/ie) {
		# no more /^\s*... because <font ..
		s/(\s*)(\w+)(\s*:[^:])/&anchor($1, $2, $3)/ie;
	    }

	    # colorize
	    #---------
	    s/(\(\s*kind\s*\=\s*\w+\s*\))/&color($c_blocks, $1)/ieg;
            s/\b($keywords1)/&color($c_blocks, $1)/ieg;
            s/\b($keywords2)\b/&color($c_blocks, $1)/ieg;
            s/\b($keywords3)\b/&color($c_anchors, $1)/ieg;
            s/\b($keywords4)(\s*\()/&color($c_intrinsics, $1, $2)/ieg;
            s/\b($keywords5)\b/&color($c_intrinsics, $1)/ieg;
            s/($keywords6)/&color($c_intrinsics, $1)/ieg;

	    print HTMLFILE;
	}

	# comment
	if ($comment) {
	    print HTMLFILE &color($c_comments, $comment), "\n";
	}
	else {
	    print HTMLFILE "\n";
	}
    }

    print HTMLFILE &footer($dir, $file, $date);
    
    close HTMLFILE;
}

#==============================================================================

sub header {

    # argument
    my($file) = @_;

    return <<EOF;
<!DOCTYPE HTML PUBLIC-//IETF//DTD HTML 3.2 Strict Level 1//EN>

<HTML>

<!--$copyright-->

<HEAD>
   <TITLE>$file</TITLE>
   <BASE TARGET = source_code>
</HEAD>

<BODY BGCOLOR = #f5f5f5>
<BODY LINK    = $rcvars{html}{c_anchors}>
<BODY ALINK   = $rcvars{html}{c_anchors}>
<BODY VLINK   = $rcvars{html}{c_anchors}>

<A NAME = toppage><CENTER><B>$file</B></CENTER></A>

<HR WIDTH = "50%" ALIGN = CENTER>

<PRE>
EOF
}

#==============================================================================

sub footer () {
    # arguments
    my ($dir, $file, $date) = @_;

    # local variables
    my @infos;

    foreach ("uses", "includes", "calls") {# FGEN03, "externals") {
	my $counter = 0;

	push(@infos, " "x(11 - length), "<I>$_</I>: ") if
	    @{ $db{$dir}{$file}->{$_} };

	my ($item, @tmp);

	foreach $item (@{ $db{$dir}{$file}->{$_} }) {
	    push(@tmp, $item);

	    if ($counter != 0 && $counter%5 == 0) { push(@infos, " "x14) }

	    if (++$counter%5 == 0) {
		push(@infos, join(", ", @tmp), "\n"); @tmp = ();
	    }
	}

	push(@infos, join(", ", @tmp), "\n") if @tmp;
    }

    return <<EOF;

<HR WIDTH = "50%" ALIGN = CENTER>
<CENTER><B>Info Section</B></CENTER>
 @infos
</PRE>

<A HREF=$file.html#toppage><B>back to top</B></A>

<HR>

ECHAM-6 (C) 2009 Max-Planck-Institut f&uuml;r Meteorologie, Hamburg<br>
$date<p>
<PRE>HTML derived from FORTRAN source by $copyright.</PRE>

</BODY>

</HTML>
EOF
}

#==============================================================================

sub color {
    # arguments
    my ($color, $text, $moretext) = @_;

    $moretext = "" unless $moretext;

    return "<FONT COLOR=$color>$text</FONT>$moretext";
}

#==============================================================================

sub hrefs {
    # argument (use|call|include)
    my ($key, $name) = my ($searchkey, $searchname, $file, $cut) = @_;

    # local variables
    my $htmlfile;

    # lowercase and strip quotes
    $searchkey  =~ tr/[A-Z]/[a-z]/;
    $searchname =~ tr/[A-Z]/[a-z]/;
    $searchname =~ s/\"|\'//g;
    #FGEN03 exchange
    if (($key =~ /use/i     && ($htmlfile = $modules    {$searchname})) ||
        ($key =~ /call/i    && ($htmlfile = $subroutines{$searchname}   ||
	                		    $entries    {$searchname}   ||
				            $interfaces {$searchname})) ||
        ($key =~ /module\s+procedure/i
	                    && ($htmlfile = $subroutines{$searchname}   ||
	                                    $functions  {$searchname})) ||
        ($key =~ /include/i && ($htmlfile = $includes   {$searchname})) ||
        ($key =~ /p/i       && ($htmlfile = $programs   {$searchname})) ||
	($key =~ /m/i       && ($htmlfile = $modules    {$searchname})) ||
	($key =~ /i/i       && ($htmlfile = $includes   {$searchname})) ||
	($key =~ /s/i       && ($htmlfile = $subroutines{$searchname})) ||
	($key =~ /f/i       && ($htmlfile = $functions  {$searchname})) ||
	($key =~ /e/i       && ($htmlfile = $entries    {$searchname})) ||
	($key =~ /w/i       && ($htmlfile = $interfaces {$searchname}))) {

	# append html extension and lowercase
	$htmlfile .= ".html";
	$htmlfile  =~ s#.*/##;

	$key = "" if $cut;

	# success
	return "$key<A HREF=$htmlfile#$searchname>$name</A>";
    }
    else { # failed => not referenced
	print "   -> not referenced: $searchkey $searchname ($file)\n";

	if ($key =~ /include|i/i) {
	    # quote and colorize include file as a string
	    return join ("", "$key",
			 &color($rcvars{html}{c_strings}, $name));
	}
	else {
	    # back without modifications
	    return "$key$name";
	}
    }
}
	
#==============================================================================

sub hrefexts {
    # argument external functions 
    #my ($name) = my ($searchname) = @_;
    my ($searchname, $file) = @_;
    my $name = $searchname;
    # local variables
    my $htmlfile;

    # lowercase
    $searchname =~ tr/[A-Z]/[a-z]/;

    if ($htmlfile  = $functions{$searchname} || $interfaces{$searchname} ||
	$file) {
	# append html extension and lowercase
	$htmlfile .= ".html";
	$htmlfile  =~ s#.*/##;

	# sucess
	return "<A HREF=$htmlfile#$searchname>$name</A>";
    }

    # failed => not referenced
    print "    not referenced: include file $searchname\n";

    # back without modifications
    return "$name";
}

#==============================================================================

sub anchor {
    # arguments: (program|module|subroutine|function|entry), name
    my ($key, $name, $rest) = @_;

    my $refname = $name;

    $refname =~ tr/[A-Z]/[a-z]/;

    $rest ||= "";

    return "$key<A NAME=$refname>$name</A>$rest";
}

#==============================================================================

sub abs2rel {
    #arguments
    my($path, $base) = @_;

    my(@parts, @common, $com, $i);

    # remove leading slash
    $path =~ s/^\///;
    $base =~ s/^\///;

    # common parts
    @parts = split("/", $path);

    foreach (split("/", $base)) {
	if ($_ eq $parts[0]) {
	    push(@common, $_); shift(@parts); ; last if (! @parts);
	}
	else { last }
    }
    
    $com  = join("/", @common);
    $path =~ s/$com[\/]?//;

    if ($base eq $com) {
	return "./$path";
    }
    else {
	$base   =~ s/^[\/]?$com//;
	@common = split("/", $base);

	for($i=0;$i<$#common;$i++) { $path = "../$path" }
	return "$path";
    }
}

#==============================================================================

__END__

# man page
# to create it on a linux system:
# $> pod2man fgen |groff -Tascii -mandoc

=head1 NAME

f2html - html source converter for fortran 77/90 source code

=head1 SYNOPSIS

f2html -B<d> dir [-B<h>] [-B<77>] [-B<f> F<fgenrc>] [dir [dir] ...]

=head1 DESCRIPTION

f2html scans files with the following extension: *.F90 *.f90 *.F *.f *.h
and converts these files to colorized html pages. The calls to subroutines,
functions, include files and modules are linked with their definition. A list
of programs, include file, subroutines, functions, entries and interfaces
definitions is also provided to ease the navigation through the soure code.

=head1 OPTIONS

=over 12

=item -B<d> dir

dir is a directory where f2html will write all html converted pages. There is
no default

=item -B<77>

This switch tells f2html that B<all> header files are in fortran 77. This seems
necessary, because there is no way to be sure if a header file is in fortran 77
or in fortran 90

=item -B<h>

Show usage and exit

=item -B<f> F<fgenrc>

F<fgenrc> is an alternate configuration file for f2html. If this option is not
specified, f2html search in this order: a user specific configuration file
F<~/.fgenrc>, then the system default configuration file
F</afs/dkrz.de/pf/m/m214030/share/fgen/fgenrc>. If this fails, that means that your installation
is broken and f2html exit with an error message

=item dir

List of search directories. If none specified, the current working directory is
taken as default. There is no limit for the number of search directories.

=back

=head1 COLORS

The colors can be customized. There is a color for comments, a color for
strings, a color for blocks (if/endif do/enddo ...), a color for intrinsics
functions and relationships operators and a color for main blocks (program,
module, subroutine, function) and anchors (include, use, calls) which can
given in the configuration file.

=head1 FILES

~/.fgenrc

/afs/dkrz.de/pf/m/m214030/share/fgen/fgenrc

=head1 BUGS

The subroutines date_time, system_clock, random_number and random seed are not
recognised as intrinsics and the warning that f2html issues can be safely
ignored. problems arise when an interface and a module procedure have the
same name

=head1 AUTHOR

Beroud Jean-Marc (ber@sma.ch). Bug reports and suggestions are welcome!

=cut
