#! /usr/bin/env perl
#_____________________________________________________________________________
#
# Luis Kornblueh, MPIM, 2012-05-31 
# Luis Kornblueh, MPIM, 2013-10-16
#                 add functions to return strings with 
#                 the revision information  
#_____________________________________________________________________________

=head1 pversion

pversion - Enabling provenance collection for VCSs

=head1 SYNOPSIS

pversion [options] [VCS working directory]

Options:

   --srcdir  <VCS working directory>
   --help    usage

=cut

=head1 OPTIONS

=over 8

=item B<--srcdir <VCS working directory>>

VCS handled top level directory which gets compiled

=item B<--help>

Prints the usage documentation and exits

=back

=cut

=head1 DESCRIPTION

This program will generate a c source file which contains persistent 
information on the repository, branch, and revision of the compiled
source code. The file is suitable for being linked into the executable.  

=cut

#_____________________________________________________________________________

use strict;
use warnings;
use constant { true => 1, false => 0 };

use Getopt::Long;
use Pod::Usage;
use File::Temp;
use File::Compare;
use File::Copy;

my $srcdir = '';
my $help = 0;

GetOptions( 'help|?' => \$help, 'srcdir=s' => \$srcdir) or pod2usage(-exitstatus => 2);

pod2usage(-exitstatus => 1) if $help;
                                                  
my $remote_url = '';
my $branch = '';
my $revision = '';

if ( -d $srcdir."/.svn" ) {
    my @svn_info = `svn info $srcdir`;
    my @remote_urls = grep(/URL/, @svn_info);
    $remote_url = $remote_urls[0];
    $remote_url =~ s/^URL:[ \t]*//;
    $remote_url =~ s/ *\n//;
    $branch = $remote_url;
    $branch =~ s/\S+(branches|trunk)\///;
    $branch =~ s/ *\n//;
    my @revisions = `svnversion $srcdir`; 
    $revision = $revisions[0];
    $revision =~ s/Revision: */r/;
    $revision =~ s/ *\n//;
} elsif ( -d $srcdir."/.git" ) {
    my @remote_urls = `git --git-dir $srcdir/.git remote -v`;
    @remote_urls = grep(/fetch/, @remote_urls);
    $remote_url = $remote_urls[0];
    $remote_url =~ s/^origin[ \t]*//;
    $remote_url =~ s/ *\(fetch\) *\n//;
    my @branches = `git --git-dir $srcdir/.git branch`;	
    @branches = grep(/^\*/, @branches); 
    $branch = $branches[0];
    $branch =~ s/\* *//;
    $branch =~ s/ *\n//;
    my @revisions = `git --git-dir $srcdir/.git --no-pager log --max-count=1`; 
    @revisions = grep(/commit/, @revisions);
    $revision = $revisions[0];
    $revision =~ s/commit *//;
    $revision =~ s/ *\n//;
} else {
    print "Unknown repository type or no working copy/repository: no support will be given.\n";
    $remote_url = "Unknown";
    $branch = "Unknown";
    $revision = "Unknown";
}

my $need_to_compare = false;
my $version_c;
my $fname;

if ( -e "version.c") {
    $need_to_compare = true;
    $version_c = File::Temp->new(UNLINK => false);
    $fname = $version_c->filename;
} else {
    open $version_c, ">", "version.c" or die "$0: open version.c: $!";
}

print $version_c "#ifdef __xlc__\n";
print $version_c "#pragma comment(user,\"$remote_url,$branch,$revision\")\n";
print $version_c "#endif\n";
print $version_c "#include <string.h>\n";
print $version_c "\n";
print $version_c "const char remote_url[] = \"$remote_url\";\n";
print $version_c "const char branch[] = \"$branch\";\n";
print $version_c "const char revision[] = \"$revision\";\n"; 
print $version_c "\n";
print $version_c "void repository_url(char *name, int *actual_len)\n";
print $version_c "{\n";
print $version_c "  if (strlen(remote_url) > *actual_len)\n";
print $version_c "    {\n";
print $version_c "      *actual_len = 0;\n";
print $version_c "    }\n";
print $version_c "  else\n";
print $version_c "    {\n";
print $version_c "      strcpy(name, remote_url);\n";
print $version_c "      *actual_len = strlen(name);\n";
print $version_c "    }\n";
print $version_c "\n";
print $version_c "  return;\n";
print $version_c "}\n";
print $version_c "\n";
print $version_c "void branch_name(char *name, int *actual_len)\n";
print $version_c "{\n";
print $version_c "  if (strlen(branch) > *actual_len)\n";
print $version_c "    {\n";
print $version_c "      *actual_len = 0;\n";
print $version_c "    }\n";
print $version_c "  else\n";
print $version_c "    {\n";
print $version_c "      strcpy(name, branch);\n";
print $version_c "      *actual_len = strlen(name);\n";
print $version_c "    }\n";
print $version_c "\n";
print $version_c "  return;\n";
print $version_c "}\n";
print $version_c "\n";
print $version_c "void revision_key(char *name, int *actual_len)\n";
print $version_c "{\n";
print $version_c "  if (strlen(revision) > *actual_len)\n";
print $version_c "    {\n";
print $version_c "      *actual_len = 0;\n";
print $version_c "    }\n";
print $version_c "  else\n";
print $version_c "    {\n";
print $version_c "      strcpy(name, revision);\n";
print $version_c "      *actual_len = strlen(name);\n";
print $version_c "    }\n";
print $version_c "\n";
print $version_c "  return;\n";
print $version_c "}\n";
print $version_c "\n";

close $version_c;

if ($need_to_compare) {
    if (compare($fname, "version.c") == 0) {
        unlink $fname;
    } else {
        move($fname, "version.c") or die "Copy failed: $!";
    }
}

exit 0;
