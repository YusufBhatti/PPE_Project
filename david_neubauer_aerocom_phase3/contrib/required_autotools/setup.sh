#! /bin/bash
#
#------------------------------------------------
# 
# Script written by Luis Kornblueh
# Passed to Sylvaine Ferrachat on 2016.04.19
#
# Usage:
# $ cd /path/to/where/your/want/to/install/these/tools
# $ /path/to/your/echam-hammoz/contrib/required_autotools/setup.sh
# --> make sure /path/to/where/your/want/to/install/these/tools/bin is prepended to the PATH
#     of users who will need this.
# --> in the same way, make sure /path/to/where/your/want/to/install/these/tools/share/man is prepended 
#     to the @MANPATH@ of all echam users of our platform.
#
#------------------------------------------------

set -eu

basedir=$(pwd)
export PATH=$basedir/bin:$PATH

srcdir=$(pwd)/src
if [[ ! -d $srcdir ]]
then
    mkdir -p $srcdir
fi

#_____________________________________________________________________________
# get and install sufficiently recent versions of tools

cd $srcdir

m4_version=$(m4 --version | awk 'NR==1{print $4}')
if [[ $m4_version < 1.4.6 ]]
then
    curl -O http://ftp.gnu.org/gnu/m4/m4-1.4.17.tar.gz
    tar xvf m4-1.4.17.tar.gz
    cd  m4-1.4.17
    ./configure --prefix=$basedir
    make install
    cd $srcdir
fi

autoconf_version=$(autoconf --version | awk 'NR==1{print $4}')
if [[ $autoconf_version < 2.69 ]]
then
    curl -O http://ftp.gnu.org/gnu/autoconf/autoconf-2.69.tar.gz
    tar xvf autoconf-2.69.tar.gz
    cd  autoconf-2.69
    ./configure --prefix=$basedir
    make install
    cd $srcdir
fi

automake_version=$(automake --version | awk 'NR==1{print $4}')
if [[ $automake_version < 1.14.1 ]]
then
    curl -O http://ftp.gnu.org/gnu/automake/automake-1.14.1.tar.gz
    tar xvf automake-1.14.1.tar.gz
    cd  automake-1.14.1
    ./configure --prefix=$basedir
    make install
    cd $srcdir
fi

libtool_version=$((libtool --version 2> /dev/null || glibtool --version) | awk 'NR==1{print $4}')
if [[ $libtool_version < 2.4.2 ]]
then
    curl -O http://ftp.gnu.org/gnu/libtool/libtool-2.4.6.tar.gz
    tar xvf libtool-2.4.6.tar.gz
    cd libtool-2.4.6
    ./configure --prefix=$basedir
    make install
    cd $srcdir
fi

svnversion=$(svn --version | awk 'NR==1{print $3}')

if [[ $svnversion < 1.6.11 ]]
then
    curl -O http://apache.openmirror.de/apr/apr-1.5.2.tar.gz
    tar xvf apr-1.5.2.tar.gz
    cd apr-1.5.2
    ./configure --prefix=$basedir --disable-libtool-lock
    make install
    cd $srcdir
 
    curl -O http://apache.openmirror.de/apr/apr-util-1.5.4.tar.gz
    tar xvf apr-util-1.5.4.tar.gz
    cd apr-util-1.5.4
    ./configure --prefix=$basedir --with-apr=$basedir
    make install
    cd $srcdir

    curl -k -O https://www.sqlite.org/2016/sqlite-amalgamation-3110100.zip
    unzip sqlite-amalgamation-3110100.zip

    curl -O http://apache.openmirror.de/subversion/subversion-1.9.3.tar.gz
    tar xvf subversion-1.9.3.tar.gz
    cd subversion-1.9.3
    mkdir sqlite-amalgamation
    cp ../sqlite-amalgamation-3110100/sqlite3.c sqlite-amalgamation/.
    ./configure --prefix=$basedir --with-apr=$basedir/bin/apr-1-config  --with-apr-util=$basedir/bin/apu-1-config --without-apxs
    make install
    cd $srcdir
fi

