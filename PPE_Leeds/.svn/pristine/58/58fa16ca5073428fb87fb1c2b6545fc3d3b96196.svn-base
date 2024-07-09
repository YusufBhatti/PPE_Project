#! /usr/bin/env python
#------------------------------------------------------------------------------

import os
import sys
import re

from optparse import OptionParser

#------------------------------------------------------------------------------
def unique(inlist):
    return list(set(inlist))
#------------------------------------------------------------------------------
def load_file(filename):
    print 'INFO: load file '+filename
    try:
        t = open(filename, 'r')
    except IOError:
        print 'ERROR: could not open '+filename+' ...'
    else:    
        lines = t.readlines()
        t.close()
    return lines
#------------------------------------------------------------------------------
def main():
    usage = '%prog [options]' 
    parser = OptionParser(usage=usage)
    parser.add_option('-e', '--exp-dir',    dest='exp_dir',    help='experiment directory')
    (options, args) = parser.parse_args()
    if options.exp_dir:
        exp_dir = options.exp_dir
    else:
        print 'ERROR: experiments directory missing ...'
        exit (1)
    files = os.listdir(exp_dir)         
    diag_pattern = re.compile('^hdiffdiag*')
    diag_files = filter(diag_pattern.search, files)

    all_lines = []
    for df in diag_files:
        f = os.path.join(exp_dir, df)
        all_lines.extend(load_file(f))

    hd = open('hdiffdiag.dat', 'w')    
    for line in sorted(unique(all_lines)):
        hd.write(line)
    hd.close()

#==============================================================================
if __name__ == "__main__":
    sys.exit(main())
