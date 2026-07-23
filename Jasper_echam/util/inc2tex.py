#! /usr/bin/env python
#-----------------------------------------------------------------------------

import re

namelist_file = open("runctl.inc")
lines = namelist_file.readlines()

for line in lines:
    if line.startswith('NAMELIST'):
        a = re.split('\/', line, 3)
        namelist = a[1]

print '\\begin{table}'
print '\\caption{\\label{'+namelist+'}NAMELIST '+namelist
print '\\begin{tabular}{p{4cm}p{10cm}'

for line in lines:
    line = line[:-1]
    if line.strip() == '': continue          # skip empty line
    if line.startswith('!'): continue        # skip lines starting with !
    if line.startswith('NAMELIST'): continue # skip lines starting with NAMELIST
    line = line.replace('&!','&')
    line = line.replace(',','')    
    line = line.replace('!','&')
    line = line.replace('_','\\_')        
    print line + ' \\\\'

print '\\end{tabular}'    
print '\\end{table}'    
