#! /usr/bin/env python
#-----------------------------------------------------------------------------
#
# Build general experiment run environment with three target directories.
#
# Note: <experiment_id> consists of the following structure:
#
#       NNNdddd               - NNN is a three character abbreviation for
#                               your name; see at:
#
# http://svn.zmaw.de/dokuwiki/doku.php?id=listofids:list_of_experimenter_ids 
#
#                             - dddd is a four digit experiment number
#
#       In case the experiment should consist of ensemble runs the
#       <experiment_id> gets the number of the respective run appended:
#
#       <experiment_id>_rrrr  - rrrr is a four digit number giving the
#                               respective ensemble member      
#
#-----------------------------------------------------------------------------
#
# 1) run script and post-processing directories
#
#    <experiment_id>/scripts
#    <experiment_id>/post
#
# 2) model work directory
#
#    <experiment_id>
#
# 3) archive (potentially other host required
#
#    standard experiments: 
#
#    <experiment_id>/<dirs>
#
#    ensemble experiments:
#
#    <experiment_id>/<experiment_id>_rrrr/<dirs>
#
#    For both sets the <dirs> directories are given by:
#
#    .../input
#    .../outdata
#    .../restart
#    .../means
#
#    And a tar file containing the full model directory used (including
#    respective binaries) and the run scripts used for production. As well,
#    it should contain a protocol on the experiment.
#    This last step we expect to be add by hand at the end of the experiment.
#
#-----------------------------------------------------------------------------
#
# Following base directories are recommended for computing hosts known to be
# regularlly used:
#
# tornado.dkrz.de:
#
# 1) $WORK/<userid>/<experiment id>/...
# 2) $WRKSHR/<experiment id>/...
# 3) gsiftp://gridftp.dkrz.de/$PRJHOME/arch/<userid>/<experiment id>/...
#
# Remark: PRJHOME is set on tornado
#
# On hurrikan.dkrz.de:
#
# hurrikan support is dropped.
#
# On vip.rzg.mpg.de:
#
# 1) $HOME/<experiment id>/...
#    $HOME/<experiment id>/...
# 2) /ptmp/<userid>/<experiment id>/...
# 3) gsiftp://gridftp.dkrz.de/$PRJHOME/arch/<userid>/<experiment id>/...
#
# Remark: PRJHOME is NOT set on vip
#
#-----------------------------------------------------------------------------

import os
import sys
import re
import pwd
import shutil
import string

from optparse import OptionParser, OptionGroup
from urlparse import urlparse, urlunparse

#-----------------------------------------------------------------------------
# get command line arguments and presets some variables

usage = '%prog [options]' 

parser = OptionParser(usage=usage)

parser.add_option('-m', '--ensemble_member', action='store', type='int',
                  dest='ensemble_member', default=0,
                  help='set current ensemble member number [default=0]')

parser.add_option('-e', '--expid', action='store', type='string',
                  dest='experiment_id', default='xxx9999',
                  help='experiment identifier NNNdddd '
                       '(NNN - your id, dddd the experiment number)')

parser.add_option('-t', '--test', action='store_false',
                  dest='dont_simulate', default=True,
                  help='test only; do not mkdir/write on any file system or file')

group = OptionGroup(parser, 'Root directories',
                    'Predefined once available for tornado and vip.')

group.add_option('-s', '--scripts', action='store', type='string',
                 dest='script_dir', default=None,
                 help='root directory for scripts and post-processing')

group.add_option('-w', '--work', action='store', type='string',
                 dest='work_dir', default=None,
                 help='root directory for model run')

group.add_option('-a', '--archive', action='store', type='string',
                 dest='archive_dir', default=None,
                 help='root directory for archiving '
                 '(could be a full URL including protocol and archiving host)')

parser.add_option_group(group)

(options, args) = parser.parse_args()

#-----------------------------------------------------------------------------
# now store the values and initialize basic directory name

experiment_id = options.experiment_id
if  options.ensemble_member == 0:
    run_number = 0

script_base  = options.script_dir
work_base    = options.work_dir
archive_base = options.archive_dir

script_dir          = os.path.join(experiment_id, 'scripts')
post_dir            = os.path.join(experiment_id, 'post')

work_dir            = os.path.join(experiment_id, 'work')

archive_input_dir   = os.path.join(experiment_id, 'input')
archive_outdata_dir = os.path.join(experiment_id, 'outdata')
archive_restart_dir = os.path.join(experiment_id, 'restart')
archive_means_dir   = os.path.join(experiment_id, 'means')

print

#-----------------------------------------------------------------------------
# this is, where we are

base_path = os.getcwd()
m  = re.search('^.*/echam-\d\.\d\.\d{2}(\-|_)*\w*', base_path)
if m is not None:
    base_dir = m.group()
else:
    print 'Could not find ECHAM identifier in current path!'
    sys.exit(1)

print 'Model base directory: ', base_dir

print

#-----------------------------------------------------------------------------
# get name of echam version

model = re.search('echam-\d\.\d\.\d{2}(\-|_)*\w*', base_path).group()

print 'Model version: ', model

print

#-----------------------------------------------------------------------------
# try to get SVN revision, filter result

cmd = 'svn info ' + base_dir
f = os.popen(cmd, 'r')
svn_info = f.readlines()
f.close()

url = '   URL could not be determined.'
revision = '   SVN revision unknown.'
for line in svn_info:
    m = re.search('^URL.*', line)
    if m is not None:
        url = m.group()
    n = re.search('^Revision.*', line)
    if n is not None:
        revision = n.group()

print 'SVN information: '
print '   ' + url
print '   ' + revision

print

#-----------------------------------------------------------------------------
# determine experiment id, effective user, and host (get hostname without
# numbers with all this clusters around) - for later use

if run_number > 0:
    expno = '%s_%3.3d' % (experiment_id, run_number)
else:
    expno = experiment_id
print 'Current experiment identifier: ', expno

user_name = pwd.getpwuid(os.getuid())[0]
print 'Current user: ', user_name

hostname = re.sub('\d', '', os.uname()[1])
print 'Current host: ', hostname

print

#-----------------------------------------------------------------------------
# now start setting directories according to command line arguments

if hostname == 'tornado':

    prj_home    = os.getenv('PRJHOME')
    if prj_home is None:
        print 'Environment variable $PRJHOME not set.'
        sys.exit(1)

    prj_dir = os.path.join(prj_home, 'arch') 

    archive_protocol = 'gsiftp'
    archive_host     = 'gridftp.dkrz.de'
    archive_root_dir =  os.path.join(prj_dir, user_name)
    archive_url = (archive_protocol, archive_host, archive_root_dir, '', '', '') 
    archive_root = urlunparse(archive_url)

    script_root = os.getenv('WORK')
    if script_root is None:
        print 'Environment variable $WORK not set.'
        sys.exit(1)
    script_root = os.path.join(script_root, user_name)
    
    work_root   = os.getenv('WRKSHR')
    if work_root is None:
        print 'Environment variable $WRKSHR not set.'
        sys.exit(1)

    print 'Setting up directories on tornado: '
    print
    print '   predefined root $WORK/<userid>/<experiment id>/...'
    print '   predefined root $WRKSHR/<experiment id>/...'
    print '   predefined root gsiftp://gridftp.dkrz.de/$PRJHOME/arch/<userid>/<experiment id>/...'     
    

elif hostname == 'vip':

    prj_home    = os.getenv('PRJHOME')
    if prj_home is None:
        print 'Environment variable $PRJHOME not set.'
        sys.exit(1)

    prj_dir = os.path.join(prj_home, 'arch') 

    archive_protocol = 'gsiftp'
    archive_host     = 'gridftp.dkrz.de'
    archive_root_dir =  os.path.join(prj_dir, user_name)
    archive_url = (archive_protocol, archive_host, archive_root_dir, '', '', '') 
    archive_root = urlunparse(archive_url)

    script_root = os.getenv('WORK')
    if script_root is None:
        print 'Environment variable $WORK not set.'
        sys.exit(1)
    script_root = os.path.join(script_root, user_name)

    work_root = '/ptmp'

    print 'Setting up directories on vip: '
    print
    print '   predefined root $HOME/<experiment id>/...'
    print '   predefined root /ptmp/<userid>/<experiment id>/...'
    print '   predefined root gsiftp://gridftp.dkrz.de/$PRJHOME/arch/<userid>/<experiment id>/...'    

else:

    if archive_base is None:
        archive_protocol = 'fileutils'
        archive_host     = ''
        archive_root     = os.path.join(base_dir, 'experiments')
    else:
        archive_url = urlparse(archive_base)
        archive_protocol = archive_url[0]
        archive_host     = archive_url[1]
        archive_root_dir = archive_url[2]
        archive_url = (archive_protocol, archive_host, archive_root_dir, '', '', '') 
        archive_root = urlunparse(archive_url)


    if script_base is None:  
        script_root = os.path.join(base_dir, 'experiments')
    else:
        script_root = script_base

    if work_base is None:  
        work_root = os.path.join(base_dir, 'experiments')        
    else:
        work_root = work_base

    print 'Setting up default directories (unknown host ' + hostname + '): '    

#-----------------------------------------------------------------------------
# now set final absolute directoy names

script_dir          = os.path.join(script_root,  script_dir)            
post_dir            = os.path.join(script_root,  post_dir)           
work_dir            = os.path.join(work_root,    work_dir)           
archive_input_dir   = os.path.join(archive_root, archive_input_dir)  
archive_outdata_dir = os.path.join(archive_root, archive_outdata_dir)
archive_restart_dir = os.path.join(archive_root, archive_restart_dir)
archive_means_dir   = os.path.join(archive_root, archive_means_dir)    

print    
print 'Final directory names:'
print '      ', script_dir    
print '      ', post_dir
print '      ', work_dir         
print '      ', archive_input_dir  
print '      ', archive_outdata_dir
print '      ', archive_restart_dir
print '      ', archive_means_dir  

print

#-----------------------------------------------------------------------------
# now create all remaining directories - start with local file systems
#
# define function for all error handling
#.............................................................................
def make_local_directory (directory):
    if not os.path.exists(directory):            
        try:
            os.makedirs(directory)            
        except IOError, e:
            print 'Could not make directory ', directory
            print e
        else:
            print 'Directory ', directory, ' created.'
    else:
        print 'Directory ', directory, ' exists.'
    print
#.............................................................................
def make_gsiftp_directory (host, directory):    
    cmd  = 'gsissh '
    # gsiftp is not part of the known URL protocols yet,
    # so we have to cheat a bit and replace gsiftp by ftp
    # before using in urlparse

    atmp = re.sub('^gsiftp', 'ftp', directory)
    args = archive_host + ' mkdir -p ' +  urlparse(atmp)[2]  

    print '*** Reduced error control for command: ' + cmd + args

    try:
        os.system(cmd + args)
    except OSError, e:
        print 'Could not make directory ', directory
        print e
    else:
        print 'Directory ', directory, ' created.'
    print
#.............................................................................    

if options.dont_simulate:

    make_local_directory(script_dir)
    make_local_directory(post_dir)    
    make_local_directory(work_dir)    


    # create archive directores, eventually on another host
                
    if archive_protocol == 'fileutils':
        
        if options.dont_simulate:
            
            make_local_directory(archive_input_dir)
            make_local_directory(archive_outdata_dir)
            make_local_directory(archive_restart_dir)
            make_local_directory(archive_means_dir)            
            
    elif archive_protocol == 'gsiftp':

        if options.dont_simulate:

            make_gsiftp_directory(archive_host, archive_input_dir)
            make_gsiftp_directory(archive_host, archive_outdata_dir)
            make_gsiftp_directory(archive_host, archive_restart_dir)
            make_gsiftp_directory(archive_host, archive_means_dir)
            
#-----------------------------------------------------------------------------
# copy job files - not the solution - need to filter files:
# read input file - process - write output files, so it is not a copy!

# HAM5RUN.START --------------------------------------------------------------

try:
    start_job = open(os.path.join(base_dir, 'run/ham5run.cluster'), 'r').readlines()
except IOError, e:
    print 'Could not open ham5run.cluster template file.'
    print e
else:
    print 'Read ham5run.cluster ...'

# parse ham5run start job

if not options.dont_simulate:
    print
    print 78 * '-'
    print

for line in range(len(start_job)):
    if re.search('^EXP=', start_job[line]):
        start_job[line] = 'EXP=' + experiment_id + '\n'
    if re.search('^RERUN=', start_job[line]):
        start_job[line] = 'RERUN=.false.\n'
    if re.search('^base_dir=', start_job[line]):
        start_job[line] = 'base_dir=' + base_dir + '\n'
    if re.search('^script_dir=', start_job[line]):
        start_job[line] = 'script_dir=' + script_dir + '\n'
    if re.search('^work_dir=', start_job[line]):
        start_job[line] = 'work_dir=' + work_dir + '\n'

    if not options.dont_simulate:
        print start_job[line],

if options.dont_simulate:
    try:
        open(os.path.join(script_dir, 'ham5run.start'),"w").writelines(start_job)
    except IOError, e:
        print 'Could not write ham5run.start.'
        print e
    else:
        print 'Write patched ham5run.start ...'

print

# HAM5RUN --------------------------------------------------------------------

try:
    run_job = open(os.path.join(base_dir, 'run/ham5run.cluster'), 'r').readlines()
except IOError, e:
    print 'Could not open ham5run.cluster template file.'    
    print e
else:
    print 'Read ham5run.cluster ...'

# parse ham5run job

if not options.dont_simulate:
    print
    print 78 * '-'
    print
    
for line in range(len(run_job)):
    if re.search('^EXP=', run_job[line]):
        run_job[line] = 'EXP=' + experiment_id + '\n'
    if re.search('^RERUN=', run_job[line]):
        run_job[line] = 'RERUN=.true.\n'
    if re.search('^base_dir=', run_job[line]):
        run_job[line] = 'base_dir=' + base_dir + '\n'
    if re.search('^script_dir=', run_job[line]):
        run_job[line] = 'script_dir=' + script_dir + '\n'
    if re.search('^work_dir=', run_job[line]):
        run_job[line] = 'work_dir=' + work_dir + '\n'

    if not options.dont_simulate:
        print run_job[line],

if options.dont_simulate:
    try:
        open(os.path.join(script_dir, 'ham5run'),"w").writelines(run_job)
    except IOError, e:
        print 'Could not write ham5run.'
        print e
    else:
        print 'Write patched ham5run ...'

print

# JOB1 -----------------------------------------------------------------------

try:
    job1 = open(os.path.join(base_dir, 'run/job1'), 'r').readlines()
except IOError, e:
    print 'Could not read job1.'
    print e
else:
    print 'Read job1 ...'

# parse job1

if not options.dont_simulate:
    print
    print 78 * '-'
    print
    
for line in range(len(job1)):
    if re.search('^EXP=', job1[line]):
        job1[line] = 'EXP=' + experiment_id + '\n'
    if re.search('^base_dir=', job1[line]):
        job1[line] = 'base_dir=' + base_dir + '\n'
    if re.search('^work_dir=', job1[line]):
        job1[line] = 'work_dir=' + work_dir + '\n'
    if re.search('^archive_restart_dir=', job1[line]):
        job1[line] = 'archive_restart_dir=' + archive_restart_dir + '\n'

    if not options.dont_simulate:
        print job1[line],

if options.dont_simulate:
    try:
        open(os.path.join(script_dir, 'job1'),"w").writelines(job1)
    except IOError, e:
        print 'Could not write job1.'
        print e
    else:
        print 'Write patched job1 ...'

print

# JOB2 -----------------------------------------------------------------------

try:
    job2 = open(os.path.join(base_dir, 'run/job2'), 'r').readlines()
except IOError, e:
    print 'Could not read job2.'
    print e
else:
    print 'Read job2 ...'

# parse job2

if not options.dont_simulate:
    print
    print 78 * '-'
    print
    
for line in range(len(job2)):
    if re.search('^EXP=', job2[line]):
        job2[line] = 'EXP=' + experiment_id + '\n'
    if re.search('^base_dir=', job2[line]):
        job2[line] = 'base_dir=' + base_dir + '\n'
    if re.search('^work_dir=', job2[line]):
        job2[line] = 'work_dir=' + work_dir + '\n'
    if re.search('^post_dir=', job2[line]):
        job2[line] = 'post_dir=' + post_dir + '\n'
    if re.search('^archive_outdata_dir=', job2[line]):
        job2[line] = 'archive_outdata_dir=' + archive_outdata_dir + '\n'
    if re.search('^archive_means_dir=', job2[line]):
        job2[line] = 'archive_means_dir=' + archive_means_dir + '\n'

    if not options.dont_simulate:
        print job2[line],

if options.dont_simulate:
    try:
        open(os.path.join(script_dir, 'job2'),"w").writelines(job2)
    except IOError, e:
        print 'Could not write job2.'
        print e
    else:
        print 'Write patched job2 ...'

print

#-----------------------------------------------------------------------------
#
sys.exit(0)
#
#-----------------------------------------------------------------------------
