#!/bin/bash

#------------------------------------------------------------------------------------------------------
# Sylvaine Ferrachat 2023-04
#
# Developed in the context of the aeroclim project (but may be useful in other contexts).
#
# Reset date (defaults to month 1st) and time (defaults to 0 UTC) in some user-defined stream (defaults to '_tracerm')
#
#------------------------------------------------------------------------------------------------------

set -e

indir="$1"  # Input dir, where to find the original files
outdir="$2" # Output dir, where to create the shifted files

suffix="${3:-_tracerm.nc}"   # Optional, stream suffix. Default: "_tracerm.nc"

new_day="${4:-01}"         # Optional, day number for the reset (Default "01")
new_time="${5:-00:00:00}"  # Optional, time for the reset (Default "00:00:00")

#-- Security: create outdir if necessary
if [ ! -s $outdir ] ; then
    mkdir -p $outdir
fi

#-- Loop over all files in indir that comply the pattern
for infile in `find $indir -type f -name "*${suffix}"` ; do
    echo "Processing: ${infile}..."

    filename=`basename $infile`
    outfile="${outdir}/${filename}"

    #-- Get the orig date from the original file
    #   (Use the 'T' character as splitter here)
    IFS='T' read date time <<< "`cdo -s showtimestamp $infile`"

    #-- Reset day to 1st of the month:
    new_date=`sed -r -e "s|^ *([0-9]{4})-([0-9]{2})-[0-9]{2}|\1-\2-${new_day}|" <<< "$date"` 

    echo "    Resetting $date --> $new_date"
    echo "    Resetting $time --> $new_time"

    #-- Apply the date/time changes
    cdo setdate,$new_date -settime,$new_time $infile $outfile

    echo
done

exit
