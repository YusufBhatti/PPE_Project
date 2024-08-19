#!/bin/bash

# Function to set svn:ignore property
set_svn_ignore() {
    local dir=$1
    local patterns=$2

    if [ -d "$dir" ]; then
        echo -e "$patterns" | svn propset svn:ignore -F - "$dir"
        echo "svn:ignore set for directory $dir"
    else
        echo "Directory $dir does not exist or is not under version control."
    fi
}

# Define ignore patterns
ignore_patterns=$(cat <<EOF
*~
*#
*.o
*.mod
/config.log
/config.status
/Makefile
/build_echam6.sh
/bin*/echam6
/config/config.h
/config/stamp-h1
/contrib/generate-scripts/*thunder*
/src/*.mod
/src/*.lst
/src/Make*
/experiments/
/lib/*.a
~
EOF
)

# Set svn:ignore for the root directory and specific paths
set_svn_ignore "." "$ignore_patterns"

# Individual directories for which the ignore patterns should be set
directories=(
    "autom4te.cache"
    "config"
    "bin"
    "contrib/generate-scripts"
    "src"
    "experiments"
    "lib"
)

# Loop through directories and set svn:ignore where applicable
for dir in "${directories[@]}"; do
    if svn info "$dir" > /dev/null 2>&1; then
        set_svn_ignore "$dir" "$ignore_patterns"
    else
        echo "Directory $dir is not under version control."
    fi
done

