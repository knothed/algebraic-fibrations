#!/bin/sh

# arg1: sage script to run
# arg2: file to write sage output to
# arg3: extra argument for sage script, optional

# relative path somehow doesn't work with sbatch, despite we are in the correct directory
# source install.sh
source ~/algebraic-fibrations/script/install.sh

exec 3>&1 4>&2 >$2 2>&1

sage $1 $3