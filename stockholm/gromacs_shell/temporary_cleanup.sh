#! /bin/bash
# finds and removes all gromacs backed up "#xxx#" files

find . -type f -name '*#*' -exec rm -i {} +