#!/bin/bash

# untar your R installation. Make sure you are using the right version!
tar -xzf R402.tar.gz
# (optional) if you have a set of packages (created in Part 1), untar them also
tar -xzf packages.tar.gz

# make sure the script will use your R installation, 
# and the working directory as its home location
export PATH=$PWD/R/bin:$PATH
export RHOME=$PWD/R
export R_LIBS=$PWD/packages

# $1 signal; $2 rank; $3 info; $4 dup
# run your script
Rscript Figure4_chtc_sample.R $1 $2 $3 $4