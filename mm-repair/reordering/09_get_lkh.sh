#!/bin/bash -l

# download LKH if not in local dir.

# sample invocation:
#   09_get_lkh.sh LKH-3.0.7 http://webhotel4.ruc.dk/~keld/research/LKH-3/LKH-3.0.7.tgz
# must be executed from the mm-repair/reordering directory 

LKH_VERSION=$1
LKH_LINK=$2

if [ ! -d "$LKH_VERSION" ]; then
  wget $LKH_LINK &&
  tar xvfz $LKH_VERSION.tgz && 
  rm -v $LKH_VERSION.tgz  &&
  cd $LKH_VERSION &&
  make &&
  cd ..
fi 