#!/bin/bash

################################################################################
# Prepare the environment for pmt_analysis. It uses UPS producs on gpvm machines
# It requres a local installtion of ROOT, Eigen and CMake otherwise
# mailto:ascarpell@bnl.gov
################################################################################

# Setup manually the EIGEN_INSTALL variable which is the most opportune for you
export EIGEN_INSTALL="/usr/include/eigen3/"

# This part of the code is effective only on a fnal.gov domain
if [ "${HOSTNAME#*.}" == "sdcc.bnl.gov" ]; then

  echo "Setup using /cvmfs/ at 'S{HOSTNAME#*.}'"

  # Setup the general icaruscode ups packages
  source /cvmfs/fermilab.opensciencegrid.org/products/larsoft/setup
  #source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh

  remove_quotes()
  {
    var=$1
    temp="${var%\"}"
    temp="${temp#\"}"
    echo "$temp"
  }

  # Find the location of the eigen heders and replace the value set before
  EIGEN_INSTALL="/cvmfs/larsoft.opensciencegrid.org/products/eigen/v3_3_5/include/eigen3/"

  # Now setup the latest version of Eigen using ups producs
  out=( "$(ups list -aK+ eigen | tail -1)" )
  line=( $out )
  version=$( remove_quotes ${line[1]} )
  setup eigen -v ${version}

  # Now setup the latest version of ROOT using ups producs
  out=( "$(ups list -aK+ root | tail -1)" )
  line=( $out )
  version=$( remove_quotes ${line[1]} )
  qualifier=$( remove_quotes ${line[3]} )
  setup root -v ${version} -q ${qualifier}

  #now setup CMAKE
  out=( "$(ups list -aK+ cmake | tail -1)" )
  line=( $out )
  version=$( remove_quotes ${line[1]} )
  setup cmake -v ${version}

fi

# Check if the EIGEN_INSTALL variable is configured
if [ "$EIGEN_INSTALL" == "" ]; then
  echo "DON'T FORGET TO SET THE EIGEN_INSTALL VARIABLE!"
  return
fi

# Check if root is configured
if [ "$ROOTSYS" == "" ]; then
  echo "DON'T FORGET TO SET A ROOT VERSION!"
  return
fi

# Get the current dir
currentdir=$(pwd)

# Add directories to the local path
export PATH="$currentdir/build/:$PATH"
export LD_LIBRARY_PATH="$currentdir/build/:$LD_LIBRARY_PATH"
export PROJECT_PATH="$currentdir"
