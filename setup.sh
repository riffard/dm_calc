#!/usr/bin/env bash

ROOT_VERSION="${ROOT_VERSION:-6}"
LZBUILD_VERSION="${LZBUILD_VERSION:-latest}"

export ROOT_VERSION
export LZBUILD_VERSION

if [ "x${BASH_ARGV[0]}" = "x" ]; then
    if [ -f setup.sh ]; then
        BACC_ROOT="$PWD"
    else
        echo ERROR: must "cd where/BACCARAT/is" before calling ". setup.sh" for this version of bash!
        BACC_ROOT=
        return 1
    fi
else
    # get param to "."
    thisdir=$(dirname ${BASH_ARGV[0]})
    BACC_ROOT=$(cd ${thisdir} > /dev/null;pwd)
fi
export BACC_ROOT


if [ "X" == "X${LZBUILD_DIR}" ] ; then
    LZBUILD_DIR=/cvmfs/lz.opensciencegrid.org/LzBuild/release-${LZBUILD_VERSION}/LzBuild
fi
source ${LZBUILD_DIR}/setup.sh


if [ "${ROOT_VERSION}" == 6 ]; then
  echo "Using ROOT 6"
  source /cvmfs/sft.cern.ch/lcg/releases/LCG_79/ROOT/6.04.02/${LZ_BUILD_CONFIG}/bin/thisroot.sh
elif [ "${ROOT_VERSION}" == 5 ]; then
  echo "Using ROOT 5"
  # ROOT5 is no longer supported, so can hardwire last valid version
  source /cvmfs/lz.opensciencegrid.org/ROOT/v5.34.32/slc6_gcc44_x86_64/root/bin/thisroot.sh
fi

