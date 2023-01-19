#!/bin/sh

#eval "$(conda shell.bash hook)"
#conda activate /eos/project/c/cms-ecal-calibration/ecal-conda/h4analysis
source /cvmfs/sft.cern.ch/lcg/views/LCG_100/x86_64-centos7-gcc8-opt/setup.sh

export LD_LIBRARY_PATH=./lib:DynamicTTree/lib/:CfgManager/lib/:$LD_LIBRARY_PATH
export ROOT_INCLUDE_PATH=./interface:DynamicTTree/interface/:CfgManager/interface/:$ROOT_INCLUDE_PATH
