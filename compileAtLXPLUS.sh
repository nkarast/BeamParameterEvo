#!/bin/bash
source /cvmfs/sft.cern.ch/lcg/external/gcc/6.2.0/x86_64-slc6-gcc62-opt/setup.sh
g++ -o calculator -O3 -Wall -Wextra -pedantic -march=native BeamParameters.cpp $(gsl-config --libs)
