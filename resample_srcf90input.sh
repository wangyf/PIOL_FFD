#!/bin/bash
bin=resample_srcinput
gfortran resample_srcinput.f90 -O3 -o $bin


echo -e "mr31\n64 64 1072\n1\n1 64 1 1 64 1 1 1072 1\n" | ./$bin