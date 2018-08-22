#!/bin/bash -e

# On the MIRA
#qsub -A GMSeismicSim -t 60 --proccount 16384 --mode c16 -n 1024 -O m0_small ./sord-mO
#qsub  -A GMSeismicSim -t 35  --proccount 16384 --mode c16 -n 1024 -O PIOL-FFD ./PIOL-FFD

qsub script.sh
