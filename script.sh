#!/bin/bash
# run on bluewater

#PBS -A baln 
#PBS -N mfinal
#PBS -l nodes=1024:ppn=32
#PBS -l walltime=4:00:00
#PBS -e stderr
#PBS -o stdout
#PBS -m abe

echo "$( date ): mfinal started" >> log

cd '/mnt/c/scratch/sciteam/wang9/NepalMomentrate/PIOL-FFD_src'

/usr/bin/time -p aprun -n 16384 -N 16 -d 2 ./PIOL-FFD

echo "$( date ): mfinal finished" >> log

