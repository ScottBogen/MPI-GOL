#!/bin/bash
find . -type f -name slurm\* -exec rm {} \;
mpicc -w -o golo gol.c
sbatch -N2 --time=00:02:00 send.sh

#while [ ! find . -maxdepth 1 -name "*slurm*"]
#do
sleep 1
#done
ls -l
#find . -name '*slurm*' -exec more {};
