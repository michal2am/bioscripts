#!/bin/bash -l

## Generic script for job launching using SLURM, compatible with Tryton and Prometheus
## michaladammichalowski@gmail.com
## 14.10.16

## useful stuff - modules ##
## module avail ?
## module load ?

## usefull stuff - grants ##
## plg-show-grants
## plg-show-grant-details ?

## usefull stuff - slurm ##
## squeue
## sbatch ?
## squeue --job ?
## squeue -u ?
## scancel ?
## sinfo ?


## job name
#SBATCH -J ADFtestjob
## nodes number
#SBATCH -N 2
## cores number (24 per node)
#SBATCH -n 48
## task number
#SBATCH --ntasks-per-node=24
## walltime
#SBATCH --time=01:00:00
## grant
#SBATCH -A testgrant
## partition
#SBATCH -p plgrid-testing
## standard output
#SBATCH --output="adf.out"
## error
#SBATCH --error="adf.err"


## przejscie do katalogu z ktorego wywolany zostal sbatch
cd $SLURM_SUBMIT_DIR

srun /bin/hostname
module load apps/adf/2014.07
adf input.adf