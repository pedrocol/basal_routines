#!/bin/bash
 
#PBS -l ncpus=24
#PBS -l mem=190GB
#PBS -l jobfs=200GB
#PBS -q normal
#PBS -P e14
#PBS -l walltime=23:00:00
#PBS -l storage=gdata/a00+scratch/a00
#PBS -l wd


module load intel-compiler/2020.3.304 #cdftools

cd /scratch/e14/pc5520/ISF/BATHY
./batinterp.exe -f namelist_gebco_BedMachineAntarctica > gebco.log
