#!/bin/env bash
#PBS -j eo -e ${HOME}/logs/metamos_mock_${LOGNAME}.log
#PBS -l mem=32GB,walltime=144:00:00
#PBS -q large
#PBS -N METAMOS_MOCK

## script for analyzing synthetic microbiome assemblies using VALET

source ${HOME}/.bashrc

root_dir=/cbcb/project2-scratch/nolson/VALET/manuscript/results/mock
read_root=${root_dir}/sim_dat/reads/mock-BV8-BC6-AO4-AB2_
metamos=/cbcb/project2-scratch/cmhill/metamos/metAMOS-1.5rc3/

## running metamos
time ${metamos}/initPipeline -q -1 ${read_root}1.fq -2 ${read_root}2.fq -d ${root_dir}/sim_dat/assemblies -i 100:500
time ${metamos}/runPipeline -v -a metavelvet,velvet,soap2,spades,masurca -X lap,n50 -k 45 -p 32 -d ${root_dir}/sim_dat/assemblies -n FunctionalAnnotation
