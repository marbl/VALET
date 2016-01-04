#!/bin/env bash
#PBS -j eo -e ${HOME}/logs/valet_synthetic_${LOGNAME}.log
#PBS -l mem=64GB,walltime=144:00:00
#PBS -q large
#PBS -N VALET_SYN

## script for analyzing synthetic microbiome assemblies using VALET

source ${HOME}/.bashrc

SERG_ASM=/cbcb/project2-scratch/sergek/metagenome/run1/Assemble/out/
#SYN_SAMPLE=/cbcb/project2-scratch/cmhill/metagenomes/SRR606249
SYN_SAMPLE=/cbcb/project2-scratch/nolson/VALET/manuscript/results/synthetic/SRS606249_subsampled_data/SRR606249_
VALET=/cbcb/project2-scratch/nolson/VALET/VALET/src/py/valet.py
RESULTDIR=/cbcb/project2-scratch/nolson/VALET/manuscript/results/synthetic/SRS606249
assemblies="${SERG_ASM}/metavelvet.55.asm.contig,${SERG_ASM}/velvet.55.asm.contig,${SERG_ASM}/soapdenovo2.55.asm.contig,${SERG_ASM}/spades.55.asm.contig"
ASM_NAME="MetaVelvet,Velvet,SOAPdenovo2,SPADES"
#time ${VALET} -a ${assemblies} -q -1  ${SYN_SAMPLE}_1.downsample.fastq -2 ${SYN_SAMPLE}_2.downsample.fastq -o ${RESULTDIR}_metamos --threads 32 --assembly-names ${ASM_NAME}
#time ${VALET} -a ${assemblies} -q -1  ${SYN_SAMPLE}_1.downsample.fastq -2 ${SYN_SAMPLE}_2.downsample.fastq -o ${RESULTDIR}_metamos_nr --threads 32 --assembly-names ${ASM_NAME} --skip-reapr
time ${VALET} -a ${assemblies} -q -1  ${SYN_SAMPLE}_1.downsample.fastq -2 ${SYN_SAMPLE}_2.downsample.fastq -o ${RESULTDIR}_metamos_5k --threads 32 --assembly-names ${ASM_NAME} -z 5000
#time ${VALET} -a ${assemblies} -q -1  ${SYN_SAMPLE}_1.downsample.fastq -2 ${SYN_SAMPLE}_2.downsample.fastq -o ${RESULTDIR}_metamos_nr_5k --threads 32 --assembly-names ${ASM_NAME} --skip-reapr -z 5000
