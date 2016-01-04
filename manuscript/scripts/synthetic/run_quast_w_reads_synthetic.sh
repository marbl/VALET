#!/bin/env bash
#PBS -j eo -e ${HOME}/logs/quast_synthetic_w_reads_${LOGNAME}.log
#PBS -l mem=64GB,walltime=144:00:00
#PBS -q large
#PBS -N QUAST_WR_SYN

## script for analyzing synthetic and mock microbiome assemblies using MetaQUAST with reads

source ${HOME}/.bashrc

QUAST_DIR=/cbcb/project2-scratch/nolson/VALET/quast-3.2
SYN_ROOT=/cbcb/project2-scratch/nolson/VALET/manuscript/results/synthetic
SYN_ASM=/cbcb/project2-scratch/sergek/metagenome/run1/Assemble/out/
SYN_READS=${SYN_ROOT}/SRS606249_subsampled_data/SRR606249_
SYN_REF=${SYN_ROOT}/SRR606249.references.fasta
RESULTDIR=${SYN_ROOT}/QUAST-w-reads

time ${QUAST_DIR}/metaquast.py --threads 32 --ambiguity-usage all --no-plots --min-contig 5000  \
        ${SYN_ASM}/metavelvet.55.asm.contig \
        -1 ${SYN_READS}_1.downsample.fastq -2 ${SYN_READS}_2.downsample.fastq \
        -R ${SYN_REF} -o ${RESULTDIR}-MetaVelvet/

time ${QUAST_DIR}/metaquast.py --threads 32 --ambiguity-usage all --no-plots --min-contig 5000  \
        ${SYN_ASM}/velvet.55.asm.contig \
        -1 ${SYN_READS}_1.downsample.fastq -2 ${SYN_READS}_2.downsample.fastq \
         -R ${SYN_REF} -o ${RESULTDIR}-Velvet/

time ${QUAST_DIR}/metaquast.py --threads 32 --ambiguity-usage all --no-plots --min-contig 5000  \
        ${SYN_ASM}/soapdenovo2.55.asm.contig \
        -1 ${SYN_READS}_1.downsample.fastq -2 ${SYN_READS}_2.downsample.fastq \
        -R ${SYN_REF} -o ${RESULTDIR}-SOAPdenovo2/

time ${QUAST_DIR}/metaquast.py --threads 32 --ambiguity-usage all --no-plots --min-contig 5000  \
        ${SYN_ASM}/spades.55.asm.contig \
        -1 ${SYN_READS}_1.downsample.fastq -2 ${SYN_READS}_2.downsample.fastq \
        -R ${SYN_REF} -o ${RESULTDIR}-SPADES/
