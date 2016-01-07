#!/bin/env bash
#PBS -j eo -e ${HOME}/logs/valet_hmp_${LOGNAME}.log
#PBS -l mem=32GB,walltime=72:00:00
#PBS -q large
#PBS -N VALET_HMP

## script for analyzing HMP vaginal microbiome metamos microbiome assemblies using VALET

source ${HOME}/.bashrc

SERG_ASM=/cbcb/project2-scratch/sergek/metagenome/run2/Assemble/out
HMP_SAMPLE=/cbcb/project2-scratch/cmhill/hmp/samples/SRS014465
VALETDIR=/cbcb/project2-scratch/nolson/VALET/VALET/src/py
RESULTDIR=/cbcb/project2-scratch/nolson/VALET/manuscript/results/hmp
assemblies="${HMP_SAMPLE}/SRS014465.scaffolds.fa,${SERG_ASM}/idba-ud.55.asm.contig,${SERG_ASM}/metavelvet.55.asm.contig,${SERG_ASM}/velvet.55.asm.contig,${SERG_ASM}/spades.55.asm.contig,${SERG_ASM}/soapdenovo2.55.asm.contig"

time ${VALETDIR}/valet.py -a ${assemblies} -q -1  ${HMP_SAMPLE}/SRS014465.denovo_duplicates_marked.trimmed.1.fastq -2 ${HMP_SAMPLE}/SRS014465.denovo_duplicates_marked.trimmed.2.fastq -o ${RESULTDIR}/SRS014465_metamos --threads 32 --assembly-names HMP,IDBA-UD,MetaVelvet,Velvet,SPAdes,SOAPdenovo2
time ${VALETDIR}/valet.py -a ${assemblies} -q -1  ${HMP_SAMPLE}/SRS014465.denovo_duplicates_marked.trimmed.1.fastq -2 ${HMP_SAMPLE}/SRS014465.denovo_duplicates_marked.trimmed.2.fastq -o ${RESULTDIR}/SRS014465_metamos_nr --threads 32 --assembly-names HMP,IDBA-UD,MetaVelvet,Velvet,SPAdes,SOAPdenovo2 --skip-reapr
time ${VALETDIR}/valet.py -a ${assemblies} -q -1  ${HMP_SAMPLE}/SRS014465.denovo_duplicates_marked.trimmed.1.fastq -2 ${HMP_SAMPLE}/SRS014465.denovo_duplicates_marked.trimmed.2.fastq -o ${RESULTDIR}/SRS014465_metamos_5k --threads 32 --assembly-names HMP,IDBA-UD,MetaVelvet,Velvet,SPAdes,SOAPdenovo2 -z 5000
time ${VALETDIR}/valet.py -a ${assemblies} -q -1  ${HMP_SAMPLE}/SRS014465.denovo_duplicates_marked.trimmed.1.fastq -2 ${HMP_SAMPLE}/SRS014465.denovo_duplicates_marked.trimmed.2.fastq -o ${RESULTDIR}/SRS014465_metamos_nr_5k --threads 32 --assembly-names HMP,IDBA-UD,MetaVelvet,Velvet,SPAdes,SOAPdenovo2 --skip-reapr -z 5000

