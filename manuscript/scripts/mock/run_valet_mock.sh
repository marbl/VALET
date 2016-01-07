#!/bin/env bash
#PBS -j eo -e ${HOME}/logs/valet_mock_${LOGNAME}.log
#PBS -l mem=32GB,walltime=144:00:00
#PBS -q large
#PBS -N VALET_MOCK

## script for analyzing synthetic microbiome assemblies using VALET

source ${HOME}/.bashrc

MOCK_ROOT=/cbcb/project2-scratch/nolson/VALET/manuscript/results/mock
MOCK_ASM=${MOCK_ROOT}/sim_dat/assemblies/Assemble/out/
MOCK_SAMPLE=${MOCK_ROOT}/sim_dat/reads/mock-BV8-BC6-AO4-AB2
VALET=/cbcb/project2-scratch/nolson/VALET/VALET/src/py/valet.py
RESULTDIR=${MOCK_ROOT}/mock
assemblies="${MOCK_ASM}/metavelvet.45.asm.contig,${MOCK_ASM}/velvet.45.asm.contig,${MOCK_ASM}/soapdenovo2.45.asm.contig,${MOCK_ASM}/masurca.45.asm.contig"
ASM_NAME="MetaVelvet,Velvet,SOAPdenovo2,MaSuRCA"
time ${VALET} -a ${assemblies} -q -1  ${MOCK_SAMPLE}_1.fq -2 ${MOCK_SAMPLE}_2.fq -o ${RESULTDIR}_metamos --threads 32 --assembly-names ${ASM_NAME}
time ${VALET} -a ${assemblies} -q -1  ${MOCK_SAMPLE}_1.fq -2 ${MOCK_SAMPLE}_2.fq -o ${RESULTDIR}_metamos_nr --threads 32 --assembly-names ${ASM_NAME} --skip-reapr
time ${VALET} -a ${assemblies} -q -1  ${MOCK_SAMPLE}_1.fq -2 ${MOCK_SAMPLE}_2.fq -o ${RESULTDIR}_metamos_5k --threads 32 --assembly-names ${ASM_NAME} -z 5000
time ${VALET} -a ${assemblies} -q -1  ${MOCK_SAMPLE}_1.fq -2 ${MOCK_SAMPLE}_2.fq -o ${RESULTDIR}_metamos_nr_5k --threads 32 --assembly-names ${ASM_NAME} --skip-reapr -z 5000
