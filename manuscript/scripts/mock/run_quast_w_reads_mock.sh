#!/bin/env bash
#PBS -j eo -e ${HOME}/logs/quast_mock_${LOGNAME}.log
#PBS -l mem=64GB,walltime=144:00:00
#PBS -q large
#PBS -N QUAST_MOCK

## script for analyzing synthetic and mock microbiome assemblies using MetaQUAST with reads

source ${HOME}/.bashrc

QUAST_DIR=/cbcb/project2-scratch/nolson/VALET/quast-3.2
MOCK_ROOT=/cbcb/project2-scratch/nolson/VALET/manuscript/results/mock
MOCK_ASM=${MOCK_ROOT}/sim_dat/assemblies/Assemble/out/
MOCK_READS=${MOCK_ROOT}/sim_dat/reads/mock-BV8-BC6-AO4-AB2
MOCK_REF=${MOCK_ROOT}/sim_dat/ref/reference.refseq.fna
RESULTDIR=${MOCK_ROOT}/QUAST-w-reads

time ${QUAST_DIR}/metaquast.py --threads 32 --ambiguity-usage all --no-plots --min-contig 5000  \
    	${MOCK_ASM}/metavelvet.45.asm.contig \
	-1 ${MOCK_READS}_1.fq -2 ${MOCK_READS}_2.fq \
	-R ${MOCK_REF} -o ${RESULTDIR}-MetaVelvet/

time ${QUAST_DIR}/metaquast.py --threads 32 --ambiguity-usage all --no-plots --min-contig 5000  \
    	${MOCK_ASM}/velvet.45.asm.contig \
	-1 ${MOCK_READS}_1.fq -2 ${MOCK_READS}_2.fq \	
	 -R ${MOCK_REF} -o ${RESULTDIR}-Velvet/

time ${QUAST_DIR}/metaquast.py --threads 32 --ambiguity-usage all --no-plots --min-contig 5000  \
	${MOCK_ASM}/soapdenovo2.45.asm.contig \
	-1 ${MOCK_READS}_1.fq -2 ${MOCK_READS}_2.fq \
	-R ${MOCK_REF} -o ${RESULTDIR}-SOAPdenovo2/

time ${QUAST_DIR}/metaquast.py --threads 32 --ambiguity-usage all --no-plots --min-contig 5000  \
    	${MOCK_ASM}/masurca.45.asm.contig \
	-1 ${MOCK_READS}_1.fq -2 ${MOCK_READS}_2.fq \
	-R ${MOCK_REF} -o ${RESULTDIR}-MaSuRCA/
