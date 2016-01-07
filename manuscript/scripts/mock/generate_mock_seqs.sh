#!/usr/bin/bash
## Generating similated datasets for VALET supplemental mock community
## 
##  --Genomes--
## Bacterioides vulgaus ATCC8482 (80X) NC_009614.1
## Bacillus cereus ATCC 10987 (60X) NC_003909.8
## Actinomyces odontolyticus ATCC 17982 (40X) NZ_DS264586.1
## Actinobacter baumanni ATCC 17978 (20X) NC_009085.1 
## 
## --Simulated Reads--
## wgsim default parameters
## number of reads determine based on the wc -c for reference genome files/200*coverage

SIMDIR=/cbcb/project2-scratch/nolson/VALET/manuscript/results/mock
GENOMEDIR=${SIMDIR}/sim_dat/ref
READSDIR=${SIMDIR}/sim_dat/reads

wgsim ${GENOMEDIR}/NC_009614.1.fa -N 2000000 -1 100 -2 100\
	${READSDIR}/NC_009614.1_1.fq ${READSDIR}/NC_009614.1_2.fq > ${READSDIR}/NC_009614.1_wgsim.out

wgsim ${GENOMEDIR}/NC_003909.8.fa -N 1600000 -1 100 -2 100\
        ${READSDIR}/NC_003909.8_1.fq ${READSDIR}/NC_003909.8_2.fq > ${READSDIR}/NC_003909.8_wgsim.out

wgsim ${GENOMEDIR}/NZ_DS264586.1.fa -N 480000 -1 100 -2 100\
        ${READSDIR}/NZ_DS264586.1_1.fq ${READSDIR}/NZ_DS264586.1_2.fq > ${READSDIR}/NZ_DS264586.1_wgsim.out

wgsim ${GENOMEDIR}/NC_009085.1.fa -N 400000 -1 100 -2 100\
        ${READSDIR}/NC_009085.1_1.fq ${READSDIR}/NC_009085.1_2.fq > ${READSDIR}/NC_009085.1_wgsim.out

cat ${READSDIR}/{NC_009614.1_1.fq,NC_003909.8_1.fq,NZ_DS264586.1_1.fq,NC_009085.1_1.fq} > ${READSDIR}/mock-BV8-BC6-AO4-AB2_1.fq
cat ${READSDIR}/{NC_009614.1_2.fq,NC_003909.8_2.fq,NZ_DS264586.1_2.fq,NC_009085.1_2.fq} > ${READSDIR}/mock-BV8-BC6-AO4-AB2_2.fq
