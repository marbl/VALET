#!/bin/env bash
#PBS -j eo -e ${HOME}/logs/syn_comp_w_reads_5k_${LOGNAME}.log
#PBS -l mem=32GB,walltime=4:00:00
#PBS -q workstation
#PBS -N COMP-QWR-SYN

## script for comparing VALET and metaQUAST with reads results
source ${HOME}/.bashrc


VALET_DIR=/cbcb/project2-scratch/nolson/VALET/VALET/src/py
ROOT_DIR=/cbcb/project2-scratch/nolson/VALET/manuscript/results/synthetic

## Compare VALET and metaQUST results

for assembler_name in MetaVelvet Velvet SOAPdenovo2 SPADES;
do
    QUAST_RESULTS_DIR=${ROOT_DIR}/QUAST-w-reads-${assembler_name}/combined_reference/contigs_reports
    VALET_RESULTS_DIR=${ROOT_DIR}/SRS606249_metamos_5k/${assembler_name}
    
    ## Directory for results
    COMP_DIR=${ROOT_DIR}/COMP_QUAST-w-reads-5k-${assembler_name}
    mkdir $COMP_DIR

    ## Filter contigs aligning to reference based on metaQUAST results
    cat ${QUAST_RESULTS_DIR}/alignments_*.tsv | \
        awk '{for (i=2; i<NF; i+=1) {print $i} }' | sort | uniq > ${COMP_DIR}/reference.alignments   

    ## Filtering VALET results to only include metaQUAST aligned contigs
    grep -f ${COMP_DIR}/reference.alignments ${VALET_RESULTS_DIR}/summary.bed \
        > ${COMP_DIR}/valet_reference.summary.bed
    grep -f ${COMP_DIR}/reference.alignments ${VALET_RESULTS_DIR}/suspicious.bed \
        > ${COMP_DIR}/valet_reference.suspicious.bed  

    ## Generating comparison results file
    ${VALET_DIR}/parse_quast.py -q  ${QUAST_RESULTS_DIR}/contigs_report*.stdout \
        -b ${COMP_DIR}/valet_reference.summary.bed \
        -m ${COMP_DIR}/quast_only_summary_features.bed \
            > ${COMP_DIR}/comparison_reference.summary.results

    ${VALET_DIR}/parse_quast.py -q  ${QUAST_RESULTS_DIR}/contigs_report*.stdout \
        -b ${COMP_DIR}/valet_reference.suspicious.bed \
        -m ${COMP_DIR}/quast_only_suspicious_features.bed \
            > ${COMP_DIR}/comparison_reference.suspicious.results
done
