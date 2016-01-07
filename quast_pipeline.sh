#! /bin/bash

METAREF_DIR="/cbcb/project2-scratch/vicky/backup/metaref"
QUAST_DIR="/cbcb/project-scratch/cmhill/tools/quast-2.3"

assembly=$1
assembly_dir=`dirname $1`

# First find the correct references using metaref.
cmd="${METAREF_DIR}/metaref.pl ${assembly} ${assembly_dir}/reference"
echo $cmd
#eval $cmd

# COPY THE REFERENCES TO THE ASSEMBLY DIR
cp /cbcb/project2-scratch/cmhill/metagenomes/omega_w200_z1000/reference.refseq.fna ${assembly_dir}/reference.refseq.fna

# Run QUAST using the provided references.
cmd="${QUAST_DIR}/metaquast.py --ambiguity-usage all --no-plots --min-contig $2  ${assembly} -R ${assembly_dir}/reference.refseq.fna -o ${assembly_dir}/quast/"
echo $cmd
eval $cmd

# Only count misassemblies on those contigs that could align.
cat ${assembly_dir}/quast/combined_quast_output/contigs_reports/alignments_filtered_assembly.tsv | awk '{for (i=2; i<NF; i+=1) {print $i} }' | sort | uniq > ${assembly_dir}/reference.alignments
grep -f ${assembly_dir}/reference.alignments ${assembly_dir}/summary.gff > ${assembly_dir}/reference.summary.gff  

cmd="/cbcb/project-scratch/cmhill/VALET/src/py/parse_quast.py -q  ${assembly_dir}/quast/combined_quast_output/contigs_reports/*.stdout -g ${assembly_dir}/reference.summary.gff -m ${assembly_dir}/missed_features > ${assembly_dir}/reference.summary.results"
echo $cmd
eval $cmd

grep -f ${assembly_dir}/reference.alignments ${assembly_dir}/suspicious.gff > ${assembly_dir}/reference.suspicious.gff      
cmd="/cbcb/project-scratch/cmhill/VALET/src/py/parse_quast.py -q  ${assembly_dir}/quast/combined_quast_output/contigs_reports/*.stdout -g ${assembly_dir}/reference.suspicious.gff -m ${assembly_dir}/missed_features > ${assembly_dir}/reference.suspicious.results"
echo $cmd
eval $cmd