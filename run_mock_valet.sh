#! /bin/bash

assembly_dir="/cbcb/project2-scratch/cmhill/metagenomes/mock/mock_metamos/Assemble/out/"
valet_dir="/cbcb/project-scratch/cmhill/VALET/"

for assembly in idba-ud-scaf idba-ud soapdenovo2 metavelvet spades ; do
    echo $assembly

    for b in 100 ; do
    for w in 100 200 ; do
    for s in `seq 2 3`; do
        #/cbcb/project2-scratch/cmhill/metagenomes/mock/mock_metamos/Assemble/out/idba-ud.55.asm.contig 
        cmd="${valet_dir}/src/py/pipeline.py -a ${assembly_dir}/${assembly}.55.asm.contig -q -1 mock_4cp_3cp_2cp_1cp.20x_1.fastq -2 mock_4cp_3cp_2cp_1cp.20x_2.fastq -o results/${assembly}_w${w}_b${b}_z1000_s${s} -w ${w} -g 10 -l 0 -X 1000 -p 32 -b ${b} -z 1000 -n fr -s ${s}"
    
        #echo $cmd
        
        cmd+="; cp /cbcb/project2-scratch/cmhill/metagenomes/mock/mock_4cp_3cp_2cp_1cp_reference.fasta results/${assembly}_w${w}_b${b}_z1000_s${s}/reference.refseq.fna" 
        #echo $cmd
        #eval $cmd

        cmd+="; ../quast_pipeline.sh results/${assembly}_w${w}_b${b}_z1000_s${s}/filtered_assembly.fasta 1000"
        echo $cmd
        qs "$cmd"
        #eval $cmd
    done
    done
    done
done