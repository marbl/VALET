#/bin/sh

sample="mock_metamos"
mate_1="mock_4cp_3cp_2cp_1cp.20x_1.fastq"
mate_2="mock_4cp_3cp_2cp_1cp.20x_2.fastq"

metamos="/cbcb/project2-scratch/cmhill/metamos/metAMOS-1.5rc3/"

#${metamos}/initPipeline -q -1 ${mate_1} -2 ${mate_2} -d ${sample}/ -i 300:700 
${metamos}/runPipeline -v -a spades,soap2,idba-ud,metavelvet -k 55 -p 32 -d ${sample} -e FindORFS -n FunctionalAnnotation