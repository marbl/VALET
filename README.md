#VALET
Pipeline for evaulating metagenomic assemblies.

## Prerequisites
VALET requires the following tools to function correctly.
* [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (Tested with version 2.2.4)
* [samtools](http://www.htslib.org) (Tested with version 1.2)
* [bedtools](http://bedtools.readthedocs.org/en/latest/) (Tested with version 2.24.0)
* Python2 (Tested with 2.7.9)
* [numpy](http://www.numpy.org) (Tested with version 1.9.2)
* [REAPR](https://www.sanger.ac.uk/resources/software/reapr/) (*OPTIONAL: REAPR usage can be disabled with --skip-reapr argument*)

Please note, if REAPR throws an error, then you need to install the following PERL libraries:
* File::Basename
* File::Copy
* File::Spec
* File::Spec::Link
* Getopt::Long
* List::Util

## Installing VALET
Once the repository has been cloned, to install the required tools run the command:
```
git clone https://github.com/cmhill/VALET.git
cd VALET/

# Let your shell know where to find the VALET pipeline.
export VALET=`pwd`/src/py/
```

## Running VALET
Included is a small set of test assemblies of *Candidatus Carsonella ruddii.* In addition to the reference genome, each assembly contains a mis-assembly: a duplication (test/c_rudii_dup.fna), a relocation (test/c_rudii_reloc.fna), and a duplication + relocation (test/c_rudii_reloc_dup.fna). Test the installation by running the following command:

```
$VALET/valet.py -a test/c_rudii_reference.fna,test/c_rudii_dup.fna,test/c_rudii_relocation.fna,test/c_rudii_reloc_dup.fna -1 test/lib1.1.fastq -2 test/lib1.2.fastq --assembly-names reference,duplication,relocation,reloc-dup
```
```
###########################################################################
PROCESSING ASSEMBLY: reference (test/c_rudii_reference.fna)
###########################################################################
---------------------------------------------------------------------------
STEP:	 FILTERING ASSEMBLY CONTIGS LESS THAN 1000 BPs
RESULTS:	 output/reference/filtered_assembly.fasta
---------------------------------------------------------------------------
STEP:	 ALIGNING READS
COMMAND:	 bowtie2-build /Users/cmhill/Work/VALET-repo/output/reference/filtered_assembly.fasta /Users/cmhill/Work/VALET-repo/output/reference/indexes/temp_JWCaAp
COMMAND:	 bowtie2 -a -x /Users/cmhill/Work/VALET-repo/output/reference/indexes/temp_JWCaAp -q -U test/lib1.1.fastq,test/lib1.2.fastq --very-sensitive -a --reorder -p 8 --un /Users/cmhill/Work/VALET-repo/output/reference/unaligned_reads/unaligned.reads -S output/reference/sam/library.sam
---------------------------------------------------------------------------
STEP:	 RUNNING SAMTOOLS
COMMAND:	 samtools view -F 0x100 -bS output/reference/sam/library.sam
COMMAND:	 samtools sort output/reference/bam/library.bam output/reference/bam/sorted_library
COMMAND:	 samtools mpileup -C50 -A -f output/reference/filtered_assembly.fasta output/reference/bam/sorted_library.bam
RESULTS:	 output/reference/coverage/mpileup_output.out
COMMAND:	 samtools index output/reference/bam/sorted_library.bam
---------------------------------------------------------------------------
STEP:	 CALCULATING CONTIG COVERAGE
RESULTS:	 output/reference/coverage/temp.cvg
---------------------------------------------------------------------------
STEP:	 PARTITIONING COVERAGE FILE
COMMAND:	 ./src/py/split_pileup.py -p output/reference/coverage/mpileup_output.out -c 8
---------------------------------------------------------------------------
STEP:	 DEPTH OF COVERAGE
COMMAND:	 ./src/py/depth_of_coverage.py -m output/reference/coverage/mpileup_output.out -w 351 -o output/reference/coverage/errors_cov.bed -g -e -c 8
COMMAND:	 bedtools sort -i output/reference/coverage/errors_cov.bed
RESULTS:	 output/reference/coverage.bed
---------------------------------------------------------------------------
STEP:	 BREAKPOINT
COMMAND:	 ./src/py/breakpoint_splitter.py -u /Users/cmhill/Work/VALET-repo/output/reference/unaligned_reads/ -o output/reference/breakpoint/split_reads/
COMMAND:	 ./src/py/breakpoint_finder.py -a output/reference/filtered_assembly.fasta -r output/reference/breakpoint/split_reads/ -b 50 -o output/reference/breakpoint/ -c output/reference/coverage/temp.cvg -p 8
COMMAND:	 bedtools sort -i output/reference/breakpoint/interesting_bins.bed
RESULTS:	 output/reference/breakpoint/../breakpoints.bed
---------------------------------------------------------------------------
STEP:	 SUMMARY
RESULTS:	 output/reference/summary.bed
RESULTS:	 output/reference/summary.tsv
###########################################################################
PROCESSING ASSEMBLY: duplicate (test/c_rudii_dup.fna)
###########################################################################
...
###########################################################################
GENERATING ASSEMBLY COMPARISON PLOTS
###########################################################################
COMMAND:	 Rscript ./src/R/compare_assemblies.R output/reference/summary.tsv,output/duplicate/summary.tsv,output/relocation/summary.tsv,output/reloc-dup/summary.tsv reference,duplicate,relocation,reloc-dup output/comparison_plots
RESULTS:	 output/comparison_plots.pdf

```

The flagged regions (potential misassemblies) are stored in two files **[OUTPUT_DIR]/[ASSEMBLY_NAME]/summary.bed** and **[OUTPUT_DIR]/[ASSEMBLY_NAME]/suspicious.bed**.
The flagged regions are stored in [BED format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) [GFF format](http://www.sanger.ac.uk/resources/software/gff/spec.html).  If multiple misassembly signatures overlap, their intersection is written to **suspicious.gff**.

```
relocref        69818   69868   Breakpoint_finder
relocref        97542   101062  Low_coverage    0       .       97543   101062  0,0,255
relocref        101082  112037  Low_coverage    0       .       101083  112037  0,0,255
relocref        112085  112756  Low_coverage    0       .       112086  112756  0,0,255
relocref        112898  117913  Low_coverage    0       .       112899  117913  0,0,255
relocref        118264  126647  Low_coverage    0       .       118265  126647  0,0,255
relocref        126737  129400  Low_coverage    0       .       126738  129400  0,0,255
relocref        129532  131981  Low_coverage    0       .       129533  131981  0,0,255
relocref        131732  131782  Breakpoint_finder
relocref        132049  139431  Low_coverage    0       .       132050  139431  0,0,255
relocref        139510  145882  Low_coverage    0       .       139511  145882  0,0,255
relocref        151010  151060  Breakpoint_finder
relocref        151061  151111  Breakpoint_finder
```

In addition, a breakdown of each contig's number of misassemblies is available in the **[OUTPUT_DIR]/[ASSEMBLY_NAME]/summary.tsv** file:

```
contig_name     contig_length   abundance       low_cov low_cov_bps     high_cov        high_cov_bps    reapr   reapr_bps       breakpoints     breakpoints_bps
relocref        196205  40      9       47419   0       0       0       0       4       204
```

Lastly, VALET produces a plot comparing how the different assemblies incur the different types of errors as they accumulate contigs. Contigs are first sorted by their abundance * length and then processed in decreasing order. A straight vertical line at x=0 indicates an assembly with no flagged regions.

![FRC plot](http://cbcb.umd.edu/~cmhill/files/frc_plot.png)


## Investigating potential misassemblies using IGV
VALET produces batch scripts to visualize the assemblies and flagged regions in IGV. Simply run:

```
/path/to/igv/igv.sh -b output/reloc-dup/IGV.batch
```

![IGV](http://cbcb.umd.edu/~cmhill/files/reloc_dup.png)


## Example usages

TODO


## Tutorial: Finding misassemblies in the Human Microbiome Project

Here we show how VALET can be used to find misassemblies in the Human Microbiome Project (http://www.hmpdacc.org/).

Before you continue, make sure the following tools are installed:
* **VALET**
* **git**
* GNU tools:
    * **wget** (can also use **curl**)
    * **tar**
* **Integrative Genomics Viewer** (http://www.broadinstitute.org/igv/home).  Any other genome browser that accepts BAM/SAM alignments and can overlay GFF files are acceptable.

### Downloading HMP sample SRS014465
```
mkdir samples
cd samples

# Download reads
wget ftp://public-ftp.hmpdacc.org/Illumina/vaginal_introitus/SRS014465.tar.bz2
tar xvjf SRS014465.tar.bz2
cd ..

# Download assembly
wget ftp://public-ftp.hmpdacc.org/HMASM/PGAs/vaginal_introitus/SRS014465.scaffolds.fa.bz2
tar xvjf SRS014465.scaffolds.fa.bz2

# Export sample directory to a path variable
export HMP_SAMPLE=`pwd`/SRS014465/

cd ..
```

### Running VALET

```
$VALET/pipeline.py -a $HMP_SAMPLE/SRS014465.scaffolds.fa \
    -1 $HMP_SAMPLE/SRS014465.denovo_duplicates_marked.trimmed.1.fastq \
    -2 $HMP_SAMPLE/SRS014465.denovo_duplicates_marked.trimmed.2.fastq \
    -o SRS014465_valet --threads 32
```

```
INFO:    Coverage file not provided, will create one.
---------------------------------------------------------------------------
STEP:    FILTERING ASSEMBLY CONTIGS LESS THAN 1000 BPs
RESULTS:         SRS014465_valet/filtered_assembly.fasta
---------------------------------------------------------------------------
STEP:    ALIGNING READS
...
```

Now you are free to explore the flagged regions!

## Options
```
Usage: valet.py [options]

Options:
  -h, --help            show this help message and exit
  -a FILE, --assembly-fasta=FILE
                        Candidate assembly files
  --assembly-names=ASSEMBLY_NAMES
                        Names for the different assemblies.
  -r FILE, --reads=FILE
                        First Read File
  -1 FIRST_MATES, --1=FIRST_MATES
                        Fastq filenames separated by commas that contain the
                        first mates.
  -2 SECOND_MATES, --2=SECOND_MATES
                        Fastq filenames separated by commas that contain the
                        second mates.
  -c COVERAGE_FILE, --coverage-file=COVERAGE_FILE
                        Assembly created per-contig coverage file
  -i CONFIG_FILE, --config-file=CONFIG_FILE
                        Config file with preset parameters.
  -o OUTPUT_DIR, --output-dir=OUTPUT_DIR
                        Output directory
  -w WINDOW_SIZE, --window-size=WINDOW_SIZE
                        Sliding window size when determining misassemblies.
  -q, --fastq           if set, input reads are fastq format (fasta by
                        default).
  -p THREADS, --threads=THREADS
                        Number of threads
  -I MIN_INSERT_SIZE, --minins=MIN_INSERT_SIZE
                        Min insert sizes for mate pairs separated by commas.
  -X MAX_INSERT_SIZE, --maxins=MAX_INSERT_SIZE
                        Max insert sizes for mate pairs separated by commas.
  -n ORIENTATION, --orientation=ORIENTATION
                        Orientation of the mates.
  -m MU, --mu=MU        average mate pair insert sizes.
  -t SIGMA, --sigma=SIGMA
                        standard deviation of mate pair insert sizes.
  -x MAX_ALIGNMENTS, --max-alignments=MAX_ALIGNMENTS
                        bowtie2 parameter to set the max number of alignments.
  -e EMAIL, --email=EMAIL
                        Email to notify when job completes
  -g MIN_COVERAGE, --min-coverage=MIN_COVERAGE
                        Minimum average coverage to run misassembly detection.
  -l COVERAGE_MULTIPLIER, --coverage-multiplier=COVERAGE_MULTIPLIER
                        When binning by coverage, the new high = high + high *
                        multiplier
  -s MIN_SUSPICIOUS_REGIONS, --min-suspicious=MIN_SUSPICIOUS_REGIONS
                        Minimum number of overlapping flagged miassemblies to
                        mark region as suspicious.
  -d SUSPICIOUS_WINDOW_SIZE, --suspicious-window-size=SUSPICIOUS_WINDOW_SIZE
                        Mark region as suspicious if multiple signatures occur
                        within this window size.
  -z MIN_CONTIG_LENGTH, --min-contig-length=MIN_CONTIG_LENGTH
                        Ignore contigs smaller than this length.
  -b IGNORE_END_DISTANCES, --ignore-ends=IGNORE_END_DISTANCES
                        Ignore flagged regions within b bps from the ends of
                        the contigs.
  -k BREAKPOINTS_BIN, --breakpoint-bin=BREAKPOINTS_BIN
                        Bin sized used to find breakpoints.
  -f ORF_FILE, --orf-file=ORF_FILE
                        gff formatted file containing orfs
  --kmer=KMER_LENGTH    kmer length used for abundance estimation
  --skip-reapr
```
