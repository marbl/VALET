#!/usr/bin/python
from __future__ import print_function
from subprocess import call
from optparse import OptionParser
from tempfile import mkstemp
import ConfigParser
import os
import random
import re
import shlex
import shutil
import subprocess
import os
import resource
import sys
from optparse import OptionParser


BASE_PATH = os.path.dirname(sys.argv[0])[:-len('src/py/')]
FILE_LIMIT = resource.getrlimit(resource.RLIMIT_NOFILE)[1]
FNULL = open('/dev/null', 'w')
COMMANDS_FILE = None


class BColors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'


def get_options():
    """ Sets the options.
    """
    parser = OptionParser()
    parser.add_option("-a", "--assembly-fasta", dest="assembly_filenames",
                      help="Candidate assembly files", metavar="FILE")
    parser.add_option("--assembly-names", dest="assembly_names",
                      help="Names for the different assemblies.")
    parser.add_option("-r", "--reads", dest="reads_filenames",
                      help="First Read File", metavar="FILE")
    parser.add_option("-1", "--1", dest="first_mates",
                      help="Fastq filenames separated by commas that contain the first mates.")
    parser.add_option("-2", "--2", dest="second_mates",
                      help="Fastq filenames separated by commas that contain the second mates.")
    parser.add_option("-c", "--coverage-file", dest="coverage_file",
                      help="Assembly created per-contig coverage file")
    parser.add_option("-i", "--config-file", dest="config_file",
                      help="Config file with preset parameters.")
    parser.add_option("-o", "--output-dir", dest="output_dir",
                      help="Output directory", default="output")
    parser.add_option("-w", "--window-size", dest="window_size",
                      help="Sliding window size when determining misassemblies.", default="501")
    parser.add_option("-q", "--fastq", dest="fastq_file",
                      default=False, action='store_true',
                      help="if set, input reads are fastq format (fasta by default).")
    parser.add_option("-p", "--threads", dest="threads",
                      help="Number of threads", default="8")
    parser.add_option("-I", "--minins", dest="min_insert_size",
                      help="Min insert sizes for mate pairs separated by commas.", default="0")
    parser.add_option("-X", "--maxins", dest="max_insert_size",
                      help="Max insert sizes for mate pairs separated by commas.", default="500")
    parser.add_option("-n", "--orientation", dest="orientation", default="fr",
                      help="Orientation of the mates.")
    parser.add_option("-m", "--mu", dest="mu", default="180",
                      help="average mate pair insert sizes.")
    parser.add_option("-t", "--sigma", dest="sigma", default="18",
                      help="standard deviation of mate pair insert sizes.")
    parser.add_option("-x", "--max-alignments", dest="max_alignments", default="10000",
                      help="bowtie2 parameter to set the max number of alignments.")
    parser.add_option("-e", "--email", dest="email",
                      help="Email to notify when job completes")
    parser.add_option("-g", "--min-coverage", dest="min_coverage", type="int", default=0,
                      help="Minimum average coverage to run misassembly detection.")
    parser.add_option("-l", "--coverage-multiplier", dest="coverage_multiplier", type=float, default=0.0,
                      help="When binning by coverage, the new high = high + high * multiplier")
    parser.add_option("-s", "--min-suspicious", dest="min_suspicious_regions", default=2, type=int,
                      help="Minimum number of overlapping flagged miassemblies to mark region as suspicious.")
    parser.add_option("-d", "--suspicious-window-size", dest="suspicious_window_size", default=2000, type=int,
                      help="Mark region as suspicious if multiple signatures occur within this window size.")
    parser.add_option('-z', "--min-contig-length", dest="min_contig_length", default=1000, type=int,
                      help="Ignore contigs smaller than this length.")
    parser.add_option('-b', "--ignore-ends", dest="ignore_end_distances", default=200, type=int,
                      help="Ignore flagged regions within b bps from the ends of the contigs.")
    parser.add_option('-k', "--breakpoint-bin", dest="breakpoints_bin", default="50", type=str,
                      help="Bin sized used to find breakpoints.")
    parser.add_option(
        '-f', "--orf-file", dest="orf_file", help="gff formatted file containing orfs")
    parser.add_option("--kmer", dest="kmer_length", help="kmer length used for abundance estimation",
                      default="15")

    (options, args) = parser.parse_args()

    # If a config file is given, read the arguments from there.
    # if options.config_file:
    #     config = ConfigParser.ConfigParser()
    #     config.read(options.config_file)

    should_err = False
    if not options.assembly_filenames:
        error("You need to provide a fasta file with -a")
        should_err = True
    if (not options.first_mates or not options.second_mates) and not options.reads_filenames:
        error("You need to provide reads")
        should_err = True
    if should_err:
        parser.print_help()
        exit(-1)

    ensure_dir(options.output_dir + '/')

    return (options, args)


def main():
    (options, args) = get_options()

    global COMMANDS_FILE
    COMMANDS_FILE = open(options.output_dir + '/commands', 'w')

    # For each assembly:
    assemblies = options.assembly_filenames.split(",")
    assembly_names = options.assembly_names.split(
        ",") if options.assembly_names else ['asm_']

    # Store the final assembly names.
    final_assembly_names = []

    for counter, assembly in enumerate(assemblies):

        # Get the current assembly's name.
        assembly_name = assembly_names[counter % len(assembly_names)] \
                if options.assembly_names else 'asm_' + str(counter)
        final_assembly_names.append(assembly_name)

        # Create the assembly output directory.
        output_dir = options.output_dir + '/' + assembly_name
        ensure_dir(output_dir)

        bold("PROCESSING ASSEMBLY: " + str(assembly_name) + ' (' + assembly + ')')

        results_filenames = []

        # Filter assembled contigs by length.
        if options.min_contig_length > 0:
            step("FILTERING ASSEMBLY CONTIGS LESS THAN " + str(options.min_contig_length) + ' BPs')
            filtered_filename = output_dir + '/filtered_assembly.fasta'
            ensure_dir(filtered_filename)
            all_contig_lengths = filter_short_contigs(assembly, options.min_contig_length, filtered_filename)
            results(filtered_filename)
            assembly = filtered_filename
            #input_fasta_saved = options.fasta_file

        # Align reads with Bowtie2
        step("ALIGNING READS")
        sam_output_location_dir = output_dir + "/sam/"
        sam_output_location = sam_output_location_dir + "library.sam"
        ensure_dir(sam_output_location)
        unaligned_dir = run_bowtie2(options, assembly, output_dir, sam_output_location)
        contig_lengths = get_contig_lengths(sam_output_location)

        # Run samtools.
        step("RUNNING SAMTOOLS")
        bam_location, sorted_bam_location, pileup_file = \
                run_samtools(options, assembly, output_dir, sam_output_location, index=True)

        # Run coverage estimation.
        step("CALCULATING CONTIG COVERAGE")
        options.coverage_file = calculate_contig_coverage(options, output_dir, pileup_file)
        #pileup_file = run_abundance_by_kmers(options, assembly, output_dir)
        results(options.coverage_file)

        contig_abundances = get_contig_abundances(options.coverage_file)

        # Calculate assembly probability.

        # If more thread, partition coverage file.
        if options.threads > 1:
            step("PARTITIONING COVERAGE FILE")
            run_split_pileup(options, pileup_file)

        # Run depth of coverage marker.
        step("DEPTH OF COVERAGE")
        results_filenames.append(run_depth_of_coverage(options, output_dir, pileup_file))

        # Run Breakpoint finder.
        outputBreakpointDir = output_dir + "/breakpoint/"
        ouputBreakpointLocation = outputBreakpointDir + "errorsDetected.csv"
        ensure_dir(outputBreakpointDir)

        step("BREAKPOINT")
        results_filenames.append(run_breakpoint_finder(options, assembly, unaligned_dir, outputBreakpointDir))

        # Run REAPR/mate-pair happiness line.

        # Generate summary files.
        step("SUMMARY")
        final_misassemblies = generate_summary_files(options, results_filenames, contig_lengths, output_dir)
        results(output_dir + "/summary.bed")

        generate_summary_table(output_dir + "/summary.tsv", all_contig_lengths, \
            contig_lengths, contig_abundances, final_misassemblies)
        results(output_dir + "/summary.tsv")

    # Generate comparison plots of all assemblies.
    bold("GENERATING ASSEMBLY COMPARISON PLOTS")
    generate_comparison_plot(options, final_assembly_names)
    results(options.output_dir + '/comparison_plots.pdf')

def filter_short_contigs(fasta_filename, min_contig_length, filtered_fasta_filename):
    """
    Filter out contigs less than a certain length.

    Args:
        fasta_filename:
        min_contig_length: Filter out contigs less than this length.
        filtered_fasta_filename: Filename where the filtered FASTA entries will be written too.

    Returns:
        A dictionary mapping contig names to lengths.
    """

    filtered_assembly_file = open(filtered_fasta_filename, 'w')
    all_contig_lengths = {}
    curr_length = 0
    with open(fasta_filename,'r') as assembly:
        for contig in contig_reader(assembly):
            curr_length = len(''.join(contig['sequence']))

            if curr_length >= min_contig_length:
                filtered_assembly_file.write(contig['name'])
                filtered_assembly_file.writelines(contig['sequence'])
                filtered_assembly_file.write('\n')

            all_contig_lengths[contig['name'].strip()[1:]] = curr_length

    filtered_assembly_file.close()

    return all_contig_lengths


def run_split_pileup(options, pileup_file):
    """ Split the pileup file into a number of chunks.

    Args:
        options: Command line options.
        pileup_file: Filename of the pileup file to split.
    """

    call_arr = [os.path.join(BASE_PATH, "src/py/split_pileup.py"), "-p", pileup_file, "-c", options.threads]
    run(call_arr)


def run_abundance_by_kmers(options, assembly_filename, output_dir):
    """ Pileup based on k-mer abundances.

    Args:
        options: Commandline options provided to VALET.
        assembly_filename: Current assembly FASTA filename.
        output_dir: Output directory of the current assembly.

    Returns:
        The filename of the resulting pileup.
    """

    coverage_filename = output_dir + '/coverage/temp_kmer.cvg'
    ensure_dir(coverage_filename)
    coverage_file = open(coverage_filename, 'w')
    abundance_error = open(output_dir + '/coverage/stderr', 'w')

    options.kmer_pileup_file = output_dir + "/coverage/kmer_pileup"
    options.coverage_file = output_dir + '/coverage/temp_kmer.cvg'

    # ./src/py/abundance_by_kmers.py -a test/test_kmer_abun.fna -r test/test_kmers_abun_lib.fastq -k 15 -t 4 -e .98 -p tmp_kmer_abun_15_30 -m 30
    call_arr = [os.path.join(BASE_PATH, "src/py/abundance_by_kmers.py"), \
            "-a", assembly_filename,\
            "-r", options.reads_filenames,\
            "-k", options.kmer_length,\
            "-t", options.threads,\
            "-e", ".98",
            "-p", options.kmer_pileup_file]
    run(call_arr, coverage_file, abundance_error)
    #call(call_arr, stdout=coverage_file)

    return options.kmer_pileup_file


def get_contig_lengths(sam_filename):
    """
    Return a dictionary of contig names => contig lengths from a SAM file.

    Args:
        sam_filename: input SAM filename.

    Returns:
        Dictionary of contig names => contig lengths.
    """

    sam_file = open(sam_filename, 'r')

    # Build dictionary of contig lengths.
    contig_lengths = {}
    pattern = re.compile('SN:(?P<contig>[\w_\|\.\-]+)\s*LN:(?P<length>\d+)')
    line = sam_file.readline()
    while line.startswith("@"):

        if line.startswith("@SQ"):
            matches = pattern.search(line)

            if len(matches.groups()) == 2:
                contig_lengths[matches.group('contig')] = int(matches.group('length'))

        line = sam_file.readline()

    return contig_lengths


def run_bowtie2(options, assembly_filename, output_dir, output_sam):
    """
    Run Bowtie2 with the given options and save the SAM file.
    """

    # Using bowtie2.
    # Create the bowtie2 index if it wasn't given as input.
    # if not assembly_index:
    if not os.path.exists(os.path.abspath(output_dir) + '/indexes'):
        os.makedirs(os.path.abspath(output_dir) + '/indexes')

    fd, index_path = mkstemp(prefix='temp_',\
            dir=(os.path.abspath(output_dir)   + '/indexes/'))
    try:
        os.mkdirs(os.path.dirname(index_path))
    except:
        pass

    fasta_file = assembly_filename

    build_bowtie2_index(os.path.abspath(index_path), os.path.abspath(fasta_file))
    assembly_index = os.path.abspath(index_path)

    unaligned_dir = os.path.abspath(output_dir) + '/unaligned_reads/'
    ensure_dir(unaligned_dir)
    unaligned_file = unaligned_dir + 'unaligned.reads'

    #input_sam_file = output_sam_file
    #read_type = " -f "
    #if options.fastq_file:
    read_type = " -q "

    bowtie2_args = ""
    bowtie2_unaligned_check_args = ""
    if options.first_mates:
        bowtie2_args = "-a -x " + assembly_index + " -1 " + options.first_mates\
                + " -2 " + options.second_mates + " -p " + options.threads\
                + " --very-sensitive -a " + " --reorder --"\
                + options.orientation + " -I " + options.min_insert_size\
                + " -X " + options.max_insert_size + " --no-mixed" #+ " --un-conc "\
                #+ unaligned_file

        bowtie2_unaligned_check_args = "-a -x " + assembly_index + read_type + " -U "\
                + options.first_mates + "," + options.second_mates + " --very-sensitive -a "\
                + " --reorder -p " + options.threads + " --un " + unaligned_file

    else:
        bowtie2_args = "-a -x " + assembly_index + read_type + " -U "\
                + options.reads_filenames + " --very-sensitive -a "\
                + " --reorder -p " + options.threads + " --un " + unaligned_file

    if not options:
        sys.stderr.write("[ERROR] No Bowtie2 options specified" + '\n')
        return

    # Using bowtie 2.
    command = "bowtie2 " + bowtie2_args + " -S " + output_sam
    run(shlex.split(command), stderr=FNULL)

    #out_cmd( FNULL.name, FNULL.name,[command])


    #args = shlex.split(command)
    #bowtie_proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=FNULL)
    #bowtie_output, err = bowtie_proc.communicate()


    if bowtie2_unaligned_check_args != "":
        command = "bowtie2 " + bowtie2_unaligned_check_args + " -S " + output_sam + "_2.sam"
        run(shlex.split(command), stderr=FNULL)

        #out_cmd( FNULL.name,  FNULL.name, [command])
        #args = shlex.split(command)
        #bowtie_proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=FNULL)
        #bowtie_output, err = bowtie_proc.communicate()

    return unaligned_dir


def build_bowtie2_index(index_filename, assembly_filename):
    """
    Build a Bowtie2 index.

    Args:
        index_filename: bowtie2 index name.
        assembly_filename: assembly filename.

    Returns:
        The path to the newly created bowtie2 index.
    """

    command = "bowtie2-build " + os.path.abspath(assembly_filename) + " " + os.path.abspath(index_filename)
    run(shlex.split(command), stdout=FNULL, stderr=FNULL)

    # Bad workaround.
    # out_cmd(FNULL.name, FNULL.name, [command])

    #bowtie2_build_proc = subprocess.Popen(command, shell = True, stdout = FNULL, stderr = FNULL)
    #bowtie_output, err = bowtie2_build_proc.communicate()
    #bowtie2_build_proc.wait()

    return index_filename


def calculate_contig_coverage(options, output_dir, pileup_file):
    """
    Calculate contig coverage.  The coverage of a contig is the mean per-bp coverage.

    Args:
        options: command line options.
        pileup_file: filename of the samtools formatted pileup file.
    Returns:
        Filename of the coverage file.
    """

    coverage_filename = output_dir + '/coverage/temp.cvg'
    coverage_file = open(coverage_filename, 'w')

    prev_contig = None
    curr_contig = None

    length = 0
    curr_coverage = 0

    for record in open(pileup_file, 'r'):
        fields = record.strip().split()

        if prev_contig != fields[0]:
            if prev_contig:
                coverage_file.write(prev_contig + '\t' + str(float(curr_coverage) / length) + '\n')

            prev_contig = fields[0]
            length = 0
            curr_coverage = 0

        curr_coverage += int(fields[3])
        length += 1
    if prev_contig:
        coverage_file.write(prev_contig + '\t' + str(float(curr_coverage) / length) + '\n')
    coverage_file.close()

    return coverage_filename


def run_samtools(options, assembly_filename, output_dir, sam_output_location, with_pileup = True, index=False):
    """ Takes a sam file and runs samtools to create bam, sorted bam, and mpileup.

    Args:
        options: Command line options.
        assembly_filename: FASTA filename of the assembly.
        output_dir: output directory for the current assembly.
        sam_output_location: SAM filename.
        with_pileup: Calculate pileup coverage with mpileup.
        index: Index the resulting BAM file.

    Returns:
        Tuple consisting of resulting BAM filename, sorted BAM filename, and pileup filename.
    """

    bam_dir = output_dir + "/bam/"
    ensure_dir(bam_dir)
    bam_location = bam_dir + "library.bam"
    sorted_bam_location = bam_dir + "sorted_library"
    bam_fp = open(bam_location, 'w+')
    error_file_location = bam_dir + "error.log"
    error_fp = open(error_file_location, 'w+')

    #warning("About to run samtools view to create bam")
    call_arr = ["samtools", "view", "-bS", sam_output_location]
    run(call_arr, bam_fp, error_fp)
    #out_cmd(bam_fp.name, error_fp.name, call_arr)
    #call(call_arr, stdout = bam_fp, stderr = error_fp)

    #warning("About to attempt to sort bam")
    call_arr = ["samtools", "sort", bam_location, sorted_bam_location]
    run(call_arr, stderr=FNULL)
    #out_cmd( "", FNULL.name, call_arr)
    #call(call_arr, stderr = FNULL)

    coverage_file_dir = output_dir + "/coverage/"
    ensure_dir(coverage_file_dir)
    pileup_file = coverage_file_dir + "mpileup_output.out"
    p_fp = open(pileup_file, 'w')

    if with_pileup:
        call_arr = ["samtools", "mpileup", "-C50", "-A", "-f", assembly_filename, sorted_bam_location + ".bam"]
        run(call_arr, p_fp, FNULL)
        #out_cmd(p_fp.name, FNULL.name, call_arr)
        results(pileup_file)
        #warning("That command outputs to file: ", pileup_file)
        #call(call_arr, stdout = p_fp, stderr = FNULL)

    if index:
        call_arr = ["samtools", "index", sorted_bam_location + ".bam"]
        run(call_arr, FNULL, FNULL)
        #out_cmd(FNULL.name, FNULL.name, call_arr)
        #call(call_arr, stdout = FNULL, stderr = FNULL)

    return (bam_location, sorted_bam_location, pileup_file)


def run_depth_of_coverage(options, output_dir, pileup_file):
    """ Run depth of coverage.

    Args:
        options: Command line options.
        output_dir: Assembly output directory.
        pileup_file: samtools pileup filename.
    Returns:
        Filename of flagged coverage regions.
    """

    coverage_error = open(output_dir + "/coverage/error.log", 'a')

    dp_fp = output_dir + "/coverage/errors_cov.bed"
    abundance_file = options.coverage_file
    #call_arr = ["src/py/depth_of_coverage.py", "-a", abundance_file, "-m", pileup_file, "-w", options.window_size, "-o", dp_fp, "-g", "-e"]
    call_arr = [os.path.join(BASE_PATH, "src/py/depth_of_coverage.py"), "-m", pileup_file, "-w", options.window_size, "-o", dp_fp, "-g", "-e", "-c", options.threads]
    run(call_arr, stdout=coverage_error, stderr=coverage_error)
    results(dp_fp)

    return dp_fp


def run_breakpoint_finder(options, assembly_filename, unaligned, breakpoint_dir):
    """Attempts to find breakpoints.

    Args:
        options: Command line options.
        assembly_filename: Assembly FASTA filename.
        unaligned: Unaligned sequences filename.
        breakpoint_dir: Output breakpoint directory.

    Returns:
        Filename of the breakpoints found.
    """

    std_err_file = open(breakpoint_dir + 'splitter_std_err.log', 'w')
    call_arr = [os.path.join(BASE_PATH,'src/py/breakpoint_splitter.py'),\
            '-u', unaligned,\
            '-o', breakpoint_dir + 'split_reads/']

    run(call_arr, stderr=std_err_file)
    #out_cmd( "", std_err_file.name, call_arr)
    #call(call_arr, stderr=std_err_file)
    std_err_file.close()

    std_err_file = open(breakpoint_dir + 'std_err.log','w')
    call_arr = [os.path.join(BASE_PATH, 'src/py/breakpoint_finder.py'),\
            '-a', assembly_filename,\
            '-r', breakpoint_dir + 'split_reads/',\
            '-b', options.breakpoints_bin, '-o', breakpoint_dir,\
            '-c', options.coverage_file,\
            '-p', options.threads]
    run(call_arr, stderr=std_err_file)
    #out_cmd( "", std_err_file.name,call_arr)
    #call(call_arr,stderr=std_err_file)
    results(breakpoint_dir + 'interesting_bins.bed')
    return breakpoint_dir + 'interesting_bins.bed'


def generate_summary_files(options, results_filenames, contig_lengths, output_dir):
    """ Generate the summary files from the individual error files.

    Args:
        options: Command line options.
        results_filenames: Filename of the flagged BED files.
        contig_lengths: Lengths of contigs.
        output_dir: Given assembly output directory.

    Returns:
        Final mis-assemblies in tuple format.
    """

    summary_file = open(output_dir + "/summary.bed", 'w')
    suspicious_file = open(output_dir + "/suspicious.bed", 'w')
    summary_table_file = open(output_dir + "/summary.tsv", 'w')
    #suspicious_table_file = open(options.output_dir + "/suspicious.tsv", 'w')

    misassemblies = []
    for results_file in results_filenames:
        if results_file:
            for line in open(results_file, 'r'):
                misassemblies.append(line.strip().split('\t'))

    # Sort misassemblies by start site.
    misassemblies.sort(key = lambda misassembly: (misassembly[0], int(misassembly[1]), int(misassembly[2])))
    final_misassemblies = []
    for misassembly in misassemblies:

        # Truncate starting/ending region if it is near the end of the contigs.
        if int(misassembly[1]) <= options.ignore_end_distances and \
            int(misassembly[2]) > options.ignore_end_distances:
          misassembly[1] = str(options.ignore_end_distances + 1)

        if int(misassembly[2]) >= (contig_lengths[misassembly[0]] - options.ignore_end_distances) and \
            int(misassembly[1]) < (contig_lengths[misassembly[0]] - options.ignore_end_distances):
          misassembly[2] = str(contig_lengths[misassembly[0]] - options.ignore_end_distances - 1)

        # Don't print a flagged region if it occurs near the ends of the contig.
        if int(misassembly[1]) > options.ignore_end_distances and \
                int(misassembly[2]) < (contig_lengths[misassembly[0]] - options.ignore_end_distances):
            summary_file.write('\t'.join(misassembly) + '\n')

            final_misassemblies.append(misassembly)

    summary_file.close()
    return final_misassemblies


def generate_summary_table(table_filename, all_contig_lengths, filtered_contig_lengths, contig_abundances, misassemblies, orf=False):
    """
    Output the misassemblies in a table format:

    contig_name  contig_length  low_cov  low_cov_bps  high_cov  high_cov_bps ...
    CONTIG1 12000   1   100 0   0 ...
    CONTIG2 100 NA  NA  NA ...
    """

    table_file = open(table_filename, 'w')

    if orf:
        table_file.write("contig_name\tcontig_length\tabundance\torf_low_cov\torf_low_cov_bps\torf_high_cov\torf_high_cov_bps\torf_reapr\torf_reapr_bps\torf_breakpoints\torf_breakpoints_bps\n")
    else:
        table_file.write("contig_name\tcontig_length\tabundance\tlow_cov\tlow_cov_bps\thigh_cov\thigh_cov_bps\treapr\treapr_bps\tbreakpoints\tbreakpoints_bps\n")

    prev_contig = None
    curr_contig = None

    # Misassembly signatures
    low_coverage = 0
    low_coverage_bps = 0
    high_coverage = 0
    high_coverage_bps = 0
    reapr = 0
    reapr_bps = 0
    breakpoints = 0
    breakpoints_bps = 0

    processed_contigs = set()

    for misassembly in misassemblies:
        """
        contig00001     REAPR   Read_orientation        88920   97033   .       .       .       Note=Warning: Bad read orientation;colour=1
        contig00001     REAPR   FCD     89074   90927   0.546142        .       .       Note=Error: FCD failure;colour=17
        contig00001     DEPTH_COV       low_coverage    90818   95238   29.500000       .       .       low=30.000000;high=70.000000;color=#7800ef
        """

        curr_contig = misassembly[0]

        if prev_contig is None:
            prev_contig = curr_contig

        if curr_contig != prev_contig:
            # Output previous contig stats.
            table_file.write(prev_contig + '\t' + str(filtered_contig_lengths[prev_contig]) + '\t' + str(contig_abundances[prev_contig]) + '\t' + \
                str(low_coverage) + '\t' + str(low_coverage_bps) + '\t' + str(high_coverage) + '\t' + \
                str(high_coverage_bps) + '\t' + str(reapr) + '\t' + str(reapr_bps) + '\t' + str(breakpoints) + '\t' + \
                str(breakpoints_bps) + '\n')

            processed_contigs.add(prev_contig)

            # Reset misassembly signature counts.
            low_coverage = 0
            low_coverage_bps = 0
            high_coverage = 0
            high_coverage_bps = 0
            reapr = 0
            reapr_bps = 0
            breakpoints = 0
            breakpoints_bps = 0

            prev_contig = curr_contig

        # Process the current contig misassembly.
        if misassembly[3] == 'REAPR':
            #if 'Warning' not in misassembly[8]:
            #    reapr += 1
            #    reapr_bps += (int(misassembly[4]) - int(misassembly[3]) + 1)
            reapr += 1
            reapr_bps += (int(misassembly[2]) - int(misassembly[1]) + 1)

        elif misassembly[3] == 'Low_coverage':
            low_coverage += 1
            low_coverage_bps += (int(misassembly[2]) - int(misassembly[1]) + 1)

        elif misassembly[3] == 'High_coverage':
            high_coverage += 1
            high_coverage_bps += (int(misassembly[2]) - int(misassembly[1]) + 1)

        elif misassembly[3] == 'Breakpoint_finder':
            breakpoints += 1
            breakpoints_bps += (int(misassembly[2]) - int(misassembly[1]) + 1)

        else:
            print("Unhandled error: " + misassembly[3])

    if prev_contig:
        # Output previous contig stats.
        table_file.write(prev_contig + '\t' + str(filtered_contig_lengths[prev_contig]) + '\t' + str(contig_abundances[prev_contig]) + '\t' + \
            str(low_coverage) + '\t' + str(low_coverage_bps) + '\t' + str(high_coverage) + '\t' + \
            str(high_coverage_bps) + '\t' + str(reapr) + '\t' + str(reapr_bps) + '\t' + str(breakpoints) + '\t' + \
            str(breakpoints_bps) + '\n')

        processed_contigs.add(prev_contig)

    # We need to add the remaining, error-free contigs.
    for contig in filtered_contig_lengths:
        if contig not in processed_contigs:
            table_file.write(contig + '\t' + str(filtered_contig_lengths[contig]) + '\t' + str(contig_abundances[contig]) + '\t' + \
                '0\t0\t0\t0\t0\t0\t0\t0\n')
            processed_contigs.add(contig)


    # Finally, add the contigs that were filtered out prior to evaluation.
    for contig in all_contig_lengths:
        if contig not in processed_contigs:
            table_file.write(contig + '\t' + str(all_contig_lengths[contig]) + '\t' + 'NA\t' + \
                'NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n')
            processed_contigs.add(contig)


def generate_comparison_plot(options, assembly_names):
    """ Generate a plot showing the cumulative assembly vs. cumulative error for all the assemblies.

    Args:
        assembly_names: Assembly names.
    """

    std_err_file = open(options.output_dir + '/error.log', 'a')

    assembly_summary_files = []
    for name in assembly_names:
        assembly_summary_files.append(options.output_dir + '/' + name + '/summary.tsv')

    call_arr = ["Rscript", os.path.join(BASE_PATH, "src/R/compare_assemblies.R"), ','.join(assembly_summary_files), ','.join(assembly_names), options.output_dir + '/comparison_plots']
    run(call_arr, stdout=std_err_file, stderr=std_err_file)


def get_contig_abundances(abundance_filename):
    """
    Return a dictionary of contig names => contig abundances from the '/coverage/temp.cvg'.

    Returns:
        A dictionary mapping from contig names to abundances.
    """

    abundance_file = open(abundance_filename, 'r')

    # Build a dictionary of contig abundances.
    contig_abundances = {}

    for line in abundance_file:
        contig_and_abundance = line.strip().split()
        contig_abundances[contig_and_abundance[0]] = int(round(float(contig_and_abundance[1])))

    return contig_abundances


def contig_reader(fasta_file):
    """ Iterator that takes a fasta file and returns entries.
    IMPORTANT NOTE: names are stripped at first space.
    """
    save_line = ""
    contig = {}
    in_contig = False
    for line in fasta_file:
        if line[0] == '>' and in_contig:
            save_line = line
            ret_contig = contig
            contig = {}
            contig['sequence'] = []
            contig['name'] = line.split()[0].strip() + "\n"
            yield ret_contig
        elif line[0] == '>':
            contig['name'] = line.split()[0].strip() + "\n"
            contig['sequence'] = []
            in_contig = True
        else:
            contig['sequence'].append(line.strip())
    yield contig


def run(command, stdout=sys.stdout, stderr=sys.stderr):
    """ Run the command.
    """

    if COMMANDS_FILE:
        if stdout != sys.stdout:
            stdout_sht = " 1>%s " % (stdout.name)
        else:
            stdout_sht = ""
        if stderr != sys.stderr:
            stderr_sht = " 2>%s " % (stderr.name)
        else:
            stderr_sht = ""
        COMMANDS_FILE.write(' '.join(command) + stdout_sht + stderr_sht + "\n")
        COMMANDS_FILE.flush()
    print(BColors.OKBLUE + "COMMAND:\t" + BColors.ENDC, ' '.join(command), file=sys.stderr)

    call(command, stdout=stdout, stderr=stderr)


def line(x):
    print ("-"*x, file=sys.stderr)


def hashmarks(x):
    print ("#"*x, file=sys.stderr)


def step(*objs):
    line(75)
    print(BColors.HEADER + "STEP:\t" + BColors.ENDC, *objs, file=sys.stderr)


def results(*objs):
    print(BColors.WARNING + "RESULTS:\t" +
          BColors.ENDC, *objs, file=sys.stderr)


def warning(*objs):
    print("INFO:\t", *objs, file=sys.stderr)


def error(*objs):
    print(BColors.WARNING + "ERROR:\t" + BColors.ENDC, *objs, file=sys.stderr)


def bold(objs):
    hashmarks(75)
    print(BColors.BOLD + BColors.HEADER + objs + BColors.ENDC, file=sys.stderr)
    hashmarks(75)


def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)
    assert os.path.exists(d)


if __name__ == '__main__':
    main()
