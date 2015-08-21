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
                      help="Output directory", default="output/")
    parser.add_option("-w", "--window-size", dest="window_size",
                      help="Sliding window size when determining misassemblies.", default="201")
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
    parser.add_option('-b', "--ignore-ends", dest="ignore_end_distances", default=0, type=int,
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

    return (options, args)


def main():
    (options, args) = get_options()

    global COMMANDS_FILE
    COMMANDS_FILE = open(options.output_dir + '/commands', 'w')

    # For each assembly:
    assemblies = options.assembly_filenames.split(",")
    assembly_names = options.assembly_names.split(
        ",") if options.assembly_names else ['asm_']

    for counter, assembly in enumerate(assemblies):

        # Get the current assembly's name.
        assembly_name = assembly_names[counter % len(assembly_names)] \
                if options.assembly_names else 'asm_' + str(counter)

        # Create the assembly output directory.
        output_dir = options.output_dir + '/' + assembly_name
        ensure_dir(output_dir)

        step("PROCESSING ASSEMBLY: " + str(assembly_name) + ' (' + assembly + ')')

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
        #step("ALIGNING READS")
        #unaligned_dir = run_bowtie2(options, sam_output_location)
        #contig_lengths = get_contig_lengths(sam_output_location)

        # Run samtools.

        # If coverage file doesn't exist, calculate it.
        if options.coverage_file is None:
            step("CALCULATING CONTIG COVERAGE")
            #options.coverage_file = calculate_contig_coverage(options, pileup_file)
            #results(options.coverage_file)
            pileup_file = run_abundance_by_kmers(options, assembly, output_dir)
            results(options.coverage_file)
        contig_abundances = get_contig_abundances(options.coverage_file)

        # Calculate assembly probability.

        # If more thread, partition coverage file.

        # Run depth of coverage marker.

        # Run Breakpoint finder.

        # Run REAPR/mate-pair happiness line.

        # Generate summary files.

    # Generate comparison plots of all assemblies.


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


def run_abundance_by_kmers(options, assembly_filename, output_dir):
    """ Pileup based on k-mer abundances."""

    coverage_filename = output_dir + '/coverage/temp_kmer.cvg'
    ensure_dir(coverage_filename)
    coverage_file = open(coverage_filename, 'w')

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
    run(call_arr, coverage_file)
    #call(call_arr, stdout=coverage_file)

    return options.kmer_pileup_file


def run_bowtie2(options = None, output_sam = 'temp.sam'):
    """
    Run Bowtie2 with the given options and save the SAM file.
    """

    # Using bowtie2.
    # Create the bowtie2 index if it wasn't given as input.
    # if not assembly_index:
    if not os.path.exists(os.path.abspath(options.output_dir) + '/indexes'):
        os.makedirs(os.path.abspath(options.output_dir) + '/indexes')

    fd, index_path = mkstemp(prefix='temp_',\
            dir=(os.path.abspath(options.output_dir)   + '/indexes/'))
    try:
        os.mkdirs(os.path.dirname(index_path))
    except:
        pass

    fasta_file = options.fasta_file

    build_bowtie2_index(os.path.abspath(index_path), os.path.abspath(fasta_file))
    assembly_index = os.path.abspath(index_path)

    unaligned_dir = os.path.abspath(options.output_dir) + '/unaligned_reads/'
    ensure_dir(unaligned_dir)
    unaligned_file = unaligned_dir + 'unaligned.reads'

    #input_sam_file = output_sam_file
    read_type = " -f "
    if options.fastq_file:
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
    command = os.path.join(base_path, "bin/bowtie2-2.2.2/bowtie2 ") + bowtie2_args + " -S " + output_sam

    out_cmd( FNULL.name, FNULL.name,[command])


    #call(command.split())
    args = shlex.split(command)
    bowtie_proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=FNULL)
    bowtie_output, err = bowtie_proc.communicate()


    if bowtie2_unaligned_check_args != "":
        command = os.path.join(base_path, "bin/bowtie2-2.2.2/bowtie2 ") + bowtie2_unaligned_check_args + " -S " + output_sam + "_2.sam"
        out_cmd( FNULL.name,  FNULL.name, [command])
        args = shlex.split(command)
        bowtie_proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=FNULL)
        bowtie_output, err = bowtie_proc.communicate()

    return unaligned_dir



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


def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)
    assert os.path.exists(d)


if __name__ == '__main__':
    main()
