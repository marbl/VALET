#! /usr/bin/env python

#from __future__ import print_function
from collections import defaultdict
from multiprocessing import Process, Lock
from optparse import OptionParser
from sys import platform as _platform
import collections
import copy
import math
import os
import SeqIO
import sys
import time
import pprint
import random
import itertools
import operator
import timeit
import subprocess
import tempfile


JOIN_PARAM = ['-t', '\t']
SORT_PARAM = []
if _platform == "darwin":
    JOIN_PARAM = []
    SORT_PARAM = []

PROG_NAME = "ABUN_BY_KMERS"

CLEAN_UP_FILES = []

# http://stackoverflow.com/questions/19570800/reverse-complement-dna
revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A','N':'N','R':'N','M':'N','Y':'N','S':'N','W':'N','K':'N'}[B] for B in x][::-1])


def setup_options():
    parser = OptionParser()
    parser.add_option("-a", "--assembly", dest="assembly_filename", help="Assembly file", metavar="FILE")
    parser.add_option("-r", "--reads", dest="reads_filename", help="input reads", metavar="FILE")
    parser.add_option("-k", "--kmer", dest="kmer_size", help="kmer size", default = 8, type="int")
    parser.add_option("-t", "--threads", dest="threads", help="number of threads", default = 8, type="int")
    parser.add_option("-e", "--error_rate", dest="error_rate", help="error rate (0% by default)", type="float", default=1.0)
    parser.add_option("-m", "--max_iterations", dest="max_iter", help="max iterations (30 by default)", type="int", default=30)
    parser.add_option("-p", "--pileup", dest="pileup_file", help="write pileup to this file.", default=None)
    parser.add_option("-s", "--save-intermediate-files", dest="save_intermediate_files", help="save all the intermediate files.", action='store_true')

    (options,args) = parser.parse_args()

    should_quit = False
    if options.assembly_filename == None or options.reads_filename == None:
        parser.error("You failed to provide the mpileup file")
    if should_quit:
        parser.help()
        exit(-1)
    return (options,args)


def build_kmer_to_contig_index(assembly_filename, kmer_size):
    """
    Return a dictionary of mappings from {kmer: {contig: [location_1, location_2, ...]} }
    """

    # Kmer inverted index
    kmer_ii = defaultdict(list)

    # Contig lengths
    contig_lengths = defaultdict(int)

    pf = SeqIO.ParseFasta(assembly_filename)
    tuple = pf.getRecord()
    kmer_revcomp = None
    kmer = None
    while tuple is not None:

        for i in xrange(0, len(tuple[1]) - kmer_size + 1):
            kmer = tuple[1][i:i + kmer_size]
            kmer_revcomp = revcompl(kmer)
            if kmer < kmer_revcomp:
                kmer_ii[kmer].append((tuple[0], i+1))
            else:
                kmer_ii[kmer_revcomp].append((tuple[0], i+1))

        contig_lengths[tuple[0]] = len(tuple[1])
        tuple = pf.getRecord()

    return kmer_ii, contig_lengths


def accumulate(l):
    """
    http://stackoverflow.com/questions/2249036/grouping-python-tuple-list
    """

    it = itertools.groupby(l, operator.itemgetter(0))
    for key, subiter in it:
        yield key, list(item[1] for item in subiter)


def assign_read_kmers_to_contigs(reads_filename, kmer_ii, kmer_size):
    """
    Given a set of reads and k-mer length, assign k-mer counts to the contigs.
    """

    contig_counts = defaultdict(int)

    pf = SeqIO.ParseFastQ(reads_filename)
    #tuple = pf.getNextReadSeq()
    kmer = None
    contig = None
    unalignable_kmers = 0
    num_reads = 0
    sum = 0
    for tuple in pf:
        #while tuple is not None:

        # For each k-mer in the read...
        for i in xrange(0, len(tuple[1]) - kmer_size + 1):

            # ... find what contigs contain it.
            kmer = tuple[1][i:i + kmer_size]
            if kmer in kmer_ii:
                # and randomly assign the count to one of the items.
                contig = random.choice(kmer_ii[kmer])[0]
                contig_counts[contig] += 1

            elif revcompl(kmer) in kmer_ii:
                contig = random.choice(kmer_ii[revcompl(kmer)])[0]
                contig_counts[contig] += 1

            else:
                unalignable_kmers += 1

        sum += len(tuple[1])
        num_reads += 1

    #    tuple = pf.getNextReadSeq()

    #print 'Unalignable k-mers:\t' + str(unalignable_kmers)
    return contig_counts, sum/num_reads


def assign_read_kmers_to_contigs_iterative(reads_filename, kmer_ii, kmer_size, contig_abundances):
    """
    Given a set of reads and k-mer length, assign k-mer counts to the contigs based on their abundances.
    """

    contig_counts = defaultdict(int)

    pf = SeqIO.ParseFastQ(reads_filename)
    #tuple = pf.getNextReadSeq()
    kmer = None
    contig = None
    contigs_containing_kmer = []
    unalignable_kmers = 0
    total_abundance = 0
    num_reads = 0
    sum = 0
    for tuple in pf:
        #while tuple is not None:

        # For each k-mer in the read...
        for i in xrange(0, len(tuple[1]) - kmer_size + 1):

            # ... find what contigs contain it.
            kmer = tuple[1][i:i + kmer_size]
            if kmer in kmer_ii or revcompl(kmer) in kmer_ii:
                if kmer not in kmer_ii:
                     kmer = revcompl(kmer)

                # and randomly assign the count to one of the items.
                contigs_containing_kmer = accumulate(kmer_ii[kmer])


                #print kmer +'\t',
                contigs_containing_kmer = list(contigs_containing_kmer)
                #print contigs_containing_kmer

                # Calculate total abundance
                for contig in contigs_containing_kmer:
                    total_abundance += contig_abundances[contig[0]]

                # Choose
                choice = random.randint(1, total_abundance)

                curr_abundance = 0
                chosen_contig_tuple = None
                for contig in contigs_containing_kmer:
                    curr_abundance += contig_abundances[contig[0]]

                    # Have we found the right contig?
                    if curr_abundance >= choice:
                        chosen_contig_tuple = contig
                        #print 'Selecting:\t',
                        #print chosen_contig_tuple
                        break

                contig_counts[chosen_contig_tuple[0]] += 1

                total_abundance = 0

            else:
                unalignable_kmers += 1

        sum += len(tuple[1])
        num_reads += 1

    return contig_counts, sum/num_reads


def build_read_kmers_index(reads_filename, kmer_ii, kmer_size):
    """
    Return two dictionaries.  One contains the counts of ambiguous k-mers, while the other
    contains the number of unique k-mers that map to a given contig.
    """

    ambiguous_kmer_counts = defaultdict(int)
    contig_counts = defaultdict(int)

    pf = SeqIO.ParseFastQ(reads_filename)
    #tuple = pf.getNextReadSeq()
    kmer = None
    contig = None
    contigs_containing_kmer = []
    unalignable_kmers = 0
    #total_abundance = 0
    num_reads = 0
    sum = 0
    for tuple in pf:

        # For each k-mer in the read...
        for i in xrange(0, len(tuple[1]) - kmer_size + 1):

            # ... find what contigs contain it.
            kmer = tuple[1][i:i + kmer_size]
            if kmer in kmer_ii or revcompl(kmer) in kmer_ii:
                if kmer not in kmer_ii:
                    kmer = revcompl(kmer)

                # and randomly assign the count to one of the items.
                contigs_containing_kmer = accumulate(kmer_ii[kmer])

                #print kmer +'\t',
                contigs_containing_kmer = list(contigs_containing_kmer)

                if len(contigs_containing_kmer) > 1:
                    ambiguous_kmer_counts[kmer] += 1
                else:
                    contig_counts[contigs_containing_kmer[0][0]] += 1
            else:
                unalignable_kmers += 1

        if num_reads % 100000 == 0:
            sys.stderr.write('Processed reads:\t' + str(num_reads) + '\r')

        sum += len(tuple[1])
        num_reads += 1

    return ambiguous_kmer_counts, contig_counts, sum/num_reads


def assign_read_kmers_to_contigs_new(kmer_ii, ambiguous_kmer_counts, unambiguous_contig_counts, contig_abundances):
    """
    Assign ambiguous read k-mers based on contig averages counts.
    """

    contig_counts = copy.deepcopy(unambiguous_contig_counts)

    contig_location_tuples = []
    total_abundance = 0
    # Cycle through all ambiguous k-mers and assign them.
    for kmer in ambiguous_kmer_counts.keys():
        # and randomly assign the count to one of the items.
        contig_location_tuples = kmer_ii[kmer]

        #print 'Kmer:\t' + kmer
        #print 'Count:\t' + str(ambiguous_kmer_counts[kmer])
        #print 'Contig_locations:'
        #pprint.pprint(contig_location_tuples)

        # and randomly assign the count to one of the items.
        #contigs_containing_kmer = accumulate(kmer_ii[kmer])

        #print kmer +'\t',
        contigs_containing_kmer = list(accumulate(contig_location_tuples))
        #print contigs_containing_kmer

        # Calculate total abundance
        for contig in contigs_containing_kmer:
            total_abundance += contig_abundances[contig[0]]

        # Assign fractional counts based on total abundances.
        for contig in contigs_containing_kmer:
            #total_abundance += contig_abundances[contig[0]]

            #print 'Assigning\t' + str(contig_abundances[contig[0]] * ambiguous_kmer_counts[kmer] / total_abundance) + '\tto\t' + contig[0]

            contig_counts[contig[0]] += (contig_abundances[contig[0]] * ambiguous_kmer_counts[kmer] / total_abundance)

        total_abundance = 0

        #for i in xrange(0, ambiguous_kmer_counts[kmer]):
        #    contig = random.choice(contig_location_tuples)[0]
        #    #print "Selecting contig:\t" + contig
        #    contig_counts[contig] += 1

    return contig_counts


def output_inverted_index(kmer_ii, output_filename = 'tmp.ii', threads = 8):
    """
    Given an inverted index from k-mers => (contig, location) tuples, output a file:
    AAAAAA  contig1:1,2,3;contig2:3,4,5
    """

    unique_output_file = open(output_filename + '.unique', 'w')
    ambiguous_output_file = open(output_filename + '.ambig', 'w')

    CLEAN_UP_FILES.append(unique_output_file.name)
    CLEAN_UP_FILES.append(ambiguous_output_file.name)

    contigs_containing_kmer = None
    for kmer in kmer_ii.keys():

        contigs_containing_kmer = list(accumulate(kmer_ii[kmer]))

        if len(contigs_containing_kmer) == 1:
            unique_output_file.write(kmer + '\t')
            unique_output_file.write(contigs_containing_kmer[0][0] + ':' + ','.join([str(i) for i in contigs_containing_kmer[0][1]]))
            unique_output_file.write('\n')
        else:
            ambiguous_output_file.write(kmer + '\t')
            #for contig_and_location in contigs_containing_kmer:
            #    #print contig_and_location
            ambiguous_output_file.write('\t'.join([contig_and_location[0] + ':' + ','.join([str(i) for i in contig_and_location[1]]) for contig_and_location in contigs_containing_kmer]))
            #    #ambiguous_output_file.write(contig_and_location[0] + ':' + ','.join([str(i) for i in contig_and_location[1]]) + '\t')

            ambiguous_output_file.write('\n')

    unique_output_file.close()
    ambiguous_output_file.close()

    try:
        call_arr = ['sort', '-T', './', '--parallel=' + str(threads), output_filename + '.unique', '-o', output_filename + '.unique.sorted']
        subprocess.check_call(call_arr)
        call_arr = ['sort', '-T', './', '--parallel=' + str(threads), output_filename + '.ambig', '-o', output_filename + '.ambig.sorted']
        subprocess.check_call(call_arr)
        # Output a combined list of inverted kmer indexes.
        call_arr = ['sort', '-m', '-T', './', '--parallel=' + str(threads), output_filename + '.unique.sorted', output_filename + '.ambig.sorted', '-o', output_filename + '.sorted']
        subprocess.check_call(call_arr)

        CLEAN_UP_FILES.append(output_filename + '.unique.sorted')
        CLEAN_UP_FILES.append(output_filename + '.ambig.sorted')
        CLEAN_UP_FILES.append(output_filename + '.sorted')
    except:
        call_arr = ['sort', '-T', './', output_filename + '.unique', '-o', output_filename + '.unique.sorted']
        subprocess.call(call_arr)

        call_arr = ['sort', '-T', './', output_filename + '.ambig', '-o', output_filename + '.ambig.sorted']
        subprocess.call(call_arr)
        # Output a combined list of inverted kmer indexes.
        call_arr = ['sort', '-m', '-T', './', output_filename + '.unique.sorted', output_filename + '.ambig.sorted', '-o', output_filename + '.sorted']
        subprocess.call(call_arr)

        CLEAN_UP_FILES.append(output_filename + '.unique.sorted')
        CLEAN_UP_FILES.append(output_filename + '.ambig.sorted')
        CLEAN_UP_FILES.append(output_filename + '.sorted')


def join_kmer_lists(kmer_ii_filename = 'tmp.ii', query_kmer_counts_filename = 'tmp.jf.sorted'):
    """
    Join the two sorted k-mer lists to remove singleton/error kmers from the query.

    AAAAAAAAAAAAATA contig00002:3042; 29
    AAAAAAAAAAAATAA contig00002:3043; 29
    AAAAAAAAAAAGAAA contig00002:2606; 41
    AAAAAAAAAAATAAT contig00002:3044; 30
    AAAAAAAAAAGAAAA contig00002:2607; 39
    AAAAAAAAAATAATT contig00002:3045; 31
    """

    # Store the unambiguous k-mer counts.
    contig_counts = defaultdict(int)

    call_arr = ['join', '-t', '\t', kmer_ii_filename +'.unique.sorted', query_kmer_counts_filename]
    #call_arr.extend(JOIN_PARAM)#, '-t', '\t']
    #print ' '.join(call_arr)
    join_proc = subprocess.Popen(call_arr, stdout = subprocess.PIPE)

    for line in join_proc.stdout:
        # AAAAAAAAAAGAAAA contig00002:2607; 39
        #print line
        tuple = line.rstrip().split('\t')

        #print tuple
        contig = tuple[1].split(':')[0]
        counts = tuple[2]
        #print 'Contig:\t' + contig + '\tcounts:\t' + counts

        contig_counts[contig] += int(counts)


    # Store the ambiguous k-mer counts.
    ambiguous_kmer_counts = defaultdict(int)

    call_arr = ['join', '-t', '\t', kmer_ii_filename +'.ambig.sorted', query_kmer_counts_filename]
    #call_arr.extend(JOIN_PARAM)#, '-t', '\t']
    #print ' '.join(call_arr)
    join_proc = subprocess.Popen(call_arr, stdout = subprocess.PIPE)

    for line in join_proc.stdout:
        tuple = line.rstrip().split('\t')

        kmer = tuple[0]
        counts = tuple[len(tuple)-1]

        ambiguous_kmer_counts[kmer] += int(counts)


    # Store the total k-mer count
    call_arr = ['join', '-t', '\t', kmer_ii_filename + '.sorted', query_kmer_counts_filename]
    #call_arr.extend(JOIN_PARAM)#, '-t', '\t', ]
    output_file = open(kmer_ii_filename + '.counts.join', 'w')
    join_proc = subprocess.call(call_arr, stdout = output_file)
    #subprocess.call(call_arr)

    return contig_counts, ambiguous_kmer_counts


def run_jellyfish(reads_filename, kmer_length = 15, threads = 8, output = 'tmp.jf'):
    """
    Run jellyfish on the reads.
    """

    call_arr = ['jellyfish', 'count', '-C', '-m', str(kmer_length), '-s', "11740184", '-o', output + '.tmp', '-t', str(threads)]
    call_arr.extend(reads_filename.split(','))
    #jellyfish count -m 25 -s 11740184 -t 8 /cbcb/project2-scratch/cmhill/metagenomes/mock/mock_4cp_3cp_2cp_1cp.20x_all.fastq --text
    subprocess.call(call_arr)

    call_arr = ['jellyfish', 'dump', '-c', '-t', output + '.tmp', '-o', output]
    subprocess.call(call_arr)

    try:
        call_arr = ['sort', '-T', './', '--parallel=' + str(threads), output, '-o', output + '.sorted']
        #sort_process = subprocess.Popen(['sort', '--parallel=' + str(threads), '-o', output + '.sorted'], stdin = jf_dump.stdout)
        #jf_dump.communicate()[0]
        #jf_dump.wait()
        #sort_process.wait()

        subprocess.check_call(call_arr)
    except:
        call_arr = ['sort', '-T', './', output, '-o', output + '.sorted']
        subprocess.check_call(call_arr)

    #print "output jf:\t" + str(output) + '.sorted'


def generate_pileup(kmer_ii, contig_abundances, pileup_filename, kmer_ii_filename = 'tmp.ii'):
    """
    Given contig abundances, generate a pileup.
    """

    """ 1. Go through the list of k-mer counts.
    AAAAAAAAAAAAAAA gi|150002608|ref|NC_009614.1| Bacteroides vulgatus ATCC 8482 chromosome, complete genome:371467,371468,371469   gi|42779081|ref|NC_003909.8| Bacillus cereus ATCC 10987, complete genome:1602347,1602348,1602349,4270406,4270407,4270408        668
    AAAAAAAAAAAAAAC gi|150002608|ref|NC_009614.1| Bacteroides vulgatus ATCC 8482 chromosome, complete genome:371470 gi|42779081|ref|NC_003909.8| Bacillus cereus ATCC 10987, complete genome:4270409        164
    AAAAAAAAAAAAAAG gi|42779081|ref|NC_003909.8| Bacillus cereus ATCC 10987, complete genome:1602346        79
    AAAAAAAAAAAAAAT gi|42779081|ref|NC_003909.8| Bacillus cereus ATCC 10987, complete genome:2966197        78
    """

    pileup_file = open(pileup_filename + '.tmp', 'w')

    tuple = None
    kmer = None
    total_abundance = 0

    for line in open(kmer_ii_filename + '.counts.join', 'r'):
        # 2. Assign the counts to the location based on contig abundances.
        tuple = line.strip().split('\t')
        kmer = tuple[0]
        contigs_and_locations = tuple[1:len(tuple)-1]
        count = int(tuple[len(tuple)-1])

        # Make first pass through to get total abundance count.
        for contig_and_locations_tuple in contigs_and_locations:
            contig_and_locations = contig_and_locations_tuple.split(':')
            locations = contig_and_locations[1].split(',')
            # Abundance of the contig is = # of times k-mer occurs within contig
            # multiplied by the abundance of the contig.
            total_abundance += contig_abundances[contig_and_locations[0]] * len(locations)

        # Make a second pass to assign the counts.
        for contig_and_locations_tuple in contigs_and_locations:
            contig_and_locations = contig_and_locations_tuple.split(':')
            locations = contig_and_locations[1].split(',')
            frac_total_abundance = (1.0 * contig_abundances[contig_and_locations[0]] * len(locations)) / total_abundance
            per_location_abundance = frac_total_abundance * count / len(locations)

            # Assign the counts
            for location in locations:
                pileup_file.write(contig_and_locations[0] + '\t' + location + '\t' + 'N' + '\t' + str(int(round(per_location_abundance))) + '\n')

        total_abundance = 0

    pileup_file.close()

    # Sort the pileup file by starting location.
    call_arr = ['sort', '-k1,1', '-k2,2n', pileup_filename + '.tmp', '-o', pileup_filename + '.tmp2']
    subprocess.call(call_arr)

    call_arr = ['rm', pileup_filename + '.tmp']
    subprocess.call(call_arr)

    # Sometimes the pileup is missing locations, add 0 to those locations.
    prev_contig = None
    output_file = open(pileup_filename, 'w')

    for record in open(pileup_filename + '.tmp2', 'r'):
        # contig00001     1       N       10
        # contig00001     2       N       8
        fields = record.strip().split()

        if prev_contig != fields[0]:
            prev_contig = fields[0]
            prev_value = int(fields[1])

        else:
            curr_value = int(fields[1])
            while (prev_value + 1) < curr_value:
                output_file.write(fields[0] + '\t' + str(prev_value + 1) + '\tN\t0\n')
                prev_value += 1

            prev_value = curr_value

        output_file.write('\t'.join(fields) + '\n')

    output_file.close()

    call_arr = ['rm', pileup_filename + '.tmp2']
    subprocess.call(call_arr)


def main():
    (options, args) = setup_options()

    start_time = timeit.default_timer()
    # Construct an index of k-mers to (contig, locations) pairs.
    kmer_ii, contig_lengths = build_kmer_to_contig_index(options.assembly_filename, options.kmer_size)
    elapsed = timeit.default_timer() - start_time
    sys.stderr.write('K-mer index built in:\t' + str(elapsed) + '\n')

    # Create a temporary file to hold the index.
    tempfile.tempdir = './'
    _, ii_filename = tempfile.mkstemp()
    #ii_file.close()
    sys.stderr.write('K-mer inverted index index:\t' + str(ii_filename) + '\n')
    output_inverted_index(kmer_ii, output_filename = ii_filename)
    #pprint.pprint(dict(kmer_ii))

    # Once we have the k-mer index built, we need to build the k-mer counts of the read data sets.
    start_time = timeit.default_timer()
    run_jellyfish(options.reads_filename, options.kmer_size, threads = options.threads, output = ii_filename + '.jf')
    elapsed = timeit.default_timer() - start_time
    #ambiguous_kmer_counts, unambiguous_contig_counts, average_read_length = build_read_kmers_index(options.reads_filename, kmer_ii, options.kmer_size)
    sys.stderr.write('Query k-mer counts built in:\t' + str(elapsed) + '\n')

    unique_contig_counts, ambiguous_kmer_counts = join_kmer_lists(ii_filename, ii_filename + '.jf.sorted')
    #pprint.pprint(dict(contig_counts))
    #pprint.pprint(dict(ambiguous_kmer_counts))

    #return
    #pprint.pprint(unambiguous_contig_counts)

    # Assign he contig counts
    contig_counts = assign_read_kmers_to_contigs_new(kmer_ii, ambiguous_kmer_counts, unique_contig_counts, unique_contig_counts)

    # Once we have k-mer index built, start processing the reads and assign read k-mers to assembly k-mers.
    #contig_counts, average_read_length = assign_read_kmers_to_contigs(options.reads_filename, kmer_ii, options.kmer_size)
    #pprint.pprint(dict(contig_counts))

    # Set the initial contig abundances.
    contig_abundances = defaultdict(int)
    average_read_length = 101

    for key in contig_counts.keys():
        #print key
        #print key + '\t' + str(contig_counts[key]) + '\t' + str(contig_lengths[key]) + '\t' + str((float(contig_counts[key]) * average_read_length) / (contig_lengths[key] * (average_read_length - options.kmer_size + 1)))

        # Old way.
        #contig_abundances[key] = float((contig_counts[key]) * average_read_length) / (contig_lengths[key] * (average_read_length - options.kmer_size + 1))

        contig_abundances[key] = float((contig_counts[key])) / (contig_lengths[key] - options.kmer_size + 1)

        #contig_abundances[key] = float((contig_counts[key]) / contig_lengths[key])

    #for contig in contig_abundances.keys():
    #    print contig + '\t' + str(contig_abundances[contig])

    #return

    max_iter = options.max_iter
    for i in range(0, max_iter):
        #print '**** ITERATION ***'

        contig_counts = assign_read_kmers_to_contigs_new(kmer_ii, ambiguous_kmer_counts, unique_contig_counts, contig_abundances)

        #contig_counts, average_read_length = assign_read_kmers_to_contigs_iterative(options.reads_filename, kmer_ii, options.kmer_size, contig_abundances)

        for key in contig_counts.keys():
            #contig_abundances[key] = float((contig_counts[key]) * average_read_length) / (contig_lengths[key] * (average_read_length - options.kmer_size + 1))
            contig_abundances[key] = float((contig_counts[key])) / (contig_lengths[key] - options.kmer_size + 1)
            #contig_abundances[key] = float(contig_counts[key]) / (contig_lengths[key])

        if i == max_iter-1:
            if options.pileup_file:
                generate_pileup(kmer_ii, contig_abundances, options.pileup_file, ii_filename)

            for contig in contig_abundances.keys():
                print contig.split()[0] + '\t' + str(contig_abundances[contig] / math.pow(options.error_rate,options.kmer_size))

    #pprint.pprint(dict(contig_abundances))
    if not options.save_intermediate_files:
        for filename in CLEAN_UP_FILES:
            os.remove(filename)

if __name__ == '__main__':
    main()
