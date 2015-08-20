#!/usr/bin/env python

#from __future__ import print_function
from multiprocessing import Process, Lock
from optparse import OptionParser
from collections import defaultdict
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


PROG_NAME = "ABUN_BY_KMERS"

# http://stackoverflow.com/questions/19570800/reverse-complement-dna
revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A','N':'N'}[B] for B in x][::-1])

def setup_options():
    parser = OptionParser()
    parser.add_option("-a", "--assembly", dest="assembly_filename", help="Assembly file", metavar="FILE")
    parser.add_option("-r", "--reads", dest="reads_filename", help="input reads", metavar="FILE")
    parser.add_option("-k", "--kmer", dest="kmer_size", help="kmer size", default = 8, type="int")

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
    while tuple is not None:

        for i in xrange(0, len(tuple[1]) - kmer_size + 1):
            kmer_ii[tuple[1][i:i + kmer_size]].append((tuple[0], i))

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


def main():
    (options, args) = setup_options()

    # Construct an index of k-mers to (contig, locations) pairs.
    kmer_ii, contig_lengths = build_kmer_to_contig_index(options.assembly_filename, options.kmer_size)
    #pprint.pprint(dict(kmer_ii))

    # Once we have k-mer index built, start processing the reads and assign read k-mers to assembly k-mers.
    contig_counts, average_read_length = assign_read_kmers_to_contigs(options.reads_filename, kmer_ii, options.kmer_size)
    #pprint.pprint(dict(contig_counts))

    # Set the initial contig abundances.
    contig_abundances = defaultdict(int)

    for key in contig_counts.keys():
        #print key
        #print key + '\t' + str(contig_counts[key]) + '\t' + str(contig_lengths[key]) + '\t' + str((float(contig_counts[key]) * average_read_length) / (contig_lengths[key] * (average_read_length - options.kmer_size + 1)))
        contig_abundances[key] = round(float((contig_counts[key]) * average_read_length) / (contig_lengths[key] * (average_read_length - options.kmer_size + 1)))

    for contig in contig_abundances.keys():
        print contig + '\t' + str(contig_abundances[contig])

    for i in range(0,10):
        print '**** ITERATION ***'

        contig_counts, average_read_length = assign_read_kmers_to_contigs_iterative(options.reads_filename, kmer_ii, options.kmer_size, contig_abundances)

        for key in contig_counts.keys():
            contig_abundances[key] = round(float((contig_counts[key]) * average_read_length) / (contig_lengths[key] * (average_read_length - options.kmer_size + 1)))

        for contig in contig_abundances.keys():
            print contig + '\t' + str(contig_abundances[contig])

    #pprint.pprint(dict(contig_abundances))

if __name__ == '__main__':
    main()