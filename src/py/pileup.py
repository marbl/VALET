#!/usr/bin/env python

from optparse import OptionParser
import os
import re
import sys


PROG_NAME = "PILEUP"


def setup_options():
    parser = OptionParser()
    parser.add_option("-s", "--sam", dest="sam_filename", help="SAM file name.", metavar="FILE")
    (options,args) = parser.parse_args()

    if options.sam_filename == None:
        parser.error("You failed to provide the SAM file")

    return (options,args)


def main():
    (options, args) = setup_options()

    sam_file = open(options.sam_filename, 'r')

    # Build dictionary of contig lengths.
    contig_lengths = {}
    pattern = re.compile('SN:(?P<contig>[\w_\|\.]+)\s*LN:(?P<length>\d+)')
    line = sam_file.readline()
    while line.startswith("@"):

        if line.startswith("@SQ"):
            matches = pattern.search(line)
            
            if len(matches.groups()) == 2:
                contig_lengths[matches.group('contig')] = int(matches.group('length'))

        line = sam_file.readline()

    # Read each SAM record and add the start and end points of each alignment.
    alignments = []
    while line:

        fields = line.strip().split()

        #Only record valid alignments.
        if fields[2] != "*":
            
            # Don't consider secondary alignments.
            if int(fields[1]) & 0x100 == 0:
                
                # Add the start position.
                alignments.append([fields[2], int(fields[3]), "START"])

                # Add the end position.
                alignments.append([fields[2], int(fields[3]) + len(fields[9]), "END"])

        line = sam_file.readline()


    # Sort alignments by contig and starting position
    alignments.sort(key = lambda alignment: (alignment[0], alignment[1]))

    """
    Example:

    contig00001     1       START
    contig00001     1       START
    contig00001     4       START
    contig00001     6       START
    contig00001     7       END
    """

    curr_contig = None
    curr_index = 0
    curr_length = -1

    start_index = 0
    end_index = 0

    for alignment in alignments:

        if curr_contig is None:
            curr_contig = alignment[0]
            curr_length = contig_lengths[curr_contig]

        # We have found a new contig, process the previous contig results.
        if alignment[0] != curr_contig:
            calc_pileup(curr_contig, curr_length, start_index, end_index, alignments)

            start_index = end_index
            curr_contig = alignment[0]
            curr_length = contig_lengths[curr_contig]

        end_index += 1

    if curr_contig is not None:
        calc_pileup(curr_contig, curr_length, start_index, end_index, alignments)


def calc_pileup(contig_name, contig_length, start_index, end_index, alignments):
    """ Print the per-bp coverage for a given contig. """

    curr_coverage = 0

    for i in xrange(1, contig_length + 1):

        while start_index < end_index and alignments[start_index][1] <= i:

            # If we encounter a START, curr_coverage increases by 1, otherwise decrease by 1.
            if alignments[start_index][2] == "START":
                curr_coverage += 1
                #print "INCREASING COVERAGE TO: " + str(curr_coverage)
            else:
                curr_coverage -= 1
                #print "DECREASING COVERAGE TO: " + str(curr_coverage)

            start_index += 1

        print contig_name + '\t' + str(i) + '\t' + str(curr_coverage)


if __name__ == '__main__':
    main()