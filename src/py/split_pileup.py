#!/usr/bin/env python

from optparse import OptionParser
import os
import re
import sys


def setup_options():
    parser = OptionParser()
    parser.add_option("-p", "--pileup-file", dest="pileup_filename", help="pileup file name", metavar="FILE")
    parser.add_option("-o", "--output", dest="output_location", help="output location for depth of coverage tool", default=None, metavar="FILE")
    parser.add_option("-c", "--chunks", dest="total_chunks", help="total chunks to split the pileup file.", type="int", default = 1)
    
    (options,args) = parser.parse_args()

    should_quit = False
    if options.pileup_filename == None:
        parser.error("You failed to provide the pileup file")
    if should_quit:
        parser.help()
        exit(-1)

    return (options,args)


def main():
   
    (options, args) = setup_options()

    # Count the number of lines in the file.
    line_count = 0
    for line in open(options.pileup_filename, 'r'):
        line_count += 1

    # Approximate lines per chunk.
    lines_per_chunk = line_count / options.total_chunks

    curr_chunk = 0
    prev_contig = None
    curr_lines_per_chunk = 0
    curr_chunk_file = open(options.pileup_filename + '.' + str(curr_chunk),'w')

    for line in open(options.pileup_filename, 'r'):
        tuple = line.strip().split()

        if prev_contig is None:
            prev_contig = tuple[0]

        # Only start a new file if we exceed the lines per contig AND we are encountering a new contig.
        if curr_lines_per_chunk > lines_per_chunk and prev_contig != tuple[0] and curr_chunk < options.total_chunks:
            curr_chunk += 1
            curr_chunk_file.close()
            curr_chunk_file = open(options.pileup_filename + '.' + str(curr_chunk),'w')
            curr_lines_per_chunk = 0


        curr_chunk_file.write(line)
        prev_contig = tuple[0]
        curr_lines_per_chunk += 1




if __name__ == '__main__':
    main()