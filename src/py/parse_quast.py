#!/usr/bin/env python

from collections import defaultdict
from collections import Counter
from optparse import OptionParser
import os
import sys
import re

def setup_options():
    parser = OptionParser()
    parser.add_option("-q", "--quast", dest="quast_filename", help="QUAST stdout filename.", metavar="FILE")
    parser.add_option("-g", "--gff", dest="gff_filename", help="GFF filename to compare quast results to.", metavar="FILE")
    parser.add_option("-m", "--missed", dest="missed_misassemblies", help="Write out missed misassemblies here", metavar="FILE", default=None)
    parser.add_option("-o", "--offset", dest="offset", type="int", default=0)
    (options,args) = parser.parse_args()

    if options.quast_filename == None or options.gff_filename == None:
        parser.error("You failed to provide the QUAST or GFF file.")

    return (options,args)


def get_misassemblies_from_quast(quast_filename):
    """
    Return misassemblies from quast output.
    """
    quast_file = open(quast_filename, 'r')

    """
    CONTIG: NODE_3_length_346423_cov_30.6391_ID_5 (346423bp)
        Top Length: 248697  Top ID: 99.97
                This contig is misassembled. 3 total aligns.
                        Real Alignment 1: 440492 477785 | 37264 1 | 37294 37264 | 99.92 | mock_Actinomyces_odontolyticus_ATCC_17982 NODE_3_length_346423_cov_30.6391_ID_5
                          Overlap between these two alignments (local misassembly). Inconsistency = -30 
                        Real Alignment 2: 211697 460461 | 266021 17325 | 248765 248697 | 99.97 | mock_Actinomyces_odontolyticus_ATCC_17982 NODE_3_length_346423_cov_30.6391_ID_5
                          Gap between these two alignments (local misassembly). Inconsistency = 0 
                        Real Alignment 3: 131295 221709 | 346423 256009 | 90415 90415 | 100.0 | mock_Actinomyces_odontolyticus_ATCC_17982 NODE_3_length_346423_cov_30.6391_ID_5
    """

    misassemblies = []
    line = quast_file.readline()

    while line:

        line = line.strip()

        if line.startswith('This contig is misassembled.'):

            line = quast_file.readline()

            prev_alignment = None
            found_misassembly = False
            type_misassembly = 'local'

            while line.isspace() is not True:
                line = line.strip()
                #print line.strip()
                
                if line.startswith('Real Alignment'):
                    if prev_alignment is None:
                        prev_alignment = line

                    if found_misassembly:
                        prev_tuple = prev_alignment.split()
                        curr_tuple = line.split()
                        #print prev_tuple
                        max_prev = max(int(prev_tuple[6]), int(prev_tuple[7]))
                        min_curr = min(int(curr_tuple[6]), int(curr_tuple[7]))

                        misassemblies.append([type_misassembly, curr_tuple[len(curr_tuple) -1], min(max_prev, min_curr), max(max_prev, min_curr)])
                        
                        found_misassembly = False
                        type_misassembly = 'local'

                if 'misassembly' in line and not 'Fake' in line:
                    found_misassembly = True
                    if 'Extensive' in line:
                        type_misassembly = 'extensive'
                else:
                    prev_alignment = line

                #if line:
                line = quast_file.readline()
                
        line = quast_file.readline()

    return misassemblies


def get_misassemblies_from_gff(gff_filename):

    misassemblies = []

    for line in open(gff_filename, 'r'):
        # NODE_100_length_16680_cov_16.76_ID_179  DEPTH_COV       Low_coverage    9253    9561    36.000000       .       .       low=36.500000;high=80.500000;color=#7800ef
        tuple = line.strip().split()

        misassemblies.append([tuple[2], tuple[0], int(tuple[3]), int(tuple[4])])

    return misassemblies


def overlap(a,b):
    return a[0] <= b[0] <= a[1] or b[0] <= a[0] <= b[1]

def getOverlap(a, b):
    return max(-1, min(a[1], b[1]) - max(a[0], b[0]))


def get_start_end_interval(ref, query, i = 0, j = -1):
    """
    Return the start and end points in the query where the ref contig matches.
    """
    #print 'start find interval'
    #print ref 
    #print i,j
    if i >= len(query):
        return i, j
    #print query[i][1] + '\tlen: '  + str(len(query))
    
    #print ref + '\t' + query[i][1]

    # If the current reference name matches the current query start, don't change anything!
    if ref == query[i][1]:

        # If the end point hasn't been set before, set it now.
        if j == -1:
            # While the ref contig equals query, increment the end position
            j = 0
            while ref >= query[j][1]:
                j += 1
                #print ref + '\t' + query[j][1]

        while j < len(query) and ref >= query[j][1]:
            #print ref + '\tvs.\t' + query[j][1]
            #print ref + '\t' + query[j][1]
            j += 1

        #print 'Not changing anything!'
        return i, j

    #i = curr_start
    #print "adjusting i"
    while i < len(query) and ref > query[i][1]:
        #print ref 
        #print query[i][1]
        #print ref + '\t' + query[i][1]
        i += 1

    if i == len(query):
        return i, i

    j = i + 1
    #print "adjusting j"
    #print 'UH: ' + query[j][1]
    # While the ref contig equals query, increment the end position
    while j < len(query) and ref >= query[j][1]:
        #print ref + '\tvs.\t' + query[j][1]
        #print ref + '\t' + query[j][1]
        j += 1

    #print 'returning ',
    #print i, j
    return i, j


def main():

    (options, args) = setup_options()

    ref_misassemblies = get_misassemblies_from_quast(options.quast_filename)
    ref_misassemblies.sort(key = lambda misassembly: (misassembly[1], int(misassembly[2])))
    #print '\n'.join(str(i) for i in ref_misassemblies)

    gff_misassemblies = get_misassemblies_from_gff(options.gff_filename)
    #print '\n'.join(str(i) for i in gff_misassemblies)
    gff_misassemblies.sort(key = lambda misassembly: (misassembly[1], int(misassembly[2])))

    ref_misassembly_count = 0
    extensive_misassembly_count = 0
    local_misassembly_count = 0
    query_misassembly_count = 0

    missed_misassemblies = None
    
    if options.missed_misassemblies:
        missed_misassemblies = open(options.missed_misassemblies, 'w')

    prev_contig = ref_misassemblies[1]
    start, end = 0, -1 #get_start_end_interval(prev_contig, gff_misassemblies, 0, -1)
    #print ref_misassembly
    #print 'start, end:\t' + str(start) + '\t' + str(end)


    gff_misassembly = None
    prev_query_misassembly = None

    for ref_misassembly in ref_misassemblies:
        found_misassembly = False

        start, end = get_start_end_interval(ref_misassembly[1], gff_misassemblies, start, end)
        #print ref_misassembly
        #print 'start, end:\t' + str(start) + '\t' + str(end)

        for index in xrange(start, end):# len(gff_misassemblies)):
            gff_misassembly = gff_misassemblies[index]
            #for gff_misassembly in gff_misassemblies:
            #print gff_misassembly
            #print 'Comparing ',
            #print ref_misassembly,
            #print '\tWITH\t',
            #print gff_misassembly

            if ref_misassembly[1] == gff_misassembly[1]:
                if getOverlap([ref_misassembly[2] - options.offset, ref_misassembly[3] + options.offset], \
                        [gff_misassembly[2], gff_misassembly[3]]) >= 0:
                    print ref_misassembly
                    print gff_misassembly
                    print '****'

                    if (gff_misassembly[1]  + '_' + str(gff_misassembly[2]) + '_' + str(gff_misassembly[3])) != prev_query_misassembly:
                        query_misassembly_count += 1
                    found_misassembly = True

                    prev_query_misassembly = gff_misassembly[1]  + '_' + str(gff_misassembly[2]) + '_' + str(gff_misassembly[3])


        if found_misassembly:
            ref_misassembly_count += 1
            if ref_misassembly[0] == 'extensive':
                extensive_misassembly_count += 1
            else:
                local_misassembly_count += 1


        else:
            if options.missed_misassemblies:
                missed_misassemblies.write('\t'.join(map(str,ref_misassembly)) + '\n')


    counter = Counter(elem[0] for elem in ref_misassemblies)
    print 'Ref misassemblies found:\t' + str(int(ref_misassembly_count)) + '\t' + str(len(ref_misassemblies)) + '\t' + str(float(ref_misassembly_count)/len(ref_misassemblies))
    print 'Extensive misassemblies found:\t' + str(float(extensive_misassembly_count)) + '\t' + str(counter['extensive']) + '\t' + str(float(extensive_misassembly_count)/counter['extensive'])
    print 'Local missed_misassemblies found:\t' + str(float(local_misassembly_count)) + '\t' + str(counter['local']) + '\t' + ("0" if counter['local'] == 0 else str(float(local_misassembly_count)/counter['local']))
    print 'Valid query misassemblies:\t' + str(float(query_misassembly_count)/len(gff_misassemblies))
    print 'False positive rate:\t' + str(len(gff_misassemblies) - int(query_misassembly_count)) + '\t' + str(len(gff_misassemblies)) + '\t' + str((len(gff_misassemblies) - float(query_misassembly_count))/len(gff_misassemblies))
                


if __name__ == '__main__':
    main()