#!/usr/bin/env python

from optparse import OptionParser
import os
import re
import sys
from collections import defaultdict
from sets import Set

PROG_NAME = "PILEUP"


def setup_options():
    parser = OptionParser()
    parser.add_option("-m", "--mis", dest="misassemblies", help="Misassemblies in gff format.", metavar="FILE")
    parser.add_option("-c", "--min-cutoff", dest="min_cutoff", help="Min-cutoff of signatures to mark suspicious.", type=int, default=2, metavar="FILE")
    
    (options,args) = parser.parse_args()

    if options.misassemblies == None:
        parser.error("You failed to provide the misassemblies file")

    return (options,args)


def main():
    (options, args) = setup_options()

    misassemblies = []

    for line in open(options.misassemblies, 'r'):
        misassemblies.append(line.strip().split('\t'))

    regions = find_suspicious_regions3(misassemblies, 2000, options.min_cutoff)

    for region in regions:
        print(region)


def find_suspicious_regions3(misassemblies, sliding_window = 2000, min_cutoff = 2):
    """
    Output any region that has multiple misassembly signature types within the sliding window.
    """

    regions =[]

    for misassembly in misassemblies:
        regions.append([misassembly[0], misassembly[3], 'START', misassembly[2]])
        regions.append([misassembly[0], misassembly[4], 'END', misassembly[2]])

    regions.sort(key = lambda region: (region[0], int(region[1])))

    """
    Example:

    relocref        36601   START   Breakpoint_finder                                    
    relocref        36801   END     Breakpoint_finder                                    
    relocref        67968   START   REAPR
    relocref        68054   START   REAPR
    relocref        69866   END     REAPR
    relocref        69867   START   REAPR
    relocref        71833   END     REAPR
    relocref        73001   START   Breakpoint_finder                                    
    relocref        73201   END     Breakpoint_finder   
    """

    # Store all the signatures, starting, and ending points within a given window.
    #start_points = deque([])
    #end_points = deque([])
    #signatures = deque([])

    signatures = []
    curr_contig = None
    count = 0
    suspicious_regions = []

    for index in xrange(0, len(misassemblies)): 

        curr_contig = misassemblies[index][0]

        count = 0
        second_index = index + 1
        signatures = [misassemblies[index][2]]
        start_point = int(misassemblies[index][3])
        end_point = int(misassemblies[index][4]) + sliding_window

        # While we are on the same contig, and still in the sliding window...
        while second_index < len(misassemblies) and \
                misassemblies[second_index][0] == curr_contig and \
                int(misassemblies[second_index][3]) < (int(misassemblies[index][4]) + sliding_window):

            if misassemblies[second_index][2] not in signatures:
                signatures.append(misassemblies[second_index][2])
                count += 1

            if int(misassemblies[second_index][4]) > end_point:
                end_point = int(misassemblies[second_index][4])
            second_index += 1

        if len(signatures) >= min_cutoff:
            suspicious_regions.append([misassemblies[index][0], '.', 'SUSPICIOUS', str(start_point), str(end_point), '.', '.', '.', 'color=#181009;' + ','.join(signatures)])
                

    # Hack to correct for overlapping suspicious regions.
    compressed_suspicious_regions = []

    prev_region = None
    for region in suspicious_regions:

        if prev_region is None:
            prev_region = region
        else:
            if prev_region[0] == region[0] and int(prev_region[4]) >= int(region[3]):
                prev_region[4] = region[4]
            else:
                compressed_suspicious_regions.append(prev_region)
                prev_region = region


    if prev_region:
        compressed_suspicious_regions.append(prev_region)

    return compressed_suspicious_regions


def find_suspicious_regions2(misassemblies, min_cutoff = 2):
    """
    Given a list of miassemblies in gff format
    """

    regions =[]

    for misassembly in misassemblies:
        regions.append([misassembly[0], misassembly[3], 'START', misassembly[2]])
        regions.append([misassembly[0], misassembly[4], 'END', misassembly[2]])

    regions.sort(key = lambda region: (region[0], int(region[1])))

    """
    Example:

    relocref        36601   START   Breakpoint_finder                                    
    relocref        36801   END     Breakpoint_finder                                    
    relocref        67968   START   REAPR
    relocref        68054   START   REAPR
    relocref        69866   END     REAPR
    relocref        69867   START   REAPR
    relocref        71833   END     REAPR
    relocref        73001   START   Breakpoint_finder                                    
    relocref        73201   END     Breakpoint_finder   
    """

    curr_contig = None
    curr_index = 0
    curr_length = -1

    start_indexes = []
    start_region = 0
    end_index = 0
    signatures = []
    recording = False

    signature_starts = defaultdict(list)

    curr_coverage = 0
    suspicious_regions = []

    for region in regions:

        if curr_contig is None:
            curr_contig = region[0]
            recording = False
            signature_starts = defaultdict(list)

        # We have found a new contig, process the previous contig results.
        if region[0] != curr_contig:

            curr_contig = region[0]
            recording = False

        if region[2] == 'START':
            curr_coverage += 1
            if region[3] not in signatures: signatures.append(region[3])
            signature_starts[region[3]].append(region[1])

            # Record start point.
            if curr_coverage == min_cutoff:
                start_region = region[1]
                recording == True

            start_indexes.append(region[1])

        else:
            
            curr_coverage -= 1

            end_index = region[1]
            if region[3] in signatures: signatures.remove(region[3]) 

            # If we were recording, and min signatures drop belows threshold,
            # then we need to output our results
            if curr_coverage < min_cutoff and recording:
                min_start = None

                suspicious_regions.append([region[0], '.', 'SUSPICIOUS', str(start_region), str(end_index), '.', '.', '.', 'color=#181009;' + ','.join(signatures)])
                signatures = []
                recording = False

        if curr_coverage >= min_cutoff:
            recording = True

    # Hack to correct for overlapping suspicious regions.
    compressed_suspicious_regions = []

    prev_region = None
    for region in suspicious_regions:

        compressed_suspicious_regions.append(region)
        """
        if prev_region is None:
            prev_region = region
        else:
            if prev_region[0] == region[0] and int(prev_region[4]) >= int(region[3]):
                print region[0]
                prev_region[4] = region[4]
            else:
                print region[0]
                compressed_suspicious_regions.append(prev_region)
                prev_region = region
        """

    if prev_region:
        compressed_suspicious_regions.append(prev_region)

    return compressed_suspicious_regions


def find_suspicious_regions(misassemblies, min_cutoff = 1):
    """
    Given a list of miassemblies in gff format
    """

    regions =[]

    for misassembly in misassemblies:
        regions.append([misassembly[0], misassembly[3], 'START', misassembly[2]])
        regions.append([misassembly[0], misassembly[4], 'END', misassembly[2]])

    regions.sort(key = lambda region: (region[0], int(region[1])))

    #for region in regions:
    #    print('\t'.join(region))

    """
    Example:

    relocref        36601   START   Breakpoint_finder                                    
    relocref        36801   END     Breakpoint_finder                                    
    relocref        67968   START   REAPR
    relocref        68054   START   REAPR
    relocref        69866   END     REAPR
    relocref        69867   START   REAPR
    relocref        71833   END     REAPR
    relocref        73001   START   Breakpoint_finder                                    
    relocref        73201   END     Breakpoint_finder   
    """

    curr_contig = None
    curr_index = 0
    curr_length = -1

    start_indexes = []
    start_region = 0
    end_index = 0
    signatures = []
    recording = False

    signature_starts = defaultdict(list)

    curr_coverage = 0

    for region in regions:

        if curr_contig is None:
            curr_contig = region[0]
            #start_indexes.append(region[1])
            recording = False
            signature_starts = defaultdict(list)

        # We have found a new contig, process the previous contig results.
        if region[0] != curr_contig:
            #calc_pileup(curr_contig, curr_length, start_index, end_index, alignments)

            #start_index = end_index
            curr_contig = region[0]
            curr_length = contig_lengths[curr_contig]
            recording = False

        if region[2] == 'START':
            curr_coverage += 1
            if region[3] not in signatures: signatures.append(region[3])
            signature_starts[region[3]].append(region[1])

            # Record start point.
            if curr_coverage == min_cutoff:
                start_region = region[1]
                recording == True


            start_indexes.append(region[1])

        else:
            
            curr_coverage -= 1

            end_index = region[1]
            if region[3] in signatures: signatures.remove(region[3]) 

            # If we were recording, and min signatures drop belows threshold,
            # then we need to output our results
            if curr_coverage < min_cutoff and recording:
                min_start = None
                mins = []
                for key in signature_starts.keys():
                    if len(signature_starts[key]) > 0:
                        mins.append(signature_starts[key][0])
                #print mins
                if min_start is None:
                    min_start = min(mins)

                print('%s\t%s\t%s\t%s\t%s\n' % (region[0], 'SUSPICIOUS', str(start_region), str(end_index), ','.join(signatures)))
                signatures = []
                recording = False

            signature_starts[region[3]].pop()

        

        if curr_coverage >= min_cutoff:
            recording = True

        #end_index += 1

    if curr_contig is not None:
        pass        

    return regions

if __name__ == '__main__':
    main()