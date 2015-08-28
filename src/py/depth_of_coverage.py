#!/usr/bin/env python

from __future__ import print_function
from multiprocessing import Process, Lock
from optparse import OptionParser
#from scipy.interpolate import UnivariateSpline
import numpy as np
import scipy.signal
import collections
import copy
import math
import os
import sys
import time


PROG_NAME = "DEPTH_COV"


def setup_options():
    parser = OptionParser()
    parser.add_option("-a", "--abundance-file", dest="abundance_file", help="Assembler provided abundance file", metavar="FILE")
    parser.add_option("-m", "--mpileup-file", dest="mp_file_loc", help="mpileup output file", metavar="FILE")
    parser.add_option("-p", "--pileup-file", dest="pileup_filename", help="pileup file name", metavar="FILE")
    parser.add_option("-o", "--output", dest="output_location", help="output location for depth of coverage tool", default=None, metavar="FILE")
    parser.add_option("-w", "--window-size", dest="window_size", help="window size to sweep over base pairs", type="int", default = 1)
    parser.add_option("-g", "--gff", dest="gff_format", default=False, action='store_true')
    parser.add_option("-e", "--empirical", dest="use_empirical", default=False, action='store_true')
    parser.add_option("-i", "--ignore", dest="ignore_ends", help="Ignore the first/last i bps of the read", type="int", default=0)
    parser.add_option("-c", "--chunks", dest="chunks", help="Number of output file chunks.", type="int", default=1)

    (options,args) = parser.parse_args()

    should_quit = False
    if options.mp_file_loc == None:
        parser.error("You failed to provide the mpileup file")
    if should_quit:
        parser.help()
        exit(-1)

    return (options,args)


def warning(*objs):
    print("WARNING: ", *objs, file=sys.stderr)


def debug(*objs):
    print("DEBUG: ", *objs, file=sys.stderr)


def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d) and d:
        print("true")
        os.makedirs(d)


def find_coverage_errors_scipy(mpileup_file, output_location, window_size, abundance_dict, write_lock):
    """
    Find coverage errors using scipy and it's rolling median implementation.
    """

    bad_cvg_file = open(output_location,'a') if output_location else sys.stdout

    prev_contig = None
    lower_hinge = 0
    upper_hinge = 0
    end_pos = 0
    region_index = -1
    flagged_regions = []

    prev_contig = None
    length = 0
    curr_coverages = []

    for record in open(mpileup_file, 'r'):
        fields = record.strip().split()

        # contig00001     1       A       1       ^~,     I
        if prev_contig is None:
            prev_contig = fields[0]
            lower_hinge = abundance_dict[prev_contig][1]
            upper_hinge = abundance_dict[prev_contig][2]
            end_pos = 0
            curr_coverages = []

        if prev_contig != fields[0]:
            # Output previous results.
            if prev_contig:
                # Calculate median coverage, and lower/upper hinges.
                #abundance_dict[prev_contig] = tukey_summary(curr_coverages, window_size)
                #coverage_meds = scipy.signal.medfilt(np.array(curr_coverages), kernel_size = window_size)
                coverages_np = np.asarray(curr_coverages)
                coverage_meds = np.convolve(coverages_np, np.ones((window_size,))/window_size, mode='full')

                i = 0
                while i < len(coverage_meds):

                    if not in_range(coverage_meds[i], lower_hinge, upper_hinge):
                        cov_type = "Low_coverage"
                        color = "#7800ef"
                        if coverage_meds[i] > upper_hinge:
                            cov_type = "High_coverage"
                            color = "#0077ee"

                        median = coverage_meds[i]
                        start_pos = i + 1
                        # Keep incrementing until we hit a region no longer outside the range and of the same coverage type.
                        end_pos = i + 1

                        i += 1
                        while i < len(coverage_meds) and not in_range(coverage_meds[i], lower_hinge, upper_hinge):

                            # Make sure the coverage error type is the same as before.
                            new_cov_type = "Low_coverage"
                            if coverage_meds[i] > upper_hinge:
                                new_cov_type = "High_coverage"

                            if new_cov_type == cov_type:
                                end_pos = i + 1
                                i += 1
                            else:
                                break

                        flagged_regions.append([prev_contig, PROG_NAME, cov_type, start_pos, end_pos, median, lower_hinge, upper_hinge, color])
                        #print(flagged_regions)

                    else:
                        i += 1

            # Reset settings.
            prev_contig = fields[0]
            lower_hinge = abundance_dict[prev_contig][1]
            upper_hinge = abundance_dict[prev_contig][2]
            end_pos = 0
            curr_coverages = []

        curr_coverages.append(int(fields[3]))
        length += 1

        ## Append the bp coverage to the window.
        #coverage_window.append(int(fields[3]))
        #end_pos += 1

        #if prev_contig != fields[0]:

        #    prev_contig = fields[0]
        #    length = 0
        #    curr_coverages = []


    # Check the last remaining contig.
    if prev_contig:
        # Calculate median coverage, and lower/upper hinges.
        #abundance_dict[prev_contig] = tukey_summary(curr_coverages, window_size)
        coverage_meds = scipy.signal.medfilt(np.array(curr_coverages), kernel_size = window_size)

        i = 0
        while i < len(coverage_meds):

            if not in_range(coverage_meds[i], lower_hinge, upper_hinge):
                cov_type = "Low_coverage"
                color = "#7800ef"
                if coverage_meds[i] > upper_hinge:
                    cov_type = "High_coverage"
                    color = "#0077ee"

                median = coverage_meds[i]
                start_pos = i + 1
                # Keep incrementing until we hit a region no longer outside the range and of the same coverage type.
                end_pos = i + 1

                i += 1
                while i < len(coverage_meds) and not in_range(coverage_meds[i], lower_hinge, upper_hinge):

                    # Make sure the coverage error type is the same as before.
                    new_cov_type = "Low_coverage"
                    if coverage_meds[i] > upper_hinge:
                        new_cov_type = "High_coverage"

                    if new_cov_type == cov_type:
                        end_pos = i + 1
                        i += 1
                    else:
                        break

                flagged_regions.append([prev_contig, PROG_NAME, cov_type, start_pos, end_pos, median, lower_hinge, upper_hinge, color])
                #print(flagged_regions)

            else:
                i += 1


    write_lock.acquire()
    for region in flagged_regions:
        bad_cvg_file.write("%s\t%s\t%s\t%d\t%d\t%f\t.\t.\tlow=%f;high=%f;color=%s\n" % (region[0], region[1], \
                    region[2], region[3], region[4], region[5], region[6], region[7], region[8]))
    write_lock.release()
    #if prev_contig:
    #    abundance_dict[prev_contig] = tukey_summary(curr_coverages, window_size)


def find_coverage_errors(mpile_file, output_location, window_size, abundance_dict, write_lock):


    bad_cvg_file = open(output_location,'a') if output_location else sys.stdout

    coverage_window = collections.deque(maxlen = window_size)
    flagged_regions = []

    prev_contig = None
    lower_hinge = 0
    upper_hinge = 0
    end_pos = 0
    region_index = -1

    for line in open(mpile_file, 'r'):
        fields = line.split()

        # contig00001     1       A       1       ^~,     I
        if prev_contig is None:
            prev_contig = fields[0]
            lower_hinge = abundance_dict[prev_contig][1]
            upper_hinge = abundance_dict[prev_contig][2]
            end_pos = 0

        if prev_contig != fields[0]:
            # Output previous results and clear deque.
            coverage_window.clear()
            prev_contig = fields[0]
            lower_hinge = abundance_dict[prev_contig][1]
            upper_hinge = abundance_dict[prev_contig][2]
            end_pos = 0

        # Append the bp coverage to the window.
        coverage_window.append(int(fields[3]))
        end_pos += 1

        # If the coverage window is full, check and compare with median.
        if len(coverage_window) >= window_size:
            # TODO: Deepcopy is inefficient.
            copy_window = sorted(copy.deepcopy(coverage_window))
            median = copy_window[len(copy_window) / 2]
            if not len(copy_window) % 2:
                median = (copy_window[len(copy_window) / 2] + copy_window[len(copy_window) / 2 - 1]) / 2.0

            if not in_range(median, lower_hinge, upper_hinge):
                cov_type = "Low_coverage"
                color = "#7800ef"
                if median > upper_hinge:
                    cov_type = "High_coverage"
                    color = "#0077ee"

                # Extend previous window?
                if len(flagged_regions) > 0 and \
                        flagged_regions[region_index][0] == prev_contig and \
                        flagged_regions[region_index][2] == cov_type and \
                        int(flagged_regions[region_index][4]) >= end_pos - len(coverage_window) + 1 and \
                        int(flagged_regions[region_index][4]) <= end_pos:
                    flagged_regions[region_index][4] = end_pos

                else:
                    flagged_regions.append(\
                        [prev_contig, PROG_NAME, cov_type, end_pos - len(coverage_window) + 1, end_pos, median, lower_hinge, upper_hinge, color])
                    region_index += 1

    write_lock.acquire()
    for region in flagged_regions:
        bad_cvg_file.write("%s\t%s\t%s\t%d\t%d\t%f\t.\t.\tlow=%f;high=%f;color=%s\n" % (region[0], region[1], \
                    region[2], region[3], region[4], region[5], region[6], region[7], region[8]))
    write_lock.release()


def main():
    (options, args) = setup_options()
    abundance_file = options.abundance_file
    mpile_file = options.mp_file_loc

    if options.output_location:
        ensure_dir(os.path.dirname(options.output_location))

    window_size = options.window_size
    abundance_dict = {}

    write_lock = Lock()

    if options.abundance_file:
        read_abundances(abundance_file, abundance_dict)
    else:
        # Calculate the coverage and window from the data itself.
        calculate_coverages(mpile_file, abundance_dict, window_size)

    processes = []

    # If we have more than one chunk, we need to append *.{chunk_num} to the file.
    if options.chunks == 1:
        p = Process(target=find_coverage_errors_scipy, args=(mpile_file, options.output_location, window_size, abundance_dict, write_lock))
        p.start()
        processes.append(p)

    else:
        for i in range(0, options.chunks):

            if os.path.isfile(mpile_file + '.' + str(i)):
                p = Process(target=find_coverage_errors_scipy, args=(mpile_file + '.' + str(i), options.output_location, window_size, abundance_dict, write_lock))
                p.start()
                processes.append(p)

    for p in processes:
        p.join()


def in_range(num, low, high):
    return num >= low and num <= high


def read_abundances(fp, a_dict):
    fp = open(fp,'r')
    for line in fp:
        (contig, avg_cov) = line.split()
        a_dict[contig] = float(avg_cov)
    fp.close()


def tukey_summary(orig_array, window_size = 1):
    """ Given an array of integers, return median, and lower/upper hinge tuple. """

    #array = []
    array = orig_array
    array = np.convolve(np.array(array), np.ones(window_size)/window_size, mode='valid')
    #array = scipy.signal.medfilt(np.array(orig_array), kernel_size = window_size)
    #spline = UnivariateSpline(range(1,len(orig_array)+1), np.array(orig_array))
    #array = spline(range(1,len(orig_array)+1))

    #s = UnivariateSpline(range(1,len(a)+1),a)
    array.sort()
    #print(len(array))

    median = array[len(array) / 2]
    if not len(array) % 2:
        median = (array[len(array) / 2] + array[len(array) / 2 - 1]) / 2.0

    first_quantile = array[len(array) / 4]
    third_quantile = array[3 * len(array) / 4]

    iqr = third_quantile - first_quantile

    lower_hinge = math.floor(first_quantile - 2 * iqr)
    upper_hinge = math.ceil(third_quantile + 2 * iqr)

    print(median, first_quantile, third_quantile, lower_hinge, upper_hinge)

    return (median, lower_hinge, upper_hinge)


def calculate_coverages(mpileup_file, abundance_dict, window_size = 1):
    """ For each contig, calculate coverage and lower/upper hinges. """

    prev_contig = None
    length = 0
    curr_coverages = []

    for record in open(mpileup_file, 'r'):
        fields = record.strip().split()

        if prev_contig != fields[0]:
            if prev_contig:
                # Calculate median coverage, and lower/upper hinges.
                abundance_dict[prev_contig] = tukey_summary(curr_coverages, window_size)

            prev_contig = fields[0]
            length = 0
            curr_coverages = []

        curr_coverages.append(int(fields[3]))
        length += 1

    if prev_contig:
        abundance_dict[prev_contig] = tukey_summary(curr_coverages, window_size)


def get_average_coverage(a_dict):
    l = a_dict.values()
    mean =  reduce(lambda x, y: x + y, l) / len(l)
    variances = map(lambda x: (x - avg)**2, l)
    var = reduce(lambda x, y: x + y, variances)/ len(variances)
    std = math.sqrt(var)
    return (mean, var, std)


if __name__ == '__main__':
    main()
