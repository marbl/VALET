#!/usr/bin/env python
from __future__ import print_function
from optparse import OptionParser
from subprocess import call
import os
import sys
import math
import time
from time import gmtime, strftime

class BreakpointFinder:

    def __init__(self):
        self.options = None
        self.assembly_file = None
        self.singleton_halves = []
        self.read_lengths = {}
        self.bin_contents = {}
        self.inverse_bin_contents = {}
        self.surviving_bins = {}
        self.fasta_file = False
        self.getOptions()
        self.base_path = os.path.dirname(sys.argv[0])[:-len("/src/py")]
        self.set_locations()
        self.contig_coverage = {}
        self.contig_lengths = {}
        self.average_read_length = 0
        self.number_of_reads = 0

    def set_locations(self):
        #self.bowtie_dir = os.path.join(self.base_path, "bin/bowtie2-2.2.2/")
        #self.bowtie_loc = self.bowtie_dir + "bowtie2"
        self.bowtie_loc = "bowtie2"
        self.bowtie_build_loc = "bowtie2-build"
        self.breakpoint_dir = self.options.output_dir

        self.bowtie_index =  self.breakpoint_dir + "bowtie-index/"
        self.index_prefix = self.bowtie_index+ "breakpoint"
        ensure_dir(self.bowtie_index)

        self.sam_output_dir = self.breakpoint_dir + "sam/"
        self.sam_output_location = self.sam_output_dir + "breakpoints.sam"
        ensure_dir(self.sam_output_dir)

        self.singleton_dir = self.breakpoint_dir + "singleton/"
        self.singleton_loc = self.singleton_dir + "singletons.csv"

        ensure_dir(self.singleton_dir)

        self.conc_dir = self.breakpoint_dir + "concordance/"
        self.conc_loc = self.conc_dir + "concordants.csv"
        ensure_dir(self.conc_dir)

        self.breakpoint_file = self.breakpoint_dir + "breakpoints.csv"
        self.sorted_breakpoint_file = self.breakpoint_dir + "sorted_breakpoints.csv"
        self.binned_breakpoint_file = self.breakpoint_dir + "binned_breakpoints.csv"
        self.collapsed_breakpoint_file = self.breakpoint_dir\
                + "collapsed_breakpoints.csv"
        self.bins_of_interest_file = self.breakpoint_dir + "interesting_bins.bed"
        self.meta_file = self.breakpoint_dir + "meta_data.data"

        self.bin_contents_file = self.breakpoint_dir + "bin_contents.csv"

        self.reciprical_file = self.breakpoint_dir + "reciprical_breakpoints.csv"

    def read_coverages(self):
        with open(self.options.coverage_file) as cov_f:
            for line in cov_f:
                s_l = line.split()
                self.contig_coverage[s_l[0]] = float(s_l[1])

    def run_bowtie_index(self):
        FNULL = open('/dev/null', 'w')
        call_arr = [self.bowtie_build_loc, self.assembly_file, self.index_prefix]
        out_cmd(call_arr)
        call(call_arr, stdout = FNULL, stderr = FNULL)

    def run_bowtie_2(self):
        for file_name in os.listdir(self.reads_dir):
            if "reads" in file_name:
                read_type = "-q"
                if self.fasta_file:
                    read_type = "-f"
                call_arr = [self.bowtie_loc, '-x', self.index_prefix, '-U',\
                         self.reads_dir + file_name,\
                         '-S', self.sam_output_dir + file_name + '.sam',\
                         '--un', self.singleton_dir + file_name + '.singletons',\
                         '--al', self.conc_dir + file_name + '.reads', read_type,\
                         '-I 50' , '-X 800', '-p', str(self.options.threads), '--mp 6']
                #out_cmd(call_arr)
                call(call_arr)
            else:
                warning("Skipping potential read file: " + file_name)

    def sort_breakpoints(self):
        call_arr = ["sort", '-T ./', "-t\t", '-k 1n,1', '-k 2n,2', '-k 5n,5', self.breakpoint_file]
        out_cmd(call_arr)
        out_file = open(self.sorted_breakpoint_file, 'w')
        warning("This call outputs to file: ", self.sorted_breakpoint_file)
        call(call_arr, stdout=out_file)
        out_file.close()

    def bin_breakpoints(self):
        warning("About to start binning breakpoints")
        with open(self.sorted_breakpoint_file,'r') as breakpoints,\
                open(self.binned_breakpoint_file,'w') as out_file:
            for contig_bundle in self.read_contig(breakpoints):
                self.w_s = 1
                self.w_e = self.w_s + self.bin_size - 1
                for match in contig_bundle:
                    match_split = match.split('\t')
                    match_start = int(match_split[1])
                    match_end = int(match_split[4]) + match_start
                    if match_start == 0:
                        warning("Match start was 0 in bin breakpoints")
                        continue
                    elif match_start >= self.w_s and match_start <= self.w_e:
                        out_file.write(match.strip()+('\t%s\n'%(self.w_s)))
                        self.add_to_bin_contents(match_split[0]+"\t"+str(self.w_s),match_split[2])
                    else:
                        self.w_s = ((match_start / (self.bin_size+1)) * (self.bin_size+1)) #+1
                        self.w_e = self.w_s + self.bin_size - 1
                        out_file.write(match.strip()+('\t%s\n'%(self.w_s)))
                        self.add_to_bin_contents(match_split[0]+"\t"+str(self.w_s), match_split[2])
        warning("Done binning contigs")

    def add_to_bin_contents(self, key, half_id):
        if key not in self.bin_contents:
            self.bin_contents[key] = [half_id]
        else:
            self.bin_contents[key].append(half_id)

        if half_id not in self.inverse_bin_contents:
            self.inverse_bin_contents[half_id] = [key]
        else:
            self.inverse_bin_contents[half_id].append(key)

    def output_bin_contents(self, only_surviving = False):
        with open(self.bin_contents_file,'w') as b_c_file:
            for bin in self.bin_contents:
                if only_surviving:
                    if bin in self.surviving_bins:
                        b_c_file.write(bin + "\n")
                        for content in self.bin_contents[bin]:
                            b_c_file.write(content + "\t")
                        b_c_file.write("\n")
                else:
                    b_c_file.write(bin + "\n")
                    for content in self.bin_contents[bin]:
                        b_c_file.write(content + "\t")
                    b_c_file.write("\n")



    def collapse_bins(self):
        warning("About to collapse bins")
        call_str = "awk '{printf \"%%s\\t%%s\\n\",$1,$6}' %s | sort -T ./| uniq -c > %s" %(self.binned_breakpoint_file, self.collapsed_breakpoint_file)
        out_cmd(call_str)
        call(call_str,shell=True)

    def trim_bins_3(self):
        avg_coverage_in_bin = {}
        should_pass = {}
        b_c_d = {}
        b_c_d_r = {}
        warning("About to assemble surviving bins.")
        ave_read_len = self.average_read_length / float(self.number_of_reads)
        warning("Beginning pass_1")
        with open(self.collapsed_breakpoint_file,'r') as pass_1:
            for line in pass_1:
                split_l = line.split()
                num_matches_in_bin = int(split_l[0])
                bin_name = split_l[1] + '\t' + split_l[2]
                #TODO: Come back here
                #if bin_name in avg_coverage_in_bin:
                if split_l[1] in self.contig_coverage:
                    avg_coverage_in_bin[bin_name] = self.contig_coverage[split_l[1]]
                else:
                    warning("Skipping breakpoint bin: %s because not found in avg_coverage dict" % (bin_name))
                    continue
                if (2*ave_read_len * num_matches_in_bin) / self.bin_size >= avg_coverage_in_bin[bin_name]/4.0:
                    self.surviving_bins[bin_name] = True
                else:
                    warning("Removing bin: %s because eq matches: %f < cutoff: %f" % (bin_name, (2* ave_read_len*num_matches_in_bin)/float(self.bin_size), avg_coverage_in_bin[bin_name]/4.0))
        warning("End pass 1")
        warning("Surviving bins assembled: %d surviving bins." % (len(self.surviving_bins)))
        warning("Ave read len: %f " % (self.average_read_length / float(self.number_of_reads)))

        warning("About to output bin contents")
        self.output_bin_contents(True)
        warning("Bin contents outputted")

        warning("About to assemble bin contents dictionary reciprocity")
        for bin_1 in self.surviving_bins:
            b_c_d_r[bin_1] = {}
            for read in self.bin_contents[bin_1]:
                sister_read = sister_name(read)
                if sister_read in self.inverse_bin_contents:
                    for bin_2 in self.inverse_bin_contents[sister_read]:
                        if bin_2 in self.surviving_bins:
                            if bin_2 in b_c_d_r[bin_1]:
                                b_c_d_r[bin_1][bin_2] += 1
                            else:
                                b_c_d_r[bin_1][bin_2] = 1

        #for bin_1 in self.bin_contents:
        #    if bin_1 in self.surviving_bins:
        #        b_c_d_r[bin_1] = {}
        #        for bin_2 in self.surviving_bins:
        #            if bin_2 in b_c_d_r and bin_1 in b_c_d_r[key_2]:
        #                b_c_d_r[bin_1][bin_2] = b_c_d_r[bin_2][bin_1]
        #                continue
        #            b_c_d_r[bin_1][bin_2] = len([half_read for half_read in self.bin_contents[bin_1] if sister_name(half_read) in self.bin_contents[bin_2]])
        warning("Assembled that dictionary size: %d" % (sys.getsizeof(b_c_d_r)))


        warning("About to output some reciprocals")
        with open(self.reciprical_file,'w') as reciprocals:
            for bin in self.surviving_bins:
                rec = find_reciprical_pair_2(b_c_d_r, bin)
                if rec:
                    reciprocals.write(bin + "\t" + rec + "\n")

        warning("Outputted reciprocals")

        warning("About to attempt to output interesting bins")
        with open(self.collapsed_breakpoint_file,'r') as pass_3,\
                open(self.bins_of_interest_file,'w') as out_file:
                    color = "#FF9900"
                    for line in pass_3:
                        split_line = line.split()
                        current_key = split_line[1] + "\t" + split_line[2]
                        if current_key in self.surviving_bins:
                            num_sisters = 0
                            if current_key not in b_c_d_r:
                                warning("Current key: %s not in b_c_d_r keys\n"\
                                        % (current_key))
                            for key in b_c_d_r[current_key]:
                                num_sisters += b_c_d_r[current_key][key]
                            #if num_sisters < 5:
                            #    warning("Skipping bin: %s because number of sisters: %d was low" % (current_key, num_sisters))
                            #    continue
                            rec = find_reciprical_pair_2(b_c_d_r, current_key)
                            if rec == None:
                                warning("No Reciprical for bin: %s so skipping" % (current_key))
                                continue
                            out_file.write("%s\t%d\t%d\tBreakpoint_finder\n"\
                                    % (split_line[1],\
                                    int(split_line[2])-1,\
                                    int(split_line[2])+self.bin_size-1))
                            # out_file.write("%s\tBreakpoint_finder"\
                            #         "\tBreakpoint_Finder_excessive_alignment"\
                            #         "\t%d\t%d\t%d\t.\t.\t"\
                            #         "singletons_aligned_in_bin=%f;color=%s"\
                            #         "number_of_sisters=%d;reciprocal=%s\n"\
                            #         %(split_line[1],\
                            #         int(split_line[2]),\
                            #         int(split_line[2])+self.bin_size-1,\
                            #         int(split_line[0]),\
                            #         float(split_line[0]),\
                            #         color,\
                            #         num_sisters,\
                            #         str(rec)))






    def trim_bins_2(self):
        avg_coverage_in_bin = {}
        should_pass = {}
        b_c_d = {}
        b_c_d_r = {}
        warning("About to assemble surviving bins.")
        ave_read_len = self.average_read_length / float(self.number_of_reads)
        with open(self.collapsed_breakpoint_file,'r') as pass_1:
            for line in pass_1:
                split_l = line.split()
                num_matches_in_bin = int(split_l[0])
                avg_coverage_in_bin[split_l[1] + "\t" + split_l[2]] = self.contig_coverage[split_l[1]]
                if (2*ave_read_len * num_matches_in_bin) / self.bin_size >= avg_coverage_in_bin[split_l[1] + "\t" + split_l[2]]/4.0:
                    self.surviving_bins.append(split_l[1] + "\t" + split_l[2])
                else:
                    warning("Removing bin: %s because eq matches: %f < cutoff: %f" % (split_l[1]+"\t"+split_l[2], (2* ave_read_len*num_matches_in_bin)/float(self.bin_size), avg_coverage_in_bin[split_l[1] + "\t" + split_l[2]]/4.0))
        warning("Surviving bins assembled: %d surviving bins." % (len(self.surviving_bins)))
        warning("Ave read len: %f " % (self.average_read_length / float(self.number_of_reads)))

        warning("About to output bin contents")
        self.output_bin_contents(True)
        warning("Bin contents outputted")

        warning("About to assemble bin contents dictionary reciprocity")
        for key in self.bin_contents.keys():
            if key in self.surviving_bins:
                b_c_d_r[key] = {}
                for key_2 in self.surviving_bins:
                    if key_2 in b_c_d_r.keys() and key in b_c_d_r[key_2].keys():
                        b_c_d_r[key][key_2] = b_c_d_r[key_2][key]
                        continue
                    b_c_d_r[key][key_2] = [val for val in self.bin_contents[key] if sister_name(val) in self.bin_contents[key_2]]
        warning("Assembled that dictionary size: %d" % (sys.getsizeof(b_c_d_r)))


        warning("About to output some reciprocals")
        with open(self.reciprical_file,'w') as reciprocals:
            for bin in self.surviving_bins:
                rec = find_reciprical_pair(b_c_d_r, bin)
                if rec:
                    reciprocals.write(bin + "\t" + rec + "\n")

        warning("Outputted reciprocals")

        warning("About to attempt to output interesting bins")
        with open(self.collapsed_breakpoint_file,'r') as pass_3,\
                open(self.bins_of_interest_file,'w') as out_file:
                    color = "#FF9900"
                    for line in pass_3:
                        split_line = line.split()
                        current_key = split_line[1] + "\t" + split_line[2]
                        if current_key in self.surviving_bins:
                            num_sisters = 0
                            if current_key not in b_c_d_r.keys():
                                warning("Current key: %s not in b_c_d_r keys\n %s"\
                                        % (current_key, str(b_c_d_r.keys)))
                            for key in b_c_d_r[current_key].keys():
                                num_sisters += len(b_c_d_r[current_key][key])
                            #if num_sisters < 5:
                            #    warning("Skipping bin: %s because number of sisters: %d was low" % (current_key, num_sisters))
                            #    continue
                            rec = find_reciprical_pair(b_c_d_r, current_key)
                            #if rec == None:
                            #    warning("No Reciprical for bin: %s so skipping" % (current_key))
                            #    continue
                            out_file.write("%s\tBreakpoint_finder"\
                                    "\tBreakpoint_Finder_excessive_alignment"\
                                    "\t%d\t%d\t%d\t.\t.\t"\
                                    "singletons_aligned_in_bin=%f;color=%s"\
                                    "number_of_sisters=%d;reciprocal=%s\n"\
                                    %(split_line[1],\
                                    int(split_line[2]),\
                                    int(split_line[2])+self.bin_size-1,\
                                    int(split_line[0]),\
                                    float(split_line[0]),\
                                    color,\
                                    num_sisters,\
                                    str(rec)))




    def trim_bins(self):
        avg_bin_size = 0
        num_bins = 0
        summed_var = 0
        std_dev = 0
        cutoff = 0
        with open(self.collapsed_breakpoint_file,'r') as pass_1:
            for line in pass_1:
                avg_bin_size += int(line.split()[0])
                num_bins += 1
            avg_bin_size = float(avg_bin_size)/num_bins
        with open(self.collapsed_breakpoint_file,'r') as pass_2:
            for line in pass_2:
                summed_var += (float(line.split()[0]) - avg_bin_size)**2
            summed_var = summed_var / num_bins
            std_dev = math.sqrt(summed_var)
        with open(self.collapsed_breakpoint_file,'r') as pass_3,\
                open(self.bins_of_interest_file,'w') as out_file:
                    cutoff = 2 * std_dev
                    color = "#FF9900"
                    for line in pass_3:
                        split_line = line.split()
                        if (float(split_line[0]) - avg_bin_size) > cutoff:
                            out_file.write("%s\tBreakpoint_finder"\
                                    "\tBreakpoint_Finder_excessive_alignment"\
                                    "\t%d\t%d\t%d\t.\t.\tavg_bin_size=%f;std_dev=%f;"\
                                    "singletons_aligned_in_bin=%f;color=%s\n"\
                                    %(split_line[1],\
                                    int(split_line[2]),\
                                    int(split_line[2])+self.bin_size-1,\
                                    int(split_line[0]),\
                                    float(avg_bin_size), \
                                    float(std_dev),\
                                    float(split_line[0]),\
                                    color))
                            self.surviving_bins.append(split_line[1] + "\t" + split_line[2])
        with open(self.meta_file, 'w') as meta_f:
            meta_f.write("Avg bin size: %f\nStd_Dev: %f\nCutoff: %f\n"\
                    % (float(avg_bin_size), float(std_dev), float(cutoff)))




    def read_contig(self, fp):
        warning("We're about to start the read_contig generator")
        line_counter = 0
        run_flag = True
        contig_bundle = []
        line_buffer = ""
        while run_flag:
            #start_time = time.time()
            current_contig = ""
            if line_buffer != "":
                current_contig = line_buffer.split('\t')[0]
                contig_bundle.append(line_buffer)
                line_buffer = ""
            while True:
                line = fp.readline()
                line_counter += 1
                if not line:
                    run_flag = False
                    warning("We've reached the end of the breakpoint file %d lines read." % ( line_counter))
                    break
                elif line.split('\t')[0] != current_contig:
                    line_buffer = line
                    break
                else:
                    contig_bundle.append(line)
            ret_bundle = contig_bundle
            contig_bundle = []
            #warning("Read contig returned in : %f" %(time.time() - start_time))
            yield ret_bundle

    def read_in_lengths(self):
        log_file = open(self.breakpoint_dir + "log.log", 'w')
        for file_name in os.listdir(self.conc_dir):
            print("Checking name %s\n" %(file_name))
            if '.reads' in file_name:
                print("Reading file: %s\n" %(file_name))
                with open(self.conc_dir + file_name, 'r') as reads_file:
                    for (read,length) in self.read_read(reads_file):
                        #print ("Read: %s Length: %s" %(read,length))
                        self.read_lengths[read] = length

    def read_read(self, fp):
        run_flag = True
        while run_flag:
            line_bundle = []
            for i in range(4):
                line_bundle.append(fp.readline())
            if line_bundle[0] == '':
                run_flag = False
            else:
                read = line_bundle[0][1:]
                length = len(line_bundle[1].strip())
                yield (read,length)

    '''
    Outputs breakpoints to file self.breakpoint_file
    Position Sequence Name Flag Length
    Tab sep
    '''
    def detect_breakpoints(self):
        with open(self.breakpoint_file, 'w') as out_file:
            for file_name in os.listdir(self.sam_output_dir):
                if "sam" in file_name:
                    with open(self.sam_output_dir+file_name, 'r') as in_file:
                        for line in in_file:
                            line_components = line.split('\t')
                            if line_components[0] == "@SQ":
                                self.contig_lengths[line_components[1].split(":")[1]] =\
                                         int(line_components[2].split(":")[1])
                            if line[0] != '@':
                                should_output = True

                                self.average_read_length += len(line_components[9])
                                self.number_of_reads += 1
                                if line_components[2] != "*":
                                    contig_length = self.contig_lengths[line_components[2]]
                                else:
                                    #warning("Not outputting due to *")
                                    should_output = False

                                if should_output and line_components[3] == '0':
                                    #warning("Not outputting because unaligned")
                                    should_output = False
                                elif should_output\
                                        and (line_components[4] < len(line_components[10]))\
                                        and ((contig_length - line_components[4]) < len(line_components[10])):
                                            should_output = False
                                            warning("Not outputting because close to ends")
                                #elif should_output and line_components[4] <= float(contig_length)/self.options.end_trim or (line_components[4] >= contig_length - (float(contig_length)/self.options.end_trim)):
                                #    warning("Not outputting because close to ends")
                                #    should_output = False
                                for component in line_components:
                                    if "NM:" in component: #Number of mismatches
                                        if should_output and  \
                                                int(component.split(":i:")[1]) >= 6:
                                                    should_output = False
                                                    #warning("Not outputting due to number of mismatches")

                                if line_components[3] != '0' and  should_output:
                                    out_file.write("%s\t%s\t%s\t%s\t%d\n" % \
                                            (line_components[2],\
                                            line_components[3],\
                                            line_components[0], line_components[1],\
                                            len(line_components[9])))

    def go(self):
        warning("Arguments used to call breakpoint finder: %s" % (" ".join(sys.argv)))
        self.read_coverages()
        self.run_bowtie_index()
        self.run_bowtie_2()
        self.detect_breakpoints()
        self.sort_breakpoints()
        self.bin_breakpoints()
        self.collapse_bins()
        self.trim_bins_3()

    def getOptions(self):
        parser = OptionParser()
        parser.add_option("-a", "--assembly-file", dest="assembly_file",\
                help="Assembly File to search", metavar="FILE")
        parser.add_option("-r", "--reads-dir", dest="reads", \
                help="Directory of Reads", metavar="PATH")
        parser.add_option("-b", "--bin-size", dest="bin_size", \
                help="Bin size", metavar="SIZE", type="int")
        parser.add_option("-o", "--output", dest="output_dir", \
                help="Output directory", metavar="DIR",\
                default="data/output/breakpoint/")
        parser.add_option("-c", "--coverage", dest="coverage_file", \
                help="Coverage File", metavar="PATH")
        parser.add_option("-t", "--trim", dest="end_trim", \
                help="Percentage of seq to trim", metavar="PERCENT", \
                default="10", type="int")
        parser.add_option("-p", "--threads", dest="threads", \
                help="Number of threads", default=10, type="int")
        parser.add_option("-f", "--fasta", dest="fasta_file", \
                help="Reads are in FASTA format file. (default is FASTQ)",
                action="store_true", default=False)
        (options, args) = parser.parse_args()
        self.options = options
        if options.reads:
            self.reads_dir = options.reads
        if options.bin_size:
            self.bin_size = options.bin_size
        else:
            self.bin_size = 500
        if options.fasta_file:
            self.fasta_file = True
        if not options.assembly_file:
            warning("Did not provide assembly file")
            parser.print_help()
            exit(1)
        else:
            self.assembly_file = options.assembly_file
        warning(self.singleton_halves)
        if not options.coverage_file:
            warning("Did not provide coverage file")
            parser.print_help()
            exit(1)

def sister_name(s):
    #s_split = s.split()
    rd = s
    if rd.endswith("/1"):
        return rd[:-2] + "/2"
    elif rd.endswith("/2"):
        return rd[:-2] + "/1"
    else:
        warning("Failed to format for sister: %s" %(s))
        return ""

def read_coverage(cov_file, contig, start,end):
    contig_lines = []
    cumulative = 0
    n = 0
    with open(cov_file,'r') as c_file:
        for line in c_file:
            if contig in line:
                contig_lines.append(line)
    warning("Contig Lines Len: %d" %(len(contig_lines)))

    for line in contig_lines:
        line_s = line.split()
        if int(line_s[1])>= start and int(line_s[1])<end:
            cumulative+= int(line_s[3])
            n+=1
    if n == 0:
        warning("While attempting to fetch coverage for contig: %s no information"\
                %(contig))
        return 0
    else:
        return float(cumulative)/float(n)


def find_reciprical_pair_2(b_c_d_r, bin):
    max_partner = []
    max_partner_key = []
    for partner in b_c_d_r[bin]:
        if len(max_partner) == 0 or b_c_d_r[bin][partner] > max_partner[0]:
            max_partner.append(b_c_d_r[bin][partner])
            max_partner_key.append(partner[:])

    for (max_partner_i, max_partner_key_i) in zip(max_partner[::-1], max_partner_key[::-1]):
        max_alternate = 0
        max_alternate_key = ""
        for alternate in b_c_d_r[max_partner_key_i]:
            if b_c_d_r[max_partner_key_i][alternate] > max_alternate:
                max_alternate = b_c_d_r[max_partner_key_i][alternate]
                max_alternate_key = alternate[:]
        if max_alternate_key == bin:
            return max_partner_key_i
        else:
            return None
    return None



def find_reciprical_pair(b_c_d_r, bin):
    max_partner = []
    max_partner_key = []
    for partner in b_c_d_r[bin].keys():
        if len(max_partner) == 0 or len(b_c_d_r[bin][partner]) > max_partner[0]:
            max_partner.append(len(b_c_d_r[bin][partner]))
            max_partner_key.append(partner[:])

    for (max_partner_i, max_partner_key_i) in zip(max_partner[::-1], max_partner_key[::-1]):
        max_alternate = 0
        max_alternate_key = ""
        for alternate in b_c_d_r[max_partner_key_i].keys():
            if len(b_c_d_r[max_partner_key_i][alternate]) > max_alternate:
                max_alternate = len(b_c_d_r[max_partner_key_i][alternate])
                max_alternate_key = alternate[:]
        if max_alternate_key == bin:
            return max_partner_key_i
        else:
            return None
    return None


def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)

def warning(*objs):
    timestamp = strftime("%Y-%m-%d %H:%M:%S", gmtime())
    timestamp = "[" + timestamp + "] "
    print(timestamp + "\tWARNING: ",*objs, file=sys.stderr)
    sys.stderr.flush()

def out_cmd(*objs):
    print("="*75, file=sys.stderr)
    print("About to exec: ", *objs, file=sys.stderr)

def main():
    '''
    splits read files for breakpoint
    '''
    finder = BreakpointFinder()
    finder.go()

if __name__=='__main__':
    main()
