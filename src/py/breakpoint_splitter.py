#!/usr/bin/env python
from __future__ import print_function
from optparse import OptionParser
import os
import sys

class ReadsSplitter:

    def __init__(self):
        self.options = None
        self.files_to_split = []

        self.getOptions()

    def go(self):
        for fn in self.files_to_split:
            self.splitFile(fn)

    def getOptions(self):
        parser = OptionParser()
        parser.add_option("-u", "--unaligned", dest="unaligned_dir", \
                help="Unaligned read directory", metavar="DIR")
        parser.add_option("-o", "--output", dest="output_dir",\
                help="Directory for output", metavar="DIR",\
                default="data/output/breakpoints/reads")

        (options, args) = parser.parse_args()
        self.options = options
        if options.unaligned_dir:
            for file_name in os.listdir(options.unaligned_dir):
                if 'unaligned' in file_name:
                    self.files_to_split.append(options.unaligned_dir + file_name)
                
    def splitFile(self, fn):
        if not os.path.isfile(fn):
            warning("%s DOES NOT EXIST" %(fn))
            exit(1)

        read_split_output_dir = self.options.output_dir
        ensure_dir(read_split_output_dir)

        read_split_output_1 = read_split_output_dir + os.path.split(fn)[1] + ".1" 
        read_split_output_2 = read_split_output_dir + os.path.split(fn)[1] + ".2"
        
        read_file = open(fn, 'r')

        r_o_1 = open(read_split_output_1, 'w')
        r_o_2 = open(read_split_output_2, 'w')

        for read in self.read_read(read_file):
            h1 = read[0].strip()
            read_contents = read[1].strip()
            h2 = read[2].strip()
            read_quality = read[3].strip()

            l = len(read_contents)
            l_1 = l/2
            l_2 = l - l_1
            #left
            h1_1 = h1 + "/1\n"
            read_contents_1 =  read_contents[0:l_1] + "\n"
            h2_1 = h2 + "/1\n"
            read_quality_1 = read_quality[0:l_1] + "\n"
            #right
            h1_2 = h1 + "/2\n"
            read_contents_2 = read_contents[l_1:]+ "\n"
            h2_2 = h2 + "/2\n"
            read_quality_2 = read_quality[l_1:] + "\n"

            r_o_1.write(h1_1)
            r_o_1.write(read_contents_1)
            r_o_1.write(h2_1)
            r_o_1.write(read_quality_1)

            r_o_2.write(h1_2)
            r_o_2.write(read_contents_2)
            r_o_2.write(h2_2)
            r_o_2.write(read_quality_2)

        r_o_1.close()
        r_o_2.close()
        read_file.close()



    def read_read(self, fp):
        while True:
            read_bundle = []
            for i in range(4):
                read_bundle.append(fp.readline())
            if not read_bundle[0]:
                break
            else:
                yield read_bundle


def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)

def warning(*objs):
    print("\tINFO: ",*objs, file=sys.stderr)

def main():
    '''
    splits read files for breakpoint
    '''
    splitter = ReadsSplitter()
    splitter.go()

if __name__=='__main__':
    main()

