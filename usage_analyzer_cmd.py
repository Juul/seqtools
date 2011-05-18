#!/usr/bin/env python2.6
import sys, math, random, copy
import time
import getopt
import json
from UsageAnalyzer import *


def usage():

    print 'Usage: ' + __file__ + " -f genbank_file.gb -m some_method"
    print "  -f <genbank_file> Specifies the genbank file to use as input [required]"
    print "  -m <method> Specifies the method to run. Explained at the end of this help section. [required]"
    print "  -d Enable debug"
    print "  -a Absolute values. Output absolute counts instead of ratios."
    print "  -o <filename> Sets an output file"
    print "  -h Print this help screen"
    print " --------------------------------"
    print " Methods (for use with the -m flag):"
    print "  amino_acid_usage: Calculate the ratio of amino acid usage for the whole genome."
    print "    Use -x stop|ambiguous|both to exclude stop codons, ambiguous codons or both from being included in the ratios."
    print "  codon_usage: Calculate the ratio of codon usage for the whole genome."
    print "  codon_usage_per_position: Calculate the ratio of codon usage for each position."
    print "  codon_usage_per_amino_acid: Calculate the ratio of codon usage for each type of amino acid."
    print "  nucleotide_usage_per_position: Calculate the ratio of nucleotide usage for each position."


    exit(-1)
    
if __name__ == "__main__":


    optlist, args = getopt.getopt(sys.argv[1:], 'ahdf:o:m:x:c:')

    debug = False
    genbank_file = None
    output_file = None
    method = None
    get_absolute_values = False
    show_help = False
    exclude = None
    cutoff = None

    for o, a in optlist:
        if o == '-d':
            debug = True
        if o == '-f':
            genbank_file = a
        if o == '-o':
            output_file = a
        if o == '-m':
            method = a
        if o == '-a':
            get_absolute_values = True
        if o == '-x':
            exclude = a
        if o == '-c':
            cutoff = int(a)
        if o == '-h':
            show_help = True

    if (genbank_file == None) or (method == None) or (show_help == True):
        usage()

    ua = UsageAnalyzer(genbank_file)

    output = ua.get_usage_with_method(method, absolute_values=get_absolute_values, exclude=exclude, cutoff=cutoff)

    if output == None:
        print "Unknown method"
        usage()

    if output_file == None:
        print str(output)
    else:
        f = open(output_file, 'w')
        f.write(json.dumps(output))
        f.close()
    exit(0)

