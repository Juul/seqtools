#!/usr/bin/env python2.6
import re
import sys, math, random, copy
import time
import getopt
from codonizer import *


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


    optlist, args = getopt.getopt(sys.argv[1:], 'ahdf:o:m:x:')

    debug = False
    genbank_file = None
    output_file = None
    method = None
    get_absolute_values = False
    show_help = False
    exclude = None

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
        if o == '-h':
            show_help = True

    if (genbank_file == None) or (method == None) or (show_help == True):
        usage()

    cod = Codonizer(genbank_file)

    if method == 'codon_usage':
        output = cod.calc_codon_usage(absolute_values=get_absolute_values)
    elif method == 'codon_usage_per_position':
        output = cod.calc_codon_usage_per_position(absolute_values=get_absolute_values)
    elif method == 'codon_usage_per_amino_acid':
        output = cod.calc_codon_usage_per_amino_acid(absolute_values=get_absolute_values)
    elif method == 'nucleotide_usage_per_position':
        output = cod.calc_nucleotide_usage_per_position(absolute_values=get_absolute_values)
    elif method == 'amino_acid_usage':
        if exclude == None:
            output = cod.calc_amino_acid_usage(absolute_values=get_absolute_values)
        elif exclude == 'stop':
            output = cod.calc_amino_acid_usage(absolute_values=get_absolute_values, exclude_stop=True)
        elif exclude == 'ambiguous':
            output = cod.calc_amino_acid_usage(absolute_values=get_absolute_values, exclude_ambiguous=True)
        elif exclude == 'both':
            output = cod.calc_amino_acid_usage(absolute_values=get_absolute_values, exclude_ambiguous=True, exclude_stop=True)
    else:
        print "Unknown method"
        usage()

    if output_file == None:
        print str(output)
    else:
        f = open(output_file, 'w')
        f.write(str(output))
        f.close()
    exit(0)

