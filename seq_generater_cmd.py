#!/usr/bin/env python2.6
import sys, math, random, copy
import time
import getopt
from SeqGenerater import *


def usage():

    print 'Usage: ' + __file__ + " [-f genbank_file] [-u usage_file] -m <method> -l <length>"
    print "  -m <method> Specifies the method to run. Explained at the end of this help section. [required]"
    print "  -f <genbank_file> Specifies the genbank file to use as input [required]"
    print "  -u <usage_file> Specifies the json usage file (as output from usage_analyzer_cmd.py)"
    print "  -l <length> Specifies the length of the sequence to generate"
    print "  -d Enable debug"
    print "  -o <filename> Sets an output file"
    print "  -h Print this help screen"
    print " --------------------------------"
    print " Methods (for use with the -m flag):"
    print "   amino_acid_seq_from_usage: Generate an amino acid sequence based on the amino acid usage."
    print "   nucleotide_seq_from_codon_usage: Generate a nucleotide sequence based on the codon usage."
    print "   nucleotide_seq_from_codon_usage_per_position: Generate a nucleotide sequence based on the per-posision codon usage."


    exit(-1)
    
if __name__ == "__main__":


    optlist, args = getopt.getopt(sys.argv[1:], 'f:u:m:l:dho:')

    debug = False
    usage_file = None
    genbank_file = None
    output_file = None
    method = None
    show_help = False
    length = None

    for o, a in optlist:
        if o == '-m':
            debug = True
        if o == '-f':
            genbank_file = a
        if o == '-u':
            usage_file = a
        if o == '-o':
            output_file = a
        if o == '-m':
            method = a
        if o == '-l':
            length = int(a)
        if o == '-h':
            show_help = True

    if (genbank_file == None) and (usage_file == None):
        print "You must specify either a genbank file or a usage file."
        usage()

    if (method == None) or (show_help == True):
        usage()

    sg = SeqGenerator()


    kwargs = {
        'genbank_file': genbank_file,
        'usage_file': usage_file,
        'length': length
        }

    output = sg.gen_seq_with_method(method, **kwargs)

    if output_file == None:
        print str(output)
    else:
        f = open(output_file, 'w')
        f.write(str(output))
        f.close()
    exit(0)

