from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
import random
import sys, copy
import constants

from UsageAnalyzer import *



class SeqGenerator:

    def get_usage_with_method(self, method, **kwargs):

        if method == 'amino_acid_seq_from_usage':
            return self.amino_acid_seq_from_usage(**kwargs)
        elif method == 'nucleotide_seq_from_codon_usage':
            return self.nucleotide_seq_from_codon_usage(**kwargs)
        else:
            raise "Unknown method: " + method


    def usage_from_genbank_or_usage_file(self, method, **kwargs):
        genbank_file = kwargs.get('genbank_file', None)
        usage_file = kwargs.get('usage_file', None)

        if usage_file != None:
            print "not implemented"
            exit(-1)
        elif genbank_file != None:
            ua = UsageAnalyzer(genbank_file)
            return ua.get_usage_with_method(method, **kwargs)
        else:
            raise "Neither usage nor genbank file specified. Cannot get usage from nothin."


    # takes a hash like: {'G': 0.3, 'A': 0.4, 'T': 0.1, 'C': 0.2}
    # and returns either 'G', 'A', 'T' or 'C' with the probabilities listed
    def biased_random_from_ratios(self, ratios):
        a = ratios.items()
        a.sort()

        r = random.random()

        prev_cumul = 0
        cumul = None

        for i in range(len(a)):
            cur = a[i]
            cur_key = a[i][0]
            cur_val = a[i][1]

            cumul = prev_cumul + cur_val

            if (r >= prev_cumul) and (r < cumul):
                return cur_key

            prev_cumul = cumul


        if r >= cumul:
            return a[len(a)-1][0]

        raise "Should never get here"



    def nucleotide_seq_from_codon_usage(self, **kwargs):

        kwargs['skip_start_codon'] = True
        kwargs['exclude_stop'] = True
        kwargs['exclude_ambiguous'] = True
        usage = self.usage_from_genbank_or_usage_file('codon_usage', **kwargs)

        length = kwargs.get('length', None)
        if length == None:
            raise "Cannot create amino acid sequence with unknown length"

        length = length / 3

        seq = ''

        for i in range(length):
            seq += self.biased_random_from_ratios(usage)

        return seq

    def amino_acid_seq_from_usage(self, **kwargs):

        kwargs['exclude'] = 'both'
        usage = self.usage_from_genbank_or_usage_file('amino_acid_usage', **kwargs)
        

        length = kwargs.get('length', None)
        if length == None:
            raise "Cannot create amino acid sequence with unknown length"

        seq = ''

        for i in range(length):
            seq += self.biased_random_from_ratios(usage)

        return seq
        
        
