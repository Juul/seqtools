##########################################################################
#                                                                        #
#  Copyright 2011 Marc Juul Christoffersen                               #
#                                                                        #
#  This file is part of SeqTools.                                        #
#                                                                        #
#  SeqTools is free software: you can redistribute it and/or modify      #
#  it under the terms of the GNU General Public License as published by  #
#  the Free Software Foundation, either version 3 of the License, or     #
#  (at your option) any later version.                                   #
#                                                                        #
#  SeqTools is distributed in the hope that it will be useful,           #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
#  GNU General Public License for more details.                          #
#                                                                        #
#  You should have received a copy of the GNU General Public License     #
#  along with Foobar.  If not, see <http://www.gnu.org/licenses/>.       #
#                                                                        #
##########################################################################

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
import random
import sys, copy
import constants
import json


from UsageAnalyzer import *



class SeqGenerator:

    def gen_seq_with_method(self, method, **kwargs):

        if method == 'amino_acid_seq_from_usage':
            return self.amino_acid_seq_from_usage(**kwargs)
        elif method == 'nucleotide_seq_from_codon_usage':
            return self.nucleotide_seq_from_codon_usage(**kwargs)
        elif method == 'nucleotide_seq_from_codon_usage_per_position':
            return self.nucleotide_seq_from_codon_usage_per_position(**kwargs)
        else:
            raise "Unknown method: " + method


    def usage_from_genbank_or_usage_file(self, method, **kwargs):
        genbank_file = kwargs.get('genbank_file', None)
        usage_file = kwargs.get('usage_file', None)

        if usage_file != None:
            usage = json.loads("\n".join(open(usage_file, 'r').readlines()))
            return usage
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
            raise "Cannot create nucleotide sequence with unknown length"

        length = length / 3

        seq = ''

        for i in range(length):
            seq += self.biased_random_from_ratios(usage)

        return seq


    def nucleotide_seq_from_codon_usage_per_position(self, **kwargs):

        kwargs['skip_start_codon'] = True
        kwargs['exclude_stop'] = True
        kwargs['exclude_ambiguous'] = True
        usage = self.usage_from_genbank_or_usage_file('codon_usage_per_position', **kwargs)

        length = kwargs.get('length', None)
        if length == None:
            raise "Cannot create nucleotide sequence with unknown length"

        length = length / 3

        seq = ''

        for i in range(length):
            seq += self.biased_random_from_ratios(usage[i])

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
        
        
