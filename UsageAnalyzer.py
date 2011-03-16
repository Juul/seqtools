from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
import sys, copy
import constants

class UsageAnalyzer:

    # Codon table 11 used for bacterial genomes
    # for yeast use 3
    # for other values see http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
    def __init__(self, genbank_file, codon_table_id=11):

        f = open(genbank_file, 'r')
        gb_generator = SeqIO.parse(f, 'genbank')
        i = 0
        for gb_record in gb_generator:
            self.gb_record = gb_record
            if i > 0:
                raise Exception("More than one record in genbank file. I don't know how to handle this.")
            i += 1
        f.close()

        self.codon_table_id = codon_table_id
        self.codon_table_amb = CodonTable.ambiguous_dna_by_id[codon_table_id]
        self.codon_table = CodonTable.unambiguous_dna_by_id[codon_table_id]
        self.coding_sequences = None


    def get_usage_with_method(self, method, **kwargs):
        exclude = kwargs.get('exclude', None)

        if method == 'codon_usage':
            return self.calc_codon_usage(**kwargs)
        elif method == 'codon_usage_per_position':
            return self.calc_codon_usage_per_position(**kwargs)
        elif method == 'codon_usage_per_amino_acid':
            return self.calc_codon_usage_per_amino_acid(**kwargs)
        elif method == 'nucleotide_usage_per_position':
            return self.calc_nucleotide_usage_per_position(**kwargs)
        elif method == 'amino_acid_usage':
            if exclude == None:
                return self.calc_amino_acid_usage(**kwargs)
            elif exclude == 'stop':
                kwargs['exclude_stop'] = True
                return self.calc_amino_acid_usage(**kwargs)
            elif exclude == 'ambiguous':
                kwargs['exclude_ambiguous'] = True
                return self.calc_amino_acid_usage(**kwargs)
            elif exclude == 'both':
                kwargs['exclude_stop'] = True
                kwargs['exclude_ambiguous'] = True
                return self.calc_amino_acid_usage(**kwargs)
        else:
            return None


    # this method does not output ambiguous amino acids, except for X
    # all otherwise ambiguous amino acids are categorized as X
    # stop codons are categorized as STOP
    def translate_single_codon(self, codon, **kwargs):

        exclude_stop = kwargs.get('exclude_stop', False)
        exclude_ambiguous = kwargs.get('exclude_ambiguous', False)

        try:
            amino = self.codon_table_amb.forward_table.get(codon)
        except:
            amino = None
        if amino == None:
            if codon in self.codon_table_amb.stop_codons:
                if exclude_stop == True:
                    return None
                else:
                    amino = 'STOP'
            else:
                if exclude_ambiguous == True:
                    return None
                else:
                    amino = 'X'

        return amino

    def get_coding_sequences(self, **kwargs):
        use_cached = kwargs.get('use_cached', False)

        if (self.coding_sequences != None) and (use_cached == True):
            return self.coding_sequences

        start_codons = kwargs.get('start_codons', self.codon_table.start_codons)
        seqs = []
        for feature in self.gb_record.features:
            if feature.type == 'CDS':
                seq = self.gb_record.seq[feature.location.nofuzzy_start:feature.location.nofuzzy_end]
                if feature.strand == -1:
                    seq = seq.reverse_complement()

                seqs.append(seq)
        self.coding_sequences = seqs
        return seqs


    def calc_codon_usage(self, **kwargs):
        seqs = self.get_coding_sequences(**kwargs)
        absolute_values = kwargs.get('absolute_values', False) # return the actual number of occurences of each codon, instead of the ratios
        skip_start_codon = kwargs.get('skip_start_codon', False)
        exclude_stop = kwargs.get('exclude_stop', False)
        exclude_ambiguous = kwargs.get('exclude_ambiguous', False)

        codon_counts = {}

        total_codons = 0

        if skip_start_codon == True:
            start_nucleotide = 3
        else:
            start_nucleotide = 0

        c = 1
        for bseq in seqs:
            seq = str(bseq)
            i = start_nucleotide
            while (i + 3) <= len(seq):
                codon = seq[i:i+3]
                amino = self.translate_single_codon(codon, **kwargs)
                if amino == None:
                    i += 3
                    continue

                total_codons += 1
                if codon in codon_counts:
                    codon_counts[codon] += 1
                else:
                    codon_counts[codon] = 1
                i += 3
                    
            c += 1

        if absolute_values == True:
            return  codon_counts

        codon_ratios = {}
        for key in codon_counts:
            codon_ratios[key] = codon_counts[key] / float(total_codons)
            
        return codon_ratios

    def calc_amino_acid_usage(self, **kwargs):

        absolute_values = kwargs.get('absolute_values', False) # return the actual number of occurences of each codon, instead of the ratios        
        exclude_stop = kwargs.get('exclude_stop', False)
        exclude_ambiguous = kwargs.get('exclude_ambiguous', False)

        amino_counts = {}
        codon_counts = self.calc_codon_usage(absolute_values=True, skip_start_codon=True)

        total = 0
        for codon in codon_counts:

            amino = self.translate_single_codon(codon, **kwargs)
            if amino == None:
                continue

            if amino not in amino_counts:
                amino_counts[amino] = codon_counts[codon]
            else:
                amino_counts[amino] += codon_counts[codon]

            total += codon_counts[codon]

        if absolute_values == True:
            return amino_counts

        for amino in amino_counts:
            amino_counts[amino] = amino_counts[amino] / float(total)

        return amino_counts

    def calc_amino_acid_usage_per_position(self, **kwargs):
        print "not implemented"
        exit(-1)



    def calc_codon_usage_per_amino_acid(self, **kwargs):
        absolute_values = kwargs.get('absolute_values', False) # return the actual number of occurences of each codon, instead of the ratios        

        amino_counts = {}
        codon_counts = self.calc_codon_usage(absolute_values=True)

        for codon in codon_counts:

            amino = self.translate_single_codon(codon)

            if amino not in amino_counts:
                amino_counts[amino] = {}
            if codon in amino_counts[amino]:
                amino_counts[amino][codon] += codon_counts[codon]
            else:
                amino_counts[amino][codon] = codon_counts[codon]
                
        if absolute_values == True:
            return amino_counts

        for amino in amino_counts:
            total = 0
            for codon in amino_counts[amino]:
                total += amino_counts[amino][codon]
            for codon in amino_counts[amino]:
                amino_counts[amino][codon] = amino_counts[amino][codon] / float(total)

        return amino_counts


    def calc_codon_usage_per_amino_acid_per_position(self, **kwargs):
        print "not implemented"
        exit(-1)

    def calc_codon_usage_per_position(self, **kwargs):

        absolute_values = kwargs.get('absolute_values', False) # return the actual number of occurences of each codon, instead of the ratios        
        min_cdss = kwargs.get('min_cdss', 0) # if there are fewer CDSs for a position than this number, then ignore (this can happen if only a few CDSs are very long)
        cdss_for_position = [] # number of CDSs for each position
        codon_usage_for_position = []
        
        seqs = self.get_coding_sequences(**kwargs)

        i = 0
        while True:
            nuc_i = i * 3
            cdss_for_position.append(0)
            codon_usage_for_position.append({})

            print "for position: " + str(i)

            c = 0
            for bseq in seqs:
                seq = str(bseq)

                codon = seq[nuc_i:nuc_i+3]
                if len(codon) < 3:
                    continue

                cdss_for_position[i] += 1
                if codon in codon_usage_for_position[i]:
                    codon_usage_for_position[i][codon] += 1
                else:
                    codon_usage_for_position[i][codon] = 1

                c += 1
                    
            if (cdss_for_position[i] < min_cdss) or (cdss_for_position[i] == 0):

                # remove the last element
                cdss_for_position.pop()
                codon_usage_for_position.pop()
                break

            i += 1


        if absolute_values == True:
            return codon_usage_for_position

        codon_ratios = []
        for i in xrange(len(codon_usage_for_position)):
            codon_ratios.append({})
            for key in codon_usage_for_position[i]:
                codon_ratios[i][key] = codon_usage_for_position[i][key] / float(cdss_for_position[i])

        return codon_ratios


    def calc_nucleotide_usage_per_position(self, **kwargs):

        absolute_values = kwargs.get('absolute_values', False) # return the actual number of occurences of each codon, instead of the ratios        
        min_cdss = kwargs.get('min_cdss', 0) # if there are fewer CDSs for a position than this number, then ignore (this can happen if only a few CDSs are very long)
        cdss_for_position = [] # number of CDSs for each position
        nuc_usage_for_position = []
        
        seqs = self.get_coding_sequences(**kwargs)

        i = 0
        while True:
            cdss_for_position.append(0)
            nuc_usage_for_position.append({})

            print "for position: " + str(i)

            c = 0
            for bseq in seqs:
                seq = str(bseq)

                nuc = seq[i:i+1]
                if len(nuc) < 1:
                    continue

                cdss_for_position[i] += 1
                if nuc in nuc_usage_for_position[i]:
                    nuc_usage_for_position[i][nuc] += 1
                else:
                    nuc_usage_for_position[i][nuc] = 1

                c += 1
                    
            if (cdss_for_position[i] < min_cdss) or (cdss_for_position[i] == 0):

                # remove the last element
                cdss_for_position.pop()
                nuc_usage_for_position.pop()
                break

            i += 1


        if absolute_values == True:
            return nuc_usage_for_position

        nuc_ratios = []
        for i in xrange(len(nuc_usage_for_position)):
            nuc_ratios.append({})
            for key in nuc_usage_for_position[i]:
                nuc_ratios[i][key] = nuc_usage_for_position[i][key] / float(cdss_for_position[i])

        return nuc_ratios
