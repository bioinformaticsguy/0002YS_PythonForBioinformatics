from random import randint

def validate_base_sequence(base_sequence, RNAflag=False):
    """Return True if the string base_sequence contains only 
    upper - or lowercase T, C, A, and G characters, 
    otherwise False"""
    seq = base_sequence.upper()
    return len(seq) == (seq.count('U' if RNAflag else 'T') 
                     + seq.count('C') 
                     + seq.count('A') 
                     + seq.count('G'))


def recognition_site(base_seq, recognition_seq):
	return base_seq.find(recognition_seq)



def gc_content(base_seq):
    """Return the percentage of G and C characters in base_seq"""
    assert validate_base_sequence(base_seq), \
            'argument has invalid characters'
    seq = base_seq.upper()
    return ((base_seq.count('G') + base_seq.count('C')) / len(base_seq))


def random_base(RNAflag = True):
	return ("UCAG" if RNAflag else "TCAG")[randint(0,3)]


def random_codon(RNAflag = False):
	return  random_base(RNAflag) + \
			random_base(RNAflag) + \
			random_base(RNAflag)



def replace_base_randomly_using_names(base_seq):
	""" Return a sequence with the base at a randomly
	selected position of base_seq replaced by a base 
	chosen randomly from the three bases that are not 
	at that position """
	position = randint(0, len(base_seq)-1) # -1 because len is one past end
	base = base_seq[position]
	bases = "TCAG"
	base.replace(base,"") #replace with empty string!
	newbase = bases[randint(0,2)]
	beginning = base_seq[0:position] # up to position
	end = base_seq[position+1:] # omitting the base at position

	return beginning + newbase + end



def replace_base_randomly_using_expression(base_seq):
	position = randint(0, len(base_seq)-1)
	return (base_seq[0:position] + "TCAG".replace(base_seq[position], "")[randint(0,2)]+ base_seq[position+1:])


def replace_base_randomly(base_seq):
	position = randint(0, len(base_seq)-1)
	bases = 'TCAG'.replace(base_seq[position], '')
	return (base_seq[0:position] + bases[randint(0,2)] + base_seq[position+1:])


def validate_base_sequence(base_sequence, RNAflag = False): 
    """Return True if the string base_sequence contains only upper
    or lowercase T (or U, if RNAflag), C, A, and G characters, 
    otherwise False"""
    return set(base_sequence) <= (RNAbases if RNAflag else DNAbases)



RNA_codon_table = {
# Second Base
# U C A G
# U
 'UUU': 'Phe', 'UCU': 'Ser', 'UAU': 'Tyr', 'UGU': 'Cys', # UxU
 'UUC': 'Phe', 'UCC': 'Ser', 'UAC': 'Tyr', 'UGC': 'Cys', # UxC
 'UUA': 'Leu', 'UCA': 'Ser', 'UAA': '---', 'UGA': '---', # UxA
 'UUG': 'Leu', 'UCG': 'Ser', 'UAG': '---', 'UGG': 'Urp', # UxG
# C
 'CUU': 'Leu', 'CCU': 'Pro', 'CAU': 'His', 'CGU': 'Arg', # CxU
 'CUC': 'Leu', 'CCC': 'Pro', 'CAC': 'His', 'CGC': 'Arg', # CxC
 'CUA': 'Leu', 'CCA': 'Pro', 'CAA': 'Gln', 'CGA': 'Arg', # CxA
 'CUG': 'Leu', 'CCG': 'Pro', 'CAG': 'Gln', 'CGG': 'Arg', # CxG
# A
 'AUU': 'Ile', 'ACU': 'Thr', 'AAU': 'Asn', 'AGU': 'Ser', # AxU
 'AUC': 'Ile', 'ACC': 'Thr', 'AAC': 'Asn', 'AGC': 'Ser', # AxC
 'AUA': 'Ile', 'ACA': 'Thr', 'AAA': 'Lys', 'AGA': 'Arg', # AxA
 'AUG': 'Met', 'ACG': 'Thr', 'AAG': 'Lys', 'AGG': 'Arg', # AxG
# G
 'GUU': 'Val', 'GCU': 'Ala', 'GAU': 'Asp', 'GGU': 'Gly', # GxU
 'GUC': 'Val', 'GCC': 'Ala', 'GAC': 'Asp', 'GGC': 'Gly', # GxC
 'GUA': 'Val', 'GCA': 'Ala', 'GAA': 'Glu', 'GGA': 'Gly', # GxA
  'GUG': 'Val', 'GCG': 'Ala', 'GAG': 'Glu', 'GGG': 'Gly' # GxG
}

# print(RNA_codon_table)

DNABases=dict((('A', 'adenine'),
      ('C', 'cytosine'),
      ('G', 'guanine'),
      ('T', 'thymine')
     ))


def read_FASTA_strings(filename):
    with open(filename) as file:
        return file.read().split(">")[1:]


def random_codons(minlength = 3, maxlength = 10, RNAflag = False):
    """ Generate a random list of codons (RNA if RNAflag, else DNA)
    between minimum and maximum length, inclusive. 
    """
    return [random_codon(RNAflag) for n in range(randint(minlength, maxlength))]


def translate_RNA_codon(codon):
    """ returns the amino acid for the given codon """
    return RNA_codon_table[codon]


def random_codons_translation(minlength = 3, maxlength = 10):
    """ Generates a random list of amino acids between minimum and 
    maximum length inclusive.
    """
    return [translate_RNA_codon(codon) for codon in random_codons(minlength, maxlength, True)]


def read_FASTA_strings(filename):
    with open(filename) as file:
        return file.read().split(">")[1:]
        
def read_FASTA_entries(filename):
    return [seq.partition("\n") for seq in read_FASTA_strings(filename)]


def read_FASTA_sequences(filename):
    return [[seq[0],
            seq[2].replace("\n", "")]
            for seq in read_FASTA_entries(filename)]


def make_indexed_sequence_dictionary(filename):
    return {info:seq for info, seq in read_FASTA_sequences(filename)}


