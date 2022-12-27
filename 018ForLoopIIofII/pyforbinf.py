from random import randint
from os.path import dirname, join


def validate_base_sequence(base_sequence, RNAflag=False):
    """Return True if the string base_sequence contains only 
    upper - or lowercase T, C, A, and G characters, 
    otherwise False"""
    seq = base_sequence.upper()
    return len(seq) == (seq.count('U' if RNAflag else 'T') 
                     + seq.count('C') 
                     + seq.count('A') 
                     + seq.count('G'))


def recognition_site(base_sequence, recognition_seq):
	return base_sequence.find(recognition_seq)


def gc_content(base_sequence):
    """Return the percentage of G and C characters in base_seq"""
    assert validate_base_sequence(base_sequence), \
            'argument has invalid characters'
    seq = base_sequence.upper()
    return ((base_sequence.count('G') + base_sequence.count('C')) / len(base_sequence))



def random_base(RNAflag = True):
	return ("UCAG" if RNAflag else "TCAG")[randint(0,3)]


def random_codon(RNAflag = False):
	return  random_base(RNAflag) + \
			random_base(RNAflag) + \
			random_base(RNAflag)

def replace_base_randomly_using_names(base_sequence):
	""" Return a sequence with the base at a randomly
	selected position of base_seq replaced by a base 
	chosen randomly from the three bases that are not 
	at that position """
	position = randint(0, len(base_sequence)-1) # -1 because len is one past end
	base = base_sequence[position]
	bases = "TCAG"
	base.replace(base,"") #replace with empty string!
	newbase = bases[randint(0,2)]
	beginning = base_sequence[0:position] # up to position
	end = base_sequence[position+1:] # omitting the base at position

	return beginning + newbase + end

def replace_base_randomly_using_expression(base_sequence):
	position = randint(0, len(base_sequence)-1)
	return (base_sequence[0:position] + "TCAG".replace(base_sequence[position], "")[randint(0,2)]+ base_sequence[position+1:])

def replace_base_randomly(base_sequence):
	position = randint(0, len(base_sequence)-1)
	bases = 'TCAG'.replace(base_sequence[position], '')
	return (base_sequence[0:position] + bases[randint(0,2)] + base_sequence[position+1:])


def validate_base_sequence_using_set(base_sequence, RNAflag = False): 
    """Return True if the string base_sequence contains only upper
    or lowercase T (or U, if RNAflag), C, A, and G characters, 
    otherwise False"""

    DNAbases = {'T', 'C', 'A', 'G'}
    RNAbases = {'U', 'C', 'A', 'G'}

    return set(base_sequence) <= (RNAbases if RNAflag else DNAbases)




# 

DNABases=dict((('A', 'adenine'),
               ('C', 'cytosine'),
               ('G', 'guanine'),
               ('T', 'thymine')))



def translate_RNA_codon(codon):
    """ returns the amino acid for the given codon """

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
    'GUG': 'Val', 'GCG': 'Ala', 'GAG': 'Glu', 'GGG': 'Gly'  # GxG
    }

    return RNA_codon_table[codon]

def random_codons(minlength = 3, maxlength = 10, RNAflag = False):
    """ Generate a random list of codons (RNA if RNAflag, else DNA)
    between minimum and maximum length, inclusive. 
    """
    return [random_codon(RNAflag) for n in range(randint(minlength, maxlength))]

def random_codons_translation(minlength = 3, maxlength = 10):
    """ Generates a random list of amino acids between minimum and 
    maximum length inclusive. Then returns the translation of each 
    codon in the form of a list.
    """
    return [translate_RNA_codon(codon) for codon in random_codons(minlength, maxlength, True)]

def read_FASTA_strings(filename):
    current_dir = dirname(__file__)
    with open(current_dir + "./" + filename) as file:
        return file.read().split(">")[1:]
        
def read_FASTA_entries(filename):
    return [seq.partition("\n") for seq in read_FASTA_strings(filename)]

def read_FASTA_sequences(filename):
    return [[seq[0],
            seq[2].replace("\n", "")]
            for seq in read_FASTA_entries(filename)]


def make_indexed_sequence_dictionary(filename):
    return {info:seq for info, seq in read_FASTA_sequences(filename)}


def getFirstSequenceFromFastaFile(filename):
    '''given the name of a FASTA file, read and return its first sequence, 
    ignoring the sequence description '''
    seq = ''

    with open(filename) as file:
        line = file.readline()
        while line and line[0] == ">":
            line = file.readline()
        while line and line[0] != '>':
            seq += line
            line = file.readline()

    return seq 


def readFastaIteration(filename):
    '''Given a fasta file returns a list of all the touples
    in which the first element of touple is the description of the fasta sequence 
    and the 2nd element is the sequence itself.'''
    sequences = []
    descrip = None
    with open(filename) as file:
        for line in file:
            if line[0] == ">":
                if descrip:
                    sequences.append((descrip, seq))
                descrip = line[1:-1].split('|')
                seq = '' 
            else:
                seq = seq + line[:-1]
        sequences.append((descrip, seq))
    
    return sequences

def read_FASTA(filename):
    with open(filename) as file:
        return [(part[0].split('|'),
                 part[2].replace('\n', ''))
                    for part in 
                         [entry.partition('\n') for entry in file.read().split('>')[1:]]]


def longest_sequence(filename):
    longest_seq = ''
    for info, seq in read_FASTA(filename):
        longest_seq = max(longest_seq, seq, key=len)
    return longest_seq


def extract_gi_id(description):
    '''Given a FASTA file description line, return its GenInfo ID if it has one'''
    if description[0] != '>':
        return None
    fields = description[1:].split('|')
    if 'gi' not in fields:
        return None
    return fields[1 + fields.index('gi')]

def get_gi_ids(filename):
    """Return a list of GenInfo IDs of all sequences in the file names filename"""
    with open(filename) as file:
        return [extract_gi_id(line) for line in file if line[0] == '>']

def get_gi_ids_from_files(filenames):
    """Return a list of the GenInfo ID's of all sequences formed in the
    files whose names are contained in the collection filenames."""
    idlst = []
    for filename in filenames:
        idlst += get_gi_ids(filename)
    return idlst

def search_FASTA_file_by_gi_id(id, filename):
    """Return the sequence with the GenInfo ID ID from the FASTA file
    named filename, reading one entry at a time until it is found"""
    id = str(id) # user might call with a number
    with open(filename) as file:
        return FASTA_search_by_gi_id(id, file) 


def FASTA_search_by_gi_id(id, file):
    for line in file:
        if (line[0] == '>' and
            str(id) == get_gi_id(line)):
            return read_FASTA_sequence(file)

def read_FASTA_sequence(file):
    seq = '' 
    for line in file:
        if not line or line[0] == '>':
            return seq
        seq += line[:-1]
   

def get_gi_id(description):
    fields = description[1:].split('|')
    if fields and 'gi' in fields:
        return fields[1 + fields.index('gi')]


def print_codon_table():
    """Print the DNA codon table in a nice, but simple, arrangement"""
    DNA_bases = ['A', 'U', 'G', 'C']
    
    for base1 in DNA_bases: # horizontal section (or "group")
        for base3 in DNA_bases: # line (or "row")
            for base2 in DNA_bases: # vertical section (or "column")    
                # the base2 loop is inside the base3 loop!
                print(base1+base2+base3,
                        translate_RNA_codon(base1+base2+base3),
                        end='     ')
            print()
        print()




if __name__ == '__main__':
    pass