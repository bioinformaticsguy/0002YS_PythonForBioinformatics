from random import randint
def random_base(RNAflag = True):
	return ("UCAG" if RNAflag else "TCAG")[randint(0,3)]

print(random_base())

def random_codon(RNAflag = False):
	return  random_base(RNAflag) + \
			random_base(RNAflag) + \
			random_base(RNAflag)


print(random_codon())

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


print(replace_base_randomly_using_names("ACAGA"))

def replace_base_randomly_using_expression(base_seq):
	position = randint(0, len(base_seq)-1)
	return (base_seq[0:position] + "TCAG".replace(base_seq[position], "")[randint(0,2)]+ base_seq[position+1:])


print(replace_base_randomly_using_expression("ACAGA"))


def replace_base_randomly(base_seq):
	position = randint(0, len(base_seq)-1)
	bases = 'TCAG'.replace(base_seq[position], '')
	return (base_seq[0:position] + bases[randint(0,2)] + base_seq[position+1:])

print(replace_base_randomly("ACAGA"))


from Bio.Seq import Seq
my_seq = Seq("AGTACACTGGT")
print(my_seq)
print(my_seq.complement())
print(my_seq.reverse_complement())






