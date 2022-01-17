print(set("ATGC"))
print({'TCAG'})
print({'TCAG', 'UCAG'})

DNABases = {'T', 'C', 'A', 'G'}
RNABases = {'U', 'C', 'A', 'G'}

DNAbases = {'T', 'C', 'A', 'G'}
RNAbases = {'U', 'C', 'A', 'G'}

print(DNABases, RNABases)


def validate_base_sequence(base_sequence, RNAflag = False): 
    """Return True if the string base_sequence contains only upper
    or lowercase T (or U, if RNAflag), C, A, and G characters, 
    otherwise False"""
    return set(base_sequence) <= (RNAbases if RNAflag else DNAbases)


print(validate_base_sequence("AGU", True))

ali@ali-Inspiron-5570:~$ python3
Python 3.8.5 (default, Jan 27 2021, 15:41:15) 
[GCC 9.3.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> range(5)
range(0, 5)
>>> set(range(5))
{0, 1, 2, 3, 4}
>>> set(range(5,10))
{5, 6, 7, 8, 9}
>>> set(range(5,10))
{5, 6, 7, 8, 9}
>>> range(5,10,2)
range(5, 10, 2)
>>> set(range(5,10,2))
{9, 5, 7}
>>> set(range(15,10,-2))
{11, 13, 15}
>>> set(range(0,-25,-5))
{0, -20, -15, -10, -5}
>>> ('TCAG', 'UCAG')
('TCAG', 'UCAG')
>>> ('TCAG',)
('TCAG',)
>>> ()
()
>>> ('TCAG')
'TCAG'
>>> tuple('TCAG')
('T', 'C', 'A', 'G')
>>> list1 = [1,2,3]
>>> list2 = [4,5]
>>> list1 + list2
[1, 2, 3, 4, 5]
>>> list1
[1, 2, 3]
>>> list1.extend(list2)
>>> list1
[1, 2, 3, 4, 5]
>>> list2
[4, 5]
>>> {'A': 'adenine', 'C': 'cytosine', 'G': 'guanine', 'T': 'thymine'}
{'A': 'adenine', 'C': 'cytosine', 'G': 'guanine', 'T': 'thymine'}
>>> dict((('A', 'adenine'),
...  ('C', 'cytosine'),
...  ('G', 'guanine'),
...  ('T', 'thymine')
...  ))
{'A': 'adenine', 'C': 'cytosine', 'G': 'guanine', 'T': 'thymine'}
>>> 
