##seq = open("seqdump.txt")
##
##for i in seq:
##    print(i)


##with open("seqdump.txt", "a") as file:
####    print(file.read(11))
####    print(file.readlines())
##    l = ["atgc", "tag", "gcgcgcgcgc"]
##    file.writelines(l)
##
##
##    


def read_FASTA_strings(filename):
    with open(filename) as file:
        return file.read().split(">")[1:]


strings = read_FASTA_strings("seqdump.txt")

print(strings[1])
