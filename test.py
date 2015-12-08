from Bio import SeqIO
from feature import *

## Read in the fucking files..
records = list(SeqIO.parse("vectors-100.gb", "genbank"))
specialTranslatedFeatures = list(SeqIO.parse("tags_epitopes.mfasta", "fasta"))
primerBindingSites = list (SeqIO.parse("common_primer.mfasta", "fasta"))



## test output
print("vector seq..")
for i in range(0,5):
     print records[i].seq

print("")
print("PBS seq..")
for i in range (0,5):
    print primerBindingSites[i].seq

print("")
print("STF seq..")
for i in range(0,5):
    print specialTranslatedFeatures[i].seq
