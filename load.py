from Bio import SeqIO
from Bio import SeqFeature
from feature import *

features = []

record_iterator = SeqIO.parse("vectors-100.gb", "genbank")
first_record = next(record_iterator)
for f in first_record.features:
    features.append(Feature(first_record.seq[f.location.start:f.location.end], f.type))
    #features.append([1, f.type, first_record.seq[f.location.start:f.location.end]])

