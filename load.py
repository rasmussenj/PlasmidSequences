from Bio import SeqIO

record_iterator = SeqIO.parse("vectors-100.gb", "genbank")
first_record = next(record_iterator)
print(first_record)
print("--------------------------")
print(next(record_iterator))