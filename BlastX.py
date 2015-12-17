from Bio import SeqIO
from Bio.Blast import NCBIWWW


#records = list(SeqIO.parse("EcoliK12.gb", "genbank"))
record = SeqIO.read("EcoliK12.gb", "genbank")

qblast_output = NCBIWWW.qblast("blastx", "nt", record.seq)

text = qblast_output.read()

print text

