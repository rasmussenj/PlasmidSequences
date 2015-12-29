from Bio import SeqIO
from Bio.Blast import NCBIWWW


#record = SeqIO.read("EcoliK12.gb", "genbank")
records = list(SeqIO.parse("vectors.gb", "genbank"))
seq = records[2].seq


qblast_output = NCBIWWW.qblast("blastx", "refseq_protein", seq)


text = qblast_output.read()

print text
