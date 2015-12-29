from Bio import SeqIO
from Bio.Blast import NCBIWWW


#record = SeqIO.read("EcoliK12.gb", "genbank")
records = list(SeqIO.parse("vectors.gb", "genbank"))
seq = records[0].seq
parser = NCBIWWW.BlastParser()

qblast_output = NCBIWWW.qblast("blastx", "refseq_protein", seq)
result = parser.parse(qblast_output)

text = qblast_output.read()

print text
