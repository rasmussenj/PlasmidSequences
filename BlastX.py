from Bio import SeqIO
from Bio.Blast import NCBIWWW


record = SeqIO.read("EcoliK12.gb", "genbank")

qblast_output = NCBIWWW.qblast("blastx", "refseq_protein", record.seq)

##
# nr = Protein = Non‐redundant GenBank CDS transla ons + RefSeq + PDB + SwissProt + PIR + PRF, excluding those in PAT, TSA, and env_nr.
# refseq_protein = Protein = Protein sequences from NCBI Reference Sequence project.
# swissport  = Protein = Last major release of the UniProtKB/SWISS‐PROT protein sequence database (no incremental updates).
# pat = Protein = Proteins from the Patent division of GenBank.
# pdb = Protein = Protein sequences from the 3‐dimensional structure records from the Protein Data Bank.
# env_nr = Protein = Protein sequences translated from the CDS annota on of metagenomic nucleo de sequences.
# tsa_nr = Protein = Protein sequences translated from CDSs annotated on transcriptome shotgun assemblies.
##

text = qblast_output.read()

print text

