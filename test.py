from Bio import SeqIO
from feature import *

## Read in the fucking files..
records = list(SeqIO.parse("vectors.gb", "genbank"))
primerBindingSites = list(SeqIO.parse("common_primer.mfasta", "fasta"))
specialTranslatedFeatures = list(SeqIO.parse("tags_epitopes.mfasta", "fasta"))

def searchPrimerBindingSites(records, primerBindingSites):
    for i in range (730, 780):
        recordsSeqToCheckFiveEnd = str(records[i].seq)[0:14]
        recordsSeqToCheckThreeEnd = (str(records[i].seq)[::-1])[0:14]

        complementRecordsSeqToCheckFiveEnd = str(records[i].seq.complement())[0:14]
        complementRecordsSeqToCheckThreeEnd = (str(records[i].seq.complement())[::-1])[0:14]

        for x in range (0, len(primerBindingSites)):
            if recordsSeqToCheckFiveEnd == str(primerBindingSites[x].seq)[0:14]:
                print "Primer Binding Sites: Match at five' end"
                print ("Record ID: " + records[i].id)

            if recordsSeqToCheckThreeEnd == str(primerBindingSites[x].seq)[0:14]:
                print "Primer Binding Sites: Match at three' end"
                print ("Record ID: " + records[i].id)

            if complementRecordsSeqToCheckFiveEnd == str(primerBindingSites[x].seq)[0:14]:
                print "Primer Binding Sites: Match at COMPLEMENT five' end"
                print ("Record ID: " + records[i].id)

            if complementRecordsSeqToCheckThreeEnd == str(primerBindingSites[x].seq)[0:14]:
                print "Primer Binding Sites: Match at COMPLEMENT TREE' end"
                print ("Record ID: " + records[i].id)


searchPrimerBindingSites(records, primerBindingSites)
