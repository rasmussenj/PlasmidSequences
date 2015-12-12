from Bio import SeqIO
from Bio.SeqUtils import six_frame_translations

records = list(SeqIO.parse("vectors.gb", "genbank"))
primerBindingSites = list(SeqIO.parse("common_primer.mfasta", "fasta"))
specialTranslatedFeatures = list(SeqIO.parse("tags_epitopes.mfasta", "fasta"))

for i in range(len(records)):
    seqRecordToCheck = str(records[i].seq)
    seqRecordToCheckComplement = str(records[i].seq.complement())

    for seq in primerBindingSites:
        primerSeq = str(seq.seq[::-1])
        length = len(seq.seq)

        if seqRecordToCheck[0:14] == primerSeq[0:14]:
            if seqRecordToCheck[0:length] == primerSeq:
                print "Record ID: " + records[i].id + " Primer Binding Sites: (complete) Match. " + "Primer ID: " + seq.id
            else:
                print "Record ID: " + records[i].id + " Primer Binding Sites: (partaial) Match. " + "Primer ID: " + seq.id

        if seqRecordToCheckComplement[0:14] == primerSeq[0:14]:
            if seqRecordToCheckComplement[0:length] == primerSeq:
                print "Record ID: " + records[i].id + " Primer Binding Sites: (complete) Match at complement. " + "Primer ID: " + seq.id
            else:
                print "Record ID: " + records[i].id + " Primer Binding Sites: (partaial) Match at complement. " + "Primer ID: " + seq.id