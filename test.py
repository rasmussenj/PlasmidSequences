from Bio import SeqIO
from Bio import SeqUtils
from Bio.SeqUtils import six_frame_translations

records = list(SeqIO.parse("vectors.gb", "genbank"))
primerBindingSites = list(SeqIO.parse("common_primer.mfasta", "fasta"))
specialTranslatedFeatures = list(SeqIO.parse("tags_epitopes.mfasta", "fasta"))


for x in range(500, 550):
    seqRecordToCheck = str(records[x].seq)
    seqRecordToCheckComplement = str(records[x].seq.complement())

    for seq in primerBindingSites:
        primerSeq = str(seq.seq).split("(")
        primerSeq = primerSeq[0][::-1]

        partialPrimerSeq = primerSeq[0:14]


        matchingPrimerPositions = SeqUtils.nt_search(seqRecordToCheck, partialPrimerSeq)

        if (len(matchingPrimerPositions) > 1):
            length = len(matchingPrimerPositions)
            for j in range(1, length):

                if primerSeq == seqRecordToCheck[matchingPrimerPositions[j] : matchingPrimerPositions[j] + len(primerSeq)]:
                    print "Complete Match at position: " + str(matchingPrimerPositions[j]) + ".." + str(matchingPrimerPositions[j] + len(primerSeq)) + " in record: " + str(records[x].id)
                else:
                    print "Partial Match at position: " + str(matchingPrimerPositions[j]) + ".." + str(matchingPrimerPositions[j] + len(partialPrimerSeq)) + " in record: " + str(records[x].id)


        if (len(matchingPrimerPositions) > 1):
            length = len(matchingPrimerPositions)
            for j in range(1, length):

                if primerSeq == seqRecordToCheckComplement[matchingPrimerPositions[j] : matchingPrimerPositions[j] + len(primerSeq)]:
                    print "Complete Match on COMPLEMENT at position: " + str(matchingPrimerPositions[j]) + ".." + str(matchingPrimerPositions[j] + len(primerSeq)) + " in record: " + str(records[x].id)
                else:
                    print "Partial Match on COMPLEMENT at position: " + str(matchingPrimerPositions[j]) + ".." + str(matchingPrimerPositions[j] + len(partialPrimerSeq)) + " in record: " + str(records[x].id)
