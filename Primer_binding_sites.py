from Bio import SeqIO
from Bio import SeqUtils
from Bio.Seq import *

records = list(SeqIO.parse("vectors.gb", "genbank"))
primerBindingSites = list(SeqIO.parse("common_primer.mfasta", "fasta"))


for x in range(450, 500):
    seqRecordToCheck = str(records[x].seq)
    seqRecordToCheckComplement = str(reverse_complement(records[x].seq))

    for seq in primerBindingSites:
        primerName = str(seq.name)
        primerSeq = str(seq.seq).split("(")
        primerSeq = primerSeq[0]

        partialPrimerSeq = primerSeq[len(primerSeq)-15::]

        matchingPrimerPositions = SeqUtils.nt_search(seqRecordToCheck, partialPrimerSeq)

        if (len(matchingPrimerPositions) > 1):
            difference = len(primerSeq) - len(partialPrimerSeq)
            length = len(matchingPrimerPositions)
            for j in range(1, length):
                if primerSeq == seqRecordToCheck[matchingPrimerPositions[j] - difference : matchingPrimerPositions[j] - difference + len(primerSeq)]:
                    print "Primer: " + primerName + " Complete Match at position: " + str(matchingPrimerPositions[j] - difference) + ".." + str(matchingPrimerPositions[j] - difference + len(primerSeq)) + " in record: " + str(records[x].id)
                else:
                    print "Primer: " + primerName + " Partial Match at position: " + str(matchingPrimerPositions[j]) + ".." + str(matchingPrimerPositions[j] + len(partialPrimerSeq)) + " in record: " + str(records[x].id)


        matchingPrimerPositions = SeqUtils.nt_search(seqRecordToCheckComplement, partialPrimerSeq)

        if (len(matchingPrimerPositions) > 1):
            difference = len(primerSeq) - len(partialPrimerSeq)
            length = len(matchingPrimerPositions)
            for j in range(1, length):
                if primerSeq == seqRecordToCheckComplement[matchingPrimerPositions[j] - difference : matchingPrimerPositions[j] - difference + len(primerSeq)]:
                    print "Primer: " + primerName + " Complete Match on COMPLEMENT at position: " + str(matchingPrimerPositions[j] - difference) + ".." + str(matchingPrimerPositions[j] - difference + len(primerSeq)) + " in record: " + str(records[x].id)
                else:
                    print "Primer: " + primerName + " Partial Match on COMPLEMENT at position: " + str(matchingPrimerPositions[j]) + ".." + str(matchingPrimerPositions[j] + len(partialPrimerSeq)) + " in record: " + str(records[x].id)

