from Bio import SeqIO, SeqUtils
from Bio.Seq import reverse_complement, translate


records = list(SeqIO.parse("vectors.gb", "genbank"))
specialTranslatedFeatures = list(SeqIO.parse("tags_epitopes.mfasta", "fasta"))

for x in range(399, 400):

    difference = len(records[x].seq) % 3

    if difference != 0:
        seqRecordToCheck = str(records[x].seq)[:-difference]
    else: seqRecordToCheck = str(records[x].seq)

    seqRecordToCheckComplement = str(reverse_complement(seqRecordToCheck))

    #Reading Frames
    firstReadingFrame = translate(seqRecordToCheck)
    secondReadingFrame = translate(seqRecordToCheck[1::] + seqRecordToCheck[0])
    thirdReadingFrame = translate(seqRecordToCheck[2::] + seqRecordToCheck[0:2])

    #Reading Frames (reverseComplement)
    firstReadingFrameComplement = translate(seqRecordToCheckComplement)
    secondReadingFrameComplement = translate(seqRecordToCheckComplement[1::] + seqRecordToCheckComplement[0])
    thirdReadingFrameComplement = translate(seqRecordToCheckComplement[2::] + seqRecordToCheckComplement[0:2])

    for feature in specialTranslatedFeatures:
        featureName = feature.name
        featureSeq = str(feature.seq)
        featureLength = len(feature.seq)

        firstFrameMatches = SeqUtils.nt_search(firstReadingFrame, featureSeq)
        secondFrameMatches = SeqUtils.nt_search(secondReadingFrame, featureSeq)
        thirdFrameMatches = SeqUtils.nt_search(thirdReadingFrame, featureSeq)

        firstFrameComplementMatches = SeqUtils.nt_search(firstReadingFrameComplement, featureSeq)
        secondFrameComplementMatches = SeqUtils.nt_search(secondReadingFrameComplement, featureSeq)
        thirdFrameComplementMatches = SeqUtils.nt_search(thirdReadingFrameComplement, featureSeq)

        if (len(firstFrameMatches) > 1):
            for j in range(1, len(firstFrameMatches)):

                print featureName + " Matches in first reading Frame on position " + str(firstFrameMatches[j]) + ".." + str(firstFrameMatches[j] + featureLength)
