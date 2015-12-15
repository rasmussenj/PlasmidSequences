import re

from Bio import SeqIO
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

        firstFrameMatches = re.finditer(featureSeq, firstReadingFrame)
        secondFrameMatches = re.finditer(featureSeq, secondReadingFrame)
        thirdFrameMatches = re.finditer(featureSeq, thirdReadingFrame)

        firstFrameComplementMatches = re.finditer(featureSeq, firstReadingFrameComplement)
        secondFrameComplementMatches = re.finditer(featureSeq, firstReadingFrameComplement)
        thirdFrameComplementMatches = re.finditer(featureSeq, firstReadingFrameComplement)


        for m in firstFrameMatches:
            print m.start(), m.end()
        for m in secondFrameMatches:
            print m.start(), m.end()
        for m in thirdFrameMatches:
            print m.start(), m.end()


        for m in firstFrameComplementMatches:
            print m.start(), m.end()
        for m in secondFrameComplementMatches:
            print m.start(), m.end()
        for m in thirdFrameComplementMatches:
            print m.start(), m.end()
