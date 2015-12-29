import re
from Bio import SeqIO
from Bio.Seq import reverse_complement, translate


records = list(SeqIO.parse("EcoliK12.gb", "genbank"))
specialTranslatedFeatures = list(SeqIO.parse("tags_epitopes.mfasta", "fasta"))

for x in range(len(records)):

    difference = len(records[x].seq) % 3

    if difference != 0:
        seqRecordToCheck = str(records[x].seq)[:-difference]
    else:
        seqRecordToCheck = str(records[x].seq)

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
        secondFrameComplementMatches = re.finditer(featureSeq, secondReadingFrameComplement)
        thirdFrameComplementMatches = re.finditer(featureSeq, thirdReadingFrameComplement)


        for m in firstFrameMatches:
            print featureName + " Matches in first reading frame at position " + str(m.start()) + ".." + str(m.end()) + " in record: " + str(records[x].id)
        for m in secondFrameMatches:
            print featureName + " Matches in second reading frame at position " + str(m.start()) + ".." + str(m.end()) + " in record: " + str(records[x].id)
        for m in thirdFrameMatches:
            print featureName + " Matches in third reading frame at position " + str(m.start()) + ".." + str(m.end()) + " in record: " + str(records[x].id)


        for m in firstFrameComplementMatches:
            print featureName + " Matches in first reading frame COMPLEMENT at position " + str(m.start()) + ".." + str(m.end()) + " in record: " + str(records[x].id)
        for m in secondFrameComplementMatches:
            print featureName + " Matches in second reading frame COMPLEMENT at position " + str(m.start()) + ".." + str(m.end()) + " in record: " + str(records[x].id)
        for m in thirdFrameComplementMatches:
            print featureName + " Matches in third reading frame COMPLEMENT at position " + str(m.start()) + ".." + str(m.end()) + " in record: " + str(records[x].id)
