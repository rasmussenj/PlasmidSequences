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
        ####
        featureLength = len(feature.seq)
        seqLength = len(seqRecordToCheck)

        firstReadingFrameCircular = firstReadingFrame + firstReadingFrame[0:featureLength-1]
        secondReadingFrameCircular = secondReadingFrame + secondReadingFrame[0:featureLength-1]
        thirdReadingFrameCircular = thirdReadingFrame + thirdReadingFrame[0:featureLength-1]

        firstReadingFrameComplementCircular = firstReadingFrameComplement + firstReadingFrameComplement[0:featureLength-1]
        secondReadingFrameComplementCircular = secondReadingFrameComplement + secondReadingFrameComplement[0:featureLength-1]
        thirdReadingFrameComplementCircular = thirdReadingFrameComplement + thirdReadingFrameComplement[0:featureLength-1]

        #Find Matches
        firstFrameMatchesCircular = re.finditer(featureSeq, firstReadingFrameCircular)
        secondFrameMatchesCircular = re.finditer(featureSeq, secondReadingFrameCircular)
        thirdFrameMatchesCircular = re.finditer(featureSeq, thirdReadingFrameCircular)

        firstFrameComplementMatchesCircular = re.finditer(featureSeq, firstReadingFrameComplementCircular)
        secondFrameComplementMatchesCircular = re.finditer(featureSeq, secondReadingFrameComplementCircular)
        thirdFrameComplementMatchesCircular = re.finditer(featureSeq, thirdReadingFrameComplementCircular)

        for m in firstFrameMatchesCircular:
            if m.end() <= seqLength:
                print featureName + " Matches in first reading frame at position " + str(m.start()) + ".." + \
                      str(m.end()) + " in record: " + str(records[x].id)
            else:
                print featureName + " Matches in first reading frame at position " + str(m.start()) + ".." + \
                      str(seqLength) + " Join " + "1.." + str(seqLength) + seqLength - m.end()

        for m in secondFrameMatchesCircular:
            if m.end() <= seqLength:
                print featureName + " Matches in second reading frame at position " + str(m.start()) + ".." + \
                      str(m.end()) + " in record: " + str(records[x].id)
            else:
                print featureName + " Matches in second reading frame at position " + str(m.start()) + ".." + \
                      str(seqLength) + " Join " + "1.." + str(seqLength) + seqLength - m.end()

        for m in thirdFrameMatchesCircular:
            if m.end() <= seqLength:
                print featureName + " Matches in third reading frame at position " + str(m.start()) + ".." + \
                      str(m.end()) + " in record: " + str(records[x].id)
            else:
                print featureName + " Matches in third reading frame at position " + str(m.start()) + ".." + \
                      str(seqLength) + " Join " + "1.." + str(seqLength) + seqLength - m.end()



        for m in firstFrameComplementMatchesCircular:
            if m.end() <= seqLength:
                print featureName + " Matches in first reading frame COMPLEMENT at position " + str(m.start()) + \
                      ".." + str(m.end()) + " in record: " + str(records[x].id)
            else:
                print featureName + " Matches in first reading frame COMPLEMENT at position " + str(m.start()) + \
                      ".." + str(seqLength) + " Join " + "1.." + str(seqLength) + seqLength - m.end()

        for m in secondFrameComplementMatchesCircular:
            if m.end() <= seqLength:
                print featureName + " Matches in second reading frame COMPLEMENT at position " + str(m.start()) + \
                      ".." + str(m.end()) + " in record: " + str(records[x].id)
            else:
                print featureName + " Matches in second reading frame COMPLEMENT at position " + str(m.start()) + \
                      ".." + str(seqLength) + " Join " + "1.." + str(seqLength) + seqLength - m.end()

        for m in thirdFrameComplementMatchesCircular:
            if m.end() <= seqLength:
                print featureName + " Matches in third reading frame COMPLEMENT at position " + str(m.start()) + \
                      ".." + str(m.end()) + " in record: " + str(records[x].id)
            else:
                print featureName + " Matches in third reading frame COMPLEMENT at position " + str(m.start()) + \
                      ".." + str(seqLength) + " Join " + "1.." + str(seqLength) + seqLength - m.end()