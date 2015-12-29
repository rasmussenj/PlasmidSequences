import re
from Bio.Seq import reverse_complement, translate

from load import *

from featuredic import *

from Bio import SeqIO
from Bio import SeqUtils
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.SeqFeature import *
import pickle
import random

file_Name = "featureObjects"

def generate():
    features_Container = getFeature()
    featureDictionary = FeatureDic(features_Container)
    featureDictionary.appendSpecialTransFeature("tags_epitopes.mfasta", "fasta")
    featureDictionary.appendPrimer("common_primer.mfasta", "fasta")
    return featureDictionary.featureDictionary

def write(featureStatistic_container):
    # open the file for writing
    fileObject = open(file_Name,'wb')
    pickle.dump(featureStatistic_container,fileObject)
    fileObject.close()


def read():
    fileObject = open(file_Name,'r')
    # load the object from the file into var b
    return pickle.load(fileObject)


def addFeatureSTF():
    if m.end() <= seqLength:
        newFeature = SeqFeature(FeatureLocation(m.start(),m.end(), strand=1), type=str(feature))
        newFeature.qualifiers['note'] = featureName
        newRecord.features.append(newFeature)
    else:
        newFeature = SeqFeature(CompoundLocation([FeatureLocation(m.start(),seqLength, strand=1), FeatureLocation(1,(seqLength - m.end()),strand=1)]), type=str(feature))
        newFeature.qualifiers['note'] = featureName
        newRecord.features.append(newFeature)


def addFeatureComplSTF():
    if m.end() <= seqLength:
        newFeature = SeqFeature(FeatureLocation(m.start(),m.end(), strand=-1), type=str(feature))
        newFeature.qualifiers['note'] = featureName
        newRecord.features.append(newFeature)
    else:
        newFeature = SeqFeature(CompoundLocation([FeatureLocation(m.start(),seqLength, strand=-1), FeatureLocation(1, (seqLength - m.end()), strand=-1)]), type=str(feature))
        newFeature.qualifiers['note'] = featureName
        newRecord.features.append(newFeature)


def writeSTF():
    global difference, seqRecordToCheck, seqRecordToCheckComplement, variation, featureName, featureSeq, seqLength, m
    difference = len(record.seq) % 3
    seqRecordToCheck = str(record.seq)
    if difference != 0:
        seqRecordToCheck = str(record.seq)[:-difference]
    else:
        seqRecordToCheck = str(record.seq)
    seqRecordToCheckComplement = str(reverse_complement(seqRecordToCheck))
    # Reading Frames
    firstReadingFrame = translate(seqRecordToCheck)
    secondReadingFrame = translate(seqRecordToCheck[1::] + seqRecordToCheck[0])
    thirdReadingFrame = translate(seqRecordToCheck[2::] + seqRecordToCheck[0:2])
    # Reading Frames (reverseComplement)
    firstReadingFrameComplement = translate(seqRecordToCheckComplement)
    secondReadingFrameComplement = translate(seqRecordToCheckComplement[1::] + seqRecordToCheckComplement[0])
    thirdReadingFrameComplement = translate(seqRecordToCheckComplement[2::] + seqRecordToCheckComplement[0:2])
    for variation in featureStatistic_container[feature]:
        featureName = variation.note
        featureSeq = str(variation.seq)
        featureLength = len(variation.seq)
        seqLength = len(seqRecordToCheck)

        firstReadingFrameCircular = firstReadingFrame + firstReadingFrame[0:featureLength - 1]
        secondReadingFrameCircular = secondReadingFrame + secondReadingFrame[0:featureLength - 1]
        thirdReadingFrameCircular = thirdReadingFrame + thirdReadingFrame[0:featureLength - 1]

        firstReadingFrameComplementCircular = firstReadingFrameComplement + firstReadingFrameComplement[
                                                                            0:featureLength - 1]
        secondReadingFrameComplementCircular = secondReadingFrameComplement + secondReadingFrameComplement[
                                                                         0:featureLength - 1]
        thirdReadingFrameComplementCircular = thirdReadingFrameComplement + thirdReadingFrameComplement[
                                                                            0:featureLength - 1]

        # Find Matches
        firstFrameMatchesCircular = re.finditer(featureSeq, firstReadingFrameCircular)
        secondFrameMatchesCircular = re.finditer(featureSeq, secondReadingFrameCircular)
        thirdFrameMatchesCircular = re.finditer(featureSeq, thirdReadingFrameCircular)

        firstFrameComplementMatchesCircular = re.finditer(featureSeq, firstReadingFrameComplementCircular)
        secondFrameComplementMatchesCircular = re.finditer(featureSeq, secondReadingFrameComplementCircular)
        thirdFrameComplementMatchesCircular = re.finditer(featureSeq, thirdReadingFrameComplementCircular)

        for m in firstFrameMatchesCircular:
            addFeatureSTF()

        for m in secondFrameMatchesCircular:
            addFeatureSTF()

        for m in thirdFrameMatchesCircular:
            addFeatureSTF()

        for m in firstFrameComplementMatchesCircular:
            addFeatureComplSTF()

        for m in secondFrameComplementMatchesCircular:
            addFeatureComplSTF()

        for m in thirdFrameComplementMatchesCircular:
            addFeatureComplSTF()


def writePBS():
    global variation, seqRecordToCheck, seqRecordToCheckComplement, difference, newFeature
    for variation in featureStatistic_container[feature]:
        primerSeq = str(variation.seq)
        primerName = variation.note

        partialPrimerSeq = primerSeq[len(primerSeq) - 15::]
        seqRecordToCheck = str(record.seq)
        seqRecordToCheckComplement = str(reverse_complement(record.seq))

        matchingPrimerPositions = SeqUtils.nt_search(seqRecordToCheck, partialPrimerSeq)

        if (len(matchingPrimerPositions) > 1):
            difference = len(primerSeq) - len(partialPrimerSeq)
            length = len(matchingPrimerPositions)
            for j in range(1, length):
                if primerSeq == seqRecordToCheck[matchingPrimerPositions[j] -
                        difference: matchingPrimerPositions[j] - difference + len(primerSeq)]:
                    newFeature = SeqFeature(FeatureLocation(matchingPrimerPositions[j] - difference,
                                                            matchingPrimerPositions[j] - difference + len(primerSeq),
                                                            strand=1), type=str(feature))
                    newFeature.qualifiers['note'] = primerName
                    newRecord.features.append(newFeature)

                else:
                    newFeature = SeqFeature(FeatureLocation(matchingPrimerPositions[j] - difference, AfterPosition(
                        matchingPrimerPositions[j] - difference + len(primerSeq)), strand=1), type=str(feature))
                    newFeature.qualifiers['note'] = primerName
                    newRecord.features.append(newFeature)

        matchingPrimerPositions = SeqUtils.nt_search(seqRecordToCheckComplement, partialPrimerSeq)

        if (len(matchingPrimerPositions) > 1):
            difference = len(primerSeq) - len(partialPrimerSeq)
            length = len(matchingPrimerPositions)
            for j in range(1, length):
                if primerSeq == seqRecordToCheckComplement[matchingPrimerPositions[j] -
                        difference: matchingPrimerPositions[j] - difference + len(primerSeq)]:
                    newFeature = SeqFeature(FeatureLocation(matchingPrimerPositions[j] - difference,
                                                            matchingPrimerPositions[j] - difference + len(primerSeq),
                                                            strand=-1), type=str(feature))
                    newFeature.qualifiers['note'] = primerName
                    newRecord.features.append(newFeature)
                else:
                    newFeature = SeqFeature(FeatureLocation(matchingPrimerPositions[j] - difference, AfterPosition(
                        matchingPrimerPositions[j] - difference + len(primerSeq)), strand=-1), type=str(feature))
                    newFeature.qualifiers['note'] = primerName
                    newRecord.features.append(newFeature)


if __name__ == "__main__":
    featureStatistic_container = generate()
    write(featureStatistic_container)
    #featureStatistic_container = read()
    # record = SeqIO.read("nanobody.fasta", "fasta")
    output_handle = open("outbput.gb", "w")
    records = list(SeqIO.parse("vectors.gb", "genbank"))
    for i in range(1,50):
        record = records[random.randint(1, len(records))]
        newRecord = SeqRecord(record.seq)

        #writing Header
        newRecord.seq.alphabet = generic_dna
        newRecord.id=record.id
        newRecord.name = record.name
        newRecord.description = record.description

        recordSeq = str(record.seq)

        for feature in featureStatistic_container:
            if feature not in ["PBS", "STF"]:
                for variation in featureStatistic_container[feature]:
                    featureSeq = str(variation.seq)
                    occurrence = SeqUtils.nt_search(recordSeq, featureSeq)
                    if (len(occurrence) > 1):
                        for i in range(1,len(occurrence)):
                            newFeature = SeqFeature(FeatureLocation(occurrence[i],occurrence[i]+len(featureSeq), strand=1), type=str(feature))
                            if variation.product != None:
                                newFeature.qualifiers['product'] = variation.product
                            if variation.gene != None:
                                newFeature.qualifiers['gene'] = variation.gene
                            if variation.bound_moiety != None:
                                newFeature.qualifiers['bound_moiety'] = variation.bound_moiety
                            if variation.mobile != None:
                                newFeature.qualifiers['mobile'] = variation.mobile
                            if variation.note != None:
                                newFeature.qualifiers['note'] = variation.note
                            newRecord.features.append(newFeature)

                    featureSeqComplement = str(variation.seq.complement())
                    occurrenceComplement = SeqUtils.nt_search(recordSeq, featureSeqComplement)
                    if (len(occurrenceComplement) > 1):
                        for i in range(1, len(occurrenceComplement)):
                            newFeature = SeqFeature(FeatureLocation(occurrenceComplement[i],occurrenceComplement[i]+len(featureSeq), strand=-1), type=str(feature))
                            if variation.product != None:
                                newFeature.qualifiers['product'] = variation.product
                            if variation.gene != None:
                                newFeature.qualifiers['gene'] = variation.gene
                            if variation.bound_moiety != None:
                                newFeature.qualifiers['bound_moiety'] = variation.bound_moiety
                            if variation.mobile != None:
                                newFeature.qualifiers['mobile'] = variation.mobile
                            if variation.note != None:
                                newFeature.qualifiers['note'] = variation.note
                            newRecord.features.append(newFeature)
            else:
                if(feature == "STF"):
                    writeSTF()

                if(feature == "PBS"):
                    writePBS()

        SeqIO.write(newRecord, output_handle, "genbank")
    output_handle.close()
