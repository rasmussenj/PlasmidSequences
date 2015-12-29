import re
from Bio.Seq import reverse_complement, translate

from load import *

from featuredic import *

from Bio import SeqIO
from Bio import SeqUtils
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import SeqFeature
import pickle

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

if __name__ == "__main__":
    featureStatistic_container = generate()
    #write(featureStatistic_container)
    #featureStatistic_container = read()
    # record = SeqIO.read("nanobody.fasta", "fasta")
    record = SeqIO.read("vectors-1.gb", "genbank")
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
            if(feature == "SFT"):
                difference = len(record.seq) % 3
                seqRecordToCheck = str(record.seq)
                if difference != 0:
                    seqRecordToCheck = str(record.seq)[:-difference]
                else:
                    seqRecordToCheck = str(record.seq)

                seqRecordToCheckComplement = str(reverse_complement(seqRecordToCheck))

                #Reading Frames
                firstReadingFrame = translate(seqRecordToCheck)
                secondReadingFrame = translate(seqRecordToCheck[1::] + seqRecordToCheck[0])
                thirdReadingFrame = translate(seqRecordToCheck[2::] + seqRecordToCheck[0:2])

                #Reading Frames (reverseComplement)
                firstReadingFrameComplement = translate(seqRecordToCheckComplement)
                secondReadingFrameComplement = translate(seqRecordToCheckComplement[1::] + seqRecordToCheckComplement[0])
                thirdReadingFrameComplement = translate(seqRecordToCheckComplement[2::] + seqRecordToCheckComplement[0:2])

                for variation in featureStatistic_container[feature]:
                    featureName = variation.note
                    featureSeq = str(variation.seq)

                    firstFrameMatches = re.finditer(featureSeq, firstReadingFrame)
                    secondFrameMatches = re.finditer(featureSeq, secondReadingFrame)
                    thirdFrameMatches = re.finditer(featureSeq, thirdReadingFrame)

                    firstFrameComplementMatches = re.finditer(featureSeq, firstReadingFrameComplement)
                    secondFrameComplementMatches = re.finditer(featureSeq, secondReadingFrameComplement)
                    thirdFrameComplementMatches = re.finditer(featureSeq, thirdReadingFrameComplement)


                    for m in firstFrameMatches:
                        print featureName + " Matches in first reading frame at position " + str(m.start()) + ".." + str(m.end()) + " in record: " + str(record.id)

                    for m in secondFrameMatches:
                        print featureName + " Matches in second reading frame at position " + str(m.start()) + ".." + str(m.end()) + " in record: " + str(record.id)
                    for m in thirdFrameMatches:
                        print featureName + " Matches in third reading frame at position " + str(m.start()) + ".." + str(m.end()) + " in record: " + str(record.id)


                    for m in firstFrameComplementMatches:
                        print featureName + " Matches in first reading frame COMPLEMENT at position " + str(m.start()) + ".." + str(m.end()) + " in record: " + str(record.id)
                    for m in secondFrameComplementMatches:
                        print featureName + " Matches in second reading frame COMPLEMENT at position " + str(m.start()) + ".." + str(m.end()) + " in record: " + str(record.id)
                    for m in thirdFrameComplementMatches:
                        print featureName + " Matches in third reading frame COMPLEMENT at position " + str(m.start()) + ".." + str(m.end()) + " in record: " + str(record.id)


            if(feature == "PBS"):
                for variation in featureStatistic_container[feature]:
                    primerSeq = str(variation.seq)
                    primerName = variation.note

                    partialPrimerSeq = primerSeq[len(primerSeq)-15::]
                    seqRecordToCheck = str(record.seq)
                    seqRecordToCheckComplement = str(reverse_complement(record.seq))


                    matchingPrimerPositions = SeqUtils.nt_search(seqRecordToCheck, partialPrimerSeq)

                    if (len(matchingPrimerPositions) > 1):
                        difference = len(primerSeq) - len(partialPrimerSeq)
                        length = len(matchingPrimerPositions)
                        for j in range(1, length):
                            if primerSeq == seqRecordToCheck[matchingPrimerPositions[j] -
                                    difference : matchingPrimerPositions[j] - difference + len(primerSeq)]:
                                print "Primer: " + primerName + " Complete Match at position: " + \
                                      str(matchingPrimerPositions[j] - difference) + ".." + \
                                      str(matchingPrimerPositions[j] - difference + len(primerSeq)) + \
                                      " in record: " + str(record.id)
                            else:
                                print "Primer: " + primerName + " Partial Match at position: " + \
                                      str(matchingPrimerPositions[j]) + ".." + \
                                      str(matchingPrimerPositions[j] + len(partialPrimerSeq)) + \
                                      " in record: " + str(record.id)


                    matchingPrimerPositions = SeqUtils.nt_search(seqRecordToCheckComplement, partialPrimerSeq)

                    if (len(matchingPrimerPositions) > 1):
                        difference = len(primerSeq) - len(partialPrimerSeq)
                        length = len(matchingPrimerPositions)
                        for j in range(1, length):
                            if primerSeq == seqRecordToCheckComplement[matchingPrimerPositions[j] -
                                    difference : matchingPrimerPositions[j] - difference + len(primerSeq)]:
                                print "Primer: " + primerName + " Complete Match on COMPLEMENT at position: " + \
                                      str(matchingPrimerPositions[j] - difference) + ".." + \
                                      str(matchingPrimerPositions[j] - difference + len(primerSeq)) + \
                                      " in record: " + str(record.id)
                            else:
                                print "Primer: " + primerName + " Partial Match on COMPLEMENT at position: " + \
                                      str(matchingPrimerPositions[j]) + ".." + \
                                      str(matchingPrimerPositions[j] + len(partialPrimerSeq)) + \
                                      " in record: " + str(record.id)



    output_handle = open("outbput.gb", "w")
    SeqIO.write(newRecord, output_handle, "genbank")
    output_handle.close()