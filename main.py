from Bio.Seq import reverse_complement

from load import *
from statistic import *
from Bio import SeqIO
from Bio import SeqUtils
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import SeqFeature
import pickle

file_Name = "featureObjects"

def generate():
    featureStatistic_container = getFeature()
    # for key in featureStatistic_container:
    #     print key, "---------------------------------new feature"
    #     for variation in featureStatistic_container[key]:
    #         if variation.count > 0 and key == 'CDS':
    #             print "-------------new variation -----------"
    #             print "count", variation.count
    #             print "note", variation.note
    #             print "gene", variation.gene
    #             print "bound_moiety", variation.bound_moiety
    #             print "mobile", variation.mobile
    #             print "product", variation.product

    featureStatistic_container = Statistic(featureStatistic_container).featureContainer

    # print "----------------------------------------------------------------- after staticstis ----------------------------"
    #
    # for key in featureStatistic_container:
    #     print key, "---------------------------------new feature"
    #     for variation in featureStatistic_container[key]:
    #         if variation.count > 0 and key == 'CDS':
    #             print "-------------new variation -----------"
    #             print "seq", variation.seq
    #             print "count", variation.count
    #             print "note", variation.note
    #             print "gene", variation.gene
    #             print "bound_moiety", variation.bound_moiety
    #             print "mobile", variation.mobile
    #             print "product", variation.product

    return featureStatistic_container

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
    featureStatistic_container
    # write(featureStatistic_container)
    featureStatistic_container
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

            featureSeqComplement = str(reverse_complement(variation.seq))
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


    output_handle = open("outbput.gb", "w")
    SeqIO.write(newRecord, output_handle, "genbank")
    output_handle.close()