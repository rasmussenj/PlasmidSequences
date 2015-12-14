from load import *
from Bio import SeqIO
from Bio import SeqUtils
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import SeqFeature
import pickle

file_Name = "featureObjects"

def generate():
    featureStatistic_container = countFeatures(getFeature())
    featureStatistic_container = Statistic(featureStatistic_container).featureContainer
    return featureStatistic_container

def write(featureStatistic_container):

    # open the file for writing
    fileObject = open(file_Name,'wb')
    pickle.dump(featureStatistic_container,fileObject)
    fileObject.close()


def read():
    ileObject = open(file_Name,'r')
    # load the object from the file into var b
    return pickle.load(fileObject)

if __name__ == "__main__":
    featureStatistic_container = generate()
    write(featureStatistic_container)

    record = SeqIO.read("nanobody.fasta", "fasta")
    newRecord = SeqRecord(record.seq)

    #writing Header
    newRecord.seq.alphabet = generic_dna
    newRecord.id=record.id
    newRecord.name = record.name
    newRecord.description = record.description

    strRecord = str(record.seq)

    for feature in featureStatistic_container:
        for variation_f in feature.varationList:

            featureSeq = str(variation_f.seq)
            featureSeqComplement = str(variation_f.seq.complement())

            occurrence = SeqUtils.nt_search(strRecord, featureSeq)

            if (len(occurrence) > 1):
                for i in range(1,len(occurrence)-1):
                    newFeature = SeqFeature(FeatureLocation(occurrence[i],occurrence[i]+len(featureSeq), strand=1), type=str(feature.name))
                    if variation_f.product != None:
                        newFeature.qualifiers['product'] = variation_f.product
                    if variation_f.gene != None:
                        newFeature.qualifiers['gene'] = variation_f.gene
                    if variation_f.bound_moiety != None:
                        newFeature.qualifiers['bound_moiety'] = variation_f.bound_moiety
                    if variation_f.mobile != None:
                        newFeature.qualifiers['mobile'] = variation_f.mobile
                    if variation_f.note != None:
                        newFeature.qualifiers['note'] = variation_f.note
                    newRecord.features.append(newFeature)

            featureSeqReverse = str(variation_f.seq.reverse_complement)
            occurrenceComplement = SeqUtils.nt_search(strRecord, featureSeqComplement)
            if (len(occurrenceComplement) > 1):
                for i in range(1, len(occurrenceComplement)-1):
                    newFeature = SeqFeature(FeatureLocation(occurrenceComplement[i],occurrenceComplement[i]+len(featureSeq), strand=-1), type=str(feature.name))
                    if variation_f.product != None:
                        newFeature.qualifiers['product'] = variation_f.product
                    if variation_f.gene != None:
                        newFeature.qualifiers['gene'] = variation_f.gene
                    if variation_f.bound_moiety != None:
                        newFeature.qualifiers['bound_moiety'] = variation_f.bound_moiety
                    if variation_f.mobile != None:
                        newFeature.qualifiers['mobile'] = variation_f.mobile
                    if variation_f.note != None:
                        newFeature.qualifiers['note'] = variation_f.note
                    newRecord.features.append(newFeature)


    output_handle = open("outbput.gb", "w")
    SeqIO.write(newRecord, output_handle, "genbank")
    output_handle.close()