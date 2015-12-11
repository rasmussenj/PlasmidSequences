import textwrap

from load import *
from Bio import SeqIO
from Bio import SeqUtils
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import SeqFeature

if __name__ == "__main__":

#    print list(getFeature(only_note))
    features_count = []
    feature_container = countFeatures(getFeature(),features_count)
    feature_container = Statistic(feature_container).featureContainer

    record = SeqIO.read("EcoliK12.gb", "genbank")
    newRecord = SeqRecord(record.seq)

    #writing Header
    newRecord.seq.alphabet = generic_dna
    newRecord.id=record.id
    newRecord.name = record.name
    newRecord.description = record.description

    strRecord = str(record.seq)

    for feature in feature_container:
        for variation_f in feature.varationList:

            featureSeq = str(variation_f.seq)
            featureSeqComplement = str(variation_f.seq.complement())

            occurrence = SeqUtils.nt_search(strRecord, featureSeq)
            #occurrenceComplement = SeqUtils.nt_search(strRecord, featureSeqComplement)

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


            # if (len(occurrenceComplement) > 1):
            #     for i in range(1, len(occurrenceComplement)-1):
            #          print "%-20s %-20s"% (str(feature.name), ("complement(" + (str(occurrenceComplement[i]) + ".." + str((occurrenceComplement[i]+len(featureSeqComplement)))) + ")"))
            #          if variation_f.product != None:
            #              print "%-20s %-20s"% ("", str("/product=\"" + variation_f.product[0] + "\""))
            #          if variation_f.gene != None:
            #              print "%-20s %-20s"% ("", str("/gene=\"" + variation_f.gene[0] + "\""))
            #          if variation_f.bound_moiety != None:
            #              print "%-20s %-20s"% ("", str("/bound_moiety=\"" + variation_f.bound_moiety[0] + "\""))
            #          if variation_f.mobile != None:
            #              print "%-20s %-20s"% ("", str("/mobile=\"" + variation_f.mobile[0] + "\""))
            #          if variation_f.note != None:
            #              print "%-20s %-20s"% ("", str("/note=\"" + variation_f.note[0] + "\""))

    output_handle = open("example.gb", "w")
    SeqIO.write(newRecord, output_handle, "genbank")
    output_handle.close()