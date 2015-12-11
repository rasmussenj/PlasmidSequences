import textwrap

from load import *
from Bio import SeqIO
from Bio import SeqUtils


if __name__ == "__main__":

#    print list(getFeature(only_note))
    features_count = []
    feature_container = countFeatures(getFeature(),features_count)
    feature_container = Statistic(feature_container).featureContainer

    ## nach dem Statistic ausgefuehrt wurde, beinhaltet der container nur noch
    #  features die oeffter als 10% vorkommen und mind. 3 mal vorkommen
    #
    #
    # print ("\n\n\n\n\nnur zum zeigen wie ihr zugreiffen koennt\n\n\n\n\n")
    # for feature in feature_container:
    #     print(feature.name)
    #     # jedes feature kann verschiedene variationen haben
    #     for variation_f in feature.varationList:
    #         print("| %-10d %10d%% %-10s "% (variation_f.count,  variation_f.present_in_percent, variation_f.seq))



record = SeqIO.read("EcoliK12.gb", "genbank")
strRecord = str(record.seq)
print strRecord[0:10]

for feature in feature_container:
    for variation_f in feature.varationList:

        featureSeq = str(variation_f.seq)
        featureSeqComplement = str(variation_f.seq.complement())

        occurrence = SeqUtils.nt_search(strRecord, featureSeq)
        occurrenceComplement = SeqUtils.nt_search(strRecord, featureSeqComplement)

        if (len(occurrence) > 1):
            # print featureSeq[0:10]
            # print(feature.name)
            # print "it occur ", (len(occurrence)-1), "times on the forward strand"

            for i in range(1,len(occurrence)-1):
                 print str(feature.name) + ": From:", occurrence[i], "To:", (occurrence[i]+len(featureSeq))

        if (len(occurrenceComplement) > 1):
            # print featureSeqComplement[0:10]
            # print(feature.name)
            # print "it occur ", (len(occurrenceComplement)-1), "times on the complement strand"

            for i in range(1, len(occurrenceComplement)-1):
                 print str(feature.name) + " In complement From:", occurrenceComplement[i], "To:", (occurrenceComplement[i]+len(featureSeqComplement))


