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


    print ("\n\n\n\n\nnur zum zeigen wie ihr zugreiffen koennt\n\n\n\n\n")
    for feature in feature_container:
        print(feature.name)
        # jedes feature kann verschiedene variationen haben
        for variation_f in feature.varationList:
            print("| %-10d %10d%% %-10s "% (variation_f.count,  variation_f.present_in_percent, variation_f.seq))



record = SeqIO.read("Ecolik12.gb", "genbank")
strRecord = str(record.seq)

for feature in feature_container:
    print(feature.name)
    featureSeq = str(feature.seq)
    startPosition = str(SeqUtils.nt_search(strRecord, featureSeq))
    endPosition = str(SeqUtils.nt_search(strRecord, featureSeq) + len(featureSeq))

    print ("Matching at position: " + startPosition + "..." + endPosition)
