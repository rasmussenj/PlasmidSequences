import textwrap

from load import *
from Bio import SeqIO
from Bio import SeqUtils


if __name__ == "__main__":

#    print list(getFeature(only_note))
    features_count = []
    feature_container = countFeatures(getFeature(),features_count)
    feature_container = Statistic(feature_container).featureContainer


record = SeqIO.read("EcoliK12.gb", "genbank")
strRecord = str(record.seq)
print strRecord[0:10]

print "%-20s %-20s"% ("FEATURES", "Location/Qualifiers")
for feature in feature_container:
    for variation_f in feature.varationList:

        featureSeq = str(variation_f.seq)
        featureSeqComplement = str(variation_f.seq.complement())

        occurrence = SeqUtils.nt_search(strRecord, featureSeq)
        occurrenceComplement = SeqUtils.nt_search(strRecord, featureSeqComplement)

        if (len(occurrence) > 1):
            for i in range(1,len(occurrence)-1):
                 print "%-20s %-20s"% (str(feature.name), (str(occurrence[i]) + ".." + str((occurrence[i]+len(featureSeq)))))
                 if variation_f.product != None:
                     print "%-20s %-20s"% ("", str("/product=\"" + variation_f.product[0] + "\""))
                 if variation_f.gene != None:
                     print "%-20s %-20s"% ("", str("/gene=\"" + variation_f.gene[0] + "\""))
                 if variation_f.bound_moiety != None:
                     print "%-20s %-20s"% ("", str("/bound_moiety=\"" + variation_f.bound_moiety[0] + "\""))
                 if variation_f.mobile != None:
                     print "%-20s %-20s"% ("", str("/mobile=\"" + variation_f.mobile[0] + "\""))
                 if variation_f.note != None:
                     print "%-20s %-20s"% ("", str("/note=\"" + variation_f.note[0] + "\""))


        if (len(occurrenceComplement) > 1):
            for i in range(1, len(occurrenceComplement)-1):
                 print "%-20s %-20s"% (str(feature.name), ("complement(" + (str(occurrenceComplement[i]) + ".." + str((occurrenceComplement[i]+len(featureSeqComplement)))) + ")"))
                 if variation_f.product != None:
                     print "%-20s %-20s"% ("", str("/product=\"" + variation_f.product[0] + "\""))
                 if variation_f.gene != None:
                     print "%-20s %-20s"% ("", str("/gene=\"" + variation_f.gene[0] + "\""))
                 if variation_f.bound_moiety != None:
                     print "%-20s %-20s"% ("", str("/bound_moiety=\"" + variation_f.bound_moiety[0] + "\""))
                 if variation_f.mobile != None:
                     print "%-20s %-20s"% ("", str("/mobile=\"" + variation_f.mobile[0] + "\""))
                 if variation_f.note != None:
                     print "%-20s %-20s"% ("", str("/note=\"" + variation_f.note[0] + "\""))

