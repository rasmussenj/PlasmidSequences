from Bio import SeqIO

from feature import *

def getFeature():
    record_iterator = SeqIO.parse("vectors-100.gb", "genbank")
    for i in range(0, 50):
        first_record = next(record_iterator)
        for f in first_record.features:
            if f.type in ["source", "promoter", "misc_feature"]:
                yield Feature(first_record.seq[f.location.start:f.location.end], f.type,  0)


def countFeatures(features, countList):
    for feature in features:
        if len(countList) == 0:
            countList.append(Feature(feature.seq, feature.name,1))
        else:
            feature_in_cout = False
            for f_count in countList:
                if f_count.name == feature.name:
                    #Check seq
                    if f_count.seq == feature.seq:
                        feature_in_cout = True
                        f_count.count += 1

            if not feature_in_cout:
                countList.append(Feature(feature.seq, feature.name, 1))

    return countList
