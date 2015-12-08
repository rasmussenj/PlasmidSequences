from Bio import SeqIO
import itertools as it
from feature import *

def getFeature(features_to_check):
    record_iterator = SeqIO.parse("vectors-100.gb", "genbank")
    for i in range(0, 50):
        first_record = next(record_iterator)
        for f in first_record.features:
            if f.type in features_to_check:
                yield Feature(first_record.seq[f.location.start:f.location.end], f.type,  0)


def countFeatures(features, countList, features_to_check):
    for feature in features:
        if len(countList) == 0:
            for f in features_to_check:
                countList.append(FeatureStat(f))

        for count_f in it.ifilter(lambda f: f.name == feature.name, countList):
            seq_in_list = False
            for variation in it.ifilter(lambda f: f.seq == feature.seq, count_f.varationList):
                variation.count += 1
                seq_in_list = True

            if not seq_in_list:
                count_f.varationList.append(FeatureStat.Types(feature.seq,1))


    return countList