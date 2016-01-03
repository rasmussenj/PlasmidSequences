import collections
from Bio import SeqIO
from feature import *


class FeatureSortDic:
    def __init__(self, featureDicContainer):
        self.featureDictionary = featureDicContainer
        self.removeSpuriousCommonFeatures()
        self.removeSpuriousAnnotations()
        self.removeMultipleQulification()

    def removeSpuriousCommonFeatures(self): # remove feature less then 3
        keysToRemove = []
        for key in self.featureDictionary:
            variationSum = 0
            for variation in self.featureDictionary[key]:
                variationSum += variation.count
            if variationSum < 3:
                keysToRemove.append(key)
        for key in keysToRemove:
            del self.featureDictionary[key]

    def get_most_commom(self, qulification):
        counter=collections.Counter(qulification)
        return counter.most_common(1)[0][0]

    def removeSpuriousAnnotations(self):
        for key in self.featureDictionary:
            variationSum = 0
            for variation in self.featureDictionary[key]:
                variationSum += variation.count
            for variation in self.featureDictionary[key]:
                present_in_percent = variation.count / variationSum * 100
                if present_in_percent > 5:
                    self.featureDictionary[key].remove(variation)

    def removeMultipleQulification(self):
        for key in self.featureDictionary:
            for variation in self.featureDictionary[key]:
                variation.note = self.get_most_commom(variation.note)
                variation.gene = self.get_most_commom(variation.gene)
                variation.bound_moiety = self.get_most_commom(variation.bound_moiety)
                variation.mobile = self.get_most_commom(variation.mobile)
                variation.product = self.get_most_commom(variation.product)

    def appendSpecialTransFeature(self, path, form):
        specialTransVariationList = []
        for record in SeqIO.parse(path, form):
            specialVariation = FeatureVariation(record.seq, 1)

            specialVariation.note = record.description
            specialTransVariationList.append(specialVariation)

        self.featureDictionary["STF"] = specialTransVariationList

    def appendPrimer(self, paht, form):
        primerVariationList = []
        for record in SeqIO.parse(paht, form):
            seq = str(record.seq).split("(")[0]
            primerVariation = FeatureVariation(seq, 1)
            primerVariation.note = record.description
            primerVariationList.append(primerVariation)

        self.featureDictionary["PBS"] = primerVariationList