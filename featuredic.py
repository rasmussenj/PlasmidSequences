import math
import collections
from Bio import SeqIO
from feature import *


class FeatureDic:
    def __init__(self, featureContainer):
        self.featureDictionary = featureContainer
        self.removeSpuriousCommonFeatures()
        self.removeSpuriousAnnotations()
        self.removeMultipleQulification()
        # self.showStaistic()

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

    def showStaistic(self):
        for f in self.featureDictionary:
            print("+---------------------------------------------------------")
            print("| %s " % f.name)
            print("| %-20s %-20s %10s %50s %50s %50s %50s %50s"%
                  ("count", "seq Start","percent", "note", "Gene", "Product",
                   "bound_moiety", "mobile"))

            for variation_f in f.variationList:
                print("|%-20d %-20s %10d %% %50s %50s %50s %50s %50s"% (variation_f.count, variation_f.seq[:10],
                                                    variation_f.present_in_percent,
                                                    variation_f.note, variation_f.gene,
                                                                  variation_f.product, variation_f.bound_moiety,
                                                                  variation_f.mobile))

        log_file = open("Log_Features_in_vector.txt", mode='w')
        log_file.write(" %s \t %s \t %s \t %s  \t %s \t %s \t %s \t %s \t %s \n"%
                       ("feature", "count", "seq Start","percent", "note", "Gene",
                        "Product", "bound_moiety", "mobile"))
        for f in self.featureDictionary:

            log_file.write("\n\n %s \n" % f.name)
            for variation_f in f.variationList:

                log_file.write("%s \t  %s \t %d %% \t %s \t %s \t %s \t %s \t %s \t %s\n"% (
                    " ",
                    variation_f.count,
                    variation_f.present_in_percent,
                    variation_f.note,
                    variation_f.gene,
                    variation_f.product,
                    variation_f.seq,
                    variation_f.bound_moiety,
                    variation_f.mobile))

        log_file.close()

    def appendPrimer(self, path, form):
        primerVariationList = []
        for record in SeqIO.parse(path, form):
            primerVariation = Variation(record.seq, 1)

            primerVariation.note = record.description
            primerVariationList.append(primerVariation)

        self.featureDictionary["PBS"] = primerVariationList
