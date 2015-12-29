import math
import collections

class Statistic:
    def __init__(self, featureContainer):
        self.featureContainer = featureContainer
        self.removeSpuriousCommonFeatures()
        self.removeSpuriousAnnotations()
        self.removeMultipleQulification()
        # self.showStaistic()

    def removeSpuriousCommonFeatures(self): # remove feature less then 3
        keysToRemove = []
        for key in self.featureContainer:
            variationSum = 0
            for variation in self.featureContainer[key]:
                variationSum += variation.count
            if variationSum < 3:
                keysToRemove.append(key)
        for key in keysToRemove:
            del self.featureContainer[key]

    def get_most_commom(self, qulification):
        counter=collections.Counter(qulification)
        return counter.most_common(1)[0][0]



    def removeSpuriousAnnotations(self):
        for key in self.featureContainer:
            variationSum = 0
            for variation in self.featureContainer[key]:
                variationSum += variation.count
            for variation in self.featureContainer[key]:
                present_in_percent = variation.count / variationSum * 100
                if present_in_percent > 5:
                    self.featureContainer[key].remove(variation)

    def removeMultipleQulification(self):
        for key in self.featureContainer:
            for variation in self.featureContainer[key]:
                variation.note = self.get_most_commom(variation.note)
                variation.gene = self.get_most_commom(variation.gene)
                variation.bound_moiety = self.get_most_commom(variation.bound_moiety)
                variation.mobile = self.get_most_commom(variation.mobile)
                variation.product = self.get_most_commom(variation.product)

    def showStaistic(self):
        for f in self.featureContainer:
            print("+---------------------------------------------------------")
            print("| %s " % f.name)
            print("| %-20s %-20s %10s %50s %50s %50s %50s %50s"%
                  ("count", "seq Start","percent", "note", "Gene", "Product",
                   "bound_moiety", "mobile"))

            for variation_f in f.varationList:
                print("|%-20d %-20s %10d %% %50s %50s %50s %50s %50s"% (variation_f.count, variation_f.seq[:10],
                                                    variation_f.present_in_percent,
                                                    variation_f.note, variation_f.gene,
                                                                  variation_f.product, variation_f.bound_moiety,
                                                                  variation_f.mobile))

        log_file = open("Log_Features_in_vector.txt", mode='w')
        log_file.write(" %s \t %s \t %s \t %s  \t %s \t %s \t %s \t %s \t %s \n"%
                       ("feature", "count", "seq Start","percent", "note", "Gene",
                        "Product", "bound_moiety", "mobile"))
        for f in self.featureContainer:

            log_file.write("\n\n %s \n" % f.name)
            for variation_f in f.varationList:

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