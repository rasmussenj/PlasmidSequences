import math


class FeatureStat:
    def __init__(self, name):
        self.varationList = []
        self.name = name

    class Types:
        def __init__(self, seq, count):
            self.seq = seq
            self.count = count
            self.present_in_percent = None
            self.product = None
            self.note = None

class Feature:
    def __init__(self, seq, name, count):
        self.seq = seq
        self.count = count
        self.name = name
        self.product = None
        self.gene = None
        self.bound_moiety = None
        self.mobile = None




class Statistic:
    def __init__(self, featureContainer):
        self.featureContainer = featureContainer

        self.removeSpuriousCommonFeatures()
        self.removeSpuriousAnnotations()
        self.showStaistic()

    def removeSpuriousCommonFeatures(self):
        tempFeatureContainer =[]

        for feature in self.featureContainer:
            ## remove all features without any variation
            if len(feature.varationList) != 0:
                ## remove feature variation if:
                #  less than three occurrences in learning file
                tempVariationList =[]
                for f_var in feature.varationList:
                    if f_var.count > 3:
                        tempVariationList.append(f_var)
                feature.varationList = tempVariationList
                tempFeatureContainer.append(feature)
        self.featureContainer = tempFeatureContainer

    def removeSpuriousAnnotations(self):
        ## remove feature variation if:
        #  present in less than 10 % of the instances
        for feature in self.featureContainer:
            countAll = math.fsum([x.count for x in feature.varationList])
            tempVariationList = []
            for f_var in feature.varationList:
                f_var.present_in_percent = f_var.count / countAll * 100
                if (f_var.present_in_percent > 10):
                    tempVariationList.append(f_var)

            feature.varationList = sorted(tempVariationList, key=lambda var: var.count)



    def showStaistic(self):
        for f in self.featureContainer:
            print("+---------------------------------------------------------")
            print("| %s " % f.name)
            print("| %-20s %-20s %10s %14s "% ("count", "seq Start","percent", "product"))

            for variation_f in f.varationList:
                print("|%-20d %-20s %10d %% %10s"% (variation_f.count, variation_f.seq[:10], variation_f.present_in_percent, variation_f.note))
