import math
import collections

class Statistic:
    def __init__(self, featureContainer):
        self.featureContainer = featureContainer
        # self.removeSpuriousCommonFeatures()
        # self.removeSpuriousAnnotations()
        # self.showStaistic()

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

    def get_most_commom(self, qulification):
        counter=collections.Counter(qulification)
        return counter.most_common(1)[0][0]

    def removeSpuriousAnnotations(self):
        ## remove feature variation if:
        #  present in less than 10 % of the instances
        for feature in self.featureContainer:
            countAll = math.fsum([x.count for x in feature.varationList])
            tempVariationList = []
            for f_var in feature.varationList:
                f_var.present_in_percent = f_var.count / countAll * 100
                if (f_var.present_in_percent > 5):
                    tempVariationList.append(f_var)

            feature.varationList = sorted(tempVariationList, key=lambda var: var.count, reverse=True)



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