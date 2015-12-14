import math
import collections

only_note = ['promoter', 'oriT', 'rep_origin', 'primer_bind', 'terminator', 'misc_signal', 'misc_recomb', 'LTR', 'enhancer']
note_and_gene = ['-35_signal', '-10_signal', 'RBS', 'polyA_signal', 'sig_peptide']
gene_and_product = ['CDS']
note_and_bound_moiety = ['protein_bind', 'misc_binding']
note_and_mobile = ['mobile_element']
gene = ['mRNA']
product = ['tRNA', 'rRNA']



class FeatureStatistic:
    def __init__(self, name):
        self.varationList = []
        self.name = name

    class Varation:
        def __init__(self, seq, count):
            self.seq = seq
            self.count = count
            self.present_in_percent = None
            self.product = []
            self.gene = []
            self.bound_moiety = []
            self.mobile = []
            self.note = []


        def isSame(self, nextFeature):
            if nextFeature.type in only_note:
                if self.note[0] == nextFeature.note:
                    self.gene.append()
                    self.bound_moiety.append()
                    self.mobile.append()
                    self.product.append()

                    return True

            if nextFeature.type in note_and_gene:
                if (self.note[0] == nextFeature.note and
                    self.gene[0] == nextFeature.gene):

                    self.bound_moiety.append(nextFeature.bound_moiety)
                    self.mobile.append(nextFeature.mobile)
                    self.product.append(nextFeature.product)
                    return True


            if nextFeature.type in gene_and_product:
                if (self.gene[0] == nextFeature.gene and
                            self.product[0] == nextFeature.product):

                    self.bound_moiety.append(nextFeature.bound_moiety)
                    self.mobile.append(nextFeature.mobile)
                    self.note.append(nextFeature.note)

                    return True

            if nextFeature.type in note_and_bound_moiety:
                if (self.note[0] == nextFeature.note and
                            self.bound_moiety[0] == nextFeature.bound_moiety):
                    self.gene.append(nextFeature.gene)

                    self.mobile.append(nextFeature.mobile)

                    self.product.append(nextFeature.product)
                    return True

            if nextFeature.type in note_and_mobile:
                if (self.note[0] == nextFeature.note and
                            self.mobile[0] == nextFeature.mobile):
                    self.gene.append(nextFeature.gene)
                    self.bound_moiety.append(nextFeature.bound_moiety)

                    self.product.append(nextFeature.product)
                    return True


            if nextFeature.type in gene:
                if (self.gene[0] == nextFeature.gene):

                    self.bound_moiety.append(nextFeature.bound_moiety)
                    self.mobile.append(nextFeature.mobile)
                    self.note.append(nextFeature.note)
                    self.product.append(nextFeature.product)
                    return True

            if nextFeature.type in product:
                if (self.product[0] == nextFeature.product):
                    self.gene.append(nextFeature.gene)
                    self.bound_moiety.append(nextFeature.bound_moiety)
                    self.mobile.append(nextFeature.mobile)
                    self.note.append(nextFeature.note)
                    return True


class Feature:
    def __init__(self, seq, name, count):
        self.seq = seq
        self.count = count
        self.name = name
        self.product = None
        self.gene = None
        self.bound_moiety = None
        self.mobile = None
        self.note = None




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