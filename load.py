from Bio import SeqIO
import itertools as it
from feature import *

only_note = ['promoter', 'oriT', 'rep_origin', 'primer_bind', 'terminator', 'misc_signal', 'misc_recomb', 'LTR', 'enhancer']
note_and_gene = ['-35_signal', '-10_signal', 'RBS', 'polyA_signal', 'sig_peptide']
gene_and_product = ['CDS']
note_and_bound_moiety = ['protein_bind', 'misc_binding']
note_and_mobile = ['mobile_element']
gene = ['mRNA']
product = ['tRNA', 'rRNA']

features_to_check_list = [only_note, note_and_gene, gene_and_product, note_and_bound_moiety, note_and_mobile, gene, product]
features_Container = {}
def getFeature():
    handle = open("vectors-100.gb", "rU")
    for record in SeqIO.parse(handle, "genbank") :
        if len(record.seq) > 1500:
            for f in record.features:
                for list in features_to_check_list:
                    if f.type in list:

                        countFeatures(f, record.seq[f.location.start:f.location.end])

    return features_Container






def testSeqLength(start, end):
    if (start + 2 > end):
        return 0
    else:
        return 1

def countFeatures(feature, seq):



        # fill the features_Container with all features as FeatureStatistic object
        if len(features_Container) == 0:
            for features_to_check in features_to_check_list:
                for f in features_to_check:
                    features_Container[f] = []
        # get the FeatureStatistic object which corresponds with the actual feature

        accVariationList = features_Container[feature.type]
        seq_in_list = False
        for variation in it.ifilter(lambda v: v.seq == seq, accVariationList):

            note = feature.qualifiers.get('note')[0] if feature.qualifiers.get('note') != None else ""
            gene = feature.qualifiers.get('gene')[0] if feature.qualifiers.get('gene') != None else ""
            bound_moiety = feature.qualifiers.get('bound_moiety')[0] if feature.qualifiers.get('bound_moiety') != None else ""
            mobile = feature.qualifiers.get('mobile')[0] if feature.qualifiers.get('mobile') != None else ""
            product = feature.qualifiers.get('product')[0] if feature.qualifiers.get('product') != None else ""

            if feature.type in only_note:
                if variation.note[0] == note:
                    variation.gene.append(gene)
                    variation.bound_moiety.append(bound_moiety)
                    variation.mobile.append(mobile)
                    variation.product.append(product)

                    variation.count += 1
                    seq_in_list = True


            elif feature.type in note_and_gene:
                if (variation.note[0] == note and gene):

                    variation.bound_moiety.append(bound_moiety)
                    variation.mobile.append(mobile)
                    variation.product.append(product)
                    variation.count += 1
                    seq_in_list = True

            elif feature.type in gene_and_product:
                if (variation.gene[0] == gene and product):

                    variation.bound_moiety.append(bound_moiety)
                    variation.mobile.append(mobile)
                    variation.note.append(note)

                    variation.count += 1
                    seq_in_list = True

            elif feature.type in note_and_bound_moiety:
                if (variation.note[0] == note and bound_moiety):
                    variation.gene.append(gene)

                    variation.mobile.append(mobile)

                    variation.product.append(product)

                    variation.count += 1
                    seq_in_list = True

            elif feature.type in note_and_mobile:
                if (variation.note[0] == note and mobile):
                    variation.gene.append(gene)
                    variation.bound_moiety.append(bound_moiety)

                    variation.product.append(product)
                    variation.count += 1
                    seq_in_list = True

            elif feature.type in gene:
                if (variation.gene[0] == gene):

                    variation.bound_moiety.append(bound_moiety)
                    variation.mobile.append(mobile)
                    variation.note.append(note)
                    variation.product.append(product)
                    variation.count += 1
                    seq_in_list = True

            elif feature.type in product:
                if (variation.product[0] == product):
                    variation.gene.append(gene)
                    variation.bound_moiety.append(bound_moiety)
                    variation.mobile.append(mobile)
                    variation.note.append(note)
                    variation.count += 1
                    seq_in_list = True

        # if seq not found in statFeature, create a new variation of the feature
        if not seq_in_list:

            new_variation = FeatureStatistic.Varation(seq,1)

            new_variation.note.append(feature.qualifiers.get('note')[0]
                                      if feature.qualifiers.get('note') != None else "")

            new_variation.bound_moiety.append(feature.qualifiers.get('bound_moiety')
                                              if feature.qualifiers.get('bound_moiety') != None else "")

            new_variation.gene.append(feature.qualifiers.get('gene')[0]
                                      if feature.qualifiers.get('gene') != None else "")

            new_variation.product.append(feature.qualifiers.get('product')[0]
                                         if feature.qualifiers.get('product') != None else "")

            new_variation.mobile.append(feature.qualifiers.get('mobile')
                                        if feature.qualifiers.get('mobile')!= None else "")
            accVariationList.append(new_variation)
