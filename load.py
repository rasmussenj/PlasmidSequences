from Bio import SeqIO
import itertools as it
from feature import *

only_note = ['promoter', 'oriT', 'rep_origin', 'primer_bind', 'terminator', 'misc_signal', 'misc_recomb', 'LTR', 'enhancer']
note_and_gene = ['-35_signal', '-10_signal', 'RBS', 'polyA_signal', 'sig_peptide']
gene_and_product = ['CDS']
note_and_bound_moiety = ['protein_bind', 'misc_binding']
note_and_mobile = ['mobile_element']
only_gene = ['mRNA']
only_product = ['tRNA', 'rRNA']

features_to_check_list = [only_note, note_and_gene, gene_and_product,
                          note_and_bound_moiety, note_and_mobile, only_gene, only_product]
features_Container = {}
# fill the features_Container with all features as FeatureStatistic object
for features_to_check in features_to_check_list:
    for f in features_to_check:
        features_Container[f] = []


################################# Functions ######################################

def getFeature():
    handle = open("vectors-5.gb", "rU")
    for record in SeqIO.parse(handle, "genbank") :
        if len(record.seq) > 1500:
            for f in record.features:
                for list in features_to_check_list:
                    if f.type in list and testSeqLength(f.location.start, f.location.end):
                        countFeatures(f, record.seq[f.location.start:f.location.end])
    return features_Container


def testSeqLength(start, end):
    if (start + 2 > end):
        return False
    else:
        return True


def countFeatures(feature, seq):
        # get the FeatureStatistic object which corresponds with the actual feature

        accVariationList = features_Container[feature.type]
        note = feature.qualifiers.get('note')[0] if feature.qualifiers.get('note') != None else None
        gene = feature.qualifiers.get('gene')[0] if feature.qualifiers.get('gene') != None else None
        bound_moiety = feature.qualifiers.get('bound_moiety')[0] if feature.qualifiers.get('bound_moiety') != None else None
        mobile = feature.qualifiers.get('mobile')[0] if feature.qualifiers.get('mobile') != None else None
        product = feature.qualifiers.get('product')[0] if feature.qualifiers.get('product') != None else None

        seq_in_list = False
        for variation in accVariationList:
            if variation.seq == seq:
                if feature.type in only_note:
                    if variation.note[0] == note:
                        seq_in_list = True
                        appendQualifier(variation.gene, gene)
                        appendQualifier(variation.bound_moiety, bound_moiety)
                        appendQualifier(variation.mobile, mobile)
                        appendQualifier(variation.product, product)
                        variation.count += 1
                        break

                elif feature.type in note_and_gene:
                    if (variation.note[0] == note and variation.gene[0] == gene):
                        seq_in_list = True
                        appendQualifier(variation.bound_moiety, bound_moiety)
                        appendQualifier(variation.mobile, mobile)
                        appendQualifier(variation.product, product)
                        variation.count += 1
                        break

                elif feature.type in gene_and_product:
                    if (variation.gene[0] == gene and variation.product[0] == product):
                        seq_in_list = True
                        appendQualifier(variation.bound_moiety, bound_moiety)
                        appendQualifier(variation.mobile, mobile)
                        appendQualifier(variation.note, note)

                        variation.count += 1
                        break

                elif feature.type in note_and_bound_moiety:
                    if (variation.note[0] == note and variation.bound_moiety[0] == bound_moiety):
                        seq_in_list = True
                        appendQualifier(variation.gene, gene)
                        appendQualifier(variation.mobile, mobile)
                        appendQualifier(variation.product, product)
                        variation.count += 1
                        break

                elif feature.type in note_and_mobile:
                    if (variation.note[0] == note and variation.mobile[0] == mobile):
                        seq_in_list = True
                        appendQualifier(variation.gene, gene)
                        appendQualifier(variation.bound_moiety, bound_moiety)
                        appendQualifier(variation.product, product)
                        variation.count += 1
                        break

                elif feature.type in only_gene:
                    if (variation.gene[0] == variation.gene):
                        seq_in_list = True
                        appendQualifier(variation.bound_moiety, bound_moiety)
                        appendQualifier(variation.product, product)
                        appendQualifier(variation.mobile, mobile)
                        appendQualifier(variation.note, note)
                        variation.count += 1
                        break

                elif feature.type in only_product:
                    if (variation.product[0] == product):
                        seq_in_list = True
                        appendQualifier(variation.bound_moiety, bound_moiety)
                        appendQualifier(variation.mobile, mobile)
                        appendQualifier(variation.note, note)
                        appendQualifier(variation.gene, gene)
                        variation.count += 1
                        break

        # if seq not found in statFeature, create a new variation of the feature
        if not seq_in_list:
            new_variation = FeatureStatistic.Varation(seq,1)
            new_variation.note.append(note)
            new_variation.bound_moiety.append(bound_moiety)
            new_variation.gene.append(gene)
            new_variation.product.append(product)
            new_variation.mobile.append(mobile)
            accVariationList.append(new_variation)


def appendQualifier(qualifier, value):
    if value is not None:
        qualifier.append(value)