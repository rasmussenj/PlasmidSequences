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

def getFeature():
    handle = open("vectors.gb", "rU")
    for record in SeqIO.parse(handle, "genbank") :
        if len(record.seq) > 1500:
            for f in record.features:
                if f.type in features_to_check_list:
                    if (testSeqLength(f.location.start,f.location.end) == 1):
                        yield f

def testSeqLength(start, end):
    if (start + 2 > end):
        return 0
    else:
        return 1

def countFeatures(features):
    features_Container = []
    for feature in features:
        # fill the features_Container with all features as FeatureStatistic object
        if len(features_Container) == 0:
            for features_to_check in features_to_check_list:
                for f in features_to_check:
                    features_Container.append(FeatureStatistic(f))
        # get the FeatureStatistic object which corresponds with the actual feature
        for statFeature in it.ifilter(lambda f: f.name == feature.type, features_Container):
            seq_in_list = False

            # search for the corresponding seq in statFeature
            for variation in it.ifilter(lambda v: v.isSame(feature), statFeature.varationList):
                variation.count += 1
                seq_in_list = True

            # if seq not found in statFeature, create a new variation of the feature
            if not seq_in_list:
                new_variation = FeatureStatistic.Varation(feature.seq,1)
                new_variation.note = {feature.note}
                new_variation.mobile = {feature.mobile}
                new_variation.bound_moiety = {feature.bound_moiety}
                new_variation.gene = {feature.gene}
                new_variation.product = {feature.product}
                statFeature.varationList.append(new_variation)

    return features_Container