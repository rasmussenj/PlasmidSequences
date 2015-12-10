from Bio import SeqIO
import itertools as it
from feature import *

def getFeature(features_to_check):
    note_and_gene = ['-35_signal', '-10_signal', 'RBS', 'polyA_signal', 'sig_peptide']
    gene_and_product = ['CDS']
    note_and_bound_moiety = ['protein_bind', 'misc_binding']
    note_and_mobile = ['mobile_element']
    gene = ['mRNA']
    product = ['tRNA', 'rRNA']
    handle = open("vectors-100.gb", "rU")
    for record in SeqIO.parse(handle, "genbank") :
        for f in record.features:

            if f.type in features_to_check:
                feature = Feature(record.seq[f.location.start:f.location.end], f.type,  0)
                feature.note = f.qualifiers.get('note')
                yield feature
            if f.type in note_and_gene:
                feature = Feature(record.seq[f.location.start:f.location.end], f.type,  0)
                feature.gene = f.qualifiers.get('gene')
                yield feature
            if f.type in gene_and_product:
                feature = Feature(record.seq[f.location.start:f.location.end], f.type,  0)
                feature.product = f.qualifiers.get('product')
                feature.gene = f.qualifiers.get('gene')
                yield feature
            if f.type in note_and_bound_moiety:
                feature = Feature(record.seq[f.location.start:f.location.end], f.type,  0)
                feature.product = f.qualifiers.get('bound_moiety')
                yield feature
            if f.type in note_and_mobile:
                feature = Feature(record.seq[f.location.start:f.location.end], f.type,  0)
                feature.product = f.qualifiers.get('moble_element')
                yield feature
            if f.type in gene:
                feature = Feature(record.seq[f.location.start:f.location.end], f.type,  0)
                feature.gene = f.qualifiers.get('gene')
                yield feature
            if f.type in product:
                feature = Feature(record.seq[f.location.start:f.location.end], f.type,  0)
                feature.product = f.qualifiers.get('product')
                yield feature



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
                newVariation = FeatureStat.Types(feature.seq,1)
                newVariation.note = feature.note
                count_f.varationList.append(newVariation)


    return countList