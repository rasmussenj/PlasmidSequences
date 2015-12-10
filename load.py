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

def getFeature():
    handle = open("vectors-100.gb", "rU")
    for record in SeqIO.parse(handle, "genbank") :
        for f in record.features:

            if f.type in only_note:
                if (testFeature(f.location.start,f.location.end) == 1):
                    feature = Feature(record.seq[f.location.start:f.location.end], f.type,  0)
                    feature.note = f.qualifiers.get('note')
                    yield feature
            if f.type in note_and_gene:
                if (testFeature(f.location.start,f.location.end) == 1):
                    feature = Feature(record.seq[f.location.start:f.location.end], f.type,  0)
                    feature.gene = f.qualifiers.get('gene')
                    yield feature
            if f.type in gene_and_product:
                if (testFeature(f.location.start,f.location.end) == 1):
                    feature = Feature(record.seq[f.location.start:f.location.end], f.type,  0)
                    feature.product = f.qualifiers.get('product')
                    feature.gene = f.qualifiers.get('gene')
                    yield feature
            if f.type in note_and_bound_moiety:
                if (testFeature(f.location.start,f.location.end) == 1):
                    feature = Feature(record.seq[f.location.start:f.location.end], f.type,  0)
                    feature.product = f.qualifiers.get('bound_moiety')
                    yield feature
            if f.type in note_and_mobile:
                if (testFeature(f.location.start,f.location.end) == 1):
                    feature = Feature(record.seq[f.location.start:f.location.end], f.type,  0)
                    feature.product = f.qualifiers.get('moble_element')
                    yield feature
            if f.type in gene:
                if (testFeature(f.location.start,f.location.end) == 1):
                    feature = Feature(record.seq[f.location.start:f.location.end], f.type,  0)
                    feature.gene = f.qualifiers.get('gene')
                    yield feature
            if f.type in product:
                if (testFeature(f.location.start,f.location.end) == 1):
                    feature = Feature(record.seq[f.location.start:f.location.end], f.type,  0)
                    feature.product = f.qualifiers.get('product')
                    yield feature

def testFeature(start, end):
    if (start +3 > end):
        return 0
    else:
        return 1

def countFeatures(features, countList):
    for feature in features:
        if len(countList) == 0:
            for q in [only_note, note_and_gene, gene_and_product, note_and_bound_moiety, note_and_mobile, gene, product]:
                for f in q:
                    countList.append(FeatureStat(f))

        for count_f in it.ifilter(lambda f: f.name == feature.name, countList):
            seq_in_list = False
            for variation in it.ifilter(lambda f: f.seq == feature.seq, count_f.varationList):
                variation.count += 1
                seq_in_list = True

            if not seq_in_list:
                newVariation = FeatureStat.Types(feature.seq,1)
                newVariation.note = feature.note
                newVariation.mobile = feature.mobile
                newVariation.bound_moiety = feature.bound_moiety
                newVariation.gene = feature.gene
                newVariation.product = feature.product
                count_f.varationList.append(newVariation)


    return countList