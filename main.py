from load import *


if __name__ == "__main__":
    only_note = ['promoter', 'oriT', 'rep_origin', 'primer_bind', 'terminator', 'misc_signal', 'misc_recomb', 'LTR', 'enhancer']
    note_and_gene = ['-35_signal', '-10_signal', 'RBS', 'polyA_signal', 'sig_peptide']
    gene_and_product = ['CDS']
    note_and_bound_moiety = ['protein_bind', 'misc_binding']
    note_and_mobile = ['mobile_element']
    gene = ['mRNA']
    product = ['tRNA', 'rRNA']
    print list(getFeature(only_note))
    features_count = []
    featureSta = countFeatures(getFeature(only_note),features_count, only_note)

    for f in featureSta:
        print(f.name)
        for variation_f in f.varationList:
            print(variation_f.count)
            print(variation_f.seq)
            print("---------------------------------------------")

        print("***********************new Feature********************************")