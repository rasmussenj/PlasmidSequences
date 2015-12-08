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
    feature_container = countFeatures(getFeature(only_note),features_count, only_note)
    feature_container = Statistic(feature_container).featureContainer

    ## nach dem Statistic ausgefuehrt wurde, beinhaltet der container nur noch
    #  features die oeffter als 10% vorkommen und mind. 3 mal vorkommen


    print ("\n\n\n\n\nnur zum zeigen wie ihr zugreiffen koennt\n\n\n\n\n")
    for feature in feature_container:
        print(feature.name)
        # jedes feature kann verschiedene variationen haben
        for variation_f in feature.varationList:
            print("| %-10d %10d%% %-10s "% (variation_f.count,  variation_f.present_in_percent, variation_f.seq))
