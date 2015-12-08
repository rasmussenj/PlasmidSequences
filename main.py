from load import *


if __name__ == "__main__":
    features_to_check = ["source", "promoter", "misc_feature"]
    print list(getFeature(features_to_check))
    features_count = []
    featureSta = countFeatures(getFeature(features_to_check),features_count, features_to_check)

    for f in featureSta:
        print(f.name)
        for variation_f in f.varationList:
            print(variation_f.count)
            print(variation_f.seq)
            print("---------------------------------------------")

        print("***********************new Feature********************************")