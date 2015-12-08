from load import *


if __name__ == "__main__":
    print list(getFeature())
    features_count = []
    print countFeatures(getFeature(),features_count)

    for f in features_count:
        print(f.name)
        print(f.seq)
        print(f.count)
