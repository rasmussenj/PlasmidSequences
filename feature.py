class FeatureStat:
    def __init__(self, name):
        self.varationList = []
        self.name = name

    class Types:
        def __init__(self, seq, count):
            self.seq = seq
            self.count = count


class Feature:
    def __init__(self, seq, name, count):
        self.seq = seq
        self.count = count
        self.name = name

