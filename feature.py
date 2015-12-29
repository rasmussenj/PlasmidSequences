class FeatureStatistic:
    def __init__(self, name):
        self.variationList = []
        self.name = name

    class Variation:
        def __init__(self, seq, count):
            self.seq = seq
            self.count = count
            self.present_in_percent = None
            self.product = []
            self.gene = []
            self.bound_moiety = []
            self.mobile = []
            self.note = []

class Feature:
    def __init__(self, seq, name, count):
        self.seq = seq
        self.count = count
        self.name = name
        self.product = None
        self.gene = None
        self.bound_moiety = None
        self.mobile = None
        self.note = None




