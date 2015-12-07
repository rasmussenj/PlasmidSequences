class Feature:
    def __init__(self, seq, name):
        self.seq = seq
        self.count = 0
        self.name = name

    def setSeq(self, seq):
        self.seq = seq

    def getSeq(self):
        return self.name

    def addOne(self):
        self.count =+1

    def getCount(self):
        return self.count

    def setName(self, name):
        self.name = name

    def getName(self):
        return self.name