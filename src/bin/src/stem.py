
class Stem:
    def __init__(self, base_pairs):
        self.base_pairs = base_pairs

    def get_start(self):
        return self.base_pairs[0][0]

    def get_end(self):
        return self.base_pairs[-1][0]
    
    def __repr__(self):
        return f"Stem(base_pairs={self.base_pairs})"
