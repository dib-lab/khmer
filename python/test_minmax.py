import khmer

MINMAXTABLE_SIZE=50

class Test_Basic:
    def __init__(self):
        self.mmt = khmer.new_minmax(MINMAXTABLE_SIZE)

    def test_min_1(self):
        mmt = self.mmt
        
        assert mmt.get_min(0) == 0
        v = mmt.add_min(0, 1)
        assert v == 1, v
        assert mmt.get_min(0) == 1

    def test_min_2(self):
        mmt = self.mmt
        
        assert mmt.get_min(0) == 0
        v = mmt.add_min(0, 2)
        assert v == 2
        
        v = mmt.add_min(0, 1)
        assert v == 1
        
        assert mmt.get_min(0) == 1
        
    def test_max_1(self):
        mmt = self.mmt
        
        assert mmt.get_max(0) == 0
        v = mmt.add_max(0, 1)
        assert v == 1, v
        assert mmt.get_max(0) == 1

    def test_max_2(self):
        mmt = self.mmt
        
        assert mmt.get_max(0) == 0
        v = mmt.add_max(0, 2)
        assert v == 2
        
        v = mmt.add_max(0, 1)
        assert v == 2
        
        assert mmt.get_max(0) == 2
