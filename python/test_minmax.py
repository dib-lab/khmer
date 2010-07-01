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

    def test_merge_1(self):
        mmt = self.mmt
        mmt2 = khmer.new_minmax(MINMAXTABLE_SIZE)

        mmt.add_min(0, 2)
        mmt.add_max(0, 5)

        mmt.merge(mmt2)

        v = mmt.get_min(0)
        assert v == 2, v
        
        v = mmt.get_max(0)
        assert v == 5, v

    def test_merge_2(self):
        mmt = self.mmt
        mmt2 = khmer.new_minmax(MINMAXTABLE_SIZE)

        mmt2.add_min(0, 2)
        mmt2.add_max(0, 5)

        mmt.merge(mmt2)

        v = mmt.get_min(0)
        assert v == 2, v
        
        v = mmt.get_max(0)
        assert v == 5, v

    def test_merge_3(self):
        mmt = self.mmt
        mmt2 = khmer.new_minmax(MINMAXTABLE_SIZE)

        mmt.add_min(0, 2)
        mmt.add_max(0, 5)
        
        mmt2.add_min(0, 3)
        mmt2.add_max(0, 4)

        mmt.merge(mmt2)

        v = mmt.get_min(0)
        assert v == 2, v
        
        v = mmt.get_max(0)
        assert v == 5, v

    def test_merge_4(self):
        mmt = self.mmt
        mmt2 = khmer.new_minmax(MINMAXTABLE_SIZE)

        mmt.add_min(0, 3)
        mmt.add_max(0, 4)
        
        mmt2.add_min(0, 2)
        mmt2.add_max(0, 5)

        mmt.merge(mmt2)

        v = mmt.get_min(0)
        assert v == 2, v
        
        v = mmt.get_max(0)
        assert v == 5, v
