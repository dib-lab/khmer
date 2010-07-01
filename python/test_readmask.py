import khmer

READTABLE_SIZE=50

class Test_Basic:
    def __init__(self):
        self.rt = khmer.new_readmask(READTABLE_SIZE)
        
    def test_set_true(self):
        rt = self.rt

        rt.set(0, True)
        assert rt.get(0)

    def test_set_false(self):
        rt = self.rt

        rt.set(0, False)
        assert not rt.get(0)

    def test_do_and_1(self):
        rt = self.rt

        rt.set(0, True)
        rt.do_and(0, True)
        assert rt.get(0)

    def test_do_and_2(self):
        rt = self.rt

        rt.set(0, True)
        rt.do_and(0, False)
        assert not rt.get(0)

    def test_do_and_3(self):
        rt = self.rt

        rt.set(0, False)
        rt.do_and(0, True)
        assert not rt.get(0)

    def test_do_and_4(self):
        rt = self.rt

        rt.set(0, False)
        rt.do_and(0, False)
        assert not rt.get(0)

    def test_merge_1(self):
        rt = self.rt
        rt2 = khmer.new_readmask(READTABLE_SIZE)

        rt.set(0, True)
        rt2.set(0, True)

        rt.merge(rt2)

        v = rt.get(0)
        assert v, v

    def test_merge_2(self):
        rt = self.rt
        rt2 = khmer.new_readmask(READTABLE_SIZE)

        rt.set(0, True)
        rt2.set(0, False)

        rt.merge(rt2)

        v = rt.get(0)
        assert not v, v

    def test_merge_3(self):
        rt = self.rt
        rt2 = khmer.new_readmask(READTABLE_SIZE)

        rt.set(0, False)
        rt2.set(0, True)

        rt.merge(rt2)

        v = rt.get(0)
        assert not v, v

    def test_merge_4(self):
        rt = self.rt
        rt2 = khmer.new_readmask(READTABLE_SIZE)

        rt.set(0, False)
        rt2.set(0, False)

        rt.merge(rt2)

        v = rt.get(0)
        assert not v, v
