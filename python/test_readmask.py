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
