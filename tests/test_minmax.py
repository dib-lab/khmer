import khmer

import khmer_tst_utils as utils

MINMAXTABLE_SIZE=50

def teardown():
    utils.cleanup()

class Test_Basic(object):
    def __init__(self):
        self.mmt = khmer.new_minmax(MINMAXTABLE_SIZE)

    def test_tablesize(self):
        assert self.mmt.tablesize() == MINMAXTABLE_SIZE

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
        
class Test_Filestuff(object):
    def __init__(self):
        self.mmt = khmer.new_minmax(MINMAXTABLE_SIZE)

    def test_saveload(self):
        filename = utils.get_temp_filename('save')

        mmt = self.mmt

        for i in range(0, MINMAXTABLE_SIZE):
            mmt.add_min(i, i)
            mmt.add_max(i, MINMAXTABLE_SIZE-i)

        print mmt
        print 'saving to:', filename
        mmt.save(filename)

        mmt2 = khmer.new_minmax(0)
        print 'mm2'
        print 'loading from:', filename
        mmt2.load(filename)

        for i in range(0, MINMAXTABLE_SIZE):
            assert mmt.get_min(i) == mmt2.get_min(i)
            assert mmt.get_min(i) == i
            
            assert mmt.get_max(i) == mmt2.get_max(i)
            assert mmt.get_max(i) == MINMAXTABLE_SIZE - i

    def test_save_no_load(self):
        filename = utils.get_temp_filename('save')

        mmt = self.mmt

        for i in range(0, MINMAXTABLE_SIZE):
            mmt.add_min(i, i)
            mmt.add_max(i, MINMAXTABLE_SIZE-i)

        mmt.save(filename)

        mmt2 = khmer.new_minmax(MINMAXTABLE_SIZE)
        # no load!

        try:
            for i in range(0, MINMAXTABLE_SIZE):
                if mmt2.get_min(i) != mmt.get_min(i):
                    raise Exception         # supposed to happen!

            assert 0
        except AssertionError:
            raise
        except Exception:
            pass

        try:
            for i in range(0, MINMAXTABLE_SIZE):
                if mmt2.get_max(i) != mmt.get_max(i):
                    raise Exception         # supposed to happen!
            assert 0
        except AssertionError:
            raise
        except Exception:
            pass
