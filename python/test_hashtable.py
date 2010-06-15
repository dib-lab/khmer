import khmer

def test_no_collision():
    kh = khmer.new_hashtable(4, 86)

    kh.count('AAAA')
    assert kh.get('AAAA') == 1

    kh.count('TTTT')
    assert kh.get('TTTT') == 1

def test_collision():
    kh = khmer.new_hashtable(4, 85)

    kh.count('AAAA')
    assert kh.get('AAAA') == 1

    kh.count('TTTT')
    assert kh.get('TTTT') == 2

def test_complete_no_collision():
    kh = khmer.new_hashtable(4, 4**4)
    kt = khmer.new_ktable(4)

    for i in range(0, kt.n_entries()):
        s = kt.reverse_hash(i)
        kh.count(s)

    for i in range(0, kt.n_entries()):
        s = kt.reverse_hash(i)
        assert kh.get(s) == 1

def test_complete_2_collision():
    kh = khmer.new_hashtable(4, 4**4 / 2)
    kt = khmer.new_ktable(4)

    for i in range(0, kt.n_entries()):
        s = kt.reverse_hash(i)
        kh.count(s)

    for i in range(0, kt.n_entries()):
        s = kt.reverse_hash(i)
        assert kh.get(s) == 2

def test_complete_4_collision():
    kh = khmer.new_hashtable(4, 4**4 / 4)
    kt = khmer.new_ktable(4)

    for i in range(0, kt.n_entries()):
        s = kt.reverse_hash(i)
        kh.count(s)

    for i in range(0, kt.n_entries()):
        s = kt.reverse_hash(i)
        assert kh.get(s) == 4
