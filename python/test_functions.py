import khmer

def test_forward_hash():
    assert khmer.forward_hash('AAAA', 4) == 0
    assert khmer.forward_hash('TTTT', 4) == 0
    assert khmer.forward_hash('CCCC', 4) == 170
    assert khmer.forward_hash('GGGG', 4) == 170
    
def test_forward_hash_no_rc():
    h = khmer.forward_hash_no_rc('AAAA', 4)
    assert h == 0, h
    
    h = khmer.forward_hash_no_rc('TTTT', 4)
    assert h == 85, h

    h = khmer.forward_hash_no_rc('CCCC', 4)
    assert h == 170, h

    h = khmer.forward_hash_no_rc('GGGG', 4)
    assert h == 255, h

def test_reverse_hash():
    s = khmer.reverse_hash(0, 4)
    assert s == "AAAA"
    
    s = khmer.reverse_hash(85, 4)
    assert s == "TTTT"
    
    s = khmer.reverse_hash(170, 4)
    assert s == "CCCC"
    
    s = khmer.reverse_hash(255, 4)
    assert s == "GGGG"
