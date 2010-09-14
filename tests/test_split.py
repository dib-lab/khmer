import khmer
from khmer.split import SplitHashtable

DNA = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"

def test_2_split():
    total = SplitHashtable(4, 4**4, 0, 1)

    half1 = SplitHashtable(4, 4**4, 0, 2)
    half2 = SplitHashtable(4, 4**4, 1, 2)

    total.consume(DNA)

    half1.consume(DNA)
    half2.consume(DNA)

    t_min, t_max = total.count(DNA)
    print t_min, t_max
    
    h1_min, h1_max = half1.count(DNA)
    print h1_min, h1_max
    
    h2_min, h2_max = half2.count(DNA)
    print h2_min, h2_max
    
    h_min = min(h1_min, h2_min)
    h_max = max(h1_max, h2_max)

    assert t_min == h_min, (t_min, h_min)
    assert t_max == h_max, (t_max, h_max)

def test_n_split():
    total = SplitHashtable(4, 4**4, 0, 1)

    x = []
    for i in range(0, 4**4):
        sub = SplitHashtable(4, 4**4, i, 4**4)
        x.append(sub)

    total.consume(DNA)

    for sub in x:
        sub.consume(DNA)

    minmaxes = [ sub.count(DNA) for sub in x ]
    sub_min = min([ m[0] for m in minmaxes ])
    sub_max = max([ m[1] for m in minmaxes ])

    t_min, t_max = total.count(DNA)

    assert t_min == sub_min
    assert t_max == sub_max

def test_n3_split():
    total = SplitHashtable(4, 4**4, 0, 1)

    x = []
    for i in range(0, 4**4 // 3):
        sub = SplitHashtable(4, 4**4, i, 4**4 // 3)
        x.append(sub)

    total.consume(DNA)

    for sub in x:
        sub.consume(DNA)

    minmaxes = [ sub.count(DNA) for sub in x ]
    sub_min = min([ m[0] for m in minmaxes ])
    sub_max = max([ m[1] for m in minmaxes ])

    t_min, t_max = total.count(DNA)

    assert t_min == sub_min
    assert t_max == sub_max
