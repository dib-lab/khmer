import khmer
import screed


test_file="test-data/simple_1.fa"
test_bsize=10
test_ksize=21

DNA = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"


class Test_KmerCardinality(object):
    
    def __init__(self):

        self.counter= khmer.KmerCardinality(test_bsize,test_ksize)
 

test_counter = Test_KmerCardinality() 


def test_b_size():
    assert 4<=test_bsize<=16

def test_k_size():
    assert 21<=test_ksize<=32

def test_alpha():
	assert test_counter.counter.alpha <= 1

def test_consume():
	assert str==type(DNA)









