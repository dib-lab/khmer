import khmer

K = 4
HT_SIZE=1000000000
N_HT =3 

ht = khmer.new_hashbits(K, HT_SIZE, N_HT)

ht.count('AAAA')
ht.count('CCCC')
#ht.consume('ATGAACCAGAGATTAGACCCCCTT')
print ht.get('CCCC')
print ht.get('CCTT')

print ht.n_occupied()
