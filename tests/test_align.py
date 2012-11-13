import khmer

def test_alignnocov():
   ch = khmer.new_counting_hash(10, 1048576, 1)
   read = "ACCTAGGTTCGACATGTACC"
   aligner = khmer.new_readaligner(ch)
   for i in range(20):
      ch.consume("AGAGGGAAAGCTAGGTTCGACAAGTCCTTGACAGAT")
   ch.consume("ACCTAGGTTCGACATGTACC")
   graphAlign, readAlign = aligner.align(read)

   assert readAlign == 'ACCTAGGTTCGACATGTACC'
   assert graphAlign  == 'ACCTAGGTTCGACATGTACC'


def test_readalign():
   ch = khmer.new_counting_hash(10, 1048576, 1)
   read = "ACCTAGGTTCGACATGTACC"
   aligner = khmer.new_readaligner(ch, 1, 20)
   for i in range(20):
      ch.consume("AGAGGGAAAGCTAGGTTCGACAAGTCCTTGACAGAT")
   ch.consume("GCTTTTAAAAAGGTTCGACAAAGGCCCGGG")
   graphAlign, readAlign = aligner.align(read)

   assert readAlign == 'ACCTAGGTTCGACATGTACC'
   assert graphAlign  == 'AGCTAGGTTCGACAAGT-CC'

def test_alignerrorregion():
   ch = khmer.new_counting_hash(10, 1048576, 1)
   read = "AAAAAGTTCGAAAAAGGCACG"
   aligner = khmer.new_readaligner(ch, 1, 20, 11)
   for i in range(20):
      ch.consume("AGAGGGAAAGCTAGGTTCGACAAGTCCTTGACAGAT")
   ch.consume("ACTATTAAAAAAGTTCGAAAAAGGCACGGG")
   graphAlign, readAlign = aligner.align(read)

   assert readAlign == ''
   assert graphAlign  == ''
