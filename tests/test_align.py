import khmer

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

test_readalign()
