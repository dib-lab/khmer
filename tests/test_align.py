import khmer

def test_readalign():
   ch = khmer.new_counting_hash(10, 1048576, 1)
   read = "ACCTAGGTTCGACATGTACC"
   aligner = khmer.new_readaligner(ch)
   ch.consume("AGAGGGAAAGCTAGGTTCGACAAGTCCTTGACAGAT")
   ch.consume("GCTTTTAAAAAGGTTCGACAAAGGCCCGGG")
   graphAlign, readAlign, score = aligner.align(read)

   assert graphAlign == 'AGCTAGGTTCGACAA-GT-CC'
   assert readAlign == 'ACCTAGGTTCGAC-ATGTACC'
   assert score == 75
