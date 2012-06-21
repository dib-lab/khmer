import khmer

def test_readalign():
   hb = khmer.new_hashbits(10, 1048576, 1)
   read = "ACCTAGGTTCGACATGTACC"
   aligner = khmer.new_readaligner(hb)
   hb.consume("AGAGGGAAAGCTAGGTTCGACAAGTCCTTGACAGAT")
   hb.consume("GCTTTTAAAAAGGTTCGACAAAGGCCCGGG")
   graphAlign, readAlign, score = aligner.align(read)

   assert graphAlign == 'AGCTAGGTTCGACAA-GT-CC'
   assert readAlign == 'ACCTAGGTTCGAC-ATGTACC'
   assert score == 75
