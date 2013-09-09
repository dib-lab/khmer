import khmer
ht = khmer.new_hashbits(32, 1, 1)

ht.add_tag('CGCAGGCTGGATTCTAGAGGCAGAGGTGAGCT')
ht.save_tagset('simple.tagset')
