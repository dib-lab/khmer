KHMER=/Users/t/dev/khmer
SCRIPTS=$(KHMER)/scripts
DATA=$(KHMER)/tests/test-data/random-20-a.fa
BASEDATA=$(shell basename $(DATA))
BASE=random

K=20

all: $(BASE).group0000.fa

clean:
	rm -f $(BASE).subset.*.pmap $(BASE).pmap.merged $(BASE).ht \
		$(BASE).tagset $(BASE).info $(BASE).group????.fa \
		$(BASEDATA).part $(BASE).dist

$(BASE).ht: load_graph

$(BASE).tagset: load_graph

load_graph: $(DATA)
	$(SCRIPTS)/load-graph.py -k $(K) -x 2e6 -k 20 -N 4 $(BASE) $(DATA)

$(BASE).subset.0.pmap: $(BASE).ht $(BASE).tagset
	$(SCRIPTS)/partition-graph.py -T 4 --subset-size 1e2 $(BASE)

$(BASE).pmap.merged: $(BASE).subset.0.pmap
	$(SCRIPTS)/merge-partitions.py -k $(K) $(BASE)
	
$(BASEDATA).part: $(BASE).pmap.merged
	$(SCRIPTS)/annotate-partitions.py -k $(K) $(BASE) $(DATA)

$(BASE).group0000.fa: $(BASEDATA).part
	$(SCRIPTS)/extract-partitions.py $(BASE) $(BASEDATA).part
