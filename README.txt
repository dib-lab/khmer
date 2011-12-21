PDBG BRANCH NOTE: To reproduce the results of the paper, "Scaling metagenome 
sequence assembly with probabilistic de Bruijn graphs," follow  the 
instructions below then go to the figuregen directory and take a 
look at the README there.

Welcome to khmer, k-mer counting, filtering and graph traversal FTW!

As of August 2011, there's a khmer mailing list at librelist.com that
you can use to get help with khmer.  To sign up, just shoot
'khmer@librelist.com' an e-mail and it will subscribe you; then send
your question/comment there.

IMPORTANT NOTE:

khmer is *pre-publication* and *research* software, so please keep in
mind that (a) the code may have undiscovered bugs in it, (b) you
should cite us, and (c) you should get in touch if you need to cite
us, as we are writing up the project.

INSTRUCTIONS:

'make all' to build.

'cd python && python setup.py test' to test (OR set PYTHONPATH, 'make test').

'make doc' to build docs.  You'll need Sphinx installed.

khmer is under the BSD license; see doc/LICENSE.txt.  Distribution,
modification and redistribution, incorporation into other software,
and pretty much everything else is allowed.

CTB 8/2011.
