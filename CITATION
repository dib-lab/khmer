..
   This file is part of khmer, https://github.com/dib-lab/khmer/, and is
   Copyright (C) 2014-2015 Michigan State University
   Copyright (C) 2015 The Regents of the University of California.
   It is licensed under the three-clause BSD license; see LICENSE.
   Contact: khmer-project@idyll.org

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are
   met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.

    * Neither the name of the Michigan State University nor the names
      of its contributors may be used to endorse or promote products
      derived from this software without specific prior written
      permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
   HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

   Contact: khmer-project@idyll.org

.. If you update this file then you may need to update the citations in
   khmer/khmer_args.py as well

*********
Citations
*********

Software Citation
=================

If you use the khmer software, you must cite:

   Crusoe et al., The khmer software package: enabling efficient nucleotide
   sequence analysis. 2015. https://doi.org/10.12688/f1000research.6924.1

.. code-block:: tex

  @article{khmer2015,
     author = "Crusoe, Michael R. and Alameldin, Hussien F. and Awad, Sherine
  and Bucher, Elmar and Caldwell, Adam and Cartwright, Reed and Charbonneau,
  Amanda and Constantinides, Bede and Edvenson, Greg and Fay, Scott and Fenton,
  Jacob and Fenzl, Thomas and Fish, Jordan and Garcia-Gutierrez, Leonor and
  Garland, Phillip and Gluck, Jonathan and GonzÃ¡lez, IvÃ¡n and Guermond, Sarah
  and Guo, Jiarong and Gupta, Aditi and Herr, Joshua R. and Howe, Adina and
  Hyer, Alex and HÃ¤rpfer, Andreas and Irber, Luiz and Kidd, Rhys and Lin, David
  and Lippi, Justin and Mansour, Tamer and McA'Nulty, Pamela and McDonald, Eric
  and Mizzi, Jessica and Murray, Kevin D. and Nahum, Joshua R. and Nanlohy,
  Kaben and Nederbragt, Alexander Johan and Ortiz-Zuazaga, Humberto and Ory,
  Jeramia and Pell, Jason and Pepe-Ranney, Charles and Russ, Zachary N and
  Schwarz, Erich and Scott, Camille and Seaman, Josiah and Sievert, Scott and
  Simpson, Jared and Skennerton, Connor T. and Spencer, James and Srinivasan,
  Ramakrishnan and Standage, Daniel and Stapleton, James A. and Stein, Joe and
  Steinman, Susan R and Taylor, Benjamin and Trimble, Will and Wiencko, Heather
  L. and Wright, Michael and Wyss, Brian and Zhang, Qingpeng and zyme, en and
  Brown, C. Titus"
     title = "The khmer software package: enabling efficient nucleotide
  sequence analysis",
     year = "2015",
     month = "08",
     publisher = "F1000",
     url = "https://doi.org/10.12688/f1000research.6924.1"
  }

If you use any of our published scientific methods you should *also*
cite the relevant paper(s) as directed below. Additionally some scripts use
the `SeqAn library <http://www.seqan.de>`_ for read parsing: the full citation
for that library is also included below.

To see a quick summary of papers for a given script just run it without using
any command line arguments.

Graph partitioning and/or compressible graph representation
===========================================================

The :program:`load-graph.py`, :program:`partition-graph.py`,
and :program:`find-knots.py` scripts are part of the compressible graph
representation and partitioning algorithms described in:

   Pell J, Hintze A, Canino-Koning R, Howe A, Tiedje JM, Brown CT.
   Scaling metagenome sequence assembly with probabilistic de Bruijn graphs
   Proc Natl Acad Sci U S A. 2012 Aug 14;109(33):13272-7.
   https://doi.org/10.1073/pnas.1121464109.
   PMID: 22847406

.. code-block:: tex

  @article{Pell2012,
      author = "Pell, Jason and Hintze, Arend and Canino-Koning, Rosangela and
  Howe, Adina and Tiedje, James M. and Brown, C. Titus",
      title = "Scaling metagenome sequence assembly with probabilistic de Bruijn
  graphs",
      volume = "109",
      number = "33",
      pages = "13272-13277",
      year = "2012",
      doi = "10.1073/pnas.1121464109",
      abstract ="Deep sequencing has enabled the investigation of a wide range of
  environmental microbial ecosystems, but the high memory requirements for de
  novo assembly of short-read shotgun sequencing data from these complex
  populations are an increasingly large practical barrier. Here we introduce a
  memory-efficient graph representation with which we can analyze the k-mer
  connectivity of metagenomic samples. The graph representation is based on a
  probabilistic data structure, a Bloom filter, that allows us to efficiently
  store assembly graphs in as little as 4 bits per k-mer, albeit inexactly. We
  show that this data structure accurately represents DNA assembly graphs in low
  memory. We apply this data structure to the problem of partitioning assembly
  graphs into components as a prelude to assembly, and show that this reduces the
  overall memory requirements for de novo assembly of metagenomes. On one soil
  metagenome assembly, this approach achieves a nearly 40-fold decrease in the
  maximum memory requirements for assembly. This probabilistic graph
  representation is a significant theoretical advance in storing assembly graphs
  and also yields immediate leverage on metagenomic assembly.",
      URL = "http://www.pnas.org/content/109/33/13272.abstract",
      eprint = "http://www.pnas.org/content/109/33/13272.full.pdf+html",
      journal = "Proceedings of the National Academy of Sciences"
  }

Digital normalization
=====================

The :program:`normalize-by-median.py` and :program:`count-median.py` scripts
are part of the digital normalization algorithm, described in:

   A Reference-Free Algorithm for Computational Normalization of
   Shotgun Sequencing Data
   Brown CT, Howe AC, Zhang Q, Pyrkosz AB, Brom TH
   arXiv:1203.4802 [q-bio.GN]
   http://arxiv.org/abs/1203.4802

.. code-block:: tex

  @unpublished{diginorm,
      author = "C. Titus Brown and Adina Howe and Qingpeng Zhang and Alexis B.
  Pyrkosz and Timothy H. Brom",
      title = "A Reference-Free Algorithm for Computational Normalization of
  Shotgun Sequencing Data",
      year = "2012",
      eprint = "arXiv:1203.4802",
      url = "http://arxiv.org/abs/1203.4802",
  }

Efficient k-mer error trimming
==============================

The :program:`script trim-low-abund.py` is described in:

   Crossing the streams: a framework for streaming analysis of short DNA
   sequencing reads
   Zhang Q, Awad S, Brown CT
   https://doi.org/10.7287/peerj.preprints.890v1

.. code-block:: tex

  @unpublished{semistream,
      author = "Qingpeng Zhang and Sherine Awad and C. Titus Brown",
      title = "Crossing the streams: a framework for streaming analysis of
          short DNA sequencing reads",
      year = "2015",
      eprint = "PeerJ Preprints 3:e1100",
      url = "https://doi.org/10.7287/peerj.preprints.890v1"
  }

K-mer counting
==============

The :program:`abundance-dist.py`, :program:`filter-abund.py`, and
:program:`load-into-counting.py` scripts implement the probabilistic k-mer
counting described in:

   These Are Not the K-mers You Are Looking For: Efficient Online K-mer
   Counting Using a Probabilistic Data Structure
   Zhang Q, Pell J, Canino-Koning R, Howe AC, Brown CT.
   https://doi.org/10.1371/journal.pone.0101271

.. code-block:: tex

  @article{khmer-counting,
      author = "Zhang, Qingpeng AND Pell, Jason AND Canino-Koning, Rosangela
  AND Howe, Adina Chuang AND Brown, C. Titus",
      journal = "PLoS ONE",
      publisher = "Public Library of Science",
      title = "These Are Not the K-mers You Are Looking For: Efficient Online
  K-mer Counting Using a Probabilistic Data Structure",
      year = "2014",
      month = "07",
      volume = "9",
      url = "https://doi.org/10.1371/journal.pone.0101271",
      pages = "e101271",
      abstract = "<p>K-mer abundance analysis is widely used for many purposes in
  nucleotide sequence analysis, including data preprocessing for de novo
  assembly, repeat detection, and sequencing coverage estimation. We present the
  khmer software package for fast and memory efficient <italic>online</italic>
  counting of k-mers in sequencing data sets. Unlike previous methods based on
  data structures such as hash tables, suffix arrays, and trie structures, khmer
  relies entirely on a simple probabilistic data structure, a Count-Min Sketch.
  The Count-Min Sketch permits online updating and retrieval of k-mer counts in
  memory which is necessary to support online k-mer analysis algorithms. On
  sparse data sets this data structure is considerably more memory efficient than
  any exact data structure. In exchange, the use of a Count-Min Sketch introduces
  a systematic overcount for k-mers; moreover, only the counts, and not the
  k-mers, are stored. Here we analyze the speed, the memory usage, and the
  miscount rate of khmer for generating k-mer frequency distributions and
  retrieving k-mer counts for individual k-mers. We also compare the performance
  of khmer to several other k-mer counting packages, including Tallymer,
  Jellyfish, BFCounter, DSK, KMC, Turtle and KAnalyze. Finally, we examine the
  effectiveness of profiling sequencing error, k-mer abundance trimming, and
  digital normalization of reads in the context of high khmer false positive
  rates. khmer is implemented in C++ wrapped in a Python interface, offers a
  tested and robust API, and is freely available under the BSD license at
  github.com/dib-lab/khmer.</p>",
      number = "7",
      doi = "10.1371/journal.pone.0101271"
  }

FASTA and FASTQ reading
=======================

Several scripts use the SeqAn library for FASTQ and FASTA reading as described
in:

   SeqAn An efficient, generic C++ library for sequence analysis
   Döring A, Weese D, Rausch T, Reinert K.
   https://doi.org/10.1186/1471-2105-9-11

.. code-block:: tex

  @Article{SeqAn,
    AUTHOR = {Doring, Andreas and Weese, David and Rausch, Tobias and Reinert,
      Knut},
    TITLE = {SeqAn An efficient, generic C++ library for sequence analysis},
    JOURNAL = {BMC Bioinformatics},
    VOLUME = {9},
    YEAR = {2008},
    NUMBER = {1},
    PAGES = {11},
    URL = {http://www.biomedcentral.com/1471-2105/9/11},
    DOI = {10.1186/1471-2105-9-11},
    PubMedID = {18184432},
    ISSN = {1471-2105},
    ABSTRACT = {BACKGROUND: The use of novel algorithmic techniques is pivotal
    to many important problems in life science. For example the sequencing of
    the human genome [1] would not have been possible without advanced assembly
    algorithms. However, owing to the high speed of technological progress and
    the urgent need for bioinformatics tools, there is a widening gap between
    state-of-the-art algorithmic techniques and the actual algorithmic
    components of tools that are in widespread use. RESULTS: To remedy this
    trend we propose the use of SeqAn, a library of efficient data types and
    algorithms for sequence analysis in computational biology. SeqAn comprises
    implementations of existing, practical state-of-the-art algorithmic
    components to provide a sound basis for algorithm testing and development.
    In this paper we describe the design and content of SeqAn and demonstrate
    its use by giving two examples. In the first example we show an application
    of SeqAn as an experimental platform by comparing different exact string
    matching algorithms. The second example is a simple version of the well-
    known MUMmer tool rewritten in SeqAn. Results indicate that our
    implementation is very efficient and versatile to use. CONCLUSION: We
    anticipate that SeqAn greatly simplifies the rapid development of new
    bioinformatics tools by providing a collection of readily usable, well-
    designed algorithmic components which are fundamental for the field of
    sequence analysis. This leverages not only the implementation of new
    algorithms, but also enables a sound analysis and comparison of existing
    algorithms.},
  }

.. vim: set filetype=rst:
