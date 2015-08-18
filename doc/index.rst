.. vim: set filetype=rst

###################################################################
The khmer software for advanced biological sequencing data analysis
###################################################################

What is khmer?
==============

The khmer software is a set of command-line tools for working with DNA
shotgun sequencing data from genomes, transcriptomes, metagenomes, and
single cells.  khmer can make *de novo* assemblies faster, and
sometimes better. khmer can also identify (and fix) problems with
shotgun data.  You can read more about khmer in `our software paper
<https://dx.doi.org/10.12688/f1000research.6924.1>`__.

khmer is free and open source software.

**To install khmer**, you will need a Linux or Mac computer, together with
Python 2.7 or Python 3.x.   See `our installation docs <user/install.html>`__
for detailed instructions.

**To use khmer**, you will generally need to work at the UNIX command
line.  See `our command line documentation <user/scripts.html>`__.

We have **additional documentation** in several places, including
`protocols for metagenome and mRNAseq assembly
<https://khmer-protocols.readthedocs.org/>`__ and `recipes for several
common research tasks <https://khmer-recipes.readthedocs.org/>`__.
You might also be interested in `papers using or citing khmer
<user/biblio.html>`__.

**To get help**, please ask questions on `the khmer mailing list
<http://lists.idyll.org/listinfo/khmer>`__ or `post an issue on
GitHub <https://github.com/dib-lab/khmer/issues/new>`__.

We welcome contributions to the khmer project!  We are friendly and
supportive of new contributors, and have a `code of conduct
<dev/CODE_OF_CONDUCT.html>`__.  Please see our `docs on getting
started on khmer development <dev/getting-started.html>`__.

Details
=======

:Authors: Michael R. Crusoe, Hussien F. Alameldin, Sherine Awad, Elmar
        Bucher, Adam Caldwell, Reed Cartwright, Amanda Charbonneau, Bede
        Constantinides, Greg Edvenson, Scott Fay, Jacob Fenton, Thomas Fenzl,
        Jordan Fish, Leonor Garcia-Gutierrez, Phillip Garland, Jonathan Gluck,
        Iván González, Sarah Guermond, Jiarong Guo, Aditi Gupta, Joshua R.
        Herr, Adina Howe, Alex Hyer, Andreas Härpfer, Luiz Irber, Rhys Kidd,
        David Lin, Justin Lippi, Tamer Mansour, Pamela McA'Nulty, Eric
        McDonald, Jessica Mizzi, Kevin D. Murray, Joshua R. Nahum, Kaben
        Nanlohy, Alexander Johan Nederbragt, Humberto Ortiz-Zuazaga, Jeramia
        Ory, Jason Pell, Charles Pepe-Ranney, Zachary N Russ, Erich Schwarz,
        Camille Scott, Josiah Seaman, Scott Sievert, Jared Simpson, Connor T.
        Skennerton, James Spencer, Ramakrishnan Srinivasan, Daniel Standage,
        James A. Stapleton, Joe Stein, Susan R Steinman, Benjamin Taylor, Will
        Trimble, Heather L. Wiencko, Michael Wright, Brian Wyss, Qingpeng
        Zhang, en zyme, C. Titus Brown

:Contact: khmer-project@idyll.org
:GitHub: https://github.com/dib-lab/khmer
:Chat: https://gitter.im/dib-lab/khmer
:License: BSD

.. toctree::
   :maxdepth: 1

   introduction
   contributors
   citations
   whats-new-2.0
   release-notes/index
   user/index
   dev/index
   roadmap
   LICENSE

There are two mailing lists dedicated to khmer, an announcements-only list and
a discussion list. To search their archives and sign-up for them, please visit
the following URLs:
    
    * Discussion: http://lists.idyll.org/listinfo/khmer

    * Announcements: http://lists.idyll.org/listinfo/khmer-announce

The archives for the khmer mailing list are available at: 
http://lists.idyll.org/pipermail/khmer/

khmer development was initially supported by AFRI Competitive Grant
no.  `2010-65205-20361
<http://ged.msu.edu/downloads/2009-usda-vertex.pdf>`__ from the USDA
NIFA, and is now funded by the National Human Genome Research
Institute of the National Institutes of Health under Award Number
`R01HG007513 <http://ged.msu.edu/downloads/2012-bigdata-nsf.pdf>`__
through May 2016, both to C. Titus Brown.  More recently, we have
received support from the Gordon and Betty Moore Foundation under
Award number GBMF4551.
