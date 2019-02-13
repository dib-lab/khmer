..
   This file is part of khmer, https://github.com/dib-lab/khmer/, and is
   Copyright (C) 2010-2015 Michigan State University
   Copyright (C) 2015-2016 The Regents of the University of California.
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
<https://doi.org/10.12688/f1000research.6924.1>`__.

khmer is free and open source software.

**To install khmer**, you will need a Linux or Mac computer, together with
Python >= 3.5.   See :doc:`our installation docs <user/install>`
for detailed instructions.

**To use khmer**, you will generally need to work at the UNIX command
line.  See :doc:`our command line documentation <user/scripts>`.

We have **additional documentation** in several places, including
`protocols for metagenome and mRNAseq assembly
<https://khmer-protocols.readthedocs.io/>`__ and `recipes for several
common research tasks <https://khmer-recipes.readthedocs.io/>`__.
You might also be interested in :doc:`papers using or citing khmer
<user/biblio>`.

**To get help**, :doc:`please follow this guide <user/getting-help>`.

We welcome contributions to the khmer project!  We are friendly and
supportive of new contributors, and have a :doc:`code of conduct
<dev/CODE_OF_CONDUCT>`.  Please see our :doc:`docs on getting
started on khmer development <dev/getting-started>`.

Details
=======

:Authors: Daniel Standage, Hussien F. Alameldin, Ali Aliyari, Sherine
        Awad, Elmar Bucher, Adam Caldwell, Reed Cartwright, Amanda Charbonneau,
        Lisa Cohen, Bede Constantinides, Michael R. Crusoe, Greg Edvenson,
        Scott Fay, Jacob Fenton, Thomas Fenzl, Jordan Fish, Leonor Garcia-
        Gutierrez, Phillip Garland, Jonathan Gluck, Iván González, Sarah
        Guermond, Jiarong Guo, Aditi Gupta, Tim Head, Joshua R. Herr, Adina
        Howe, Alex Hyer, Andreas Härpfer, Luiz Irber, Shannon EK Joslin, Rhys
        Kidd, Nicole Kingsley, David Lin, Justin Lippi, Tamer Mansour, Pamela
        McA'Nulty, Eric McDonald, Jessica Mizzi, Kevin D. Murray, Joshua R.
        Nahum, Kaben Nanlohy, Russell Neches, Alexander Johan Nederbragt,
        Humberto Ortiz-Zuazaga, Jeramia Ory, Jason Pell, Charles Pepe-Ranney,
        Zachary N Russ, Erich Schwarz, Camille Scott, Josiah Seaman, Ryan
        Shean, Scott Sievert, Jared Simpson, Connor T. Skennerton, James
        Spencer, Ramakrishnan Srinivasan, James A. Stapleton, Joe Stein, Sascha
        Steinbiss, Susan R Steinman, Cait Sydney, Benjamin Taylor, Will
        Trimble, Heather L. Wiencko, Michael Wright, Brian Wyss, Qingpeng
        Zhang, en zyme, C. Titus Brown


:Contact: khmer-project@idyll.org
:GitHub: https://github.com/dib-lab/khmer
:License: BSD

.. toctree::
   :maxdepth: 1

   introduction
   contributors
   citations
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
