SeqAn Changelog
---------------

This file summarizes the changes to the SeqAn library and apps.


Release 2.0.0
~~~~~~~~~~~~~

Major release with many new features and applications. 
Note, the majority of the modules are backward compatible to the previous version.
Some modules, e.g. I/O-modules, have some adapted easier-to-use or unified interfaces.

Library Updates
^^^^^^^^^^^^^^^

- Faster and easier-to-use modules for basic and formatted file I/O:
    - ``stream``
    - ``seq_io``
    - ``bam_io``
    - ``vcf_io``
    - ``gff_io``
- Faster data structures:
    - FMIndex (up to 4X).
    - Packed Strings.
- New alignment modules:
    - X-Drop extension for alignments (``align_extend``)
    - Sequence-profile alignments (``align_profile``)
- New AminoAcid-Dna translation module (``translation``)
- The motif finding module (``find_module``) has been removed.

Infrastructure Updates
^^^^^^^^^^^^^^^^^^^^^^

- The repository has been migrated to GitHub (https://github.com/seqan/seqan).
- Continuous integration builds happen on TravisCI.
- The manual has been migrated to sphinx (http://seqan.readthedocs.org).
- The ``core`` and ``extras`` subfolders have been removed.

New Apps
^^^^^^^^

- ANISE and BASIL
    - Methods for the detection and assembly of inserted sequence in High-Throughput Sequencing Data.

- BS Tools
    - Bisulfite read mapping and SNP and methylation level calling.

- Fiona
    - A parallel and automatic strategy for read error correction.

- Gustaf
    - Generic mUlti-SpliT Alignment Finder.

- Mason 2
    - A read simulator.

- NGS ROI
    - Region of Interest Analysis for NGS Data.

- Samcat
    - Concatenate and convert SAM/BAM files (faster than samtools).

- Seqcons 2
    - Compute consensus from sequences sequences with and without approximate alignment information.

- Yara
    - Yet another read aligner (replaces Masai).


Release 1.4.2
~~~~~~~~~~~~~

Documentation-only release backward compatible with 1.4.1.


Release 1.4.1
~~~~~~~~~~~~~

This minor release should be backward compatible with 1.4. It contains small fixes and many demos for improving the API documentation. Some file format functionality has been added.

Highlights
^^^^^^^^^^

- Many new demos and improved API documentation throughout the library.
- New file format support and tutorials for this functionality: VCF I/O, BED I/O, and improvements to GFF and GTF I/O.

Selected Bug Fixes
^^^^^^^^^^^^^^^^^^

- ``gff_io.h`` does not contain corrupt includes any more
- Gapped X-drop seed extension now works with score matrices such as BLOSUM60.
- SAM writer code now writes ``255`` for invalid ``MAPQ`` and ``0`` for invalid/unapplicable ``TLEN`` instead of ``*``.
- Fix in Postorder ParentLinks VSTree Iterator.
- ``SEQAN_PATH_TO_ROOT()`` can now be used in demo programs.
- Removing duplicate definition of ``SEQAN_ENABLE_TESTING`` in build system.
- Write support for ``char *`` for ``BamTagsDict``.
- Fix in ``StringEnumerator``.
- Fix writing out of file extension when writing KNIME plugins.

Release 1.4
~~~~~~~~~~~

Highlights
^^^^^^^^^^

- New read mappers applications Masai and RazerS 3.
- Extended and more robust I/O functionality in ``stream``, ``seq_io``, ``bam_io``, and ``gff_io``.
- Module arg_parse creates improved command line help and supports workflow engine integration.
    - Also see https://github.com/genericworkflownodes
- Greatly improved alignment module with better performance and interfaces.
- Greatly improved build system, ``find_package(SeqAn)`` for your CMake build systems.

New Apps
^^^^^^^^

- ALF
    - Alignment free sequence comparison.

- Breakpoint Calculator
    - Breakpoint computation for genomic alignments.

- Masai
    - Fast index-based read mapper.

- RazerS 3
    - Fast filtration-based, parallel read mapper.

- SnpStore
    - SNP and small indel calling.

Major App Updates
^^^^^^^^^^^^^^^^^

- All applications now use the ArgumentParser and have better CLI help.

- Rabema
    - Rewritten from scratch, includes BAM support.
    - Greatly lowered memory requirements.

- SeqCons
    - Fixing input bugs, supports SAM I/O now.

- Stellar
    - Major update improving running time, including bug fixes, and
      allowing for various alphabet types.

- MicroRazerS
    - Adding support for SAM output.

Major Library Updates
^^^^^^^^^^^^^^^^^^^^^

- Modules ``seq_io``, ``bam_io``, ``gff_io`` with I/O functionality.
- FM Index in module ``index``.
- Rewritten ``align`` module with better performance, more consistent interfaces.
- Split alignment module ``align_split``.
- Metaprogramming: introducing ``EnableIf``, ``DisableIf``, ``EnableIf2``, and ``DisableIf2`` metafunctions
- Module ``alignment_free`` for alignment free sequence comparison.
- Module ``journaled_set`` for managing many similar sequences.
- Faster open addressing q-gram index.
- generic support for memory mapped files via FileMapping class
- Adding module ``parallel`` with atomic operations in C++98.
- Greatly improved FragmentStore documentation.
- Adding ``position()``, ``operator-()``, ``operator[]`` with proxy functionality and relation operators to journaled string iterator.
- Pigeonhole-based filter algorithm.
- Parallel repeat finding.
- Clang support, C++11 support

Major Library Bug Fixes
^^^^^^^^^^^^^^^^^^^^^^^

- Fixing repeat finding on Dna5Q.
- Fixing insert size computation in store_all.h
- Fixing memory initialization problem in ``appendValue()`` for Block String.
- Default constructor of Iter modified, such that data_container and data_position are initialized.
- Fixed error loading Fasta on Windows.
- Fixed wrong StringSet size types, allow to easily subclass Alloc strings
- Now supports SAM files with missing read sequences
- Fixing SeqAn code for C++11
- FragmentStore fixes.

Miscellaneous
^^^^^^^^^^^^^

- Experimental support added platforms for ICC and PGI compilers.
- Experimental support for CUDA.
- Build System
    - Large updates to build system.
    - Includes ``FindSeqAn.cmake`` for easily using SeqAn in your own CMake build system.
    - Packaging now based on CPack
- Xcode plugin for MacPorts LLVM/Clang in Xcode 3 and 4
- Improved code generator ``skel.py``.
- Many minor bug fixes
- Cleaned code base
- Added test cases (e.g. Stellar)
- Improved documentation and added examples (Mason, Rabema, RazerS, etc.)
- Improving coding style compliance of Array String implementation.
- Various tool improvements (e.g. RazerS 3)
- Performance improvements.
