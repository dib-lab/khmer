.. vim: set filetype=rst

Known Issues
============

Some users have reported that normalize-by-median.py will utilize more
memory than it was configured for. This is being investigated in
https://github.com/dib-lab/khmer/issues/266

BSD users will need to download khmer and install using `make CC="cc -fPIC"
install`. https://github.com/dib-lab/khmer/issues/719

Some very badly formatted FASTA/FASTQ files will silently be accepted without
error. https://github.com/dib-lab/khmer/issues/884
