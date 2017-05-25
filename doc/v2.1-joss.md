---
title: 'khmer release v2.1: software for biological sequence analysis'
tags:
    - bioinformatics
    - sequence analysis
    - Bloom filter
    - Count-Min sketch
    - de Bruijn graph
    - assembly
    - graph traversal
    - streaming
    - quality control
authors:
    - name: Daniel Standage
      orcid: 0000-0003-0342-8531
      affiliation: 1
    - name: Ali Aliyari
      orcid: 0000-0003-0925-4886
      affiliation: 2
    - name: Lisa J. Cohen
      orcid: 0000-0002-3600-7218
      affiliation:
          - 1
          - 3
    - name: Michael R. Crusoe
      orcid: 0000-0002-2961-9670
      affiliation: 4
    - name: Tim Head
      orcid: 0000-0003-0931-3698
      affiliation: 5
    - name: Luiz Irber
      orcid: 0000-0003-4371-9659
      affiliation:
          - 1
          - 6
    - name: Shannon EK Joslin
      orcid: 0000-0001-5470-1193
      affiliation: 2
    - name: N. B. Kingsley
      orcid: 0000-0001-8574-8592
      affiliation: 2
    - name: Kevin D. Murray
      orcid: 0000-0002-2466-1917
      affiliation: 7
    - name: Russell Neches
      orcid: 0000-0002-2055-8381
      affiliation: 8
    - name: Camille Scott
      orcid: 0000-0001-8822-8779
      affiliation:
          - 1
          - 6
    - name: Ryan Shean
      affiliation: 9
    - name: Sascha Steinbiss
      orcid: 0000-0002-2151-0574
      affiliation: 10
    - name: Cait Sydney
      affiliation: 11
    - name: C. Titus Brown
      orcid: 0000-0001-6001-2677
      affiliation:
          - 1
          - 12
affiliations:
    - name: Lab for Data Intensive Biology; School of Veterinary Medicine; University of California, Davis
      index: 1
    - name: Integrative Genetics and Genomics Graduate Group; University of California, Davis
      index: 2
    - name: Molecular, Cellular, and Integrative Physiology Graduate Group; University of California, Davis
      index: 3
    - name: Common Workflow Language Project
      index: 4
    - name: Wild Tree Tech
      index: 5
    - name: Computer Science Graduate Group; University of California, Davis
      index: 6
    - name: ARC Centre of Excellence in Plant Energy Biology; Australian National University
      index: 7
    - name: Microbiology Graduate Group; Univerity of California, Davis
      index: 8
    - name: University of Washington
      index: 9
    - name: Debian Project
      index: 10
    - name: Google, Inc.
      index: 11
    - name: Department of Population Health and Reproduction; School of Veterinary Medicine; University of California, Davis
      index: 12
date: 2017-05-25
---

# Summary

The khmer software is a set of command-line tools built around a Python library designed for analysis of large DNA sequence collections.
Functionality in khmer has primarily been motivated by scaling issues with (meta)genome and (meta)transcriptome assembly (Crusoe *et al.*, 2015).
khmer provides convenient access to several k-mer based operations on DNA sequence collections, such as abundance filtering, error trimming, assembly graph partitioning, and most notably, abundance normalization of reads (Brown *et al.*, 2012) and streaming error trimming of reads (Zhang, Awad, and Brown, 2015).
All of these operations utilize khmer's implementation of two primary data structures, the Bloom filter and the Count-Min Sketch, for efficient probabalistic storage of k-mer presence or k-mer abundance, respectively (Pell *et al.*, 2012; Zhang *et al.*, 2014).

Release version 2.1 of the khmer software includes several new features that extend its utility to a wider set of sequence processing and analysis problems.
These include the following:
support for variable-coverage trimming of sequence reads;
support for k > 32 using the non-reversible hash function MurmurHash3 (Appleby, 2010);
a new optional Count-Min Sketch implementation providing increased storage efficiency;
support for assembly directly from a k-mer graph;
a script for computing a compact de Bruijn graph from a k-mer graph;
and several examples of khmer's Python and C++ APIs for those interested in using and extending the library.

# References

**Appleby A** (2010) SMHasher. *GitHub*, https://github.com/aappleby/smhasher, accessed May 22, 2017.

**Brown CT, Howe A, Zhang Q, Pyrkosz AB, Brom TH** (2012) A Reference-Free Algorithm for Computational Normalization of Shotgun Sequencing Data. *arXiv*, arXiv:1203.4802, https://arxiv.org/abs/1203.4802.

**Crusoe MR, Alameldin HF, Awad S, Bucher E, *et al.*** (2015) The khmer software package: enabling efficient nucleotide sequence analysis. *F1000 Research*, doi:10.12688/f1000research.6924.1.

**Pell J, Hintze A, Canino-Koning R, Howe A, Tiedje JM, Brown CT** (2012). Scaling metagenome sequence assembly with probabilistic de Bruijn graphs. *Proc Natl Acad Sci U S A*, 109(33):13272-7, doi:10.1073/pnas.1121464109.

**Zhang Q, Awad S, Brown CT** (2015) Crossing the streams: a framework for streaming analysis of short DNA sequencing reads. *PeerJ PrePrints*, 3:e890v1, doi:10.7287/peerj.preprints.890v1.

**Zhang Q, Pell J, Canino-Koning R, Howe A, Brown CT** (2014) These Are Not the K-mers You Are Looking For: Efficient Online K-mer Counting Using a Probabilistic Data Structure. *PLoS ONE*, e101271, doi:10.1371/journal.pone.0101271.
