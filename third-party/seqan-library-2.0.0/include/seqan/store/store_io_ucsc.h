// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_INCLUDE_STORE_STORE_IO_UCSC_H
#define SEQAN_INCLUDE_STORE_STORE_IO_UCSC_H

namespace SEQAN_NAMESPACE_MAIN {

//////////////////////////////////////////////////////////////////////////////
// add a UCSC record to FragmentStore

template <typename TSpec, typename TConfig>
inline void
_storeAnnotationRecord(
    FragmentStore<TSpec, TConfig> & fragStore,
    UcscRecord const & record,
    UcscKnownGene)
{
    typedef FragmentStore<TSpec, TConfig>               TFragmentStore;
    typedef typename TFragmentStore::TAnnotationStore   TAnnotationStore;
    typedef typename Value<TAnnotationStore>::Type      TAnnotation;
    typedef typename TAnnotation::TId                   TId;

    SEQAN_ASSERT_EQ(length(fragStore.annotationStore), length(fragStore.annotationNameStore));

    // add transcript and CDS
    TId transId = TAnnotation::INVALID_ID;
    _storeAppendAnnotationName(fragStore, transId, record.transName, (TId) TFragmentStore::ANNO_MRNA);
    TId cdsId = length(fragStore.annotationStore);
    appendName(fragStore.annotationNameStore, record.proteinName, fragStore.annotationNameStoreCache);

    resize(fragStore.annotationStore, cdsId + 1 + length(record.exonBegin), Generous());
    resize(fragStore.annotationNameStore, cdsId + 1 + length(record.exonBegin), Generous());

    // add contig name
    TId contigId = TAnnotation::INVALID_ID;
    _storeAppendContig(fragStore, contigId, record.contigName);

    // update transcript data
    TAnnotation & transcript = fragStore.annotationStore[transId];
    TId geneId = transcript.parentId;
    if (geneId == TAnnotation::INVALID_ID)
        geneId = 0;
    transcript.parentId = geneId;
    transcript.contigId = contigId;
    transcript.typeId = TFragmentStore::ANNO_MRNA;

    // add CDS entry
    TAnnotation & cds = fragStore.annotationStore[cdsId];
    cds.parentId = transId;
    cds.contigId = contigId;
    cds.typeId = TFragmentStore::ANNO_CDS;
    cds.beginPos = record.cdsBegin;
    cds.endPos = record.cdsEnd;
    _adjustParent(transcript, cds);

    // add exons
    for (unsigned i = 0; i < length(record.exonBegin); ++i)
    {
        TAnnotation & exon = fragStore.annotationStore[cdsId + 1 + i];
        exon.parentId = transId;
        exon.contigId = contigId;
        exon.typeId = TFragmentStore::ANNO_EXON;
        exon.beginPos = record.exonBegin[i];
        exon.endPos = record.exonEnds[i];
        _adjustParent(transcript, exon);
    }
    if (geneId != 0)
        _adjustParent(fragStore.annotationStore[geneId], transcript);
}

template <typename TSpec, typename TConfig>
inline void
_storeAnnotationRecord(
    FragmentStore<TSpec, TConfig> & fragStore,
    UcscRecord const & record,
    UcscKnownIsoforms const &)
{
    typedef FragmentStore<TSpec, TConfig>               TFragmentStore;
    typedef typename TFragmentStore::TAnnotationStore   TAnnotationStore;
    typedef typename Value<TAnnotationStore>::Type      TAnnotation;
    typedef typename TAnnotation::TId                   TId;

    SEQAN_ASSERT_EQ(length(fragStore.annotationStore), length(fragStore.annotationNameStore));

    TId geneId = TAnnotation::INVALID_ID;
    TId transId = TAnnotation::INVALID_ID;

    // add transcript and CDS
    _storeAppendAnnotationName(fragStore, geneId, record.transName, (TId) TFragmentStore::ANNO_GENE);
    _storeAppendAnnotationName(fragStore, transId, record.contigName, (TId) TFragmentStore::ANNO_MRNA);

    // set parent link locus->root
    TAnnotation & locus = fragStore.annotationStore[geneId];
    locus.parentId = 0;
    locus.typeId = TFragmentStore::ANNO_GENE;

    // set parent link transcript->locus
    TAnnotation & transcript = fragStore.annotationStore[transId];
    transcript.parentId = geneId;
    transcript.typeId = TFragmentStore::ANNO_MRNA;

    _adjustParent(locus, transcript);
}

// support for dynamically chosen file formats
template <typename TSpec, typename TConfig>
inline void
_storeAnnotationRecord(
    FragmentStore<TSpec, TConfig> & /*store*/,
    UcscRecord const & /*record*/,
    TagSelector<> const & /*format*/)
{
    SEQAN_FAIL("AnnotationStore: File format not specified.");
}

template <typename TSpec, typename TConfig, typename TTagList>
inline void
_storeAnnotationRecord(
    FragmentStore<TSpec, TConfig> & store,
    UcscRecord const & record,
    TagSelector<TTagList> const & format)
{
    typedef typename TTagList::Type TFormat;

    if (isEqual(format, TFormat()))
        _storeAnnotationRecord(store, record, TFormat());
    else
        _storeAnnotationRecord(store, record, static_cast<typename TagSelector<TTagList>::Base const &>(format));
}



//////////////////////////////////////////////////////////////////////////////
// read a whole UCSC stream into FragmentStore

template <typename TFSSpec, typename TConfig, typename TSpec>
inline void
readRecords(FragmentStore<TFSSpec, TConfig> & fragStore,
            FormattedFile<Ucsc, Input, TSpec> & file)
{
    typedef FormattedFile<Ucsc, Input, TSpec> TFile;
    typename DirectionIterator<TFile, Input>::Type iter = directionIterator(file, Input());

    if (atEnd(iter))
        return;

    // get first character from the stream
    UcscRecord record;
    UcscIOContext ctx;

    refresh(fragStore.contigNameStoreCache);
    refresh(fragStore.annotationNameStoreCache);
    refresh(fragStore.annotationTypeStoreCache);

    while (!atEnd(iter))
    {
        readRecord(record, ctx, iter, format(file));
        _storeAnnotationRecord(fragStore, record, format(file));
    }
    _storeClearAnnoBackLinks(fragStore.annotationStore);
    _storeCreateAnnoBackLinks(fragStore.annotationStore);
    _storeRemoveTempAnnoNames(fragStore);
}

//////////////////////////////////////////////////////////////////////////////
// extract FragmentStore annotation into a UCSC record

template <typename TRecord, typename TSpec, typename TConfig, typename TAnnotation, typename TId>
inline bool
_fillAnnotationRecord(
    TRecord & record,
    FragmentStore<TSpec, TConfig> & fragStore,
    TAnnotation & annotation,
    TId id,
    UcscKnownGene)
{
    typedef FragmentStore<TSpec, TConfig> TFragmentStore;

    if (annotation.typeId != TFragmentStore::ANNO_MRNA)
        return false;

    record.transName = getAnnoUniqueName(fragStore, id);
    if (annotation.contigId < length(fragStore.contigNameStore))
        record.contigName = fragStore.contigNameStore[annotation.contigId];
    else
        clear(record.contigName);

    clear(record.proteinName);
    clear(record.exonBegin);
    clear(record.exonEnds);

    record.annotationBeginPos = annotation.beginPos;
    record.annotationEndPos = annotation.endPos;

    TId lastChildId = annotation.lastChildId;
    TId i = lastChildId;
    do
    {
        i = fragStore.annotationStore[i].nextSiblingId;
        TAnnotation & anno = fragStore.annotationStore[i];
        if (anno.typeId == TFragmentStore::ANNO_CDS)
        {
            if (i < length(fragStore.annotationNameStore))
                record.proteinName = fragStore.annotationNameStore[i];
            record.cdsBegin = anno.beginPos;
            record.cdsEnd = anno.endPos;
        }
        if (anno.typeId == TFragmentStore::ANNO_EXON)
        {
            appendValue(record.exonBegin, anno.beginPos, Generous());
            appendValue(record.exonEnds, anno.endPos, Generous());
        }
    }
    while (i != lastChildId);
    return true;
}

template <typename TRecord, typename TSpec, typename TConfig, typename TAnnotation, typename TId>
inline bool
_fillAnnotationRecord(
    TRecord & record,
    FragmentStore<TSpec, TConfig> & fragStore,
    TAnnotation & annotation,
    TId id,
    UcscKnownIsoforms)
{
    typedef FragmentStore<TSpec, TConfig> TFragmentStore;

    if (annotation.typeId != TFragmentStore::ANNO_MRNA)
        return false;

    if (annotation.parentId == TAnnotation::INVALID_ID || annotation.parentId == 0)
        return false;

    record.transName = getAnnoUniqueName(fragStore, annotation.parentId);
    record.contigName = getAnnoUniqueName(fragStore, id);
    return true;
}

//////////////////////////////////////////////////////////////////////////////
// write FragmentStore to a stream in UCSC format

template <typename TSpec, typename TFSSpec, typename TFSConfig, typename TFormat>
inline void
writeRecords(FormattedFile<Ucsc, Output, TSpec> & ucscFile,
             FragmentStore<TFSSpec, TFSConfig> & store,
             TFormat const & format)
{
    typedef FragmentStore<TFSSpec, TFSConfig>                       TFragmentStore;
    typedef typename TFragmentStore::TAnnotationStore               TAnnotationStore;
    typedef typename Value<TAnnotationStore>::Type                  TAnnotation;
    typedef typename Iterator<TAnnotationStore, Standard>::Type     TAnnoIter;
    typedef typename Id<TAnnotation>::Type                          TId;

    typename DirectionIterator<FormattedFile<Ucsc, Output, TSpec>, Output>::Type iter = directionIterator(ucscFile, Output());
    UcscRecord record;

    TAnnoIter it = begin(store.annotationStore, Standard());
    TAnnoIter itEnd = end(store.annotationStore, Standard());

    for (TId id = 0; it != itEnd; ++it, ++id)
        if (_fillAnnotationRecord(record, store, *it, id, format))
            writeRecord(iter, record, format);
}

template <typename TSpec, typename TFSSpec, typename TFSConfig>
inline void
writeRecords(FormattedFile<Ucsc, Output, TSpec> & ucscFile,
             FragmentStore<TFSSpec, TFSConfig> & store)
{
    writeRecords(ucscFile, store, format(ucscFile));
}

} // namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_INCLUDE_STORE_STORE_IO_UCSC_H
