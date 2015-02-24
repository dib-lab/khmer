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

#ifndef SEQAN_HEADER_STORE_IO_GFF_H
#define SEQAN_HEADER_STORE_IO_GFF_H

namespace SEQAN_NAMESPACE_MAIN {

//////////////////////////////////////////////////////////////////////////////
// Read Gff
//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec = void>
struct IOContextGff_
{
    typedef typename TFragmentStore::TAnnotationStore   TAnnotationStore;
    typedef typename Value<TAnnotationStore>::Type      TAnnotation;
    typedef typename TAnnotation::TId                   TId;

    CharString contigName;
    CharString typeName;
    CharString annotationName;
    CharString parentKey;
    CharString parentName;

    CharString _key;
    CharString _value;
    StringSet<CharString> keys;
    StringSet<CharString> values;

    CharString gtfGeneId;
    CharString gtfGeneName;
    CharString gtfTranscriptName;       // transcipt_id is stored in parentName

    TId annotationId;
    TAnnotation annotation;
};

template <typename TFragmentStore, typename TSpec>
inline void clear(IOContextGff_<TFragmentStore, TSpec> & ctx)
{
    typedef typename TFragmentStore::TAnnotationStore   TAnnotationStore;
    typedef typename Value<TAnnotationStore>::Type      TAnnotation;

    clear(ctx.contigName);
    clear(ctx.typeName);
    clear(ctx.annotationName);
    clear(ctx.parentKey);
    clear(ctx.parentName);
    clear(ctx._key);
    clear(ctx._value);
    clear(ctx.gtfGeneId);
    clear(ctx.gtfGeneName);
    clear(ctx.gtfTranscriptName);
    clear(ctx.keys);
    clear(ctx.values);
    ctx.annotationId = TAnnotation::INVALID_ID;
    clear(ctx.annotation.values);
}

//////////////////////////////////////////////////////////////////////////////
// _readOneAnnotation
//
// reads in one annotation line from a Gff file

template <typename TFragmentStore, typename TSpec>
inline void
_readOneAnnotation(
    IOContextGff_<TFragmentStore, TSpec> & ctx,
    GffRecord const & record)
{
//IOREV _nodoc_ _hasCRef_
    typedef typename TFragmentStore::TContigPos         TContigPos;

    clear(ctx);

    // read column 1: contig name
    ctx.contigName = record.ref;

    // skip column 2
    // read column 3: type
    ctx.typeName = record.type;

    // read column 4 and 5: begin and endposition
    ctx.annotation.beginPos = record.beginPos;
    ctx.annotation.endPos = record.endPos;

    // skip column 6
    // read column 7: orientation
    if (record.strand == '-')
    {
        TContigPos tmp = ctx.annotation.beginPos;
        ctx.annotation.beginPos = ctx.annotation.endPos;
        ctx.annotation.endPos = tmp;
    }

    // skip column 8
    // read column 9: name
    for (unsigned i = 0; i < length(record.tagNames); ++i)
    {
        ctx._key = record.tagNames[i];
        ctx._value = record.tagValues[i];
        if (ctx._key == "ID")
        {
            ctx.annotationName = ctx._value;
        }
        else if (!empty(ctx._key) && !empty(ctx._value))
        {
            appendValue(ctx.keys, ctx._key);
            appendValue(ctx.values, ctx._value);
        }

        if (ctx._key == "Parent" || ctx._key == "ParentID" || ctx._key == "transcript_id")
        {
            ctx.parentKey = ctx._key;
            ctx.parentName = ctx._value;
        }
        else if (ctx._key == "transcript_name")
        {
            ctx.gtfTranscriptName = ctx._value;
        }
        else if (ctx._key == "gene_id")
        {
            ctx.gtfGeneId = ctx._value;
        }
        else if (ctx._key == "gene_name")
        {
            ctx.gtfGeneName = ctx._value;
        }

        clear(ctx._key);
        clear(ctx._value);
    }
}

template <typename TAnnotation>
inline void
_adjustParent(
    TAnnotation & parent,
    TAnnotation const & child)
{
    if (child.contigId == TAnnotation::INVALID_ID || child.beginPos == TAnnotation::INVALID_POS || child.endPos == TAnnotation::INVALID_POS)
        return;

    parent.contigId = child.contigId;

    // Has parent an invalid begin and end position?
    if ((parent.beginPos == TAnnotation::INVALID_POS) && (parent.endPos == TAnnotation::INVALID_POS))
    {
        parent.beginPos = child.beginPos;
        parent.endPos = child.endPos;
        return;
    }

    if ((parent.beginPos == TAnnotation::INVALID_POS) || (parent.endPos == TAnnotation::INVALID_POS))
        return;

    typename TAnnotation::TPos childBegin, childEnd;
    if (child.beginPos < child.endPos)
    {
        childBegin = child.beginPos;
        childEnd = child.endPos;
    }
    else
    {
        childBegin = child.endPos;
        childEnd = child.beginPos;
    }

    // Keep parent's orientation and maximize begin and end using child's boundaries.
    if (parent.beginPos < parent.endPos)
    {
        if (parent.beginPos == TAnnotation::INVALID_POS || parent.beginPos > childBegin)
            parent.beginPos = childBegin;
        if (parent.endPos == TAnnotation::INVALID_POS || parent.endPos < childEnd)
            parent.endPos = childEnd;
    }
    else
    {
        if (parent.endPos == TAnnotation::INVALID_POS || parent.endPos > childBegin)
            parent.endPos = childBegin;
        if (parent.beginPos == TAnnotation::INVALID_POS || parent.beginPos < childEnd)
            parent.beginPos = childEnd;
    }
}

template <typename TFragmentStore, typename TSpec>
inline void
_storeOneAnnotation(
    TFragmentStore & fragStore,
    IOContextGff_<TFragmentStore, TSpec> & ctx)
{
    typedef typename TFragmentStore::TAnnotationStore   TAnnotationStore;
    typedef typename Value<TAnnotationStore>::Type      TAnnotation;
    typedef typename TAnnotation::TId                   TId;

    SEQAN_ASSERT_EQ(length(fragStore.annotationStore), length(fragStore.annotationNameStore));

    // for lines in Gtf format get/add the parent gene first
    TId geneId = TAnnotation::INVALID_ID;
    if (!empty(ctx.gtfGeneId))
        _storeAppendAnnotationName(fragStore, geneId, ctx.gtfGeneId, (TId) TFragmentStore::ANNO_GENE);

    // if we have a parent transcript, get/add the parent transcript then
    if (!empty(ctx.parentName))
    {
// From now, we support gtf files with genes/transcripts having the same name.
//
//        // if gene and transcript names are equal (like in some strange gtf files)
//        // try to make the transcript name unique
//        if (ctx.gtfGeneId == ctx.parentName)
//            append(ctx.parentName, "_1");

        if (ctx.parentKey == "transcript_id")
            // type is implicitly given (mRNA)
            _storeAppendAnnotationName(fragStore, ctx.annotation.parentId, ctx.parentName, (TId) TFragmentStore::ANNO_MRNA);
        else
            // type is unknown
            _storeAppendAnnotationName(fragStore, ctx.annotation.parentId, ctx.parentName);
    }
    else
        ctx.annotation.parentId = 0;    // if we have no parent, we are a child of the root

    // add contig and type name
    _storeAppendContig(fragStore, ctx.annotation.contigId, ctx.contigName);
    _storeAppendType(fragStore, ctx.annotation.typeId, ctx.typeName);

    // add annotation name of the current line
    _storeAppendAnnotationName(fragStore, ctx.annotationId, ctx.annotationName, ctx.annotation.typeId);

    for (unsigned i = 0; i < length(ctx.keys); ++i)
    {
        // don't store gene_name as key/value pair unless it is a gene
        if (ctx.keys[i] == "gene_name" && ctx.annotation.typeId != TFragmentStore::ANNO_GENE)
            continue;

        // don't store transcript_name as key/value pair unless it is a transcript
        if (ctx.keys[i] == "transcript_name" && ctx.annotation.typeId != TFragmentStore::ANNO_MRNA)
            continue;

        // don't store Parent, transcript_id or gene_id as key/value pair (the are used to link annotations)
        if (ctx.keys[i] != ctx.parentKey && ctx.keys[i] != "gene_id")
            annotationAssignValueByKey(fragStore, ctx.annotation, ctx.keys[i], ctx.values[i]);
    }

    fragStore.annotationStore[ctx.annotationId] = ctx.annotation;

    TAnnotation & parent = fragStore.annotationStore[ctx.annotation.parentId];
    if (ctx.annotation.parentId != 0 && parent.parentId == TAnnotation::INVALID_ID)
        parent.parentId = 0;    // if our parent has no parent, it becomes a child of the root

    if (geneId != TAnnotation::INVALID_ID)
    {
        // link and adjust our gtf ancestors
        TAnnotation & gene = fragStore.annotationStore[geneId];
//        TAnnotation & transcript = fragStore.annotationStore[ctx.annotation.parentId];

        gene.parentId = 0;
        gene.typeId = TFragmentStore::ANNO_GENE;
        _adjustParent(gene, ctx.annotation);

        if (!empty(ctx.gtfGeneName))
            annotationAssignValueByKey(fragStore, gene, "gene_name", ctx.gtfGeneName);

        parent.parentId = geneId;
        parent.typeId = TFragmentStore::ANNO_MRNA;
        _adjustParent(parent, ctx.annotation);
        if (!empty(ctx.gtfTranscriptName))
            annotationAssignValueByKey(fragStore, parent, "transcript_name", ctx.gtfTranscriptName);
    }
}

template <typename TSpec, typename TConfig>
inline void
readRecords(FragmentStore<TSpec, TConfig> & fragStore,
            GffFileIn & gffFile)
{
    typedef FragmentStore<TSpec, TConfig> TFragmentStore;

    if (atEnd(gffFile))
        return;

    refresh(fragStore.contigNameStoreCache);
    refresh(fragStore.annotationNameStoreCache);
    refresh(fragStore.annotationTypeStoreCache);

    GffRecord record;
    IOContextGff_<TFragmentStore> ctx;

    while (!atEnd(gffFile))
    {
        readRecord(record, gffFile);
        _readOneAnnotation(ctx, record);
        _storeOneAnnotation(fragStore, ctx);
    }
    _storeClearAnnoBackLinks(fragStore.annotationStore);
    _storeCreateAnnoBackLinks(fragStore.annotationStore);
    _storeRemoveTempAnnoNames(fragStore);
}

//////////////////////////////////////////////////////////////////////////////
// Write Gff
//////////////////////////////////////////////////////////////////////////////

// This function write the information that are equal for gff and gtf files.
template <typename TSpec, typename TConfig, typename TAnnotation, typename TId>
inline void
_writeCommonGffGtfInfo(
    GffRecord & record,
    FragmentStore<TSpec, TConfig> & store,
    TAnnotation & annotation,
    TId /*id*/)
{
    typedef FragmentStore<TSpec, TConfig>       TFragmentStore;
    typedef typename TFragmentStore::TContigPos TContigPos;

    clear(record);

    // write column 1: contig name
    if (annotation.contigId < length(store.contigNameStore))
        if (length(store.contigNameStore[annotation.contigId]) > 0u)
            record.ref = store.contigNameStore[annotation.contigId];

    // skip column 2: source
    record.source = ".";

    // write column 3: type
    if (annotation.typeId < length(store.annotationTypeStore))
        if (length(store.annotationTypeStore[annotation.typeId]) > 0u)
            record.type = store.annotationTypeStore[annotation.typeId];

    TContigPos beginPos = annotation.beginPos;
    TContigPos endPos = annotation.endPos;
    char orientation = '+';
    if (endPos < beginPos)
    {
        TContigPos tmp = beginPos;
        beginPos = endPos;
        endPos = tmp;
        orientation = '-';
    }

    // write column 4: begin position
    if (beginPos != TAnnotation::INVALID_POS)
        record.beginPos = beginPos;

    // write column 5: end position
    if (endPos != TAnnotation::INVALID_POS)
        record.endPos = endPos;

    // skip column 6: score

    // write column 7: orientation
    record.strand = orientation;
}

template <typename TRecord, typename TSpec, typename TConfig, typename TAnnotation, typename TId>
inline bool
_fillAnnotationRecord(
    TRecord & record,
    FragmentStore<TSpec, TConfig> & store,
    TAnnotation & annotation,
    TId id,
    Gff)
{
    if (id == 0)
        return false;

    _writeCommonGffGtfInfo(record, store, annotation, id);

    // write column 9: group
    // write column 9.1: annotation id
    if (id < length(store.annotationNameStore) && !empty(getAnnoName(store, id)))
    {
        appendValue(record.tagNames, "ID");
        appendValue(record.tagValues, getAnnoName(store, id));
    }
    else if (annotation.lastChildId != TAnnotation::INVALID_ID)
    {
        appendValue(record.tagNames, "ID");
        appendValue(record.tagValues, getAnnoUniqueName(store, id));
    }

    // write column 9.2: parent id
    if (store.annotationStore[annotation.parentId].typeId > 1)  // ignore root/deleted nodes
    {
        appendValue(record.tagNames, "Parent");
        appendValue(record.tagValues, getAnnoUniqueName(store, annotation.parentId));
    }

    // write column 9.3-...: key, value pairs
    for (unsigned keyId = 0; keyId < length(annotation.values); ++keyId)
        if (!empty(annotation.values[keyId]))
        {
            appendValue(record.tagNames, store.annotationKeyStore[keyId]);
            appendValue(record.tagValues, annotation.values[keyId]);
        }

    return true;
}

template <typename TRecord, typename TSpec, typename TConfig, typename TAnnotation, typename TId>
inline bool
_fillAnnotationRecord(
    TRecord & record,
    FragmentStore<TSpec, TConfig> & store,
    TAnnotation & annotation,
    TId id,
    Gtf)
{
    typedef FragmentStore<TSpec, TConfig> TFragmentStore;

    if (annotation.typeId <= TFragmentStore::ANNO_MRNA)
        return false;

    _writeCommonGffGtfInfo(record, store, annotation, id);

    // write column 9: group

    // step up until we reach a transcript
    TId transcriptId = annotation.parentId;
    while (transcriptId < length(store.annotationStore) && store.annotationStore[transcriptId].typeId != TFragmentStore::ANNO_MRNA)
        transcriptId = store.annotationStore[transcriptId].parentId;

    // step up until we reach a gene
    TId geneId = transcriptId;
    while (geneId < length(store.annotationStore) && store.annotationStore[geneId].typeId != TFragmentStore::ANNO_GENE)
        geneId = store.annotationStore[geneId].parentId;

    typename Id<TAnnotation>::Type valueId;
    if (geneId < length(store.annotationStore) &&
        (valueId = annotationGetValueIdByKey(store, store.annotationStore[geneId], "gene_name")) != TAnnotation::INVALID_ID)
    {
        appendValue(record.tagNames, "gene_name");
        appendValue(record.tagValues, store.annotationStore[geneId].values[valueId]);
    }
    if (transcriptId < length(store.annotationStore) &&
        (valueId = annotationGetValueIdByKey(store, store.annotationStore[transcriptId], "transcript_name")) != TAnnotation::INVALID_ID)
    {
        appendValue(record.tagNames, "transcript_name");
        appendValue(record.tagValues, store.annotationStore[transcriptId].values[valueId]);
    }

    if (id < length(store.annotationNameStore) && !empty(getAnnoName(store, id)))
    {
        appendValue(record.tagNames, "ID");
        appendValue(record.tagValues, getAnnoName(store, id));
    }

    // write key, value pairs
    for (unsigned keyId = 0; keyId < length(annotation.values); ++keyId)
        if (!empty(annotation.values[keyId]))
        {
            appendValue(record.tagNames, store.annotationKeyStore[keyId]);
            appendValue(record.tagValues, annotation.values[keyId]);
        }

    // The GTF format version 2.2 requires the keys gene_id and transcript_id to be the last keys of line
    // read http://mblab.wustl.edu/GTF22.html and http://www.bioperl.org/wiki/GTF

    if (geneId < length(store.annotationStore))
    {
        appendValue(record.tagNames, "gene_id");
        appendValue(record.tagValues, getAnnoUniqueName(store, geneId));
    }

    if (transcriptId < length(store.annotationStore))
    {
        appendValue(record.tagNames, "transcript_id");
        appendValue(record.tagValues, getAnnoUniqueName(store, transcriptId));
    }
    return true;
}

// support for dynamically chosen file formats
template <typename TRecord, typename TSpec, typename TConfig, typename TAnnotation, typename TId>
inline bool
_fillAnnotationRecord(
    TRecord & /*record*/,
    FragmentStore<TSpec, TConfig> & /*store*/,
    TAnnotation & /*annotation*/,
    TId /*id*/,
    TagSelector<> const & /*format*/)
{
    SEQAN_FAIL("AnnotationStore: File format not specified.");
    return false;
}

template <typename TRecord, typename TSpec, typename TConfig, typename TAnnotation, typename TId, typename TTagList>
inline bool
_fillAnnotationRecord(
    TRecord & record,
    FragmentStore<TSpec, TConfig> & store,
    TAnnotation & annotation,
    TId id,
    TagSelector<TTagList> const & format)
{
    typedef typename TTagList::Type TFormat;

    if (isEqual(format, TFormat()))
        return _fillAnnotationRecord(record, store, annotation, id, TFormat());
    else
        return _fillAnnotationRecord(record, store, annotation, id, static_cast<typename TagSelector<TTagList>::Base const &>(format));
}


template <typename TTargetStream, typename TSpec, typename TConfig, typename TFormat>
inline void
_writeGffGtf(
    TTargetStream & target,
    FragmentStore<TSpec, TConfig> & store,
    TFormat const &format)
{
    typedef FragmentStore<TSpec, TConfig>                           TFragmentStore;
    typedef typename TFragmentStore::TAnnotationStore               TAnnotationStore;
    typedef typename Value<TAnnotationStore>::Type                  TAnnotation;
    typedef typename Iterator<TAnnotationStore, Standard>::Type     TAnnoIter;
    typedef typename Id<TAnnotation>::Type                          TId;

    TAnnoIter it = begin(store.annotationStore, Standard());
    TAnnoIter itEnd = end(store.annotationStore, Standard());

    GffRecord record;
    typename DirectionIterator<TTargetStream, Output>::Type iter = directionIterator(target, Output());

    for (TId id = 0; it != itEnd; ++it, ++id)
    {
        if (_fillAnnotationRecord(record, store, *it, id, format))
            writeRecord(iter, record, format);
    }
}

template <typename TSpec, typename TFSSpec, typename TFSConfig>
inline void
writeRecords(FormattedFile<Gff, Output, TSpec> & gffFile,
             FragmentStore<TFSSpec, TFSConfig> & store)
{
    _writeGffGtf(gffFile, store, format(gffFile));
}

} // namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
