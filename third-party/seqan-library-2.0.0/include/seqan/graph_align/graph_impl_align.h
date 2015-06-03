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
// Author: Tobias Rausch <rausch@embl.de>
// Author: Anne-Katrin Emde <anne-katrin.emde@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_GRAPH_ALIGN_GRAPH_IMPL_ALIGN_H_
#define SEQAN_INCLUDE_SEQAN_GRAPH_ALIGN_GRAPH_IMPL_ALIGN_H_

namespace seqan {

//////////////////////////////////////////////////////////////////////////////
// Alignment Graph
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////






//////////////////////////////////////////////////////////////////////////////
// Alignment Graph Output Tags
//////////////////////////////////////////////////////////////////////////////

/*!
 * @defgroup AlignmentGraphFormatTags Alignment Graph Formats
 * @brief File formats to write an @link AlignmentGraph @endlink.
 *
 * @tag AlignmentGraphFormatTags#MsfFormat
 * @headerfile <seqan/graph_align.h>
 * @brief MSF Format to write an @link AlignmentGraph @endlink.
 *
 * @signature typedef Tag<MsfFormat_> const MsfFormat;
 *
 * @tag AlignmentGraphFormatTags#FastaFormat
 * @headerfile <seqan/graph_align.h>
 * @brief Fasta Format to write an @link AlignmentGraph @endlink.
 *
 * @signature typedef Tag<FastaFormat_> const FastaFormat;
 *
 * @tag AlignmentGraphFormatTags#CgVizFormat
 * @headerfile <seqan/graph_align.h>
 * @brief CgViz Format to write an @link AlignmentGraph @endlink.
 *
 * @signature typedef Tag<CgVizFormat_> const CgVizFormat;
 */

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

struct MsfFormat_;
typedef Tag<MsfFormat_> const MsfFormat;

//////////////////////////////////////////////////////////////////////////////

struct FastaFormat_;
typedef Tag<FastaFormat_> const FastaFormat;

//////////////////////////////////////////////////////////////////////////////

struct CgVizFormat_;
typedef Tag<CgVizFormat_> const CgVizFormat;



//////////////////////////////////////////////////////////////////////////////
// Default Alignment Graph
//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo = unsigned int, typename TSpec = Default>
struct Alignment;





//////////////////////////////////////////////////////////////////////////////
// Metafunctions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec>
struct EdgeType<Graph<Alignment<TStringSet, TCargo, TSpec> > const> {
    typedef typename EdgeType<Graph<Undirected<TCargo, TSpec> > const>::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec>
struct EdgeType<Graph<Alignment<TStringSet, TCargo, TSpec> > > {
    typedef typename EdgeType<Graph<Undirected<TCargo, TSpec> > >::Type Type;
};


//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec>
struct Host<Graph<Alignment<TStringSet, TCargo, TSpec> > > {
    typedef TStringSet Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec>
struct Host<Graph<Alignment<TStringSet, TCargo, TSpec> > const> {
    typedef TStringSet const Type;
};






//////////////////////////////////////////////////////////////////////////////
// Actual Alignment Graph Class
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template<typename TId = unsigned int, typename TSize = unsigned int>
class FragmentInfo {
public:
    TId data_seq_id;
    TSize data_begin;
    TSize data_length;

    FragmentInfo() :
        data_seq_id(0),
        data_begin(0),
        data_length(0)
    {
    }

    FragmentInfo(TId id, TSize beg, TSize len) :
        data_seq_id(id),
        data_begin(beg),
        data_length(len)
    {
    }

};

//////////////////////////////////////////////////////////////////////////////

/*!
 * @class AlignmentGraph
 * @extends Graph
 * @headerfile <seqan/graph_align.h>
 * @brief Alignment graph type.
 *
 * <img src="alignmentGraph.png" tite="An alignment graph with 3 sequences." />
 *
 * @signature template <typename TStringSet[, typename TCargo[, typename TSpec]]>
 *            class Graph<Alignment<TStringSet, TCargo, TSpec> >;
 *
 * @tparam TStringSet The type of the @link StringSet @endlink containing the sequence information.  <b>Must be a
 *                    @link DependentStringSet @endlink.</b>
 * @tparam TCargo     The cargo type that can be attached to edges.  Default: <tt>unsigned</tt>.
 * @tparam TSpec      The specializing type.  Default: <tt>Default</tt>.  Use <tt>WithoutEdgeId</tt> here to omit
 *                    edge ids.  NB: if edges do not store ids then external property maps do not work.
 *
 * @section Example
 *
 * @include demos/graph_align/graph_align.cpp
 *
 * @code{.console}
 * Score = -8
 * Alignment matrix:
 *       0     .    :
 *         ATCGAATGCGGA
 *         |     |||  |
 *         ACTCGTTGC--A
 * @endcode
 */

template<typename TString, typename TSpecial, typename TCargo, typename TSpec>
class Graph<Alignment<StringSet<TString, Dependent<TSpecial> >, TCargo, TSpec> >
{
    public:
        typedef typename Id<Graph>::Type TIdType_;
        typedef typename VertexDescriptor<Graph>::Type TVertexDescriptor_;
        typedef typename Size<Graph>::Type TSize_;
        typedef std::pair<TIdType_, TSize_> TKey_;
        typedef std::map<TKey_, TVertexDescriptor_> TPosToVertexMap_;
        typedef FragmentInfo<TIdType_, TSize_> TFragmentInfo_;

        // Alignment graph
        Graph<Undirected<TCargo, TSpec> > data_align;

        // Sequences
        Holder<StringSet<TString, Dependent<TSpecial> > > data_sequence;

        // Alignment specific members
        String<TFragmentInfo_> data_fragment;

        // STL Map to retrieve a vertex given SeqId, Position
        TPosToVertexMap_ data_pvMap;


        Graph() {
        }


        template <typename TDefault>
        Graph(StringSet<TString, Dependent<TDefault> > const& sSet) {
            SEQAN_CHECKPOINT
            data_sequence = sSet;

            // Cover all sequences with nil vertices
            TVertexDescriptor_ nilVertex = getNil<TVertexDescriptor_>();
            TSize_ lenSet = length(sSet);
            for(TSize_ k=0; k<lenSet;++k)
                data_pvMap.insert(std::make_pair(TKey_(positionToId(sSet,k), length(sSet[k])), nilVertex));
        }

        template <typename TDefault>
        Graph(StringSet<TString, Owner<TDefault> > const& sSet) {
            SEQAN_CHECKPOINT
            StringSet<TString, Dependent<> > depStr(sSet);
            data_sequence = depStr;

            // Cover all sequences with nil vertices
            TVertexDescriptor_ nilVertex = getNil<TVertexDescriptor_>();
            TSize_ lenSet = length(sSet);
            for(TSize_ k=0; k<lenSet;++k)
                data_pvMap.insert(std::make_pair(TKey_(positionToId(const_cast<StringSet<TString, Owner<TDefault> >&>(sSet),k), length(sSet[k])), nilVertex));
        }


        ~Graph() {
            SEQAN_CHECKPOINT
            clear(*this);
        }

        Graph(Graph const & _other)
        {
            SEQAN_CHECKPOINT
            _copyGraph(_other, *this);
        }

        Graph const& operator = (Graph const & _other) {
            SEQAN_CHECKPOINT
            if (this == &_other) return *this;
            clear(*this);
            _copyGraph(_other, *this);
            return *this;
        }
};

//////////////////////////////////////////////////////////////////////////////
// INTERNAL FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec>
inline String<typename EdgeType<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type*>&
_getVertexString(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g) {
    SEQAN_CHECKPOINT
    typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;
    return const_cast<String<TEdgeStump*>&>(g.data_align.data_vertex);
}

/////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec>
inline typename VertexIdHandler<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type&
_getVertexIdManager(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g) {
    SEQAN_CHECKPOINT
    typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
    typedef typename VertexIdHandler<TGraph>::Type TVertexIdManager;
    return const_cast<TVertexIdManager&>(g.data_align.data_id_managerV);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec>
inline typename EdgeIdHandler<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type&
_getEdgeIdManager(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g) {
    SEQAN_CHECKPOINT
    typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
    typedef typename EdgeIdHandler<TGraph>::Type TEdgeIdManager;
    return const_cast<TEdgeIdManager&>(g.data_align.data_id_managerE);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec>
inline void
_copyGraph(Graph<Alignment<TStringSet, TCargo, TSpec> > const& source,
           Graph<Alignment<TStringSet, TCargo, TSpec> >& dest,
           bool)
{
    SEQAN_CHECKPOINT
    clear(dest);
    dest.data_align = source.data_align;
    dest.data_sequence = source.data_sequence;
    dest.data_fragment = source.data_fragment;
    dest.data_pvMap = source.data_pvMap;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec>
inline void
_copyGraph(Graph<Alignment<TStringSet, TCargo, TSpec> > const& source,
           Graph<Alignment<TStringSet, TCargo, TSpec> >& dest)
{
    _copyGraph(source, dest, false); // Never transpose, underlying graph is undirected
}


//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec>
inline void
transpose(Graph<Alignment<TStringSet, TCargo, TSpec> > const& source,
          Graph<Alignment<TStringSet, TCargo, TSpec> >& dest)
{
    SEQAN_CHECKPOINT
    // Alignment graph, no transpose just copy
    _copyGraph(source, dest, false);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec>
inline void
transpose(Graph<Alignment<TStringSet, TCargo, TSpec> >&)
{
    SEQAN_CHECKPOINT
    // Nothing to do in an alignment graph
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec>
inline typename Size<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type
numEdges(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g)
{
    SEQAN_CHECKPOINT
    return numEdges(g.data_align);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec>
inline typename Size<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type
numVertices(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g)
{
    SEQAN_CHECKPOINT
    return numVertices(g.data_align);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec>
inline bool
empty(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g)
{
    SEQAN_CHECKPOINT
    return empty(g.data_align);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec>
inline void
clearEdges(Graph<Alignment<TStringSet, TCargo, TSpec> >& g)
{
    SEQAN_CHECKPOINT
    clearEdges(g.data_align);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec>
inline void
clearVertices(Graph<Alignment<TStringSet, TCargo, TSpec> >& g)
{
    SEQAN_CHECKPOINT
    typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
    typedef typename Size<TStringSet>::Type TSize;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename TGraph::TKey_ TKey;

    clear(g.data_fragment);
    g.data_pvMap.clear();
    clearVertices(g.data_align);


    // Don't forget to cover the sequences with nil vertices again
    TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();
    if(!empty(value(g.data_sequence))) {
        TSize lenSet = length(stringSet(g));
        for(TSize k=0; k<lenSet;++k)
            g.data_pvMap.insert(std::make_pair(TKey(positionToId(stringSet(g),k), length(stringSet(g)[k])), nilVertex));
    }
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec>
inline void
clear(Graph<Alignment<TStringSet, TCargo, TSpec> >& g)
{
    SEQAN_CHECKPOINT
    // Only clear also removes the sequences
    clear(value(g.data_sequence));
    clearVertices(g);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename Size<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type
outDegree(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
          TVertexDescriptor const vertex)
{
    SEQAN_CHECKPOINT
    return outDegree(g.data_align, vertex);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename Size<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type
inDegree(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
         TVertexDescriptor const vertex)
{
    SEQAN_CHECKPOINT
    return inDegree(g.data_align, vertex);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename Size<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type
degree(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
       TVertexDescriptor const vertex)
{
    SEQAN_CHECKPOINT
    return degree(g.data_align, vertex);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TId, typename TPos, typename TLength>
inline typename VertexDescriptor<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type
addVertex(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
          TId id,
          TPos begin,
          TLength len)
{
    SEQAN_CHECKPOINT
    typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename Id<TGraph>::Type TIdType;
    typedef typename Size<TGraph>::Type TSize;
    typedef typename TGraph::TFragmentInfo_ TFragmentInfo;
    typedef typename TGraph::TKey_ TKey;
    typedef typename TGraph::TPosToVertexMap_ TPosToVertexMap;
    //typedef typename TPosToVertexMap::key_type TKey;

    //for(TPosToVertexMap::const_iterator p = g.data_pvMap.begin(); p != g.data_pvMap.end(); ++p) {
    //    std::cout << p->first.first << ',' << p->first.second << ':' << p->second << std::endl;
    //}

    TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();

    // Store the new fragment
    typename TPosToVertexMap::iterator interval = g.data_pvMap.lower_bound(TKey((TIdType)id, (TSize)begin + (TSize)len));
    // Segment does not belong to Sequence anymore
    SEQAN_ASSERT(interval != g.data_pvMap.end());
    // Segment end must be assigned to nil so far
    SEQAN_ASSERT(interval->second == nilVertex);
    // Segment must belong to the whole old interval
    SEQAN_ASSERT(*interval == *g.data_pvMap.upper_bound(TKey((TIdType)id, (TSize)begin)));

    // Insert new vertex
    TVertexDescriptor vd = addVertex(g.data_align);
    if (length(g.data_fragment) <= vd) resize(g.data_fragment, vd + 1, Generous());
    assignProperty(g.data_fragment, vd, TFragmentInfo(id, begin, len));

    // Update position to vertex map
    // Does the end of the new fragment coincides with the end of the interval?
    if ( (TSize) begin + len == (TSize) interval->first.second) {
        // Does the beginning of the new fragment coincides with the beginning of the interval?
        if ((begin == 0) ||
            (g.data_pvMap.find(TKey((TIdType)id, (TSize)begin)) != g.data_pvMap.end())) {
            // Replace interval
            interval->second = vd;
        } else {
            // Split interval once
            g.data_pvMap.insert(std::make_pair(TKey(interval->first.first,(TSize)begin), interval->second));
            g.data_pvMap.erase(interval);
            g.data_pvMap.insert(std::make_pair(TKey((TIdType)id,(TSize)begin+len), vd));
        }
    } else {
        // Does the beginning of the new fragment coincides with the beginning of the interval?
        if ((begin == 0) ||
            (g.data_pvMap.find(TKey((TIdType)id, (TSize)begin)) != g.data_pvMap.end())) {
            // Split interval once
            // Just insert here because we store interval ends
            g.data_pvMap.insert(std::make_pair(TKey((TIdType)id,(TSize)begin+len), vd));
        } else {
            // Split interval twice
            TSize tmp = interval->first.second;
            g.data_pvMap.insert(std::make_pair(TKey(interval->first.first,begin), interval->second));
            g.data_pvMap.erase(interval);
            g.data_pvMap.insert(std::make_pair(TKey((TIdType)id,(TSize)begin+len), vd));
            g.data_pvMap.insert(std::make_pair(TKey((TIdType)id,tmp), nilVertex));
        }
    }
    return vd;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TVD>
inline void
removeVertex(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
             TVD const v)
{
    typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename TGraph::TKey_ TKey;
    typedef typename TGraph::TPosToVertexMap_ TPosToVertexMap;
    typename TPosToVertexMap::iterator interval = (g.data_pvMap.lower_bound(TKey(sequenceId(g,v), fragmentBegin(g,v) + fragmentLength(g,v))));

    // Clear the interval
    interval->second = getNil<TVertexDescriptor>();
    typename TPosToVertexMap::iterator interval_iter = interval;
    if (interval_iter != g.data_pvMap.begin()) {
        --interval_iter;
        if ((interval_iter->second == getNil<TVertexDescriptor>()) &&
            (interval_iter->first.first == interval->first.first)) {
            g.data_pvMap.erase(interval_iter);
        }
    }
    interval_iter = interval;
    ++interval_iter;
    if (interval_iter != g.data_pvMap.end()) {
        if ((interval_iter->second == getNil<TVertexDescriptor>()) &&
            (interval_iter->first.first == interval->first.first)) {
            g.data_pvMap.erase(interval);
        }
    }

    // Remove the vertex
    removeVertex(g.data_align,v);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename EdgeDescriptor<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type
addEdge(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
        TVertexDescriptor const source,
        TVertexDescriptor const target)
{
    SEQAN_CHECKPOINT
    return addEdge(g.data_align, source, target);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TVertexDescriptor, typename TCargo2>
inline typename EdgeDescriptor<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type
addEdge(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
        TVertexDescriptor const source,
        TVertexDescriptor const target,
        TCargo2 const cargo)
{
    SEQAN_CHECKPOINT
    return addEdge(g.data_align, source, target, cargo);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void
removeEdge(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
           TVertexDescriptor const source,
           TVertexDescriptor const target)
{
    SEQAN_CHECKPOINT
    removeEdge(g.data_align, source, target);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TEdgeDescriptor>
inline void
removeEdge(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
           TEdgeDescriptor const edge)
{
    SEQAN_CHECKPOINT
    removeEdge(g.data_align, sourceVertex(g.data_align,edge), targetVertex(g.data_align,edge));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void
removeOutEdges(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
               TVertexDescriptor const v)
{
    SEQAN_CHECKPOINT
    removeOutEdges(g.data_align, v);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void
removeInEdges(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
              TVertexDescriptor const v)
{
    SEQAN_CHECKPOINT
    removeInEdges(g.data_align,v);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TEdgeDescriptor>
inline typename VertexDescriptor<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type
targetVertex(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
             TEdgeDescriptor const edge)
{
    SEQAN_CHECKPOINT
    return targetVertex(g.data_align, edge);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TEdgeDescriptor>
inline typename VertexDescriptor<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type
sourceVertex(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
             TEdgeDescriptor const edge)
{
    SEQAN_CHECKPOINT
    return sourceVertex(g.data_align, edge);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TMatrix>
inline void
getAdjacencyMatrix(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
                   TMatrix& mat)
{
    SEQAN_CHECKPOINT
    getAdjacencyMatrix(g.data_align, mat);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename EdgeDescriptor<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type
findEdge(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
         TVertexDescriptor const v,
         TVertexDescriptor const w)
{
    SEQAN_CHECKPOINT
    return findEdge(g.data_align, v, w);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TStringSet, typename TCargo, typename TSpec>
inline void
write(TFile & target,
      Graph<Alignment<TStringSet, TCargo, TSpec> > const& g)
{
//IOREV _nodoc_ specialization not documented (at least not obviously)
    typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
    typedef typename Size<TGraph>::Type TSize;
    typedef typename TGraph::TFragmentInfo_ TSegment;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;
    typedef typename Iterator<String<TEdgeStump*> const, Standard>::Type TIterConst;

    String<char> align;
    if (!convertAlignment(g, align)) {
        write(target,"Adjacency list:\n");
        TIterConst it = begin(g.data_align.data_vertex, Standard());
        TIterConst itEnd = end(g.data_align.data_vertex, Standard());
        TVertexDescriptor pos = 0;
        for(;it!=itEnd; ++it, ++pos) {
            if (!idInUse(g.data_align.data_id_managerV, pos)) continue;
            TVertexDescriptor sourceV = pos;
            appendNumber(target, (int)sourceV);
            TSegment seg = getProperty(g.data_fragment, sourceV);
            write(target," (SeqId:");
            appendNumber(target, (int)seg.data_seq_id);
            write(target," ,Begin:");
            appendNumber(target, (int)seg.data_begin);
            write(target," ,Length:");
            appendNumber(target, (int)seg.data_length);
            write(target,") -> ");
            TEdgeStump* current = *it;
            while(current!=0) {
                TVertexDescriptor adjV = getTarget(current);
                if (adjV != sourceV) {
                    appendNumber(target, (int)adjV);
                    writeValue(target, ',');
                    current=getNextS(current);
                } else {
                    adjV = getSource(current);
                    appendNumber(target, (int)adjV);
                    writeValue(target, ',');
                    current=getNextT(current);
                }
            }
            writeValue(target, '\n');
        }
        write(target,"Edge list:\n");
        it = begin(g.data_align.data_vertex, Standard());
        pos = 0;
        for(;it!=itEnd; ++it, ++pos) {
            TVertexDescriptor sourceV = pos;
            TEdgeStump* current = *it;
            while(current!=0) {
                TVertexDescriptor targetV = getTarget(current);
                if (sourceV != targetV) {
                    write(target,"Source: ");
                    appendNumber(target, (int)sourceV);
                    write(target, ",Target: ");
                    appendNumber(target, (int)targetV);
                    write(target, " (Id: ");
                    appendNumber(target, (int)_getId(current));
                    write(target, ")\n");
                    current=getNextS(current);
                } else {
                    current=getNextT(current);
                }
            }
        }
    } else {
        TSize nseq = length(stringSet(g));
        TSize colLen = length(align) / nseq;

        TSize baseCount=0;
        TSize leftSpace=6;
        TSize xPos = 0;
        write(target,"Alignment matrix:\n");
        while (xPos < colLen) {
            TSize windowSize = 50;
            if ((xPos + windowSize)>colLen) windowSize = colLen - xPos;

            // Print header line
            TSize offset=0;
            // Larger numbers need to be further left
            if (baseCount != 0) offset = (unsigned int) floor(log((double)baseCount) / log((double)10));
            for(TSize j = 0;j<leftSpace-offset;++j) writeValue(target, ' ');
            appendNumber(target, (int)baseCount);
            baseCount+=windowSize;
            writeValue(target, ' ');
            for(TSize col = 1;col<=windowSize;++col) {
                if ((col % 10)==0) writeValue(target, ':');
                else if ((col % 5)==0) writeValue(target, '.');
                else writeValue(target, ' ');
            }
            write(target, " \n");

            // Print sequences
            for(TSize row=0;row<2*nseq-1;++row) {
                for(TSize col = 0;col<leftSpace+2;++col)
                    writeValue(target, ' ');
                if ((row % 2)==0) {
                    for(TSize col = xPos;col<xPos+windowSize;++col)
                        writeValue(target, align[(row/2)*colLen+col]);
                } else {
                    for(TSize col = xPos;col<xPos+windowSize;++col) {
                        if ((align[((row-1)/2)*colLen + col] != gapValue<char>()) &&
                            (align[((row+1)/2)*colLen + col] != gapValue<char>()) &&
                            (align[((row-1)/2)*colLen + col] == align[((row+1)/2)*colLen + col]))
                            writeValue(target, '|');
                        else
                            writeValue(target, ' ');
                    }
                }
                writeValue(target, '\n');
            }
            writeValue(target, '\n');
            xPos+=windowSize;
        }
        writeValue(target, '\n');
    }
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TSpec, typename TNames>
inline void
write(TFile & file,
      Graph<TSpec> const& g,
      TNames const& names,
      FastaFormat)
{
//IOREV _nodoc_ specialization not documented
    SEQAN_CHECKPOINT
    typedef Graph<TSpec> TGraph;
    typedef typename Size<TGraph>::Type TSize;

    typename DirectionIterator<TFile, Output>::Type target = directionIterator(file, Output());

    String<char> align;
    if (convertAlignment(g, align)) {
        TSize nseq = length(stringSet(g));
        TSize colLen = length(align) / nseq;
        typedef typename Iterator<String<char>, Standard>::Type TIter;
        TIter it = begin(align, Standard());
        for(TSize i = 0; i<nseq; ++i) {
            writeValue(target, '>');
            write(target,names[i]);
            writeValue(target, '\n');
            TSize col = 0;
            while(col < colLen) {
                TSize max = ((colLen - col) < 60) ? colLen - col : 60;
                for(TSize finger = 0; finger<max; ++finger, ++col, ++it)
                    writeValue(target, *it);
                writeValue(target, '\n');
            }
        }
    }
}



//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TSpec, typename TNames>
inline void
write(TFile & file,
      Graph<TSpec> const& g,
      TNames const& names,
      MsfFormat)
{
//IOREV _nodoc_ specialization not documented
    SEQAN_CHECKPOINT
    typedef Graph<TSpec> TGraph;
    typedef typename Size<TGraph>::Type TSize;

    typename DirectionIterator<TFile, Output>::Type target = directionIterator(file, Output());

    String<char> align;
    if (convertAlignment(g, align)) {
        TSize nseq = length(stringSet(g));
        TSize colLen = length(align) / nseq;

        write(target,"PileUp\n");
        writeValue(target, '\n');
        write(target," MSF: ");
        appendNumber(target, (unsigned)colLen);
        write(target," Type: P Check: 0 ..\n\n");
        TSize offset = 0;
        for(TSize i = 0; i<nseq; ++i) {
            write(target," Name: ");
            write(target,names[i]);
            write(target," oo  Len:  ");
            TSize len = length(names[i]);
            if (len > offset) offset = len;
            appendNumber(target, (unsigned)colLen);
            write(target," Check: 0");
            write(target," Weight: 1.00");
            writeValue(target, '\n');
        }
        offset += 5;
        writeValue(target, '\n');
        write(target, "//\n");
        writeValue(target, '\n');
        writeValue(target, '\n');
        TSize col = 0;
        while(col < colLen) {
            TSize max = 0;
            for(TSize i = 0; i<nseq; ++i) {
                max = ((colLen - col) < 50) ? colLen - col : 50;
                write(target,names[i]);
                for(TSize j = 0; j<offset - length(names[i]); ++j)
                    writeValue(target, ' ');
                for(TSize finger = col; finger<col+max; ++finger) {
                    if ((finger - col) % 10 == 0) writeValue(target, ' ');
                    if (align[i*colLen + finger] == '-') writeValue(target, '.');
                    else writeValue(target, align[i*colLen + finger]);
                }
                writeValue(target, '\n');
            }
            col += max;
            writeValue(target, '\n');
            writeValue(target, '\n');
        }
    }
}


//////////////////////////////////////////////////////////////////////////////
template <typename TFile, typename TStringSet, typename TSpec, typename TEdge>
inline void
_writeCargo(TFile & file,
             Graph<Alignment<TStringSet, void, TSpec> > const&,
             TEdge const&)
{
    appendNumber(file, 0);
}

//////////////////////////////////////////////////////////////////////////////
template <typename TFile, typename TStringSet, typename TCargo, typename TSpec, typename TEdge>
inline void
_writeCargo(TFile & file,
             Graph<Alignment<TStringSet, TCargo, TSpec> > const&,
             TEdge const& edge)
{
    appendNumber(file, (int)getCargo(edge));
}

//////////////////////////////////////////////////////////////////////////////


template <typename TFile, typename TStringSet, typename TCargo, typename TSpec, typename TNames>
inline void
write(TFile & file,
      Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
      TNames const& names,
      CgVizFormat)
{
//IOREV _nodoc_ specialization not documented
    typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
    typedef typename Size<TGraph>::Type TSize;
    typedef typename Id<TGraph>::Type TId;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;
    typedef typename Iterator<String<TEdgeStump*> const, Rooted>::Type TIterConst;

    typename DirectionIterator<TFile, Output>::Type target = directionIterator(file, Output());

    TStringSet& str = stringSet(g);
    TSize nseq = length(str);

    write(target, "{DATA Data\n");
    write(target, "\t[__GLOBAL__] tracks=");
    appendNumber(target, nseq);
    writeValue(target, '\n');
    for(TSize i = 0; i<nseq; ++i) {
        write(target, "\tfasta_id=\"");
        write(target,names[i]);
        write(target, "\" sequence=\"");
        write(target,str[i]);
        write(target, "\" track=");
        appendNumber(target, i);
        write(target, " type=\"DNA\": 0 ");
        appendNumber(target, length(str[i])- 1);
        writeValue(target, '\n');
    }
    write(target, "}\n");
    for(TSize i = 0; i<nseq; ++i) {
        write(target, "{DATA ");
        appendNumber(target, i);
        //write(target,names[i]);
        write(target, "-seqlen\n\t[__GLOBAL__]\n\tlength=");
        appendNumber(target, length(str[i]));
        write(target, ":\t0 ");
        appendNumber(target, length(str[i])- 1);
        write(target, "\n}\n");
    }
    for(TSize i=0; i<nseq; ++i) {
        for(TSize j=i+1; j<nseq; ++j) {
            write(target, "{DATA ");
            appendNumber(target, i);
            //write(target,names[i]);
            write(target, "-vs-");
            appendNumber(target, j);
            //write(target,names[j]);
            writeValue(target, '\n');
            write(target, "\t[__GLOBAL__]\n");
            for(TIterConst it = begin(g.data_align.data_vertex);!atEnd(it);goNext(it)) {
                TVertexDescriptor sourceV = position(it);
                TId id1 = sequenceId(g, sourceV);
                if ((positionToId(str, id1) != i) &&
                    (positionToId(str, id1) != j)) continue;
                TEdgeStump* current = getValue(it);
                while(current!=0) {
                    TVertexDescriptor targetV = getTarget(current);
                    TId id2 = sequenceId(g, targetV);
                    if (sourceV != targetV) {
                        if ((positionToId(str, id2) != i) &&
                            (positionToId(str, id2) != j)) {
                                current=getNextS(current);
                                continue;
                        }
                        write(target, "\t");
                        write(target, "source=");
                        appendNumber(target, (int)sourceV);
                        writeValue(target, ' ');
                        write(target, "target=");
                        appendNumber(target, (int)targetV);
                        writeValue(target, ' ');
                        write(target, "edgeId=");
                        appendNumber(target, (int)_getId(current));
                        writeValue(target, ' ');
                        write(target, "cargo=");
                        _writeCargo(target,g,current);
                        writeValue(target, ' ');
                        write(target, "label=");
                        write(target,label(g,sourceV));
                        writeValue(target, ' ');
                        write(target, "labelOpp=");
                        write(target,label(g,targetV));
                        writeValue(target, ':');
                        writeValue(target, '\t');
                        appendNumber(target, (int)fragmentBegin(g, sourceV));
                        writeValue(target, ' ');
                        appendNumber(target, (int)fragmentBegin(g, targetV));
                        writeValue(target, ' ');
                        appendNumber(target, (int)(fragmentBegin(g, sourceV) + fragmentLength(g, sourceV)));
                        writeValue(target, ' ');
                        appendNumber(target, (int)(fragmentBegin(g, targetV) + fragmentLength(g, targetV)));
                        writeValue(target, '\n');
                        current=getNextS(current);
                    } else {
                        current=getNextT(current);
                    }
                }
            }
            write(target, "}\n");
        }
    }
}


//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn AlignmentGraph#assignStringSet
 * @brief Assigns a new StringSet to an AlignmentGraph.
 *
 * @signature void assignStringSet(g, stringSet);
 *
 * @param[in,out] g         An AlignmentGraph.
 * @param[in]     stringSet The @link StringSet @endlink to assign to.
 *
 * @see AlignmentGraph#getStringSet
 * @see AlignmentGraph#stringSet
 */

template<typename TString, typename TDefault, typename TCargo, typename TSpec, typename TDefault2>
inline void
assignStringSet(Graph<Alignment<StringSet<TString, Dependent<TDefault> >, TCargo, TSpec> >& g,
                StringSet<TString, Dependent<TDefault2> > const& sStr)
{
    SEQAN_CHECKPOINT
    typedef Graph<Alignment<StringSet<TString, Dependent<TDefault> >, TCargo, TSpec> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename TGraph::TKey_ TKey;
    typedef typename Size<TGraph>::Type TSize;

    TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();
    clear(g);
    g.data_sequence = (StringSet<TString, Dependent<TDefault> >) sStr;
    TSize lenSet = length(sStr);
    for(TSize k=0; k<lenSet;++k)
        g.data_pvMap.insert(std::make_pair(TKey(positionToId(sStr,k), length(sStr[k])), nilVertex));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TDefault, typename TCargo, typename TSpec, typename TDefault2>
inline void
assignStringSet(Graph<Alignment<StringSet<TString, Dependent<TDefault> >, TCargo, TSpec> >& g,
                StringSet<TString, Owner<TDefault2> > const& sStr)
{
    SEQAN_CHECKPOINT
    StringSet<TString, Dependent<> > depStr(sStr);
    assignStringSet(g, depStr);
}


//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn AlignmentGraph#getStringSet
 * @brief Gets the StringSet of an AlignmentGraph.
 *
 * @signature THost getStringSet(g);
 *
 * @param[in,out] g An AlignmentGraph.
 *
 * @return THost A const reference to the StringSet's host.
 *
 * @see AlignmentGraph#assignStringSet
 * @see AlignmentGraph#stringSet
 */

template<typename TStringSet, typename TCargo, typename TSpec>
inline typename Host<Graph<Alignment<TStringSet, TCargo, TSpec> > const>::Type&
getStringSet(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g)
{
    SEQAN_CHECKPOINT
    return value(g.data_sequence);
}


//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn AlignmentGraph#stringSet
 * @brief Return reference to the StringSet of an AlignmentGraph.
 *
 * @signature THostRef stringSet(g);
 *
 * @param[in,out] g An AlignmentGraph.
 *
 * @return THost A reference to the StringSet's host.
 *
 * @see AlignmentGraph#assignStringSet
 * @see AlignmentGraph#getStringSet
 */

template <typename T>
struct StringSetType;

template <typename TStringSet, typename TCargo, typename TSpec>
struct StringSetType<Graph<Alignment<TStringSet, TCargo, TSpec> > >
{
    typedef typename Host<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type & Type;
};


template<typename TStringSet, typename TCargo, typename TSpec>
inline typename StringSetType<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type
stringSet(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g)
{
    SEQAN_CHECKPOINT
    return const_cast<TStringSet&>(value(g.data_sequence));
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn AlignmentGraph#label
 * @brief Get the label associated with this vertex descriptor or the sequence that is associated with a fragment.
 *
 * @signature TInfix label(ag, v);
 *
 * @param[in] ag The AlignmentGraph to query.
 * @param[in] v  The vertex descriptor to query.
 *
 * @return TInfix An infix representing the sequence label.
 */

template<typename TStringSet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename Infix<typename Value<TStringSet>::Type>::Type
label(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
      TVertexDescriptor const v)
{
    SEQAN_CHECKPOINT
    typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
    typedef typename TGraph::TFragmentInfo_ TSegment;
    TSegment seg = getProperty(g.data_fragment, v);
    //std::cout << seg.data_seq_id << ",";
    //std::cout << seg.data_begin << ",";
    //std::cout << seg.data_length << ",";
    //std::cout << getValueById(value(g.data_sequence), seg.data_seq_id) << std::endl;
    //std::cout << infix(getValueById(value(g.data_sequence), seg.data_seq_id), seg.data_begin, seg.data_begin + seg.data_length) << std::endl;
    return infix(getValueById(value(g.data_sequence), seg.data_seq_id), seg.data_begin, seg.data_begin + seg.data_length);
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn AlignmentGraph#sequenceId
 * @brief Gets the sequence id that is associated with this vertex descriptor or with a sequence of a fragment.
 *
 * @signature TId sequenceId(ag, v);
 *
 * @param[in] ag The AlignmentGraph.
 * @param[in] v  Vertex descriptor of the vertex to retrieve sequence id for.
 *
 * @return TId The sequence id.
 */

template<typename TStringSet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename Id<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type&
sequenceId(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
           TVertexDescriptor const v)
{
    SEQAN_CHECKPOINT
    return const_cast<typename Id<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type&>(g.data_fragment[v].data_seq_id);
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn AlignmentGraph#fragmentBegin
 * @brief Gets the begin position for this vertex descriptor in the sequence.
 *
 * @signature TPos fragmentBegin(ag, v);
 *
 * @param[in] ag The AlignmentGraph to query.
 * @param[in] v  The vertex descriptor of the vertex.
 *
 * @return TPos The begin position.
 */

template<typename TStringSet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename Position<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type&
fragmentBegin(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
              TVertexDescriptor const v)
{
    SEQAN_CHECKPOINT
    return const_cast<typename Position<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type&>(g.data_fragment[v].data_begin);
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn AlignmentGraph#fragmentLength
 * @brief Gets the length of the fragment for this vertex descriptor in the sequence.
 *
 * @signature TSize fragmentLength(ag, v);
 *
 * @param[in] ag The AlignmentGraph to query.
 * @param[in] v  The vertex descriptor of the vertex.
 *
 * @return TPos The fragment size.
 */

template<typename TStringSet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename Size<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type&
fragmentLength(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
               TVertexDescriptor const v)
{
    SEQAN_CHECKPOINT
    return const_cast<typename Size<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type&>(g.data_fragment[v].data_length);
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn AlignmentGraph#findVertex
 * @brief Finds a vertx given a sequence id and a position.
 *
 * @signature TVertexDescriptor findVertex(ag, id, pos);
 *
 * @param[in] ag  The AlignmentGraph to search in.
 * @param[in] id  The sequence id to search for.
 * @param[in] pos The position to search for.
 *
 * @return TVertexDescriptor The vertex covering the given position on the specified sequence,
 * <tt>getNil&lt;TVertexDescriptor&gt;()</tt> if none could be found.
 */

template<typename TStringSet, typename TCargo, typename TSpec, typename TSeqId, typename TPos>
inline typename VertexDescriptor<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type
findVertex(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
           TSeqId id,
           TPos pos)
{
    SEQAN_CHECKPOINT
    typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
    typedef typename TGraph::TKey_ TKey;
    typedef typename Size<TGraph>::Type TSize;
    typedef typename Id<TGraph>::Type TIdType;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;

    return (pos >= (TPos) length(getValueById(stringSet(g),id))) ? getNil<TVertexDescriptor>() : g.data_pvMap.upper_bound(TKey((TIdType)id, (TSize)pos))->second;
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn AlignmentGraph#getProjectedPosition
 * @brief Projects a position of one sequence taking part in a pairwise match onto the other sequence.
 *
 * @signature getProjectedPosition(ag, seqId, pos, seqId2, pos2);
 *
 * @param[in]  ag     The @link AlignmentGraph @endlink to query.
 * @param[in]  seqId  The sequenceId.
 * @param[in]  pos    The position to project from.
 * @param[out] seqId2 The resulting id of the sequence that was projected onto.
 * @param[out] pos2   The reuslting position of the sequence was was projected onto.
 */

template<typename TStringSet, typename TCargo, typename TSpec, typename TSeqId, typename TPosition, typename TSeqId2, typename TPosition2>
inline void
getProjectedPosition(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
                     TSeqId const id1,
                     TPosition const pos1,
                     TSeqId2& id2,
                     TPosition2& pos2)
{
    SEQAN_ASSERT(length(stringSet(g)) == 2);

    typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;
    TStringSet& str = stringSet(g);
    TVertexDescriptor sV = findVertex(g, id1, pos1);

    // Case 1: No projection possible
    if (sV == getNil<TVertexDescriptor>()) {
        if ( (TSeqId) positionToId(str, 0) == id1) id2 = (TSeqId2) positionToId(str,1);
        else id2 = (TSeqId2) positionToId(str,0);
        pos2 = 0;
        return;
    }

    // Case 2: Projection is possible
    TEdgeStump* current = getValue(g.data_align.data_vertex, sV);
    if(current != (TEdgeStump*) 0) {
        TVertexDescriptor tV = target(current);
        if (tV == sV) tV = source(current);
        pos2 = (TPosition2) (fragmentBegin(g,tV) + (pos1 - fragmentBegin(g, sV)));
        id2 = (TSeqId2) sequenceId(g, tV);
        return;
    } else {
        // If no out-going edge, get the preceding or following vertex
        if (fragmentBegin(g, sV) == 0) {
            getProjectedPosition(g, id1, fragmentBegin(g,sV) + fragmentLength(g, sV), id2, pos2);
            return;
        } else {
            getProjectedPosition(g, id1, fragmentBegin(g,sV) - 1, id2, pos2);
            return;
        }
    }
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TValue, typename TCargo, typename TSpec, typename TSeqId, typename TPosition, typename TSeqId2, typename TPosition2>
inline void
getProjectedPosition(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
                     TValue seg_num,
                     TSeqId const id1,
                     TPosition const pos1,
                     TSeqId2& id2,
                     TPosition2& pos2)
{
    SEQAN_ASSERT(length(stringSet(g)) == 2);

    typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;
    TStringSet& str = stringSet(g);
    TVertexDescriptor sV = findVertex(g, id1, pos1);

    // Case 1: No projection possible
    if (sV == getNil<TVertexDescriptor>()) {
        if ( (TSeqId) positionToId(str, 0) == id1) id2 = (TSeqId2) positionToId(str,1);
        else id2 = (TSeqId2) positionToId(str,0);
        pos2 = 0;
        return;
    }

    // Case 2: Projection is possible
    TEdgeStump* current = getValue(g.data_align.data_vertex, sV);
    if(current != (TEdgeStump*) 0) {
        TVertexDescriptor tV = target(current);
        if (seg_num == 0) tV = source(current); // segnum 0 is defined as the sourceVertex, segnum 1 is targetVertex
        pos2 = (TPosition2) (fragmentBegin(g,tV) + (pos1 - fragmentBegin(g, sV)));
        id2 = (TSeqId2) sequenceId(g, tV);
        return;
    } else {
        // If no out-going edge, get the preceding or following vertex
        if (fragmentBegin(g, sV) == 0) {
            getProjectedPosition(g, seg_num, id1, fragmentBegin(g,sV) + fragmentLength(g, sV), id2, pos2);
            return;
        } else {
            getProjectedPosition(g, seg_num, id1, fragmentBegin(g,sV) - 1, id2, pos2);
            return;
        }
    }
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn AlignmentGraph#getFirstCoveredPosition
 * @brief Finds the first position in a sequence that is not assigned to a nil vertex.
 *
 * @signature TPos getFirstCoveredPosition(ag, id);
 *
 * @param[in] ag An @link AlignmentGraph @endlink.
 * @param[in] id A sequence id.
 *
 * @return TPos A sequence position.
 *
 * @see AlignmentGraph#getLastCoveredPosition
 */

template<typename TStringSet, typename TCargo, typename TSpec, typename TSeqId>
inline typename Position<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type
getFirstCoveredPosition(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
                        TSeqId const id)
{
    typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename TGraph::TPosToVertexMap_ TPosToVertexMap;


    typename TPosToVertexMap::const_iterator it = g.data_pvMap.upper_bound(std::make_pair(id, 0));

    // Case 1: id is not covered
    if (it == g.data_pvMap.end()) return length(getValueById(stringSet(g), id));

    // Case 2: We found a nil vertex, go one forward
    if (it->second == getNil<TVertexDescriptor>()) {
        ++it;
        if (it == g.data_pvMap.end()) return length(getValueById(stringSet(g), id));
    }

    // Now we have the right vertex, return the beginning if the sequence id still fits
    if (it->first.first != id) return length(getValueById(stringSet(g), id));
    else return fragmentBegin(g, it->second);
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn AlignmentGraph#getLastCoveredPosition
 * @brief Finds the last position in a sequence that is not assigned to a nil vertex.
 *
 * @signature TPos getLastCoveredPosition(ag, id);
 *
 * @param[in] ag An @link AlignmentGraph @endlink.
 * @param[in] id A sequence id.
 *
 * @return TPos A sequence position.
 *
 * @see AlignmentGraph#getFirstCoveredPosition
 */

template<typename TStringSet, typename TCargo, typename TSpec, typename TSeqId>
inline typename Position<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type
getLastCoveredPosition(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
                       TSeqId id)
{
    typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename TGraph::TPosToVertexMap_ TPosToVertexMap;

    TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();
    typename TPosToVertexMap::const_iterator it = g.data_pvMap.lower_bound(std::make_pair(id, length(getValueById(stringSet(g), id))));

    // Case 1: No last position (all nil)
    if ((it == g.data_pvMap.begin()) && (it->second == nilVertex)) return 0;

    // Case 2: Found a nil position but there is a vertex before
    if (it->second == nilVertex) {
        --it;
    }

    // If the sequence id still matches return the position behind the last position belonging to this vertex
    if (it->first.first != id) return 0;
    return fragmentBegin(g, it->second) + fragmentLength(g, it->second);
}


//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn AlignmentGraph#convertAlignment
 * @brief Converts an alignment graph into an alignment matrix.
 *
 * @signature bool convertAlignment(g, component, order, compLength);
 * @signature bool convertAlignment(g, matrix);
 *
 * @param[in]  g          AlignmentGraph to convert.  An alignment graph.
 * @param[out] component  Vertex to component mapping.
 * @param[out] order      The order of the component graph when sorting topologically.
 * @param[out] compLength Component sizes.
 * @param[out] matrix     A string that represents an alignment matrix.
 *
 * @return bool true iff the alignment graph is a valid alignment and false otherwise.
 *
 * @section Remarks
 *
 * The variant with <tt>component</tt> and <tt>order</tt> computes a topological sorting of connected components.
 */

template<typename TStringSet, typename TCargo, typename TSpec, typename TComponentMap, typename TOrderMap, typename TComponentLength>
inline bool
convertAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
                 TComponentMap& component,
                 TOrderMap& order,
                 TComponentLength& compLength)
{
    SEQAN_CHECKPOINT
    typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename Id<TGraph>::Type TIdType;
    typedef typename Size<TGraph>::Type TSize;
    typedef typename Value<TComponentMap>::Type TComponent;
    typedef typename TGraph::TPosToVertexMap_ TPosToVertexMap;
    typedef typename TComponentLength::mapped_type TMappedType;
    typedef typename TComponentLength::key_type TKeyType;
    TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();

    // Check for empty graph
    if (empty(g)) return false;

    // Connected Components
    TSize numComponents = connectedComponents(component, g);

    // Make a directed graph to represent the ordering of the components
    // Note: Multiple vertices might have the same component
    Graph<Directed<void, WithoutEdgeId> > componentGraph;
    reserve(_getVertexString(componentGraph), numComponents);
    for(TSize i = 0; i<numComponents;++i) addVertex(componentGraph);

    TSize nseq = length(value(g.data_sequence));
    String<std::set<TComponent> > componentsPerSeq;
    typedef String<String<TComponent> > TOrderedComponents;
    TOrderedComponents orderedComponentsPerSeq;
    resize(componentsPerSeq, nseq);
    resize(orderedComponentsPerSeq, nseq);
    typename TPosToVertexMap::const_iterator it1 = g.data_pvMap.begin();
    typename TPosToVertexMap::const_iterator it1End = g.data_pvMap.end();
    for(;it1!=it1End;++it1) {
        // If sections are not assigned to a vertex -> no alignment
        if (it1->second == nilVertex) return false;

        // Remember the sequence that component belongs to
        TSize currentSeq = idToPosition(value(g.data_sequence), it1->first.first);

        // Append component
        TComponent c = getProperty(component, it1->second);
        if ((value(componentsPerSeq,currentSeq)).empty()) {
            String<TComponent> tmp;
            value(orderedComponentsPerSeq, currentSeq) = tmp;
        }
        appendValue(value(orderedComponentsPerSeq, currentSeq), c, Generous());
        // If two components appear twice in the same sequence -> no alignment
        if (!((value(componentsPerSeq,currentSeq)).insert(c)).second) return false;
    }
    clear(componentsPerSeq);

    // Draw edges for the components within a sequence
    typedef typename Iterator<TOrderedComponents>::Type TIterTOrderedComponents;
    TIterTOrderedComponents itBegin = begin(orderedComponentsPerSeq);
    TIterTOrderedComponents itEnd = end(orderedComponentsPerSeq);
    for(;itBegin != itEnd; ++itBegin) {
        TSize n = length(*itBegin);
        for(TSize i = 0; i<n-1; ++i) {
            addEdge(componentGraph, value((*itBegin), i), value((*itBegin), i+1));
        }
    }

    // Make a topological sort of the component graph
    topologicalSort(order, componentGraph);

    //// Debug code
    //std::cout << "Topological sort: " << std::endl;
    //for(TSize i = 0; i<length(order);++i) {
    //    std::cout << order[i] << ',';
    //}
    //std::cout << std::endl;

    // Walk through all sequences and check the component order
    unsigned int compIndex = 0;
    unsigned int compIndexLen = length(order);
    typename TPosToVertexMap::const_iterator it = g.data_pvMap.begin();
    TIdType currentSeq = it->first.first;
    for(; it != g.data_pvMap.end(); ++it) {
        if (it->first.first != currentSeq) {
            compIndex = 0;
            currentSeq = it->first.first;
        }
        TKeyType c = (TKeyType) getProperty(component, it->second);
        if(!compLength.insert(std::make_pair(c, (TMappedType) fragmentLength(g, it->second))).second) {
            compLength[c] = _max(compLength[c], (TMappedType) fragmentLength(g, it->second));
        }
        while ((compIndex < compIndexLen) && (order[compIndex] != c)) ++compIndex;
        // Crossing components -> no alignment
        if (compIndex >= compIndexLen) return false;
        // Next component
        ++compIndex;
    }

    return true;
}


//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TMatrix>
inline bool
convertAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
                 TMatrix& mat)
{
    SEQAN_CHECKPOINT
    typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
    typedef typename Value<TMatrix>::Type TValue;
    typedef typename Size<TGraph>::Type TSize;
    typedef typename Id<TGraph>::Type TIdType;
    typedef typename TGraph::TPosToVertexMap_ TPosToVertexMap;
    typedef std::map<unsigned int, unsigned int> TComponentLength;

    // Strongly Connected Components, topological sort, and length of each component
    String<unsigned int> component;
    String<unsigned int> order;
    TComponentLength compLength;

    if (!convertAlignment(g, component, order, compLength)) return false;

    // Create the matrix
    TSize len = 0;
    TSize nseq = length(stringSet(g));
    for(TComponentLength::iterator cIt=compLength.begin(); cIt != compLength.end(); ++cIt) len+=cIt->second;
    char gapChar = gapValue<char>();
    resize(mat, len * nseq, gapChar);

    // Fill the matrix
    TSize row = 0;
    TSize col = 0;
    typename TPosToVertexMap::const_iterator it = g.data_pvMap.begin();
    unsigned int compIndex = 0;
    unsigned int compIndexLen = length(order);
    TIdType currentSeq = it->first.first;
    for(; it != g.data_pvMap.end(); ++it) {
        if (it->first.first != currentSeq) {
            SEQAN_ASSERT(col <= len);
            //std::cout << std::endl;
            ++row;col=0;
            compIndex = 0;
            currentSeq = it->first.first;
        }
        unsigned int c = getProperty(component, it->second);
        while ((compIndex < compIndexLen) && (order[compIndex] != c)) {
            for(TSize i=0;i<compLength[order[compIndex]];++i, ++col)
                mat[row*len + col] = gapValue<char>();
            ++compIndex;
        }
        typedef typename Value<TStringSet>::Type TStringSetStr;
        TStringSetStr str = label(g,it->second);
        typedef typename Iterator<TStringSetStr, Standard>::Type TStringSetStrIter;
        TStringSetStrIter itStr = begin(str, Standard());
        TStringSetStrIter itStrEnd = end(str, Standard());
        for(;itStr != itStrEnd; goNext(itStr), ++col)
            mat[row*len + col] = (TValue) (*itStr);
        for(TSize i = length(str); i < compLength[order[compIndex]]; ++i, ++col)
            mat[row*len + col] = gapValue<char>();
        ++compIndex;
    }
    SEQAN_ASSERT(row + 1 == nseq);
    //std::cout << std::endl;

    return true;
}



//////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TSpec1, typename TStringSet, typename TCargo, typename TSpec2>
inline bool
convertAlignment(String<TValue, TSpec1> const& mat,
                 Graph<Alignment<TStringSet, TCargo, TSpec2> >& g)
{
    typedef String<TValue, TSpec1> TMatrix;
    typedef Graph<Alignment<TStringSet, TCargo, TSpec2> > TGraph;
    typedef typename Size<TGraph>::Type TSize;
    typedef typename Iterator<TMatrix, Standard>::Type TMatIter;
    clearVertices(g);
    TValue gapChar = gapValue<TValue>();
    TSize nseq = length(stringSet(g));
    TSize alignLen = length(mat) / nseq;
    String<Fragment<> > matches;
    for(TSize seq1 = 0; seq1<nseq; ++seq1) {
        for(TSize seq2 = seq1 + 1; seq2<nseq; ++seq2) {
            TMatIter seq1It = begin(mat);
            seq1It += seq1 * alignLen;
            TMatIter seq2It = begin(mat);
            seq2It += seq2 * alignLen;
            TSize alignPos = 0;
            TSize length = 0;
            TSize offset1 = 0;
            TSize offset2 = 0;
            for(TSize col = 0; col<alignLen; ++col, ++seq1It, ++seq2It, ++alignPos) {
                if ((*seq1It == gapChar) || (*seq2It == gapChar)) {
                    if (length) {
                        appendValue(matches, Fragment<>(seq1, alignPos - offset1 - length, seq2, alignPos - offset2 - length, length));
                        length = 0;
                    }
                    if (*seq1It == gapChar) ++offset1;
                    if (*seq2It == gapChar) ++offset2;
                } else ++length;
            }
            if (length) appendValue(matches, Fragment<>(seq1, alignPos - offset1 - length, seq2, alignPos - offset2 - length, length));
        }
    }
    //_debugMatches(stringSet(g), matches);
    matchRefinement(matches,stringSet(g),g);
    return true;
}



//////////////////////////////////////////////////////////////////////////////


template<typename TStringSet, typename TCargo, typename TSpec>
inline void
rebuildGraph(Graph<Alignment<TStringSet, TCargo, TSpec> >& g)
{
    typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
    typedef typename Size<TGraph>::Type TSize;

    // Initialization
    typedef Fragment<> TFragment;
    typedef String<TFragment> TFragmentString;
    TFragmentString matches;
    TSize nseq = length(stringSet(g));

    // Collect all character pairs
    typedef std::pair<TSize, TSize> TResiduePair;
    typedef std::set<TResiduePair> TResiduePairSet;
    String<TResiduePairSet> resPair;
    resize(resPair, nseq * nseq);
    TEdgeIterator itE(g);
    for(;!atEnd(itE);++itE) {
        TVertexDescriptor sV = sourceVertex(itE);
        TVertexDescriptor tV = targetVertex(itE);
        TSize seq1 = idToPosition(stringSet(g), sequenceId(g, sV));
        TSize seq2 = idToPosition(stringSet(g), sequenceId(g, tV));
        TSize index = 0;
        TSize pos1 = 0;
        TSize pos2 = 0;
        if (seq1 < seq2) {
            index = seq1 * nseq + seq2;
            pos1 = fragmentBegin(g, sV);
            pos2 = fragmentBegin(g, tV);
        } else {
            index = seq2 * nseq + seq1;
            pos1 = fragmentBegin(g, tV);
            pos2 = fragmentBegin(g, sV);
        }
        for(TSize i = 0; i<fragmentLength(g, sV); ++i) {
            resPair[index].insert(std::make_pair(pos1 + i, pos2 + i));
        }
    }

    // Rebuild the graph with maximal segments
    for(TSize i = 0; i<length(resPair); ++i) {
        if (resPair[i].empty()) continue;
        TSize seq1 = i / nseq;
        TSize seq2 = i % nseq;
        typename TResiduePairSet::const_iterator pos = resPair[i].begin();
        typename TResiduePairSet::const_iterator posEnd = resPair[i].end();
        TSize startMatch1 = pos->first;
        TSize startMatch2 = pos->second;
        TSize len = 1;
        ++pos;
        while(pos != posEnd) {
            if ((startMatch1 + len == pos->first) && (startMatch2 + len == pos->second)) ++len;
            else {
                appendValue(matches, TFragment(seq1, startMatch1, seq2, startMatch2, len), Generous());
                startMatch1 = pos->first;
                startMatch2 = pos->second;
                len = 1;
            }
            ++pos;
        }
        appendValue(matches, TFragment(seq1, startMatch1, seq2, startMatch2, len), Generous());
    }
    clearVertices(g);
    matchRefinement(matches,stringSet(g),g);
}





//////////////////////////////////////////////////////////////////////////////
// Heaviest Common Subsequence adaptation
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TSize2, typename TSpec2, typename TPositions, typename TSize, typename TVertexDescriptor, typename TString>
inline void
_heaviestCommonSubsequence(Graph<Alignment<TStringSet, TCargo, TSpec> > const&,
                            String<TSize2, TSpec2> const& /*slotToPos*/,
                            TPositions const&,
                            TSize const,
                            TSize const,
                            TVertexDescriptor const,
                            TString const&,
                            TString const&,
                            Nothing&)
{
    SEQAN_CHECKPOINT
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TSize2, typename TSpec2, typename TPositions, typename TSize, typename TString, typename TOutString>
inline void
_heaviestCommonSubsequence(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
                            String<TSize2, TSpec2> const& slotToPos,
                            TPositions const& positions,
                            TSize const m,
                            TSize const n,
                            TString const& str1,
                            TString const& str2,
                            TOutString& align)
{
    SEQAN_CHECKPOINT
    typedef typename Value<TString>::Type TVertexSet;
    typedef typename Iterator<TString const, Standard>::Type TStringIter;
    typedef typename Iterator<TString, Standard>::Type TSIter;
    typedef typename Iterator<TVertexSet const, Standard>::Type TVertexSetIter;

    // Create the alignment sequence
    TSize numMatches = length(positions);
    TSize alignLength = numMatches + (n - numMatches) + (m - numMatches);
    clear(align);
    resize(align, alignLength, TVertexSet(), Exact() );
    TSIter pointerAlign = begin(align, Standard());
    TSIter pointerAlignEnd = end(align, Standard());
    TStringIter pointerStr1 = begin(str1, Standard());
    TSize posStr1 = 0;
    TSize posStr2 = 0;
    TStringIter pointerStr2 = begin(str2, Standard());
    int p = length(positions) - 1;
    while(pointerAlign != pointerAlignEnd) {
        TSize i = m;
        TSize j = n;
        if (p>=0) {
            i = (TSize) (slotToPos[positions[p]] / (TSize) n);   // Get the index in str1
            j = n - 1 - (TSize) (slotToPos[positions[p]] % (TSize) n); // Get the index in str2
        };

        // In what order do we insert gaps? -> Only important at the beginning and at the end, not between matches
        bool firstI = true;
        if ((i != posStr1) && (j != posStr2))
        {
            if ((posStr1 == 0) && (posStr2 == 0)) {
                TStringIter tmpPointerStr1 = pointerStr1;
                TStringIter tmpPointerStr2 = pointerStr2;
                TSize tmpPosStr1 = posStr1;
                TSize tmpPosStr2 = posStr2;
                TSize len1 = 0;
                TSize len2 = 0;
                for(;i != tmpPosStr1; ++tmpPosStr1, ++tmpPointerStr1) len1 += fragmentLength(g, value(*tmpPointerStr1, 0));
                for(;j != tmpPosStr2; ++tmpPosStr2, ++tmpPointerStr2) len2 += fragmentLength(g, value(*tmpPointerStr2, 0));
                if (len1 > len2) firstI = false;
            } else if ((i == m) && (i == n)) {
                TStringIter tmpPointerStr1 = pointerStr1;
                TStringIter tmpPointerStr2 = pointerStr2;
                TSize tmpPosStr1 = posStr1;
                TSize tmpPosStr2 = posStr2;
                TSize len1 = 0;
                TSize len2 = 0;
                for(;i != tmpPosStr1; ++tmpPosStr1, ++tmpPointerStr1) len1 += fragmentLength(g, value(*tmpPointerStr1, 0));
                for(;j != tmpPosStr2; ++tmpPosStr2, ++tmpPointerStr2) len2 += fragmentLength(g, value(*tmpPointerStr2, 0));
                if (len1 < len2) firstI = false;
            }
        }
        if (firstI) {
            // Gaps in seq 2
            while (i != posStr1) {
                TVertexSetIter itV = begin(*pointerStr1, Standard());
                TVertexSetIter itVEnd = end(*pointerStr1, Standard());
                for(;itV != itVEnd;++itV) appendValue(*pointerAlign, *itV, Generous());
                ++pointerAlign;
                ++pointerStr1; ++posStr1;
            }
            // Gaps in seq 1
            while (j != posStr2) {
                TVertexSetIter itV = begin(*pointerStr2, Standard());
                TVertexSetIter itVEnd = end(*pointerStr2, Standard());
                for(;itV != itVEnd;++itV) appendValue(*pointerAlign, *itV, Generous());
                ++pointerAlign;
                ++pointerStr2; ++posStr2;
            }
        } else {
            // Gaps in seq 1
            while (j != posStr2) {
                TVertexSetIter itV = begin(*pointerStr2, Standard());
                TVertexSetIter itVEnd = end(*pointerStr2, Standard());
                for(;itV != itVEnd;++itV) appendValue(*pointerAlign, *itV, Generous());
                ++pointerAlign;
                ++pointerStr2; ++posStr2;
            }
            // Gaps in seq 2
            while (i != posStr1) {
                TVertexSetIter itV = begin(*pointerStr1, Standard());
                TVertexSetIter itVEnd = end(*pointerStr1, Standard());
                for(;itV != itVEnd;++itV) appendValue(*pointerAlign, *itV, Generous());
                ++pointerAlign;
                ++pointerStr1; ++posStr1;
            }
        }

        // Matches
        if (p>=0) {
            TVertexSetIter itV = begin(*pointerStr1, Standard());
            TVertexSetIter itVEnd = end(*pointerStr1, Standard());
            for(;itV != itVEnd;++itV) appendValue(*pointerAlign, *itV, Generous());
            TVertexSetIter itV2 = begin(*pointerStr2, Standard());
            TVertexSetIter itVEnd2 = end(*pointerStr2, Standard());
            for(;itV2 != itVEnd2;++itV2) appendValue(*pointerAlign, *itV2, Generous());
            ++pointerAlign;
            ++pointerStr1; ++posStr1;
            ++pointerStr2; ++posStr2;
            --p;
        }
    }
}




//////////////////////////////////////////////////////////////////////////////

// TODO(holtgrew): Where does this belong? This is also there in graph_algorithms.h...

template<typename TStringSet, typename TCargo, typename TSpec, typename TString, typename TOutString>
inline TCargo
heaviestCommonSubsequence(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
                          TString const& str1,
                          TString const& str2,
                          TOutString& align)
{
    SEQAN_CHECKPOINT
    typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
    typedef typename Size<TGraph>::Type TSize;
    typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
    typedef typename Value<TString>::Type TVertexSet;

    TSize m = length(str1);  // How many sets of vertex descriptors in seq1
    TSize n = length(str2);  // How many sets of vertex descriptors in seq2

    // Size of the sequences
    // Note for profile alignments every member of the sequence is a String!!! of vertex descriptors

    // Fill the vertex to position map for str1
    // Remember for each vertex descriptor the position in the sequence
    typedef String<TSize> TMapVertexPos;
    TMapVertexPos map;
    resize(map, getIdUpperBound(_getVertexIdManager(g)), MaxValue<TSize>::VALUE);
    typedef typename Iterator<TString const, Standard>::Type TStringIterConst;
    typedef typename Iterator<TVertexSet const, Standard>::Type TVertexSetIterConst;
    TStringIterConst itStr1 = begin(str1, Standard());
    TStringIterConst itStrEnd1 = end(str1, Standard());
    TSize pos = 0;
    TVertexSetIterConst itV;
    TVertexSetIterConst itVEnd;
    for(;itStr1 != itStrEnd1;++itStr1, ++pos) {
        itV = begin(*itStr1, Standard());
        itVEnd = end(*itStr1, Standard());
        for(;itV != itVEnd;++itV) map[*itV] = pos;
    }

    // We could create the full graph -> too expensive
    // Remember which edges are actually present
    typedef String<TSize> TOccupiedPositions;
    typedef typename Iterator<TOccupiedPositions, Standard>::Type TOccIter;
    TOccupiedPositions occupiedPositions;
    TStringIterConst itStr2 = begin(str2, Standard());
    TStringIterConst itStrEnd2 = end(str2, Standard());
    TSize posItStr2 = 0;
    TSize pPos = 0;
    for(;itStr2 != itStrEnd2;++itStr2, ++posItStr2) {
        itV = begin(*itStr2, Standard());
        itVEnd = end(*itStr2, Standard());
        for(;itV != itVEnd;++itV) {
            TOutEdgeIterator itOut(g, *itV);
            for(;!atEnd(itOut); ++itOut) {
                // Target vertex must be in the map
                pPos = map[targetVertex(itOut)];
                if (pPos != MaxValue<TSize>::VALUE)
                    appendValue(occupiedPositions, pPos * n + (TSize) (n - posItStr2 - 1), Generous());
            }
        }
    }
    std::sort(begin(occupiedPositions, Standard()), end(occupiedPositions, Standard()));
    // Get all occupied positions
    typedef String<TSize> TSlotToPos;
    typedef typename Iterator<TSlotToPos, Standard>::Type TSlotToPosIter;
    TSlotToPos slotToPos;
    TSize counter = 0;
    TSize oldVal = MaxValue<TSize>::VALUE;
    TOccIter occIt = begin(occupiedPositions, Standard());
    TOccIter occItEnd = end(occupiedPositions, Standard());
    for(;occIt != occItEnd; ++occIt) {
        if (oldVal != *occIt) {
            appendValue(slotToPos, *occIt, Generous());
            oldVal = *occIt;
            ++counter;
        }
    }
    clear(occupiedPositions);

    // Walk through str2 and fill in the weights of the actual edges
    typedef String<TCargo> TWeights;
    TWeights weights;
    resize(weights, length(slotToPos), 0);
    itStr2 = begin(str2, Standard());
    posItStr2 = 0;
    for(;itStr2 != itStrEnd2;++itStr2, ++posItStr2) {
        itV = begin(*itStr2, Standard());
        itVEnd = end(*itStr2, Standard());
        for(;itV != itVEnd;++itV) {
            TOutEdgeIterator itOut(g, *itV);
            for(;!atEnd(itOut); ++itOut) {
                // Target vertex must be in the map
                pPos = map[targetVertex(itOut)];
                if ( pPos != MaxValue<TSize>::VALUE)
                    weights[std::distance(begin(slotToPos, Standard()), std::lower_bound(begin(slotToPos, Standard()), end(slotToPos, Standard()), pPos * n + (TSize) (n - posItStr2 - 1)))] += (TCargo) cargo(*itOut);
            }
        }
    }
    clear(map);

    // Now the tough part: Find the right number for a given position
    typedef String<TSize> TSequenceString;
    typedef typename Iterator<TSequenceString, Standard>::Type TSeqIter;
    TSequenceString seq;
    resize(seq, length(slotToPos));
    TSeqIter itSeq = begin(seq, Standard());
    TSlotToPosIter itSlotPos = begin(slotToPos, Standard());
    TSlotToPosIter itSlotPosEnd = end(slotToPos, Standard());
    for(;itSlotPos != itSlotPosEnd; ++itSlotPos, ++itSeq)
        *itSeq = n - 1 - (*itSlotPos % n);

    // Calculate the heaviest increasing subsequence
    String<TSize> positions;
    TCargo score = (TCargo) heaviestIncreasingSubsequence(seq, weights, positions);

    // Retrieve the alignment sequence
    _heaviestCommonSubsequence(g, slotToPos, positions, m, n, str1, str2, align);

    return score;

}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TString>
inline TCargo
heaviestCommonSubsequence(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
                          TString const& str1,
                          TString const& str2)
{
    SEQAN_CHECKPOINT
    Nothing noth;
    return heaviestCommonSubsequence(g, str1, str2, noth);
}

template <typename TStringSet, typename TCargo, typename TSpec, typename TSequenceH, typename TSequenceV, typename TId, typename TPos, typename TTraceValue>
inline void
_alignTracePrint(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
                 TSequenceH const &,
                 TSequenceV const &,
                 TId const id1,
                 TPos const pos1,
                 TId const id2,
                 TPos const pos2,
                 TPos const segLen,
                 TTraceValue const tv)
{
    // TraceBack values
    TTraceValue Diagonal = 0; TTraceValue Horizontal = 1; TTraceValue Vertical = 2;

    if (segLen == 0) return;

    if (tv == Horizontal)
        addVertex(g, id1, pos1, segLen);
    else if (tv == Vertical)
        addVertex(g, id2, pos2, segLen);
    else if (tv == Diagonal)
        addEdge(g, addVertex(g, id1, pos1, segLen), addVertex(g, id2, pos2, segLen));
}

// ----------------------------------------------------------------------------
// Function _alignTracePrint()                            [String<Fragment<> >]
// ----------------------------------------------------------------------------

// TODO(holtgrew): Are TSequenceH and TSequenceV used *anywhere*.

template <typename TFragment, typename TSequenceH, typename TSequenceV, typename TId, typename TPos, typename TTraceValue>
inline void
_alignTracePrint(String<TFragment>& matches,
                 TSequenceH const &,
                 TSequenceV const &,
                 TId const id1,
                 TPos const pos1,
                 TId const id2,
                 TPos const pos2,
                 TPos const seqLen,
                 TTraceValue const tv)
{
    // Only the diagonal case
    if ((seqLen) && (tv == 0))
        appendValue(matches, TFragment(id1, pos1, id2, pos2, seqLen), Generous());
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_GRAPH_ALIGN_GRAPH_IMPL_ALIGN_H_
