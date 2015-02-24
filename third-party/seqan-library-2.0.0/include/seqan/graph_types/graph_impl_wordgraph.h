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

#ifndef SEQAN_HEADER_GRAPH_IMPL_WORDGRAPH_H
#define SEQAN_HEADER_GRAPH_IMPL_WORDGRAPH_H

namespace SEQAN_NAMESPACE_MAIN
{

template <typename TSpec = Default>
struct WordGraph;

//////////////////////////////////////////////////////////////////////////////
// Graph - WordGraph
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/*!
 * @class WordGraph
 * @extends Automaton
 * @headerfile <seqan/graph_type.h>
 * @brief A special automaton that stores words instead of single characters along its edges.
 *
 * @signature template <[typename TAlphabet[, typename Spec]]>
 *            class Graph<Automaton<TAlphabet, String<TAlphabet>, WordGraph<TSpec> >;
 *
 * @tparam TAlphabet The alphabet type.
 * @tparam TSpec     The specializing types.
 */

template<typename TAlphabet, typename TSpec>
class Graph<Automaton<TAlphabet, String<TAlphabet>, WordGraph<TSpec> > >
{
    public:
        typedef typename VertexIdHandler<Graph>::Type TVertexIdManager_;
        typedef typename EdgeIdHandler<Graph>::Type TEdgeIdManager_;
        typedef typename VertexDescriptor<Graph>::Type TVertexDescriptor_;
        typedef typename EdgeType<Graph>::Type TEdge_;

        String<AutomatonEdgeArray<TEdge_, TAlphabet> > data_vertex;        // List of tables
        TVertexIdManager_ data_id_managerV;
        TEdgeIdManager_ data_id_managerE;
        TVertexDescriptor_ data_root;


//____________________________________________________________________________


        Graph() : data_root(0) {
            SEQAN_CHECKPOINT
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
            _copyGraph(_other, *this);
            return *this;
        }
};


//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TSpec, typename TVertexDescriptor>
inline typename EdgeDescriptor<Graph<Automaton<TAlphabet, String<TAlphabet>, WordGraph<TSpec> > > >::Type
addEdge(Graph<Automaton<TAlphabet, String<TAlphabet>, WordGraph<TSpec> > >& g,
        TVertexDescriptor const source,
        TVertexDescriptor const target,
        String<TAlphabet> const & label)
{
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT(idInUse(g.data_id_managerV, source));
    SEQAN_ASSERT(idInUse(g.data_id_managerV, target));

    typedef Graph<Automaton<TAlphabet, String<TAlphabet>, WordGraph<TSpec> > > TGraph;
    typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    typedef typename Id<TGraph>::Type TId;

    TAlphabet firstChar = getValue(label, 0);
    TEdgeDescriptor e = findEdge(g, source, firstChar);
    TId id = obtainId(g.data_id_managerE);
    _assignId(e, id);
    assignTarget(e, target);
    String<TAlphabet> suf(suffix(label,1));
    assignCargo(e, suf);
    return e;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TSpec, typename TVertexDescriptor, typename TChars>
inline typename EdgeDescriptor<Graph<Automaton<TAlphabet, String<TAlphabet>, WordGraph<TSpec> > > >::Type
addEdge(Graph<Automaton<TAlphabet, String<TAlphabet>, WordGraph<TSpec> > >& g,
        TVertexDescriptor const source,
        TVertexDescriptor const target,
        TChars const* chars)
{
    SEQAN_CHECKPOINT
    return addEdge(g,source,target,String<TAlphabet>(chars));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TSpec, typename TVertexDescriptor, typename TLabel, typename TEdgeCargo>
inline typename EdgeDescriptor<Graph<Automaton<TAlphabet, String<TAlphabet>, WordGraph<TSpec> > > >::Type
addEdge(Graph<Automaton<TAlphabet, String<TAlphabet>, WordGraph<TSpec> > >& /*g*/,
        TVertexDescriptor const /*source*/,
        TVertexDescriptor const /*target*/,
        TLabel const /*label*/,
        TEdgeCargo const /*cargo*/)
{
    // No additional cargo allowed. Cargo is used for the words in the graph.
    // Use external property map.
    SEQAN_ASSERT_FAIL("No additional cargo allowed. Cargo is used for the words in the graph. Use external property map.");
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TSpec, typename TVertexDescriptor>
inline void
removeEdge(Graph<Automaton<TAlphabet, String<TAlphabet>, WordGraph<TSpec> > >& g,
        TVertexDescriptor const source,
        TVertexDescriptor const target,
        String<TAlphabet> const& label)
{
    SEQAN_CHECKPOINT;
    (void)target;  // In case it is compiled without assertions.
    SEQAN_ASSERT(idInUse(g.data_id_managerV, source));
    SEQAN_ASSERT(idInUse(g.data_id_managerV, target));

    TAlphabet firstChar = getValue(label, 0);
    removeEdge(g, findEdge(g,source, firstChar));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TAlphabet, typename TCargo, typename TSpec>
inline void
write(TFile & target,
      Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> > > const& g)
{
//IOREV _nodoc_
    typedef Graph<Automaton<TAlphabet, String<TAlphabet>, WordGraph<TSpec> > > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename EdgeType<TGraph>::Type TEdge;
    typedef typename Size<TAlphabet>::Type TSize;
    TSize table_length = ValueSize<TAlphabet>::VALUE;
    TVertexDescriptor nilVal = getNil<TVertexDescriptor>();

    write(target, "WordGraph - Directed:\n");
    typedef typename Iterator<String<AutomatonEdgeArray<TEdge, TAlphabet> > const, Rooted>::Type TIterConst;
    for(TIterConst it = begin(g.data_vertex);!atEnd(it);goNext(it)) {
        if (!idInUse(g.data_id_managerV, position(it))) continue;
        TVertexDescriptor sourceVertex = position(it);
        for(TSize i=0;i<table_length;++i) {
            TEdge const* ed = &g.data_vertex[sourceVertex].data_edge[i];
            if (getTarget(ed) ==  nilVal) continue;
            appendNumber(target, (int)sourceVertex);
            write(target, "->");
            appendNumber(target, (int)getTarget(ed));
            writeValue(target, ' ');
            writeValue(target, ' ');
            write(target, "Label: ");
            writeValue(target, TAlphabet(i));
            write(target, CharString(getCargo(ed)));
            writeValue(target, '\n');
        }
    }
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn WordGraph#getSuccessor
 * @brief Gets the successor for a given vertex and an edge label.
 *
 * @signature TVertexDescriptor getSuccessor(a, v, str);
 *
 * @param[in] a   The WordGraph to query for its successor.
 * @param[in] v   The descriptor fo the vertex to get the successor for.
 * @param[in] str The label.
 *
 * @return TVertexDescriptor A vertex descriptor or nil if successor is not defined.
 *
 * @see WordGraph#parseString
 * @see Graph#getNil
 */

template<typename TAlphabet, typename TSpec, typename TVertexDescriptor, typename TCharacters>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, String<TAlphabet>, WordGraph<TSpec> > > >::Type
getSuccessor(Graph<Automaton<TAlphabet, String<TAlphabet>, WordGraph<TSpec> > > const& g,
             TVertexDescriptor vertex,
             TCharacters const& chars)
{
    SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex));
    typedef Graph<Automaton<TAlphabet, String<TAlphabet>, WordGraph<TSpec> > > TGraph;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;
    TEdgeStump* ed = findEdge(g, vertex, getValue(chars, 0));
    if (getCargo(ed) == suffix(chars, 1)) {
        return getTarget(ed);
    } else {
        return getNil<TVertexDescriptor>();
    }
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TSpec, typename TVertexDescriptor, typename TCharacters>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, String<TAlphabet>, WordGraph<TSpec> > > >::Type
getSuccessor(Graph<Automaton<TAlphabet, String<TAlphabet>, WordGraph<TSpec> > > const& g,
             TVertexDescriptor vertex,
             TCharacters const* chars)
{
    SEQAN_CHECKPOINT
    return getSuccessor(g,vertex,String<TAlphabet>(chars));
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn WordGraph#parseString
 * @brief Parses a string one character at a time and moves accordingly in the WordGraph.
 *
 * @signature TVertexDescriptor parseString(a, v, beginIt, endIt);
 * @signature TVertexDescriptor parseString(a, v, str);
 *
 * @param[in]     a       An WordGraph.
 * @param[in]     v       The descriptor of the vertex to start at.
 * @param[in]     str     The @link ContainerConcept @endlink to parse.
 * @param[in,out] beginIt Begin iterator to sequence to parse.  Set to the first character that could not be parsed
 *                        or to the value of endIt if all of the string was parsed.
 * @param[in]     endIt   End iterator to sequence to parse.
 *
 * @return TVertexDescriptor The vertex descriptor of the state that was reached after parsing.
 *
 * The parsing stops before @link WordGraph#getSuccessor @endlink reaches the nil state or if the complete sequence is
 * read.
 */

template<typename TAlphabet, typename TSpec, typename TVertexDescriptor, typename TIterator>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, String<TAlphabet>, WordGraph<TSpec> > > >::Type
parseString(Graph<Automaton<TAlphabet, String<TAlphabet>, WordGraph<TSpec> > > const& g,
            TVertexDescriptor const vertex,
            TIterator beginIt,
            TIterator endIt)
{
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex));
    typedef Graph<Automaton<TAlphabet, String<TAlphabet>, WordGraph<TSpec> > > TGraph;
    typedef typename Size<TGraph>::Type TSize;
    TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
    TVertexDescriptor succ = vertex;
    while (beginIt!=endIt) {
        String<TAlphabet> label(*beginIt);
        TSize range = 1;
        TVertexDescriptor tmp = getSuccessor(g,succ,label);
        while ((tmp == nilVal) &&
                (beginIt+range != endIt))
        {
            appendValue(label, *(beginIt + range));
            tmp = getSuccessor(g,succ,label);
            ++range;
        }
        if (tmp == nilVal) break;
        succ = tmp;
        beginIt = beginIt+range;
    }
    return succ;
}

/*!
 * @fn WordGraph#canParseString
 * @brief Test whether an WordGraph can parse a string completely.
 *
 * @signature bool canParseString(a[, v], str);
 *
 * @param[in] a   The WordGraph to use for parsing.
 * @param[in] v   Optionally, the descriptor of the vertex to start at.  Defaults to the root.
 * @param[in] str The string to parse.
 *
 * @return bool <tt>true</tt> if the WordGraph parses <tt>str</tt> , starting at <tt>v</tt>, completely and
 *              <tt>false</tt> otherwise.
 *
 * @section Remarks
 *
 * This has not implemented yet.
 */

// TODO(holtgrew): Not implemented yet!


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
