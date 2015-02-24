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

#ifndef SEQAN_HEADER_GRAPH_IMPL_HMM_H
#define SEQAN_HEADER_GRAPH_IMPL_HMM_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Graph - HMM
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/*!
 * @class HmmGraph
 * @extends Graph
 * @headerfile <seqan/graph_types.h>
 * @brief A grap representing a Hidden Markov Model (HMM).
 *
 * An HMM is a directed graph with edges labeled with transition probabilities and emission profiles for each vertex.
 * Vertices correspond to states in an HMM.
 *
 * @signature template <[typename TAlphabet[, typename TCargo[, typename TSpec]]]>
 *            class Graph<Hmm<TAlphabet, TCargo, TSpec> >;
 *
 * @tparam TAlphabet The alphabet of the HMM.  Default: @link Dna @endlink.
 * @tparam TCargo    The edge cargos.  Default: <tt>void</tt>.
 * @tparam TSpec     The specialization type.  Default: <tt>Default</tt>.  Use <tt>WithoutEdgeIds</tt> here to omit
 *                   edge ids.  NB: if edges do not store ids then external property maps do not work.
 */

template<typename TAlphabet, typename TCargo, typename TSpec>
class Graph<Hmm<TAlphabet, TCargo, TSpec> >
{
public:
    typedef typename VertexDescriptor<Graph>::Type TVertexDescriptor_;

    // HMM Model
    Graph<Directed<TCargo, TSpec> > data_model;

    // Emission probabilities
    String<TCargo> data_emission;

    // Silent state map
    String<bool> data_silent;

    // Begin and end state
    TVertexDescriptor_ data_begin;
    TVertexDescriptor_ data_end;

    Graph() : data_begin(), data_end()
    {}

    ~Graph()
    {
        clear(*this);
    }

    Graph(Graph const & _other)
    {
        _copyGraph(_other, *this);
    }

    Graph const& operator = (Graph const & _other)
    {
        if (this == &_other)
            return *this;
        clear(*this);
        _copyGraph(_other, *this);
        return *this;
    }
};

//////////////////////////////////////////////////////////////////////////////
// INTERNAL FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
inline String<typename EdgeType<Graph<Hmm<TAlphabet, TCargo, TSpec> > >::Type*>&
_getVertexString(Graph<Hmm<TAlphabet, TCargo, TSpec> > const& g) {
    SEQAN_CHECKPOINT
    typedef Graph<Hmm<TAlphabet, TCargo, TSpec> > TGraph;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;
    return const_cast<String<TEdgeStump*>&>(g.data_model.data_vertex);
}

/////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
inline typename VertexIdHandler<Graph<Hmm<TAlphabet, TCargo, TSpec> > >::Type&
_getVertexIdManager(Graph<Hmm<TAlphabet, TCargo, TSpec> > const& g) {
    SEQAN_CHECKPOINT
    typedef Graph<Hmm<TAlphabet, TCargo, TSpec> > TGraph;
    typedef typename VertexIdHandler<TGraph>::Type TVertexIdManager;
    return const_cast<TVertexIdManager&>(g.data_model.data_id_managerV);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
inline typename EdgeIdHandler<Graph<Hmm<TAlphabet, TCargo, TSpec> > >::Type&
_getEdgeIdManager(Graph<Hmm<TAlphabet, TCargo, TSpec> > const& g) {
    SEQAN_CHECKPOINT
    typedef Graph<Hmm<TAlphabet, TCargo, TSpec> > TGraph;
    typedef typename EdgeIdHandler<TGraph>::Type TEdgeIdManager;
    return const_cast<TEdgeIdManager&>(g.data_model.data_id_managerE);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
inline void
_copyGraph(Graph<Hmm<TAlphabet, TCargo, TSpec> > const& source,
           Graph<Hmm<TAlphabet, TCargo, TSpec> >& dest,
           bool transp)
{
    SEQAN_CHECKPOINT
    clear(dest);
    if (transp) {
        transpose(source.data_model, dest.data_model);
        dest.data_emission = source.data_emission;
        dest.data_silent = source.data_silent;
        dest.data_begin = source.data_end;
        dest.data_end = source.data_begin;
    } else {
        dest.data_model = source.data_model;
        dest.data_emission = source.data_emission;
        dest.data_silent = source.data_silent;
        dest.data_begin = source.data_begin;
        dest.data_end = source.data_end;
    }
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
inline void
_copyGraph(Graph<Hmm<TAlphabet, TCargo, TSpec> > const& source,
           Graph<Hmm<TAlphabet, TCargo, TSpec> >& dest)
{
    _copyGraph(source, dest, false);
}

//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
inline void
transpose(Graph<Hmm<TAlphabet, TCargo, TSpec> > const& source,
          Graph<Hmm<TAlphabet, TCargo, TSpec> >& dest)
{
    SEQAN_CHECKPOINT
    _copyGraph(source, dest, true);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
inline void
transpose(Graph<Hmm<TAlphabet, TCargo, TSpec> >& g)
{
    SEQAN_CHECKPOINT
    Graph<Hmm<TAlphabet, TCargo, TSpec> > dest;
    _copyGraph(g, dest, true);
    g = dest;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
inline typename Size<Graph<Hmm<TAlphabet, TCargo, TSpec> > >::Type
numEdges(Graph<Hmm<TAlphabet, TCargo, TSpec> > const& g)
{
    SEQAN_CHECKPOINT
    return numEdges(g.data_model);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
inline typename Size<Graph<Hmm<TAlphabet, TCargo, TSpec> > >::Type
numVertices(Graph<Hmm<TAlphabet, TCargo, TSpec> > const& g)
{
    SEQAN_CHECKPOINT
    return numVertices(g.data_model);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
inline bool
empty(Graph<Hmm<TAlphabet, TCargo, TSpec> > const& g)
{
    SEQAN_CHECKPOINT
    return empty(g.data_model);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
inline void
clearEdges(Graph<Hmm<TAlphabet, TCargo, TSpec> >& g)
{
    SEQAN_CHECKPOINT
    clearEdges(g.data_model);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
inline void
clearVertices(Graph<Hmm<TAlphabet, TCargo, TSpec> >& g)
{
    g.data_begin = 0;
    g.data_end = 0;
    clear(g.data_emission);
    clear(g.data_silent);
    clearVertices(g.data_model);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
inline void
clear(Graph<Hmm<TAlphabet, TCargo, TSpec> >& g)
{
    SEQAN_CHECKPOINT
    clearVertices(g);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename Size<Graph<Hmm<TAlphabet, TCargo, TSpec> > >::Type
outDegree(Graph<Hmm<TAlphabet, TCargo, TSpec> > const& g,
          TVertexDescriptor const vertex)
{
    SEQAN_CHECKPOINT
    return outDegree(g.data_model, vertex);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename Size<Graph<Hmm<TAlphabet, TCargo, TSpec> > >::Type
inDegree(Graph<Hmm<TAlphabet, TCargo, TSpec> > const& g,
         TVertexDescriptor const vertex)
{
    SEQAN_CHECKPOINT
    return inDegree(g.data_model, vertex);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename Size<Graph<Hmm<TAlphabet, TCargo, TSpec> > >::Type
degree(Graph<Hmm<TAlphabet, TCargo, TSpec> > const& g,
       TVertexDescriptor const vertex)
{
    SEQAN_CHECKPOINT
    return degree(g.data_model, vertex);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
inline typename VertexDescriptor<Graph<Hmm<TAlphabet, TCargo, TSpec> > >::Type
addVertex(Graph<Hmm<TAlphabet, TCargo, TSpec> >& g,
          bool silent)
{
    SEQAN_CHECKPOINT
    typedef Graph<Hmm<TAlphabet, TCargo, TSpec> > TGraph;
    typedef typename Size<TAlphabet>::Type TSize;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename Iterator<String<TCargo>, Standard>::Type TEmisIter;
    TSize alph_size = ValueSize<TAlphabet>::VALUE;

    TVertexDescriptor vd = addVertex(g.data_model);
    if (length(g.data_emission) / alph_size <= vd) resize(g.data_emission, (vd + 1) * alph_size, Generous());
    if (length(g.data_silent) <= vd) resize(g.data_silent, (vd + 1), Generous());
    g.data_silent[vd] = silent;
    TEmisIter it = begin(g.data_emission, Standard());
    it += vd * alph_size;
    for(TSize counter = 0; counter < alph_size; ++counter, ++it) *it = (TCargo) 0.0;
    return vd;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
inline typename VertexDescriptor<Graph<Hmm<TAlphabet, TCargo, TSpec> > >::Type
addVertex(Graph<Hmm<TAlphabet, TCargo, TSpec> >& g)
{
    SEQAN_CHECKPOINT
    return addVertex(g, false);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TEmission>
inline typename VertexDescriptor<Graph<Hmm<TAlphabet, TCargo, TSpec> > >::Type
addVertex(Graph<Hmm<TAlphabet, TCargo, TSpec> > & g,
          String<TEmission> const & emis)
{
    typedef Graph<Hmm<TAlphabet, TCargo, TSpec> > TGraph;
    typedef typename Size<TAlphabet>::Type TSize;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename Iterator<String<TCargo>, Standard>::Type TEmisIter;
    typedef typename Iterator<String<TEmission> const, Standard>::Type TInputIter;
    TSize alph_size = ValueSize<TAlphabet>::VALUE;

    SEQAN_ASSERT_EQ(alph_size, length(emis));

    TVertexDescriptor vd = addVertex(g);
    TEmisIter it = begin(g.data_emission, Standard());
    it += vd * alph_size;
    TInputIter itIn = begin(emis, Standard());
    for (TSize counter = 0; counter < alph_size; ++counter, ++itIn, ++it)
        *it = *itIn;
    return vd;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TEmission>
inline typename VertexDescriptor<Graph<Hmm<TAlphabet, TCargo, TSpec> > >::Type
addVertex(Graph<Hmm<TAlphabet, TCargo, TSpec> >& g,
          TEmission const& emis,
          bool silent)
{
    SEQAN_CHECKPOINT
    typedef Graph<Hmm<TAlphabet, TCargo, TSpec> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    TVertexDescriptor vd = addVertex(g, emis);
    g.data_silent[vd] = silent;
    return vd;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void
removeVertex(Graph<Hmm<TAlphabet, TCargo, TSpec> >& g,
             TVertexDescriptor const v)
{
    SEQAN_CHECKPOINT
    // Remove the vertex
    removeVertex(g.data_model,v);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename EdgeDescriptor<Graph<Hmm<TAlphabet, TCargo, TSpec> > >::Type
addEdge(Graph<Hmm<TAlphabet, TCargo, TSpec> >& g,
        TVertexDescriptor const source,
        TVertexDescriptor const target)
{
    SEQAN_CHECKPOINT
    return addEdge(g, source, target, (TCargo) 0.0);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor, typename TCargo2>
inline typename EdgeDescriptor<Graph<Hmm<TAlphabet, TCargo, TSpec> > >::Type
addEdge(Graph<Hmm<TAlphabet, TCargo, TSpec> >& g,
        TVertexDescriptor const source,
        TVertexDescriptor const target,
        TCargo2 const cargo)
{
    SEQAN_CHECKPOINT
    return addEdge(g.data_model, source, target, (TCargo) cargo);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void
removeEdge(Graph<Hmm<TAlphabet, TCargo, TSpec> >& g,
           TVertexDescriptor const source,
           TVertexDescriptor const target)
{
    SEQAN_CHECKPOINT
    removeEdge(g.data_model, source, target);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TEdgeDescriptor>
inline void
removeEdge(Graph<Hmm<TAlphabet, TCargo, TSpec> >& g,
           TEdgeDescriptor const edge)
{
    SEQAN_CHECKPOINT
    removeEdge(g.data_model, sourceVertex(g.data_model,edge), targetVertex(g.data_model,edge));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void
removeOutEdges(Graph<Hmm<TAlphabet, TCargo, TSpec> >& g,
               TVertexDescriptor const v)
{
    SEQAN_CHECKPOINT
    removeOutEdges(g.data_model, v);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void
removeInEdges(Graph<Hmm<TAlphabet, TCargo, TSpec> >& g,
              TVertexDescriptor const v)
{
    SEQAN_CHECKPOINT
    removeInEdges(g.data_model,v);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TEdgeDescriptor>
inline typename VertexDescriptor<Graph<Hmm<TAlphabet, TCargo, TSpec> > >::Type
targetVertex(Graph<Hmm<TAlphabet, TCargo, TSpec> > const& g,
             TEdgeDescriptor const edge)
{
    SEQAN_CHECKPOINT
    return targetVertex(g.data_model, edge);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TEdgeDescriptor>
inline typename VertexDescriptor<Graph<Hmm<TAlphabet, TCargo, TSpec> > >::Type
sourceVertex(Graph<Hmm<TAlphabet, TCargo, TSpec> > const& g,
             TEdgeDescriptor const edge)
{
    SEQAN_CHECKPOINT
    return sourceVertex(g.data_model, edge);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TMatrix>
inline void
getAdjacencyMatrix(Graph<Hmm<TAlphabet, TCargo, TSpec> > const& g,
                   TMatrix& mat)
{
    getAdjacencyMatrix(g.data_model, mat);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename EdgeDescriptor<Graph<Hmm<TAlphabet, TCargo, TSpec> > >::Type
findEdge(Graph<Hmm<TAlphabet, TCargo, TSpec> > const& g,
         TVertexDescriptor const v,
         TVertexDescriptor const w)
{
    SEQAN_CHECKPOINT
    return findEdge(g.data_model, v, w);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TAlphabet, typename TCargo, typename TSpec>
inline void
write(TFile & target,
      Graph<Hmm<TAlphabet, TCargo, TSpec> > const& g)
{
//IOREV _nodoc_
    typedef Graph<Hmm<TAlphabet, TCargo, TSpec> > TGraph;
    typedef typename Size<TAlphabet>::Type TSize;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;
    typedef typename Iterator<String<TEdgeStump*>, Standard>::Type TIterConst;
    TSize alph_size = ValueSize<TAlphabet>::VALUE;

    typedef typename Iterator<String<TCargo> const, Standard>::Type TEmisIter;

    // Alphabet
    write(target, "Alphabet:\n");
    writeValue(target, '{');
    for(TSize counter = 0; counter<alph_size-1;++counter) {
        writeValue(target,TAlphabet(counter));
        writeValue(target, ',');
    }
    writeValue(target,TAlphabet(alph_size-1));
    writeValue(target, '}');
    writeValue(target, '\n');

    // States
    write(target, "States:\n");
    writeValue(target, '{');
    TIterConst it = begin(_getVertexString(g), Standard());
    TIterConst itEnd = end(_getVertexString(g), Standard());
    bool first = true;
    for(TSize pos = 0;it!=itEnd;goNext(it), ++pos) {
        if (!idInUse(_getVertexIdManager(g), pos)) continue;
        if (!first) writeValue(target, ',');
        else first = false;
        appendNumber(target, (int)pos);
        if (isSilent(g, pos)) write(target, " (Silent)");
    }
    writeValue(target, '}');
    writeValue(target, '\n');

    // Begin and end state
    write(target, "Begin state: ");
    appendNumber(target, (int)getBeginState(g));
    writeValue(target, '\n');
    write(target, "End state: ");
    appendNumber(target, (int)getEndState(g));
    writeValue(target, '\n');

    // Transition probabilities
    write(target, "Transition probabilities:\n");
    itEnd = end(_getVertexString(g));
    it = begin(_getVertexString(g));
    for(TSize pos = 0;it!=itEnd;goNext(it), ++pos) {
        if (!idInUse(_getVertexIdManager(g), pos)) continue;
        TEdgeStump* current = getValue(it);
        appendNumber(target, (int)pos);
        write(target, " -> ");
        first = true;
        while(current!=0) {
            if (!first) writeValue(target, ',');
            else first = false;
            appendNumber(target, (int)getTarget(current));
            write(target, " (");
            write(target, (double)cargo(current));
            write(target, ") ");
            current=getNextT(current);
        }
        writeValue(target, '\n');
    }

    // Emission probabilities
    write(target, "Emission probabilities:\n");
    TEmisIter itEmis = begin(g.data_emission, Standard());
    itEnd = end(_getVertexString(g), Standard());
    it = begin(_getVertexString(g), Standard());
    first = true;
    for(TSize pos = 0;it!=itEnd;++it, ++pos) {
        if (!idInUse(_getVertexIdManager(g), pos)) continue;
        if (isSilent(g, pos)) continue;
        if (!first) writeValue(target, '\n');
        else first = false;
        appendNumber(target, (int)pos);
        write(target, ": ");
        bool my_first = true;
        itEmis = begin(g.data_emission, Standard());
        itEmis += pos * alph_size;
        for(TSize counter = 0; counter < alph_size; ++itEmis, ++counter) {
            if (!my_first) writeValue(target, ',');
            else my_first = false;
            writeValue(target, TAlphabet(counter));
            write(target, " (");
            appendNumber(target, (double)*itEmis);
            write(target, ") ");

        }
    }

}


//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn HmmGraph#assignBeginState
 * @headerfile <seqan/graph_types.h>
 * @brief Assigns a begin state.
 *
 * @signature void assignBeginState(hmm, v);
 *
 * @param[in,out] hmm The @link HmmGraph @endlink to assign the begin state of.
 * @param[in]     v   The descriptor of the vertex to assign as begin state.
 */

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void
assignBeginState(Graph<Hmm<TAlphabet, TCargo, TSpec> >& g,
                 TVertexDescriptor const vertex)
{
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT(idInUse(_getVertexIdManager(g), vertex));

    g.data_begin = vertex;
    g.data_silent[vertex] = true;
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn HmmGraph#assignEndState
 * @headerfile <seqan/graph_types.h>
 * @brief Assigns a end state.
 *
 * @signature void assignEndState(hmm, v);
 *
 * @param[in,out] hmm The @link HmmGraph @endlink to assign the end state of.
 * @param[in]     v   The descriptor of the vertex to assign as end state.
 */

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void
assignEndState(Graph<Hmm<TAlphabet, TCargo, TSpec> >& g,
               TVertexDescriptor const vertex)
{
    SEQAN_CHECKPOINT
    SEQAN_ASSERT(idInUse(_getVertexIdManager(g), vertex));

    g.data_end = vertex;
    g.data_silent[vertex] = true;
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn HmmGraph#beginState
 * @brief Return a reference to the begin state vertex descriptor.
 *
 * @signature TVertexDescriptor beginState(hmm);
 *
 * @param[in] hmm The @link HmmGraph @endlink to get the vertex descriptor of.
 *
 * @return TVertexDescriptor A reference to the begin state vertex descriptor.
 */

template<typename TAlphabet, typename TCargo, typename TSpec>
inline typename VertexDescriptor<Graph<Hmm<TAlphabet, TCargo, TSpec> > >::Type&
beginState(Graph<Hmm<TAlphabet, TCargo, TSpec> >& g)
{
    SEQAN_CHECKPOINT
    return (g.data_begin);
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn HmmGraph#endState
 * @brief Return a reference to the end state vertex descriptor.
 *
 * @signature TVertexDescriptor endState(hmm);
 *
 * @param[in] hmm The @link HmmGraph @endlink to get the vertex descriptor of.
 *
 * @return TVertexDescriptor A reference to the end state vertex descriptor.
 */

template<typename TAlphabet, typename TCargo, typename TSpec>
inline typename VertexDescriptor<Graph<Hmm<TAlphabet, TCargo, TSpec> > >::Type&
endState(Graph<Hmm<TAlphabet, TCargo, TSpec> >& g)
{
    SEQAN_CHECKPOINT
    return (g.data_end);
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn HmmGraph#getBeginState
 * @brief Return the begin state vertex descriptor.
 *
 * @signature TVertexDescriptor getBeginState(hmm);
 *
 * @param[in] hmm The @link HmmGraph @endlink to get the vertex descriptor of.
 *
 * @return TVertexDescriptor The begin state vertex descriptor.
 */

template<typename TAlphabet, typename TCargo, typename TSpec>
inline typename VertexDescriptor<Graph<Hmm<TAlphabet, TCargo, TSpec> > >::Type
getBeginState(Graph<Hmm<TAlphabet, TCargo, TSpec> > const& g)
{
    SEQAN_CHECKPOINT
    return (g.data_begin);
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn HmmGraph#getEndState
 * @brief Return the end state vertex descriptor.
 *
 * @signature TVertexDescriptor getEndState(hmm);
 *
 * @param[in] hmm The @link HmmGraph @endlink to get the vertex descriptor of.
 *
 * @return TVertexDescriptor The end state vertex descriptor.
 */

template<typename TAlphabet, typename TCargo, typename TSpec>
inline typename VertexDescriptor<Graph<Hmm<TAlphabet, TCargo, TSpec> > >::Type
getEndState(Graph<Hmm<TAlphabet, TCargo, TSpec> > const& g)
{
    SEQAN_CHECKPOINT
    return (g.data_end);
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn HmmGraph#getTransitionProbability
 * @brief Returns the transition probability, stored as the cargo.
 *
 * @signature TCargo getTransitionProbability(hmm, state1, state2);
 * @signature TCargo getTransitionProbability(hmm, e);
 *
 * @param[in] hmm    The @link HmmGraph @endlink to get transition probability for.
 * @param[in] state1 Return transition probability from start state1.
 * @param[in] state2 Return transition probability to target state2.
 * @param[in] e      An edge descriptor.
 *
 * @return TCargo Return the transition probability.
 */

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline TCargo
getTransitionProbability(Graph<Hmm<TAlphabet, TCargo, TSpec> > const& g,
                         TVertexDescriptor const state1,
                         TVertexDescriptor const state2)
{
    SEQAN_CHECKPOINT
    typedef Graph<Hmm<TAlphabet, TCargo, TSpec> > const TGraph;
    typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    TEdgeDescriptor e = findEdge(g, state1, state2);
    if (e == 0) return 0.0;
    else return cargo(e);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TEdgeDescriptor>
inline TCargo
getTransitionProbability(Graph<Hmm<TAlphabet, TCargo, TSpec> > const&,
                         TEdgeDescriptor const e)
{
    SEQAN_CHECKPOINT
    return getCargo(e);
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn HmmGraph#transitionProbability
 * @brief Returns a reference to transition probability, stored as the cargo.
 *
 * @signature TCargo getTransitionProbability(hmm, state1, state2);
 * @signature TCargo getTransitionProbability(hmm, e);
 *
 * @param[in] hmm    The @link HmmGraph @endlink to get transition probability for.
 * @param[in] state1 Return transition probability from start state1.
 * @param[in] state2 Return transition probability to target state2.
 * @param[in] e      An edge descriptor.
 *
 * @return TCargo Returns a reference to the transition probability.
 */

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline TCargo&
transitionProbability(Graph<Hmm<TAlphabet, TCargo, TSpec> >& g,
                      TVertexDescriptor const state1,
                      TVertexDescriptor const state2)
{
    SEQAN_CHECKPOINT
    typedef Graph<Hmm<TAlphabet, TCargo, TSpec> > TGraph;
    typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    TEdgeDescriptor e = findEdge(g, state1, state2);
    return cargo(e);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TEdgeDescriptor>
inline TCargo&
transitionProbability(Graph<Hmm<TAlphabet, TCargo, TSpec> >&,
                      TEdgeDescriptor e)
{
    SEQAN_CHECKPOINT
    return cargo(e);
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn HmmGraph#assignTransitionProbability
 * @brief Assigns a new transition probability to an existing edge.
 *
 * @signature void assignTransitionProbability(hmm, state1, state2, prob);
 * @signature void assignTransitionProbability(hmm, e, prob);
 *
 * @param[in] hmm    The @link HmmGraph @endlink to get transition probability for.
 * @param[in] state1 Assign transition probability from start state1.
 * @param[in] state2 Assign probability to target state2.
 * @param[in] e      Assign probability to edge with this edge descriptor.
 */

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor, typename TTransProb>
inline void
assignTransitionProbability(Graph<Hmm<TAlphabet, TCargo, TSpec> >& g,
                            TVertexDescriptor const state1,
                            TVertexDescriptor const state2,
                            TTransProb const t)
{
    SEQAN_CHECKPOINT
    typedef Graph<Hmm<TAlphabet, TCargo, TSpec> > TGraph;
    typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    TEdgeDescriptor e = findEdge(g, state1, state2);
    cargo(e) = t;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TEdgeDescriptor, typename TTransProb>
inline void
assignTransitionProbability(Graph<Hmm<TAlphabet, TCargo, TSpec> >&,
                            TEdgeDescriptor e,
                            TTransProb const t)
{
    SEQAN_CHECKPOINT
    cargo(e) = t;
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn HmmGraph#getEmissionsProbability
 * @brief Returns the emissions probability.
 *
 * @signature TCargo getEmissionProbability(hmm, state, symbol);
 *
 * @param[in] hmm    The @link HmmGraph @endlink to query.
 * @param[in] state  A vertex descriptor with the source state.
 * @param[in] symbol The symbol.  Type: <tt>TAlphabet</tt>.
 *
 * @return TCargo The emission probability.
 */

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline TCargo
getEmissionProbability(Graph<Hmm<TAlphabet, TCargo, TSpec> > const& g,
                       TVertexDescriptor const state,
                       TAlphabet const symbol)
{
    typedef typename Size<TAlphabet>::Type TSize;
    return g.data_emission[state * (TSize) ValueSize<TAlphabet>::VALUE + ordValue(symbol)];
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn HmmGraph#emissionProbability
 * @brief Returns a reference to the emissions probability.
 *
 * @signature TCargo emissionProbability(hmm, state, symbol);
 *
 * @param[in] hmm    The @link HmmGraph @endlink to query.
 * @param[in] state  A vertex descriptor with the source state.
 * @param[in] symbol The symbol.  Type: <tt>TAlphabet</tt>.
 *
 * @return TCargo Reference to the emission probability.
 */

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline TCargo&
emissionProbability(Graph<Hmm<TAlphabet, TCargo, TSpec> >& g,
                    TVertexDescriptor const state,
                    TAlphabet const symbol)
{
    typedef typename Size<TAlphabet>::Type TSize;
    return g.data_emission[state * (TSize) ValueSize<TAlphabet>::VALUE + ordValue(symbol)];
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn HmmGraph#assignEmissionsProbability
 * @brief Assigns a new emission probability.
 *
 * @signature void assignEmissionProbability(hmm, state, symbol, eProb);
 *
 * @param[in,out] hmm    The @link HmmGraph @endlink to query.
 * @param[in]     state  A vertex descriptor with the source state.
 * @param[in]     symbol The symbol.  Type: <tt>TAlphabet</tt>.
 * @param[in]     eProb  Emission probability to assign.
 */

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor, typename TEmisProb>
inline void
assignEmissionProbability(Graph<Hmm<TAlphabet, TCargo, TSpec> >& g,
                          TVertexDescriptor const state,
                          TAlphabet const symbol,
                          TEmisProb const eProb)
{
    typedef typename Size<TAlphabet>::Type TSize;
    g.data_emission[state * (TSize) ValueSize<TAlphabet>::VALUE + ordValue(symbol)] = (TCargo) eProb;
}



//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn HmmGraph#assignSilentStatus
 * @brief Assigns a silent status to a state.
 *
 * @signature void assignSilentStatus(hmm, v, silent);
 *
 * @param[in] hmm    The HmmGraph to assign the silent status in.
 * @param[in] v      Vertex descriptor of the state to assign silent status.
 * @param[in] silent A bool with the silent status.
 */

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void
assignSilentStatus(Graph<Hmm<TAlphabet, TCargo, TSpec> >& g,
                   TVertexDescriptor const vertex,
                   bool const silent)
{
    SEQAN_CHECKPOINT
    SEQAN_ASSERT(idInUse(_getVertexIdManager(g), vertex));
    g.data_silent[vertex] = silent;
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn HmmGraph#silentStatus
 * @brief Return reference to silent status flag.
 *
 * @signature TBoolRef silentStatus(hmm, v);
 *
 * @param[in] hmm    The HmmGraph to assign the silent status in.
 * @param[in] v      Vertex descriptor of the state to assign silent status.
 *
 * @return TBoolRef Reference to bool silent status flag.
 */

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline bool&
silentStatus(Graph<Hmm<TAlphabet, TCargo, TSpec> >& g,
             TVertexDescriptor const vertex)
{
    SEQAN_CHECKPOINT
    return g.data_silent[vertex];
}


//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn HmmGraph#isSilent
 * @brief Indicates whether a state is silent or not.
 *
 * @signature bool isSilent(hmm, v);
 *
 * @param[in] hmm The HmmGraph to query.
 * @param[in] v   Vertex descriptor of the state to query.
 *
 * @return bool true if the state is silent and false if not.
 */

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline bool
isSilent(Graph<Hmm<TAlphabet, TCargo, TSpec> > const& g,
         TVertexDescriptor const vertex)
{
    SEQAN_CHECKPOINT
    return g.data_silent[vertex];
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
