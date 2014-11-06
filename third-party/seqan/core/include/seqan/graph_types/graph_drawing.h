// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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


#ifndef SEQAN_HEADER_GRAPH_DRAWING_H
#define SEQAN_HEADER_GRAPH_DRAWING_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Graph - Drawing
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// WRITING
//////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////
// INTERNAL FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template <typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor, typename TAttributes>
inline void 
_markRootVertex(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g,
				TVertexDescriptor const& v,
				TAttributes& str)
{
	SEQAN_CHECKPOINT
	if (isRoot(g,v)) {
		append(str, ", shape = doublecircle");
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec, typename TVertexDescriptor, typename TAttributes>
inline void 
_markRootVertex(Graph<Directed<TCargo, TSpec> > const&,
				TVertexDescriptor const&,
				TAttributes&)
{
	SEQAN_CHECKPOINT
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec, typename TVertexDescriptor, typename TAttributes>
inline void 
_markRootVertex(Graph<Undirected<TCargo, TSpec> > const&,
				TVertexDescriptor const&,
				TAttributes&)
{
	SEQAN_CHECKPOINT
}


//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec, typename TVertexDescriptor, typename TAttributes>
inline void 
_markRootVertex(Graph<Tree<TCargo, TSpec> > const& g,
				TVertexDescriptor const& v,
				TAttributes& str)
{
	SEQAN_CHECKPOINT
	if (isRoot(g,v)) {
		append(str, ", shape = doublecircle");
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TPosition, typename TNodeMap>
inline void
_createTrieNodeAttributes(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g,
						  String<String<TPosition> > pos,
						  TNodeMap& nodeMap)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	resizeVertexMap(g, nodeMap);
	typedef typename Iterator<TGraph, VertexIterator>::Type TConstIter;
	TConstIter it(g);
	for(;!atEnd(it);++it) {
		String<char> tmp;
		std::stringstream s;
		s << *it;
		String<TPosition> endPositions = getProperty(pos,*it);
		if (!empty(endPositions)) {
			s <<  " {";
			append(tmp, "shape = box, ");
			typename Iterator<String<TPosition>, Rooted>::Type itP = begin(endPositions);
			typename Iterator<String<TPosition>, Rooted>::Type beginP = itP;
			for(;!atEnd(itP);goNext(itP)) {
				if (beginP != itP) s << ", ";
				s << *itP;
			}
			s << "}";
		}
		
		append(tmp, "label = \"");
		append(tmp, s.str().c_str());
		append(tmp, "\"");
		_markRootVertex(g, *it, tmp);
		assignProperty(nodeMap, *it, tmp);
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TNodeAttributes>
inline void
_createNodeAttributes(Graph<TSpec> const& g,
					  TNodeAttributes& nodeMap)
{
    typedef Graph<TSpec> TGraph;
	resizeVertexMap(g, nodeMap);

	typedef typename Iterator<TGraph, VertexIterator>::Type TConstIter;
	TConstIter it(g);
	for(;!atEnd(it);++it) {
		std::ostringstream outs;
		outs << "label = \"";
		outs << *it;
		outs << "\"";
		String<char> tmp;
		append(tmp, outs.str().c_str());
		_markRootVertex(g, *it, tmp);
		assignProperty(nodeMap, *it, tmp);
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TNodeAttributes, typename TNameMap>
inline void
_createNodeAttributes(Graph<TSpec> const& g,
					  TNodeAttributes& nodeMap,
					  TNameMap const& nameMap)
{
    typedef Graph<TSpec> TGraph;
	resizeVertexMap(g, nodeMap);

	typedef typename Iterator<TGraph, VertexIterator>::Type TConstIter;
	TConstIter it(g);
	for(;!atEnd(it);++it) {
		std::ostringstream outs; 
		outs << "label = \"";
		outs << getProperty(nameMap,*it);
		outs << "\"";
		String<char> tmp;
		append(tmp, outs.str().c_str());
		_markRootVertex(g, *it, tmp);
		assignProperty(nodeMap, *it, tmp);
	}
}

//////////////////////////////////////////////////////////////////////////////
template<typename TSpec, typename TEdgeAttributes>
inline void
_createEmptyEdgeAttributes(Graph<TSpec> const& g,
						   TEdgeAttributes& edgeMap)
{
	typedef Graph<TSpec> TGraph;
	resizeEdgeMap(g, edgeMap);

	typedef typename Iterator<TGraph, EdgeIterator>::Type TConstEdIter;
	TConstEdIter itEd(g);
	for(;!atEnd(itEd);++itEd) {
		assignProperty(edgeMap, *itEd, String<char>(""));
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec, typename TEdgeAttributes>
inline void
_createEdgeAttributes(Graph<Directed<TCargo, TSpec> > const& g,
					  TEdgeAttributes& edgeMap)
{
	_createEmptyEdgeAttributes(g,edgeMap);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec, typename TEdgeAttributes>
inline void
_createEdgeAttributes(Graph<Undirected<TCargo, TSpec> > const& g,
					  TEdgeAttributes& edgeMap)
{
	SEQAN_CHECKPOINT
	_createEmptyEdgeAttributes(g,edgeMap);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSpec, typename TEdgeAttributes>
inline void
_createEdgeAttributes(Graph<Tree<void, TSpec> > const& g,
					  TEdgeAttributes& edgeMap)
{
	SEQAN_CHECKPOINT
	_createEmptyEdgeAttributes(g,edgeMap);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec, typename TEdgeAttributes>
inline void
_createEdgeAttributes(Graph<Tree<TCargo, TSpec> > const& g,
					  TEdgeAttributes& edgeMap)
{
	typedef Graph<Tree<TCargo, TSpec> > TGraph;
	resizeEdgeMap(g, edgeMap);

	typedef typename Iterator<TGraph, EdgeIterator>::Type TConstEdIter;
	TConstEdIter itEd(g);
	for(;!atEnd(itEd);++itEd) {
		std::ostringstream outs; 
		outs << "label = \"";
		outs << (TCargo) getCargo(*itEd);
		outs << "\"";
		append(property(edgeMap, *itEd), outs.str().c_str());
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TAlphabet, typename TCargo, typename TSpec, typename TEdgeAttributes>
inline void
_createEdgeAttributes(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g,
					  TEdgeAttributes& edgeMap)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	resizeEdgeMap(g, edgeMap);

	typedef typename Iterator<TGraph, EdgeIterator>::Type TConstEdIter;
	TConstEdIter itEd(g);
	for(;!atEnd(itEd);++itEd) {
		String<char> tmp("label = \"");
		append(tmp, label(itEd));
		append(tmp, "\"");
		assignProperty(edgeMap, *itEd, tmp);
	}
}


//////////////////////////////////////////////////////////////////////////////

template <typename TAlphabet, typename TCargo, typename TSpec, typename TEdgeAttributes>
inline void
_createEdgeAttributes(Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> > > const& g,
					  TEdgeAttributes& edgeMap)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> > > TGraph;
	resizeEdgeMap(g, edgeMap);

	typedef typename Iterator<TGraph, EdgeIterator>::Type TConstEdIter;
	TConstEdIter itEd(g);
	for(;!atEnd(itEd);++itEd) {
		String<TAlphabet> labelTmp = getCargo(*itEd);
		String<char> str;
		resize(str,length(labelTmp)+1);
		value(str,0) = label(itEd);
		typename Iterator<String<TAlphabet>, Rooted>::Type it = begin(labelTmp);
		for(;!atEnd(it);++it) {
			char c = convert<char>(getValue(it));
			value(str,position(it) + 1) = c;
		}
		String<char> tmp("label = \"");
		append(tmp, str);
		append(tmp, "\"");
		assignProperty(edgeMap, *itEd, tmp);
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TCargo, typename TSpec>
inline void
_writeGraphFooter(TFile &,
				  Graph<Directed<TCargo, TSpec> > const&,
				  DotDrawing)
{
//IOREV
	SEQAN_CHECKPOINT
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TCargo, typename TSpec>
inline void
_writeGraphFooter(TFile &,
				  Graph<Undirected<TCargo, TSpec> > const&,
				  DotDrawing)
{
//IOREV
	SEQAN_CHECKPOINT
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TCargo, typename TSpec>
inline void
_writeGraphFooter(TFile &,
				  Graph<Tree<TCargo, TSpec> > const&,
				  DotDrawing)
{
//IOREV
	SEQAN_CHECKPOINT
}


//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TAlphabet, typename TCargo, typename TSpec>
inline void
_writeGraphFooter(TFile &,
				  Graph<Automaton<TAlphabet, TCargo, TSpec> > const&,
				  DotDrawing)
{
//IOREV
	SEQAN_CHECKPOINT
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TAlphabet, typename TCargo, typename TSpec>
inline void
_writeGraphType(TFile & file,
				Graph<Automaton<TAlphabet, TCargo, TSpec> > const&,
				DotDrawing)
{
//IOREV
	streamPut(file, "digraph");
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TCargo, typename TSpec>
inline void
_writeGraphType(TFile & file,
				Graph<Directed<TCargo, TSpec> > const&,
				DotDrawing)
{
//IOREV
	streamPut(file, "digraph");
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TCargo, typename TSpec>
inline void
_writeGraphType(TFile & file,
				Graph<Undirected<TCargo, TSpec> > const&,
				DotDrawing)
{
//IOREV
	streamPut(file, "graph");
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TCargo, typename TSpec>
inline void
_writeGraphType(TFile & file,
				Graph<Tree<TCargo, TSpec> > const&,
				DotDrawing)
{
//IOREV
	streamPut(file, "digraph");
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TAlphabet, typename TCargo, typename TSpec>
inline void
_writeEdgeType(TFile & file,
			   Graph<Automaton<TAlphabet, TCargo, TSpec> > const&,
			   DotDrawing)
{
//IOREV
	streamPut(file, " -> ");
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TCargo, typename TSpec>
inline void
_writeEdgeType(TFile & file,
			   Graph<Directed<TCargo, TSpec> > const&,
			   DotDrawing)
{
//IOREV
	streamPut(file, " -> ");
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TCargo, typename TSpec>
inline void
_writeEdgeType(TFile & file,
			   Graph<Undirected<TCargo, TSpec> > const&,
			   DotDrawing)
{
//IOREV
	streamPut(file, " -- ");
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TCargo, typename TSpec>
inline void
_writeEdgeType(TFile & file,
			   Graph<Tree<TCargo, TSpec> > const&,
			   DotDrawing)
{
//IOREV
	streamPut(file, " -> ");
}

//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

/**
.Function.write
..class:Class.Graph
..signature:write(file, graph, nodeMap, edgeMap, tag)
..param.graph:The graph to write out.
...type:Class.Graph
..param.nodeMap:A mapping from vertex descriptor to vertex label.
..param.edgeMap:A mapping from edge descriptor to edge label.
..param.tag:A tag to select the output format.
...type:Tag.DotDrawing
..include:seqan/graph_types.h
 */
template <typename TFile, typename TSpec, typename TNodeAttributes, typename TEdgeAttributes>
void
write(TFile & file,
	  Graph<TSpec> const& g,
	  TNodeAttributes const& nodeMap,
	  TEdgeAttributes const& edgeMap,
	  DotDrawing)
{
//IOREV _doc_ _batchreading_
	typedef Graph<TSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;

	_writeGraphType(file,g,DotDrawing());
	streamPut(file, " G {\n");
	streamPut(file, '\n');
	streamPut(file, "/* Graph Attributes */\n");
	streamPut(file, "graph [rankdir = LR];\n");
	streamPut(file, '\n');
	streamPut(file, "/* Node Attributes */\n");
	streamPut(file, "node [shape = rectangle, fillcolor = white, style = filled, fontname = \"Times-Italic\"];\n");
	streamPut(file, '\n');
	streamPut(file, "/* Edge Attributes */\n");
	streamPut(file, "edge [fontname = \"Times-Italic\", arrowsize = 0.75, fontsize = 16];\n");
	streamPut(file, '\n');

	streamPut(file, "/* Nodes */\n");
	typedef typename Iterator<TGraph, VertexIterator>::Type TConstIter;
	TConstIter it(g);
	for(;!atEnd(it);++it) {
		streamPut(file, (int)*it);
		streamPut(file, " [");
		streamPut(file, getProperty(nodeMap, *it));
		streamPut(file, "];\n");
	}
	streamPut(file, '\n');

	streamPut(file, "/* Edges */\n");
	typedef typename Iterator<TGraph, EdgeIterator>::Type TConstEdIter;
	TConstEdIter itEd(g);
	for(;!atEnd(itEd);++itEd) {
		TVertexDescriptor sc = sourceVertex(itEd);
		TVertexDescriptor tr = targetVertex(itEd);
		streamPut(file, (int)sc);
		_writeEdgeType(file, g, DotDrawing());
		streamPut(file, (int)tr);
		streamPut(file, " [");
		streamPut(file, getProperty(edgeMap, *itEd));
		streamPut(file, "];\n");
	}
	streamPut(file, '\n');

	_writeGraphFooter(file,g,DotDrawing());

	streamPut(file, "}\n");
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.write
..signature:write(file, graph, nodeMap, tag)
..include:seqan/graph_types.h
 */
template <typename TFile, typename TSpec, typename TNodeAttributes>
inline void
write(TFile & file,
	  Graph<TSpec> const& g, 
	  TNodeAttributes const& nodeMap,
	  DotDrawing) 
{
//IOREV _doc_ _batchreading_
	String<String<char> > edgeMap;
	_createEdgeAttributes(g,edgeMap);
	write(file,g,nodeMap,edgeMap,DotDrawing());
}


//////////////////////////////////////////////////////////////////////////////

/**
.Function.write
..signature:write(file, graph, tag)
..include:seqan/graph_types.h
 */
template <typename TFile, typename TSpec>
inline void
write(TFile & file,
	  Graph<TSpec> const& g, 
	  DotDrawing) 
{
//IOREV _doc_ _batchreading_
	String<String<char> > nodeMap;
	_createNodeAttributes(g,nodeMap);
	String<String<char> > edgeMap;
	_createEdgeAttributes(g,edgeMap);
	write(file,g,nodeMap,edgeMap,DotDrawing());
}





//////////////////////////////////////////////////////////////////////////////
// READING
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TStatement, typename TNodeAttributes, typename TEdgeAttributes, typename TNodeIdMap>
inline void
_addNode(Graph<TSpec>& g,
		 TStatement& node_id,
		 TStatement& attr_list,
		 TNodeAttributes& nodeMap,
		 TEdgeAttributes&,			  
		 TNodeIdMap& nodeIdMap)
{
	typedef Graph<TSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;

	if (nodeIdMap.find(node_id) == nodeIdMap.end()) {
		TVertexDescriptor _id = addVertex(g);
		nodeIdMap.insert(std::make_pair(node_id, _id));
		resizeVertexMap(g, nodeMap);
		assignProperty(nodeMap, _id, attr_list);
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TVertexDescriptor, typename TNodeAttributes, typename TEdgeAttributes, typename TStatement>
inline void
_addEdge(Graph<Directed<TCargo, TSpec> >& g,
		 TVertexDescriptor sourceV,
		 TVertexDescriptor targetV,
		 TNodeAttributes&,
		 TEdgeAttributes& edgeMap,
		 TStatement& attr_list)
{
	typedef Graph<Directed<TCargo, TSpec> > TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	TEdgeDescriptor e = addEdge(g, sourceV, targetV);
	resizeEdgeMap(g, edgeMap);
	assignProperty(edgeMap, e, attr_list);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TVertexDescriptor, typename TNodeAttributes, typename TEdgeAttributes, typename TStatement>
inline void
_addEdge(Graph<Undirected<TCargo, TSpec> >& g,
		 TVertexDescriptor sourceV,
		 TVertexDescriptor targetV,
		 TNodeAttributes&,
		 TEdgeAttributes& edgeMap,
		 TStatement& attr_list)
{
	typedef Graph<Undirected<TCargo, TSpec> > TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	TEdgeDescriptor e = addEdge(g, sourceV, targetV);
	resizeEdgeMap(g, edgeMap);
	assignProperty(edgeMap, e, attr_list);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TVertexDescriptor, typename TNodeAttributes, typename TEdgeAttributes, typename TStatement>
inline void
_addEdge(Graph<Tree<TCargo, TSpec> >& g,
		 TVertexDescriptor sourceV,
		 TVertexDescriptor targetV,
		 TNodeAttributes&,
		 TEdgeAttributes& edgeMap,
		 TStatement& attr_list)
{
	typedef Graph<Tree<TCargo, TSpec> > TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	TEdgeDescriptor e = addEdge(g, sourceV, targetV);
	resizeEdgeMap(g, edgeMap);
	assignProperty(edgeMap, e, attr_list);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TString>
inline typename Alphabet<Graph<Automaton<TAlphabet, TCargo, TSpec> > >::Type
_getInternalLabel(Graph<Automaton<TAlphabet, TCargo, TSpec> >&,
				  TString& str)
{
	return str[0];
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TString>
inline String<TAlphabet>
_getInternalLabel(Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> > >&,
				  TString& str)
{
	return str;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor, typename TNodeAttributes, typename TEdgeAttributes, typename TStatement>
inline void
_addEdge(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g,
		 TVertexDescriptor sourceV,
		 TVertexDescriptor targetV,
		 TNodeAttributes&,
		 TEdgeAttributes& edgeMap,
		 TStatement& attr_list)
{
	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;

	// We need the label
	typedef typename Value<TStatement>::Type TValue;
	typedef typename Iterator<TStatement>::Type TIter;
	typedef typename Position<TIter>::Type TPos;
	
	String<TValue> label;
	TIter it = begin(attr_list);
	bool found = false;
	for(;!atEnd(it);goNext(it)) {
		TPos pos = position(it);
		if (*it == ',') {
			found = false;
		} else if (found) {
			append(label, *it);
		} else if ((pos + 5 < length(attr_list)) &&
			(infix(attr_list, it, it + 5) == "label")) 
		{
				found = true;
				it += 5;
		}
	}
	TEdgeDescriptor e = addEdge(g, sourceV, targetV, _getInternalLabel(g, label));
	resizeEdgeMap(g, edgeMap);
	assignProperty(edgeMap, e, attr_list);
}


//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TStatement, typename TNodeAttributes, typename TEdgeAttributes, typename TNodeIdMap>
inline void
_addEdge(Graph<TSpec>& g,
		 TStatement& left_node_id,
		 TStatement& right_node_id,
		 TStatement& attr_list,
		 TNodeAttributes& nodeMap,
		 TEdgeAttributes& edgeMap,
		 TNodeIdMap& nodeIdMap)
{
	typedef Graph<TSpec> TGraph;
	typedef typename Value<TStatement>::Type TValue;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef std::map<String<TValue>, TVertexDescriptor> TMap;

	TVertexDescriptor sourceV;
	TVertexDescriptor targetV;

	typename TMap::iterator pos;
	pos = nodeIdMap.find(left_node_id);
	if (pos == nodeIdMap.end()) return;
	else sourceV = pos->second;

	pos = nodeIdMap.find(right_node_id);
	if (pos == nodeIdMap.end()) return;
	else targetV = pos->second;

	_addEdge(g, sourceV, targetV, nodeMap, edgeMap, attr_list);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TStatement, typename TNodeAttributes, typename TEdgeAttributes, typename TNodeIdMap>
inline void
_processNodeStatement(Graph<TSpec>& g,
					  TStatement& stmt,
					  TNodeAttributes& nodeMap,
					  TEdgeAttributes& edgeMap,
					  TNodeIdMap& nodeIdMap) 
{
	typedef typename Value<TStatement>::Type TValue;
	typedef typename Iterator<TStatement>::Type TIter;
	
	String<TValue> node_id;
	String<TValue> attr_list;  // Multiple attribute lists are ignored
	bool inAttr = false;
	TIter it = begin(stmt);
	for(;!atEnd(it);goNext(it)) {
		if (*it == '[') {
			inAttr = true;
			continue;
		} else if (*it == ']') {
			// Finished
			break;
		} else if ((*it == ' ') ||
			(*it == '"')) {
			continue;
		}
		if (inAttr) {
			append(attr_list, *it);
		} else {
			append(node_id, *it);
		}
	}
	_addNode(g, node_id, attr_list, nodeMap, edgeMap, nodeIdMap);
}


//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TStatement, typename TNodeAttributes, typename TEdgeAttributes, typename TPosition, typename TNodeIdMap>
inline void
_processEdgeStatement(Graph<TSpec>& g,
					  TStatement& stmt,
					  TNodeAttributes& nodeMap,
					  TEdgeAttributes& edgeMap,
					  TPosition pos,
					  TNodeIdMap& nodeIdMap) 
{
	typedef typename Value<TStatement>::Type TValue;
	typedef typename Iterator<TStatement>::Type TIter;
	
	String<TValue> left_node_id;
	String<TValue> right_node_id;
	String<TValue> attr_list;  // Multiple attribute lists are ignored
	bool inAttr = false;
	TIter it = begin(stmt);
	unsigned int localPos = 0;
	for(;!atEnd(it);goNext(it), ++localPos) {
		if (*it == '[') {
			inAttr = true;
			continue;
		} else if (*it == ']') {
			// Finished
			break;
		} else if ((*it == ' ') ||
			(*it == '"')) {
			continue;
		}
		if (inAttr) {
			append(attr_list, *it);
		} else if (localPos < pos) {
			append(left_node_id, *it);
		} else if (localPos > pos+1) {
			append(right_node_id, *it);
		}
	}
	//std::cout << left_node_id << "," << right_node_id << "," << std::endl;
	_addEdge(g, left_node_id, right_node_id, attr_list, nodeMap, edgeMap, nodeIdMap);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TStatement, typename TNodeAttributes, typename TEdgeAttributes, typename TNodeIdMap>
inline void
_processStatement(Graph<TSpec>& g,
				  TStatement& stmt,
				  TNodeAttributes& nodeMap,
				  TEdgeAttributes& edgeMap,
				  TNodeIdMap& nodeIdMap) 
{
	// Clear everything up to the last line
	typedef typename Value<TStatement>::Type TValue;
	typedef typename Iterator<TStatement>::Type TIter;

	// Exclude header and empty lines
	TIter it = begin(stmt);
	String<TValue> _id;
	for(;!atEnd(it);goNext(it)) {
	  if ((*it != '\t') && (*it != ' ') && (*it != '\n') && (*it != '\r')) {
	    append(_id, *it);
	  } else {
	    // Exclude any graph, subgraph, node and edge processing attributes
	    if ((_id == "graph") || (_id == "node") || (_id == "edge") || (_id == "subgraph") || (length(_id)<1)) {
	      clear(stmt);
	      return;
	    } else break; 
	  }
	}

	// Process Edges
	it = begin(stmt);
	clear(_id);
	_id = "00";
	unsigned int pos = 0;
	for(;!atEnd(it);goNext(it), ++pos) {
	  _id[pos % 2] = *it;
	  if ((_id == "--") || (_id == "->")) {
	    //std::cout << stmt << std::endl;
	    _processEdgeStatement(g, stmt, nodeMap, edgeMap, pos - 1, nodeIdMap);
	    clear(stmt);
	    return;
	  }
	}

	// Process nodes
	//std::cout << stmt << std::endl;
	_processNodeStatement(g, stmt, nodeMap, edgeMap, nodeIdMap);
	clear(stmt);
}


//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TSpec, typename TNodeAttributes, typename TEdgeAttributes>
void read(TFile & file,
		  Graph<TSpec>& g,
		  TNodeAttributes& nodeMap,
		  TEdgeAttributes& edgeMap,
		  DotDrawing)
{
    // TODO(holtgrew): Could be adapted to use RecordReader for parsing.
	typedef Graph<TSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Value<TFile>::Type TValue;
	typedef std::map<String<TValue>, TVertexDescriptor> TMap;
	TMap nodeIdMap;

	TValue c;
	String<TValue> stmt;
	while (!streamEOF(file)) {
        streamReadChar(c, file);

		if (c == ';') _processStatement(g,stmt, nodeMap, edgeMap, nodeIdMap);
		else if ((c == '\n') ||
				(c == '\r')) {
					clear(stmt);
		}
		else append(stmt,c);
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TSpec>
void read(TFile & file,
		  Graph<TSpec>& g,
		  DotDrawing) 
{
//IOREV _batchreading_ 
	String<String<char> > nodeMap;
	String<String<char> > edgeMap;
	read(file,g,nodeMap,edgeMap,DotDrawing());
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
