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

#ifndef SEQAN_HEADER_INDEX_ESA_DRAWING_H
#define SEQAN_HEADER_INDEX_ESA_DRAWING_H

namespace SEQAN_NAMESPACE_MAIN
{

template <typename TFile, typename TText, typename TESASpec>
void writeRecords(
    TFile & file,
    Index<TText, IndexEsa<TESASpec> > & stree,
    DotDrawing)
{
//IOREV _nodoc_
    SEQAN_CHECKPOINT
    typedef Index<TText, IndexEsa<TESASpec> > TIndex;

    typename DirectionIterator<TFile, Output>::Type iter = directionIterator(file, Output());

    write(iter, "digraph G {\n");
    writeValue(iter, '\n');
    write(iter, "/* Graph Attributes */\n");
    write(iter, "graph [rankdir = LR];\n");
    writeValue(iter, '\n');
    write(iter, "/* Node Attributes */\n");
    write(iter, "node [shape = ellipse, fillcolor = lightgrey, style = filled, fontname = \"Times-Italic\"];\n");
    writeValue(iter, '\n');
    write(iter, "/* Edge Attributes */\n");
    write(iter, "edge [fontname = \"Times-Italic\", arrowsize = 0.75, fontsize = 16];\n");
    writeValue(iter, '\n');

    write(iter, "/* Edges */\n");
    typedef typename Iterator<TIndex, TopDown<ParentLinks<Preorder> > >::Type TIterator;
    typedef typename Iterator<TIndex, TopDown<> >::Type TIteratorSimple;
    TIterator it(stree);

    for(;!atEnd(it);++it)
    {
        // dump node
        write(iter, "\"[");
         appendNumber(iter, range(it).i1);
        writeValue(iter, ':');
        appendNumber(iter, range(it).i2);
        write(iter, ")\"");
        if (!isRightTerminal(it))
            write(iter, " [style = dashed]");
        write(iter, ";\n");

        // dump edge from parent (if not root)
        if (!isRoot(it))
        {
            TIteratorSimple src(container(it), nodeUp(it));

            write(iter, "\"[");
            appendNumber(iter, range(src).i1);
            writeValue(iter, ':');
            appendNumber(iter, range(src).i2);
            write(iter, ")\"");

            write(iter, " -> ");

            write(iter, "\"[");
            appendNumber(iter, range(it).i1);
            writeValue(iter, ':');
            appendNumber(iter, range(it).i2);
            write(iter, ")\"");

            write(iter, " [label = \"");
            write(iter, parentEdgeLabel(it));
            write(iter, "\"];\n");
        }
    }
    writeValue(iter, '\n');

    write(iter, "}\n");
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
