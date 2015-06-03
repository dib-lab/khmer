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

#ifndef SEQAN_HEADER_MISC_SVG_H
#define SEQAN_HEADER_MISC_SVG_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// SVG File Context
//////////////////////////////////////////////////////////////////////////////

struct SVGFile
{
//IOREV _nodoc_ no documentation whatsoever; uses custom IO, fully independent (maybe postpone adaption)
    std::fstream file;

    Pair<int> cursor;
    Pair<int> size;
    String<CharString> style;

    int dpSequence;
    int dpMatrix;
    int dpTrace;

    int text;

    int readForward;
    int readReverse;
    int readText;
    int rulerTextTicks;
    int rulerTextLabel;

    friend inline void svgWriteHeader(SVGFile &svg);
    friend inline void svgWriteFooter(SVGFile &svg);

    friend inline bool open(SVGFile &svg, char const * fileName);
    friend inline bool close(SVGFile &svg);

    SVGFile():
        cursor(0,0),
        size(0,0)
    {
        resetStyles();
    }

    SVGFile(char const *fileName):
        cursor(0,0),
        size(0,0)
    {
        resetStyles();
        open(*this, fileName);
    }

    ~SVGFile()
    {
        close(*this);
    }

private:

    inline void resetStyles()
    {
        clear(style);
        readText = dpSequence = length(style);
        appendValue(style, "style=\"text-anchor:middle; font-family:Bitstream Vera Sans,Verdana,sans-serif;\" font-size=\"18px\"");
        text = length(style);
        appendValue(style, "style=\"text-anchor:middle; font-family:Courier New,Verdana,sans-serif; font-weight:bold;\" font-size=\"20px\"");
        rulerTextTicks = length(style);
        appendValue(style, "style=\"text-anchor:middle; font-family:Verdana,sans-serif;\" font-size=\"6px\"");
        rulerTextLabel = length(style);
        appendValue(style, "style=\"text-anchor:left; font-family:Verdana,sans-serif; font-weight:bold;\" font-size=\"8px\"");

        dpMatrix = length(style);
        appendValue(style, "style=\"stroke:lightgray;stroke-width:1;\" marker-end=\"url(#startMarkerNormal)\"");
        dpTrace = length(style);
        appendValue(style, "style=\"stroke:black;stroke-width:2;\" marker-end=\"url(#startMarkerBold)\"");

        readForward = length(style);
        appendValue(style, "style=\"stroke:darkblue;stroke-width:3;\"");
        appendValue(style, "style=\"stroke:lightblue;stroke-width:1;\"");
        appendValue(style, "style=\"stroke:darkblue;stroke-width:3;\" marker-end=\"url(#startMarkerForward)\"");
        appendValue(style, "style=\"stroke:lightblue;stroke-width:1;\" marker-end=\"url(#startMarkerForward)\"");
        readReverse = length(style);
        appendValue(style, "style=\"stroke:darkred;stroke-width:3;\"");
        appendValue(style, "style=\"stroke:salmon;stroke-width:1;\"");
        appendValue(style, "style=\"stroke:darkred;stroke-width:3;\" marker-end=\"url(#startMarkerReverse)\"");
        appendValue(style, "style=\"stroke:salmon;stroke-width:1;\" marker-end=\"url(#startMarkerReverse)\"");
    }
};

template <>
struct DirectionIterator<SVGFile, Output>
{
    typedef SVGFile* Type;
};

inline void svgResize(SVGFile &svg, int width, int height)
{
    svg.size.i1 = width;
    svg.size.i2 = height;
}

inline void svgWriteHeader(SVGFile &svg)
{
    svg.file << "<?xml version=\"1.0\"?>" << std::endl;
    svg.file << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"" << std::endl;
    svg.file << "  \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">" << std::endl;
    svg.file << std::endl;
    svg.file << "<svg version=\"1.1\"" << std::endl;
    svg.file << "    xmlns=\"http://www.w3.org/2000/svg\"" << std::endl;
    svg.file << "    xmlns:xlink=\"http://www.w3.org/1999/xlink\"" << std::endl;
    if (svg.size.i1 != 0)
    {
        svg.file << "    width=\"" << svg.size.i1 * 20 << "px\" height=\"" << svg.size.i2 * 20 << "px\"" << std::endl;
    }
    svg.file << ">" << std::endl;
    svg.file << "<defs>" << std::endl;

    // trace back arrow markers
    svg.file << "    <g id=\"arrowMarker\" stroke-linecap=\"round\">" << std::endl;
    svg.file << "        <line x1=\"-6\" y1=\"1.5\" x2=\"0\" y2=\"0\" />" << std::endl;
    svg.file << "        <line x1=\"-6\" y1=\"-1.5\" x2=\"0\" y2=\"0\" />" << std::endl;
    svg.file << "    </g>" << std::endl;
    svg.file << "    <marker id=\"startMarkerNormal\" markerWidth=\"6\" markerHeight=\"3\" style=\"overflow:visible\" orient=\"auto\" markerUnits=\"userSpaceOnUse\">" << std::endl;
    svg.file << "        <use xlink:href=\"#arrowMarker\" stroke=\"lightgray\" />" << std::endl;
    svg.file << "    </marker>" << std::endl;
    svg.file << "    <marker id=\"startMarkerBold\" markerWidth=\"6\" markerHeight=\"4\" style=\"overflow:visible\" orient=\"auto\" markerUnits=\"userSpaceOnUse\">" << std::endl;
    svg.file << "        <use xlink:href=\"#arrowMarker\" stroke-width=\"2\" stroke=\"black\" />" << std::endl;
    svg.file << "    </marker>" << std::endl;

    // read mapping arrow markers
    svg.file << "    <marker id=\"startMarkerForward\" markerWidth=\"10\" markerHeight=\"4\" style=\"overflow:visible\" orient=\"auto\" markerUnits=\"userSpaceOnUse\">" << std::endl;
    svg.file << "        <polyline points=\"-8,0 -9,-6 3,0 -9,6 -8,0\" fill=\"darkblue\" />" << std::endl;
    svg.file << "    </marker>" << std::endl;
    svg.file << "    <marker id=\"startMarkerReverse\" markerWidth=\"10\" markerHeight=\"4\" style=\"overflow:visible\" orient=\"auto\" markerUnits=\"userSpaceOnUse\">" << std::endl;
    svg.file << "        <polyline points=\"-8,0 -9,-6 3,0 -9,6 -8,0\" fill=\"darkred\" />" << std::endl;
    svg.file << "    </marker>" << std::endl;

    svg.file << " </defs>" << std::endl;
    svg.file << std::endl;
}

inline void svgWriteFooter(SVGFile &svg)
{
    svg.file << std::endl;
    svg.file << "</svg>" << std::endl;
}


inline bool open(SVGFile &svg, char const * fileName)
{
    svg.file.open(fileName, std::ios_base::out);
    if (!svg.file.is_open()) return false;
    svgWriteHeader(svg);
    return true;
}

inline bool close(SVGFile &svg)
{
    if (svg.file.is_open())
    {
        svgWriteFooter(svg);
        svg.file.close();
    }
    return true;
}

SVGFile *
directionIterator(SVGFile & svg, Output)
{
    return &svg;
}


template <typename TCharacter>
inline void
writeValue(SVGFile *svg, TCharacter character)
{
    char x = static_cast<char>(character);
    if (x == '\n')
    {
        ++svg->cursor.i2;
        svg->cursor.i1 = 0;
    }
    else if (x == '\t')
    {
        svg->cursor.i1 = (svg->cursor.i1 & ~7) + 8;
    }
    else
    {
        if (x != ' ')
        {
            svg->file << "<g transform=\"translate(" << svg->cursor.i1*20+10 << ',' << svg->cursor.i2*20+10
                      << ")\"><text y=\"0.3em\" " << svg->style[svg->text] << '>';
            svg->file << x;
            svg->file << "</text></g>" << std::endl;
        }
        ++svg->cursor.i1;
    }
}

template <typename TContigGaps, typename TContigName>
inline void _printContig(
    SVGFile &svg,
    AlignedReadLayout &,
    TContigGaps &contigGaps,
    TContigName const &contigName)
{
    typedef typename Iterator<TContigGaps, Standard>::Type TContigIterator;

    TContigIterator cit = begin(contigGaps, Standard());
    TContigIterator citEnd = end(contigGaps, Standard());
    for (__int64 ofs = 0; cit != citEnd; ++cit, ++ofs)
    {
        if (!isGap(cit))
        {
            if (ofs == 0)
            {
                svg.file << "<g transform=\"translate(" << ofs*20+2 << ',' << svg.cursor.i2*20+10 << ")\"><text y=\"0.3em\" " << svg.style[svg.rulerTextLabel] << '>';
                svg.file << contigName << "</text></g>" << std::endl;
            }

            __int64 seqPos = cit.current.seqPos + 1;
            if (seqPos % 5 == 0)
            {
                if (seqPos % 10 == 0)
                {
                    if (ofs >= 5)
                    {
                        svg.file << "<g transform=\"translate(" << ofs*20+10 << ',' << svg.cursor.i2*20+10 << ")\"><text y=\"0.3em\" " << svg.style[svg.rulerTextTicks] << '>';
                        svg.file << seqPos << "</text></g>" << std::endl;
                    }

                    svg.file << "<line x1=\"" << ofs*20+10 << "\" x2=\"" << ofs*20+10 << "\" ";
                    svg.file << "y1=\"" << svg.cursor.i2*20+12 << "\" y2=\"" << svg.cursor.i2*20+15 << "\" ";
                } else {
                    svg.file << "<line x1=\"" << ofs*20+10 << "\" x2=\"" << ofs*20+10 << "\" ";
                    svg.file << "y1=\"" << svg.cursor.i2*20+12 << "\" y2=\"" << svg.cursor.i2*20+15 << "\" ";
                }
                svg.file << "stroke-width=\"1\" stroke=\"gray\" />" << std::endl;
            }
        }
    }
    writeValue(&svg, '\n');

    int savedStyle = svg.text;
    svg.text = svg.readText;
    svg << contigGaps;
    svg.text = savedStyle;
}

template <typename TContigGaps, typename TReadGaps, typename TAlignedRead, typename TLine>
inline void _printRead(
    SVGFile &svg,
    AlignedReadLayout &layout,
    TContigGaps &contigGaps,
    TReadGaps &readGaps,
    TAlignedRead &alignedRead,
    TLine line)
{
    typedef typename Iterator<TContigGaps, Standard>::Type TContigIterator;
    typedef typename Iterator<TReadGaps, Standard>::Type TIterator;

    __int64 xEnd = svg.cursor.i1 * 20;
    __int64 x = 0;

    int style, arrow = 0;
    const char *first;
    const char *second;
    __int64 ofs;

    if (alignedRead.beginPos < alignedRead.endPos)
    {
        ofs = alignedRead.beginPos - svg.cursor.i1 + beginPosition(readGaps);
        first = "<line x1=\"";
        second = "\" x2=\"";
        style = svg.readForward;
    } else {
        ofs = alignedRead.endPos - svg.cursor.i1 + beginPosition(readGaps);
        first = "<line x2=\"";
        second = "\" x1=\"";
        style = svg.readReverse;
    }

    if (beginPosition(readGaps) == 0)
    {
        if (alignedRead.beginPos < alignedRead.endPos)
            x = xEnd + 5;
        else {
            x = xEnd + 10;
            arrow = 2;
        }
    }
    line = svg.cursor.i2 * 20 + 10;

    if (length(layout.mateCoords) <= alignedRead.pairMatchId)
        resize(layout.mateCoords, alignedRead.pairMatchId + 1, Pair<int>(-1,-1));
    else
    {
        if (layout.mateCoords[alignedRead.pairMatchId].i2 != -1)
        {
            Pair<__int64> a((alignedRead.beginPos - ofs) * 20, line);
            Pair<__int64> b(layout.mateCoords[alignedRead.pairMatchId]);
            if (a.i1 < b.i1)
            {
                Pair<__int64> tmp = a;
                a = b;
                b = tmp;
            }
            __int64 dx = (b.i1 - a.i1);
            __int64 dy = (b.i2 - a.i2);

            svg.file << "<path d=\"M " << a.i1 << ',' << a.i2;
            svg.file << " C " << a.i1+dy/10 << ',' << a.i2-dx/10;
            svg.file << ' ' << b.i1+dy/10 << ',' << b.i2-dx/10;
            svg.file << ' ' << b.i1 << ',' << b.i2;
            svg.file << "\" stroke-width=\"2\" stroke=\"black\" stroke-opacity=\"0.2\" fill=\"none\"/>";
        }
        else
            layout.mateCoords[alignedRead.pairMatchId] = Pair<int>((alignedRead.beginPos - ofs) * 20, line);
    }


    TContigIterator cit = begin(contigGaps, Standard()) + (_min(alignedRead.beginPos, alignedRead.endPos) + beginPosition(readGaps));
    TIterator it = begin(readGaps, Standard());
    TIterator itEnd = end(readGaps, Standard());
    int lastWasGap = -1;
    int inGap;

    for (; it != itEnd; ++it, ++cit, xEnd += 20, ++svg.cursor.i1)
    {
        inGap = isGap(it);
        if (lastWasGap != inGap || inGap != static_cast<int>(isGap(cit)) || (!inGap && convert<Dna5>(*cit) != convert<Dna5>(*it)))
        {
            if (x < xEnd && lastWasGap != -1)
            {
                svg.file << first << x << "\" y1=\"" << line << second << xEnd;
                svg.file << "\" y2=\"" << line << "\" " << svg.style[style + arrow + lastWasGap] << " />" << std::endl;
                arrow = 0;
            }
            lastWasGap = inGap;
            x = xEnd;
            if (!inGap && convert<Dna5>(*cit) != convert<Dna5>(*it))
            {
                svg.file << "<g transform=\"translate(" << xEnd + 10 << ',' << line << ")\"><text y=\"0.3em\" " << svg.style[svg.readText] << '>';
                writeValue(&svg, convert<char>(*it));
                svg.file << "</text></g>" << std::endl;
                x += 20;
                arrow = 0;
            }
        }
    }
    if (x < xEnd && lastWasGap != -1)
    {
        if (static_cast<int>(_unclippedLength(readGaps)) == static_cast<int>(endPosition(readGaps)))
        {
            if (alignedRead.beginPos < alignedRead.endPos)
            {
                arrow = 2;
                xEnd -= 10;
            }
            else
            {
                xEnd -= 5;
            }
        }
        svg.file << first << x << "\" y1=\"" << line << second << xEnd;
        svg.file << "\" y2=\"" << line << "\" " << svg.style[style + arrow + lastWasGap] << " />" << std::endl;
    }
}


//////////////////////////////////////////////////////////////////////////////
// stream operators
//////////////////////////////////////////////////////////////////////////////
/*
template <typename TSource>
inline SVGFile &
operator << (SVGFile & target,
             TSource  source)
{
    typename DirectionIterator<TStream, Output>::Type it = directionIterator(target, Output());
    write(it, source);
    return target;
}
*/

//inline SVGFile &
//operator << (SVGFile & target, char source)
//{
//  typename DirectionIterator<SVGFile, Output>::Type it = directionIterator(target, Output());
//  write(it, source);
//    return target;
//}
//
//inline SVGFile &
//operator << (SVGFile & target, char const *source)
//{
//  typename DirectionIterator<SVGFile, Output>::Type it = directionIterator(target, Output());
//  write(it, source);
//    return target;
//}


template <typename TStringSet, typename TTrace, typename TIndexPair>
void
_alignNeedlemanWunschMatrix(SVGFile& svg,
                              TStringSet const& str,
                              TTrace const& trace,
                              TIndexPair const&)
{
SEQAN_CHECKPOINT
    typedef typename Size<TStringSet>::Type TSize;
    //typedef typename Value<TTrace>::Type TTraceValue;

    // TraceBack values
    const int Diagonal = 0;
    const int Horizontal = 1;
    const int Vertical = 2;

    // Initialization
    TSize numCols = length(str[0]);
    TSize numRows = length(str[1]);

    // Print trace matrix
    for (TSize pos0 = 0; pos0 < numCols; ++pos0)
    {
        for (TSize pos1 = 0; pos1 < numRows; ++pos1)
        {
            int tv =(int)trace[pos0*numRows + pos1];
                if (tv & (1 << Diagonal))
                    svg.file << "<line x2=\"" << pos0*20+10 << "\" y2=\"" << pos1*20+10 << "\"" << " x1=\"" << (pos0+1)*20+10 << "\" y1=\"" << (pos1+1)*20+10 << "\" " << svg.style[svg.dpMatrix] << " />" << std::endl;

                if (tv & (1 << Horizontal))
                    svg.file << "<line x2=\"" << pos0*20+10 << "\" y2=\"" << (pos1+1)*20+10 << "\"" << " x1=\"" << (pos0+1)*20+10 << "\" y1=\"" << (pos1+1)*20+10 << "\" " << svg.style[svg.dpMatrix] << " />" << std::endl;

                if (tv & (1 << Vertical))
                    svg.file << "<line x2=\"" << (pos0+1)*20+10 << "\" y2=\"" << pos1*20+10 << "\"" << " x1=\"" << (pos0+1)*20+10 << "\" y1=\"" << (pos1+1)*20+10 << "\" " << svg.style[svg.dpMatrix] << " />" << std::endl;
        }
    }

    // Print sequences
    for (TSize pos0 = 0; pos0 < numCols; ++pos0)
        svg.file << "<g transform=\"translate(" << pos0*20+20 << ",5)\"><text y=\"0.3em\" " << svg.style[svg.dpSequence] << '>' << str[0][pos0] << "</text></g>" << std::endl;

    for (TSize pos1 = 0; pos1 < numRows; ++pos1)
        svg.file << "<g transform=\"translate(0," << pos1*20+25 << ")\"><text y=\"0.3em\" " << svg.style[svg.dpSequence] << '>' << str[1][pos1] << "</text></g>" << std::endl;
}

template <typename TStringSet, typename TId, typename TPos, typename TTraceValue>
inline void
_alignTracePrint(SVGFile& svg,
                   TStringSet const&,
                   TId const,
                   TPos pos1,
                   TId const,
                   TPos pos2,
                   TPos const segLen,
                   TTraceValue const tv)
{
    // TraceBack values
    TTraceValue Diagonal = 0; TTraceValue Horizontal = 1; TTraceValue Vertical = 2;

    svg.file << "<line x2=\"" << pos1*20+10 << "\" y2=\"" << pos2*20+10 << "\"";
    if (tv == Diagonal) {
        pos1 += segLen; pos2 += segLen;
    } else if (tv == Horizontal) {
        pos1 += segLen;
    } else if (tv == Vertical) {
        pos2 += segLen;
    }
    svg.file << " x1=\"" << pos1*20+10 << "\" y1=\"" << pos2*20+10;
    svg.file << "\" " << svg.style[svg.dpTrace] << " />" << std::endl;
}

//////////////////////////////////////////////////////////////////////////////



} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
