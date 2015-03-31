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

#ifndef SEQAN_HEADER_GRAPH_UTILITY_PARSING_H
#define SEQAN_HEADER_GRAPH_UTILITY_PARSING_H


#include <fstream>

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// File reading
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

// This is unused and a deletion candidate.
template<typename TPath, typename TStringSet, typename TNames>
inline unsigned int
_loadSequences(TPath const& in_path,
			   TStringSet& origStrSet,
			   TNames& names)
{
//IOREV _nodoc_ uses custom IO, with ifstream, should be adapted to File()
	typedef typename Size<TStringSet>::Type TSize;
	
	// Count sequences and read names
	TSize seqCount = 0;
	std::ifstream file;
	std::stringstream input;
	input << in_path;
	file.open(input.str().c_str(), std::ios_base::in | std::ios_base::binary);
	if (!file.is_open()) return 0;
	while (!file.eof()) {
		String<char> id;
		readID(file, id, Fasta());
		appendValue(names, id);
		goNext(file, Fasta());
		++seqCount;
	}

	// Load sequences
	file.clear();
	file.seekg(0, std::ios_base::beg);
	resize(origStrSet, seqCount);
	TSize count = 0;
	for(TSize i = 0; (i < seqCount) && !file.eof(); ++i) 	{
		read(file, origStrSet[i], Fasta());
		count += length(origStrSet[i]);
	}
    file.close();

	return count;
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
