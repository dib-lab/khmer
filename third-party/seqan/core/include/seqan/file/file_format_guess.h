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

#ifndef SEQAN_HEADER_FILE_GUESS_H
#define SEQAN_HEADER_FILE_GUESS_H

/* IOREV
 * _tested_
 * _doc_
 *
 * only works if seeking is possible
 * 
 * the whole FileFormat class seems useless, see file_format.h
 * this could aswell return the specific tag
 *
 * description in code comments not plausible
 */

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// guessFileFormat
//////////////////////////////////////////////////////////////////////////////

//guessFileFormat braucht auch data, weil die FileFormat-Klasse von TData
//abhaengig, und das ist so, weil sonst die Kombination von Templates mit
//virtuellen Funktionen nicht funktionieren wuerde.
/**
.Function.guessFileFormat
..cat:Input/Output
..summary:Tries to determine the format of a file.
..signature:guessFileFormat(file, data)
..param.file: An input file.
..param.data: The target container.
...remarks:This container is not modified by this function.
..returns:A file format object instance that represents the determined file format.
...type:Class.FileFormat
..remarks:The $data$-argument is used here as a tag to determine the type of the target.
..see:Function.FileFormat#read
..see:Tag.File Format
..include:seqan/file.h
*/
template <typename TFile, typename TData, typename TMeta>
inline FileFormat<TFile, TData, TMeta, void> 
guessFileFormat(TFile & file,
				TData & data)
{
//IOREV _doc_ see head of file_format_guess.h
SEQAN_CHECKPOINT
	typename Position<TFile>::Type old_pos = _streamTellG(file);
	typename Value<TFile>::Type c;

	_streamSeekG(file, 0); /// move to beginning of file
	c = _streamGet(file);
		
	if (c=='>') 
	{
		_streamSeekG(file, old_pos);
		return getFileFormatInstance<TFile, TData, Fasta, TMeta>();
	}
	
	if (c=='L')
	{
		_streamSeekG(file, old_pos);
		return getFileFormatInstance<TFile, TData, Genbank, TMeta>();
	}

	if (c=='I')
	{
		_streamSeekG(file, old_pos);
		return getFileFormatInstance<TFile, TData, Embl, TMeta>();
	}

	else
	{
		_streamSeekG(file, old_pos);
		return getFileFormatInstance<TFile, TData, Raw, TMeta>();
	}
}

//////////////////////////////////////////////////////////////////////////////

/* DOCH NICHT:
template <typename TTarget, typename TSource>
inline void
read(TTarget & target,
	 TSource & source)
{
SEQAN_CHECKPOINT
	read(target, source, guessFileFormat(target, source));
}
*/

//////////////////////////////////////////////////////////////////////////////
} //namespace SEQAN_NAMESPACE_MAIN

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef SEQAN_HEADER_...
