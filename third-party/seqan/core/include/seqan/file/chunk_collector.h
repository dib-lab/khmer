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

#ifndef SEQAN_HEADER_CHUNK_COLLECTOR_H
#define SEQAN_HEADER_CHUNK_COLLECTOR_H

/* IOREV
 * _tested_
 * _nodoc_
 * 
 * this class is only used file_format_raw
 * it is tested by test_file which tests file_format_raw
 * it is sparsely documented, but marked as Internal anyway (does it need doku
 * then?)
 * TODO not fully understood implications of assign and replace
 */



namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////

/**
.Internal.ChunkCollector_:
..cat:Classes
..summary:Reads piecewise from stream, collects pieces (chunks) in a vector.
..signature:ChunkCollector_<Host>
..param.Host:Type of host object that is used as allocator.
*/

//////////////////////////////////////////////////////////////////////////////

template <typename T>
struct ChunkLength
{
//IOREV
	enum { VALUE = 1024 };
};

//////////////////////////////////////////////////////////////////////////////
// StreamChunkCollector_ class: collects content of a stream in chunks
//////////////////////////////////////////////////////////////////////////////

template <typename THost>
class ChunkCollector_
{
//IOREV
public:
	THost * data_host;
	typename Size<THost>::Type data_length;

	typedef ::std::vector<typename Value<THost>::Type *, ToStdAllocator<THost, typename Value<THost>::Type *> > Chunk_Holder;
	Chunk_Holder data_chunks; 

public:
	static int const CHUNK_LENGTH = ChunkLength<ChunkCollector_>::VALUE;

public:
	ChunkCollector_(THost & _host):
		data_host(& _host),
		data_length(0),
		data_chunks(typename Chunk_Holder::allocator_type(_host))
	{
	}

	~ChunkCollector_()
	{
		clear(*this);
	}

};

template <typename THost>
inline void
clear(ChunkCollector_<THost> & me)
{
//IOREV
   typedef ChunkCollector_<THost> TChunkCollector;
   typedef typename TChunkCollector::Chunk_Holder Chunk_Holder;
      
	typename Chunk_Holder::iterator it = me.data_chunks.begin();
	typename Chunk_Holder::iterator it_end = me.data_chunks.end();

	for (; it != it_end; ++it)
	{
		deallocate(me.data_host, *it, TChunkCollector::CHUNK_LENGTH);
	}

	me.data_chunks.clear();
	me.data_length = 0;
}

template <typename THost>
inline typename Size<THost>::Type
length(ChunkCollector_<THost> const & me)
{
//IOREV 
	return me.data_length;
}

template <typename THost>
inline void
_setLength(ChunkCollector_<THost> & me, typename Size<THost>::Type new_length)
{
//IOREV
	me.data_length = new_length;
}

template <typename THost>
inline int
chunkCount(ChunkCollector_<THost> const & me)
{
//IOREV
	return me.data_chunks.size();
}

template <typename THost>
inline typename Value<THost>::Type *
getChunk(ChunkCollector_<THost> const & me, int chunk_number)
{
//IOREV
	return me.data_chunks[chunk_number];
}

template <typename THost>
inline typename Value<THost>::Type *
createChunk(ChunkCollector_<THost> & me)
{
//IOREV
   typedef ChunkCollector_<THost> TChunkCollector;
	typename Value<THost>::Type * new_chunk;
	allocate(me.data_host, new_chunk, TChunkCollector::CHUNK_LENGTH);
	me.data_chunks.push_back(new_chunk);
	return new_chunk;
}
	
//////////////////////////////////////////////////////////////////////////////

template <typename THost>
struct Host<ChunkCollector_<THost> >
{
//IOREV
	typedef THost Type;
};

template <typename THost>
struct Host<ChunkCollector_<THost> const >
{
//IOREV
	typedef THost Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename THost>
struct Value<ChunkCollector_<THost> >
{
//IOREV
	typedef typename Value<THost>::Type Type;
};

template <typename THost>
struct Value<ChunkCollector_<THost> const >
{
//IOREV
	typedef typename Value<THost>::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename THost>
struct GetValue<ChunkCollector_<THost> >
{
//IOREV
	typedef typename GetValue<THost>::Type Type;
};

template <typename THost>
struct GetValue<ChunkCollector_<THost> const >
{
//IOREV
	typedef typename GetValue<THost>::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename THost>
struct Size<ChunkCollector_<THost> >
{
//IOREV
	typedef typename Size<THost>::Type Type;
};

template <typename THost>
struct Size<ChunkCollector_<THost> const >
{
//IOREV
	typedef typename Size<THost>::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

struct AssignStreamToChunkCollector_
{
//IOREV
	template <typename THost, typename TSource>
	static inline void 
	assign_(ChunkCollector_<THost> & target,
		TSource & source)
	{
		clear(target);

		while (!_streamEOF(source))
		{
			typename Value<THost>::Type * chunk = createChunk(target);
			typename Size<THost>::Type count = _streamRead(chunk, source, ChunkLength< ChunkCollector_<THost> >::VALUE);
			_setLength(target, length(target) + count);
		}
	}

	template <typename THost, typename TSource>
	static inline void 
	assign_(ChunkCollector_<THost> & target,
		TSource & source,
		typename Size< ChunkCollector_<THost> >::Type limit)
	{
		clear(target);

		while (!_streamEOF(source))
		{
			typename Value<THost>::Type * chunk = createChunk(target);
			typename Size<THost>::Type count = _streamRead(chunk, source, ChunkLength< ChunkCollector_<THost> >::VALUE);
			_setLength(target, length(target) + count);

			if (length(target) >= limit)
			{
				_setLength(target, limit);
				break;
			}
		}
	}
};

//////////////////////////////////////////////////////////////////////////////

template <typename THost, typename TSource>
inline void 
assign(ChunkCollector_<THost> & target,
	   TSource & source)
{
//IOREV
	AssignStreamToChunkCollector_::assign_(target, source);
}
template <typename THost, typename TSource>
inline void 
assign(ChunkCollector_<THost> & target,
	   TSource const & source)
{
//IOREV
	AssignStreamToChunkCollector_::assign_(target, source);
}

template <typename THost, typename TSource, typename TSize>
inline void 
assign(ChunkCollector_<THost> & target,
	   TSource & source,
	   TSize limit)
{
//IOREV
	AssignStreamToChunkCollector_::assign_(target, source, limit);
}
template <typename THost, typename TSource, typename TSize>
inline void 
assign(ChunkCollector_<THost> & target,
	   TSource const & source,
	   TSize limit)
{
//IOREV
	AssignStreamToChunkCollector_::assign_(target, source, limit);
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

template <typename TExpand>
struct AssignChunkCollectorToString_
{
//IOREV
	template <typename TTarget, typename TSource>
	static void assign_(
		TTarget & target, 
		TSource & source)
	{
		typename Size<TTarget>::Type part_length = _clearSpace(target, length(source), TExpand());

		int i_end = chunkCount(source);
		typename Value<TTarget>::Type * pos = begin(target);
		for (int i = 0; i < i_end; ++i)
		{
			bool is_last_chunk = ( part_length <= ChunkLength<TSource>::VALUE);
			typename Size<TTarget>::Type chunk_length = (is_last_chunk) ? part_length : ChunkLength<TSource>::VALUE;
			typename Value<TSource>::Type * chunk = getChunk(source, i);
			
			arrayConstructCopy(chunk, chunk + chunk_length, pos);
			if (is_last_chunk) break;
			pos += chunk_length;
		}
	}

	template <typename TTarget, typename TSource>
	static void assign_(
		TTarget & target, 
		TSource & source,
		typename Size<TTarget>::Type limit)
	{
		typename Size<TTarget>::Type part_length = _clearSpace(target, length(source), limit, TExpand());

		int i_end = chunkCount(source);
		typename Value<TTarget>::Type * pos = begin(target);
		for (int i = 0; i < i_end; ++i)
		{
			bool is_last_chunk = ( part_length <= ChunkLength<TSource>::VALUE);
			typename Size<TTarget>::Type chunk_length = (is_last_chunk) ? part_length : ChunkLength<TSource>::VALUE;
			typename Value<TSource>::Type * chunk = getChunk(source, i);
			
			arrayConstructCopy(chunk, chunk + chunk_length, pos);
			if (is_last_chunk) break;
			pos += chunk_length;
		}
	}
};

//////////////////////////////////////////////////////////////////////////////

template <typename TTargetValue, typename TTargetSpec, typename TSourceHost, typename TExpand>
inline void 
assign(String<TTargetValue, TTargetSpec> & target,
	   ChunkCollector_<TSourceHost> const & source,
	   Tag<TExpand> /*tag*/)
{
//IOREV
	AssignChunkCollectorToString_<Tag<TExpand> >::assign_(target, source);
}
template <typename TTargetValue, typename TTargetSpec, typename TSourceHost, typename TExpand>
inline void 
assign(String<TTargetValue, TTargetSpec> & target,
	   ChunkCollector_<TSourceHost> const & source,
	   typename Size< String<TTargetValue, TTargetSpec> >::Type limit,
	   Tag<TExpand> /*tag*/)
{
//IOREV
	AssignChunkCollectorToString_<Tag<TExpand> >::assign_(target, source, limit);
}


//////////////////////////////////////////////////////////////////////////////

template <typename TExpand>
struct AppendChunkCollectorToString_
{
//IOREV
	template <typename TTarget, typename TSource>
	static void append_(
		TTarget & target, 
		TSource & source)
	{
		typedef typename Size<TTarget>::Type TSize;
		TSize target_length_old = length(target);
		TSize part_length = _clearSpace(target, length(source), target_length_old, target_length_old, TExpand());

		int i_end = chunkCount(source);
		typename Value<TTarget>::Type * pos = begin(target) + target_length_old; //begin(target) was possibly changed by _clearSpace
		for (int i = 0; i < i_end; ++i)
		{
			bool is_last_chunk = ( part_length <= ChunkLength<TSource>::VALUE);
			typename Size<TTarget>::Type chunk_length = (is_last_chunk) ? part_length : (TSize) ChunkLength<TSource>::VALUE;
			typename Value<TSource>::Type * chunk = getChunk(source, i);
			
			arrayConstructCopy(chunk, chunk + chunk_length, pos);
			if (is_last_chunk) break;
			pos += chunk_length;
		}
	}

	template <typename TTarget, typename TSource>
	static void append_(
		TTarget & target, 
		TSource & source,
		typename Size<TTarget>::Type limit)
	{
		typedef typename Size<TTarget>::Type TSize;
		TSize target_length_old = length(target);
		TSize part_length = _clearSpace(target, length(source), target_length_old, target_length_old, limit, TExpand());

		int i_end = chunkCount(source);
		typename Value<TTarget>::Type * pos = begin(target) + target_length_old; //begin(target) was possibly changed by _clearSpace
		for (int i = 0; i < i_end; ++i)
		{
			bool is_last_chunk = ( part_length <= ChunkLength<TSource>::VALUE);
			typename Size<TTarget>::Type chunk_length = (is_last_chunk) ? part_length : (TSize) ChunkLength<TSource>::VALUE;
			typename Value<TSource>::Type * chunk = getChunk(source, i);
			
			arrayConstructCopy(chunk, chunk + chunk_length, pos);
			if (is_last_chunk) break;
			pos += chunk_length;
		}
	}
};

//////////////////////////////////////////////////////////////////////////////

template <typename TTargetValue, typename TTargetSpec, typename TSourceHost, typename TExpand>
inline void 
append(String<TTargetValue, TTargetSpec> & target,
	   ChunkCollector_<TSourceHost> const & source,
	   Tag<TExpand> )
{
//IOREV
	AppendChunkCollectorToString_<Tag<TExpand> >::append_(target, source);
}
template <typename TTargetValue, typename TTargetSpec, typename TSourceHost, typename TExpand>
inline void 
append(String<TTargetValue, TTargetSpec> & target,
	   ChunkCollector_<TSourceHost> const & source,
	   typename Size< String<TTargetValue, TTargetSpec> >::Type limit,
	   Tag<TExpand> )
{
//IOREV
	AppendChunkCollectorToString_<Tag<TExpand> >::append_(target, source, limit);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TExpand>
struct ReplaceChunkCollectorToString_
{
//IOREV
	template <typename TTarget, typename TSource>
	static void replace_(
		TTarget & target,
		typename Size<TTarget>::Type pos_begin,
		typename Size<TTarget>::Type pos_end,
		TSource & source)
	{
		typename Size<TTarget>::Type part_length = _clearSpace(target, length(source), pos_begin, pos_end, TExpand());

		int i_end = chunkCount(source);
		typename Value<TTarget>::Type * pos = begin(target) + pos_begin;
		for (int i = 0; i < i_end; ++i)
		{
			bool is_last_chunk = ( part_length <= ChunkLength<TSource>::VALUE);
			typename Size<TTarget>::Type chunk_length = (is_last_chunk) ? part_length : ChunkLength<TSource>::VALUE;
			typename Value<TSource>::Type * chunk = getChunk(source, i);
			
			arrayConstructCopy(chunk, chunk + chunk_length, pos);
			if (is_last_chunk) break;
			pos += chunk_length;
		}
	}

	template <typename TTarget, typename TSource>
	static void replace_(
		TTarget & target, 
		typename Size<TTarget>::Type pos_begin,
		typename Size<TTarget>::Type pos_end,
		TSource & source,
		typename Size<TTarget>::Type limit)
	{
		typename Size<TTarget>::Type part_length = _clearSpace(target, length(source), pos_begin, pos_end, limit, TExpand());

		int i_end = chunkCount(source);
		typename Value<TTarget>::Type * pos = begin(target) + pos_begin;
		for (int i = 0; i < i_end; ++i)
		{
			bool is_last_chunk = ( part_length <= ChunkLength<TSource>::VALUE);
			typename Size<TTarget>::Type chunk_length = (is_last_chunk) ? part_length : ChunkLength<TSource>::VALUE;
			typename Value<TSource>::Type * chunk = getChunk(source, i);
			
			arrayConstructCopy(chunk, chunk + chunk_length, pos);
			if (is_last_chunk) break;
			pos += chunk_length;
		}
	}
};

//////////////////////////////////////////////////////////////////////////////

template <typename TTargetValue, typename TTargetSpec, typename TSourceHost, typename TExpand>
inline void 
replace(String<TTargetValue, TTargetSpec> & target,
		typename Size< String<TTargetValue, TTargetSpec> >::Type pos_begin,
		typename Size< String<TTargetValue, TTargetSpec> >::Type pos_end,
	   ChunkCollector_<TSourceHost> const & source,
	   Tag<TExpand> /*tag*/)
{
//IOREV
	ReplaceChunkCollectorToString_<Tag<TExpand> >::replace_(target, pos_begin, pos_end, source);
}

template <typename TTargetValue, typename TTargetSpec, typename TSourceHost, typename TExpand>
inline void 
replace(String<TTargetValue, TTargetSpec> & target,
		typename Size< String<TTargetValue, TTargetSpec> >::Type pos_begin,
		typename Size< String<TTargetValue, TTargetSpec> >::Type pos_end,
	   ChunkCollector_<TSourceHost> const & source,
	   typename Size< String<TTargetValue, TTargetSpec> >::Type limit,
	   Tag<TExpand> /*tag*/)
{
//IOREV
	ReplaceChunkCollectorToString_<Tag<TExpand> >::replace_(target, pos_begin, pos_end, source, limit);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TTargetValue, typename TSourceHost, typename TExpand>
inline void 
replace(TTargetValue * target,
		size_t pos_begin,
		size_t pos_end,
		ChunkCollector_<TSourceHost> const & source,
		Tag<TExpand> /*tag*/)
{
//IOREV
	ReplaceChunkCollectorToString_<Tag<TExpand> >::replace_(target, pos_begin, pos_end, source);
}

template <typename TTargetValue, typename TSourceHost, typename TExpand>
inline void 
replace(TTargetValue * target,
		size_t pos_begin,
		size_t pos_end,
		ChunkCollector_<TSourceHost> const & source,
		size_t limit,
		Tag<TExpand> /*tag*/)
{
//IOREV
	ReplaceChunkCollectorToString_<Tag<TExpand> >::replace_(target, pos_begin, pos_end, source, limit);
}
//____________________________________________________________________________
/*
template <typename TTargetValue, typename TSourceHost, typename TExpand>
inline void 
replace(TTargetValue * target,
		size_t pos_begin,
		size_t pos_end,
	   ChunkCollector_<TSourceHost> const & source,
	   Tag<TExpand> tag)
{
	ReplaceChunkCollectorToString_<Tag<TExpand> >::replace_(target, pos_begin, pos_end, source);
}

template <typename TTargetValue, typename TTargetSpec, typename TSourceHost, typename TExpand>
inline void 
replace(String<TTargetValue, TTargetSpec> & target,
		typename Size< String<TTargetValue, TTargetSpec> >::Type pos_begin,
		typename Size< String<TTargetValue, TTargetSpec> >::Type pos_end,
	   ChunkCollector_<TSourceHost> const & source,
	   typename Size< String<TTargetValue, TTargetSpec> >::Type limit,
	   Tag<TExpand> tag)
{
	ReplaceChunkCollectorToString_<Tag<TExpand> >::replace_(target, pos_begin, pos_end, source, limit);
}
*/
//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
