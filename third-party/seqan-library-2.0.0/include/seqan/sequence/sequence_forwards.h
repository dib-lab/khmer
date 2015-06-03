// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
// Copyright (c) 2013 NVIDIA Corporation
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Manual forwards for the sequence module.
// ==========================================================================

#ifndef SEQAN_HEADER_SEQUENCE_FORWARDS_H
#define SEQAN_HEADER_SEQUENCE_FORWARDS_H

#if !defined(_MSC_VER) || _MSC_VER <= 1600

//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{

// TODO(holtgrew): Still required since we dropped support for VS2003?
// Workaround (copied from generated forwards) for VS 2003.
#if defined(_MSC_VER) && (_MSC_VER < 1400)
template <unsigned int SPACE > struct Block;        // "include/seqan/sequence\string_stack.h"(48)
template <typename THostspec > struct Packed;           // "include/seqan/sequence\string_packed.h"(33)
template <typename TValue, typename TSpec > class String;           // "include/seqan/sequence\string_base.h"(54)
template <typename TString, typename TSpec > class StringSet;           // "include/seqan/sequence\sequence_multiple.h"(98)

template <typename TValue, typename THostspec, typename TTag> inline typename Iterator<String<TValue, Packed<THostspec> >, Tag<TTag> const>::Type end(String<TValue, Packed<THostspec> > & me, Tag<TTag> const tag_);           // "include/seqan/sequence\string_packed.h"(470)
template <typename TValue, typename THostspec, typename TTag> inline typename Iterator<String<TValue, Packed<THostspec> > const, Tag<TTag> const>::Type end(String<TValue, Packed<THostspec> > const & me, Tag<TTag> const tag_);           // "include/seqan/sequence\string_packed.h"(478)
template <typename TValue, unsigned int SPACE, typename TSpec> inline typename Iterator<String<TValue, Block<SPACE> >, Tag<TSpec> const >::Type end(String<TValue, Block<SPACE> > & me, Tag<TSpec> const);           // "include/seqan/sequence\string_stack.h"(209)
template <typename TValue, unsigned int SPACE, typename TSpec> inline typename Iterator<String<TValue, Block<SPACE> > const, Tag<TSpec> const>::Type end(String<TValue, Block<SPACE> > const & me, Tag<TSpec> const);        // "include/seqan/sequence\string_stack.h"(217)
template <typename TString, typename TSpec, typename TTag> inline typename Iterator< StringSet< TString, TSpec >, Tag<TTag> const>::Type end(StringSet< TString, TSpec > & me, Tag<TTag> const tag);        // "include/seqan/sequence\sequence_multiple.h"(1398)
template <typename TString, typename TSpec, typename TTag> inline typename Iterator< StringSet< TString, TSpec > const, Tag<TTag> const>::Type end(StringSet< TString, TSpec > const & me, Tag<TTag> const tag);        // "include/seqan/sequence\sequence_multiple.h"(1405)
#endif  // defined(_MSC_VER) && (_MSC_VER < 1400)

// ==========================================================================
// Adaption Forwards
// ==========================================================================

// TODO(holtgrew): I wonder whether everything below will still be necessary after auto-sequence feature removal? See note below.

// NOTE(holtgrew): My guess / understanding why we need forwards here.
//
// The problem with needing forwards here appears to be that we have default
// implementations of metafunctions Reference<>, Length<> etc.
//
// Consider the setup for a template function A() using getValue(): A tries to
// use getValue().  If there was no default implementation of Reference<> then
// the function getValue() would not be instantiated at this point.  Since
// there is such a default implementation, however, getValue() gets instantiated,
// returns a Reference<std::string> == (std::string &) and this is where the
// compiler balks.
//
// I think that instantiation would get deferred in the case of Reference<> not
// being defined, but I am not sure.  I need to do more research about this.

// --------------------------------------------------------------------------
// Forwards from sequence_interface.h.
// --------------------------------------------------------------------------

template <typename T> struct AllowsFastRandomAccess;
template <typename T> struct DefaultOverflowExplicit;
template <typename T> struct DefaultOverflowImplicit;
template <typename T> struct IsContiguous;
template <typename T> struct IsSequence;
struct TagExact_;
struct TagGenerous_;
struct TagInsist_;
struct TagLimit_;
typedef Tag<TagExact_> Exact;
typedef Tag<TagGenerous_> Generous;
typedef Tag<TagInsist_> Insist;
typedef Tag<TagLimit_> Limit;
typedef Tag<TagInsist_> Tight;
template <typename T> SEQAN_HOST_DEVICE inline typename Iterator<T, Standard>::Type _beginDefault(T & me, Standard);
template <typename T> SEQAN_HOST_DEVICE inline typename Iterator<T const, Standard>::Type _beginDefault(T const & me, Standard);
template <typename T> SEQAN_HOST_DEVICE inline typename Iterator<T, Rooted>::Type _beginDefault(T & me, Rooted);
template <typename T> SEQAN_HOST_DEVICE inline typename Iterator<T const, Rooted>::Type _beginDefault(T const & me, Rooted);
template <typename T, typename TSize, typename TExpand> inline typename Size<T>::Type _capacityReturned(T & me, TSize, Tag<TExpand>);
template <typename T, typename TSize> inline typename Size<T>::Type _capacityReturned(T &, TSize new_capacity, Insist const & );
template <typename T, typename TSize> inline TSize _computeSizeForCapacity(T const & , TSize capacity);
template <typename T> SEQAN_HOST_DEVICE inline typename Iterator<T, Standard>::Type _endDefault(T & me, Standard);
template <typename T> SEQAN_HOST_DEVICE inline typename Iterator<T const, Standard>::Type _endDefault(T const & me, Standard);
template <typename T> SEQAN_HOST_DEVICE inline typename Iterator<T, Rooted>::Type _endDefault(T & me, Rooted);
template <typename T> SEQAN_HOST_DEVICE inline typename Iterator<T const, Rooted>::Type _endDefault(T const & me, Rooted);
template <typename TTarget, typename TSource> inline void append(TTarget & target, TSource & source);
template <typename TTarget, typename TSource> inline void append(TTarget const & target, TSource & source);
template <typename TTarget, typename TSource> inline void append(TTarget & target, TSource const & source);
template <typename TTarget, typename TSource> inline void append(TTarget const & target, TSource const & source);
template <typename TTarget, typename TSource> inline void append(TTarget & target, TSource & source, typename Size<TTarget>::Type limit);
template <typename TTarget, typename TSource> inline void append(TTarget const & target, TSource & source, typename Size<TTarget>::Type limit);
template <typename TTarget, typename TSource> inline void append(TTarget & target, TSource const & source, typename Size<TTarget>::Type limit);
template <typename TTarget, typename TSource> inline void append(TTarget const & target, TSource const & source, typename Size<TTarget>::Type limit);
template <typename T, typename TValue> inline void appendValue(T SEQAN_FORWARD_ARG me, TValue SEQAN_FORWARD_CARG _value);
#ifndef SEQAN_CXX11_STANDARD
template <typename T, typename TValue> inline void appendValue(T const & me, TValue const & _value);
#endif
template <typename TTarget, typename TSource> inline void assign(TTarget & target, TSource & source, typename Size<TTarget>::Type limit);
template <typename TTarget, typename TSource> inline void assign(TTarget const & target, TSource & source, typename Size<TTarget>::Type limit);
template <typename TTarget, typename TSource> inline void assign(TTarget & target, TSource const & source, typename Size<TTarget>::Type limit);
template <typename TTarget, typename TSource> inline void assign(TTarget const & target, TSource const & source, typename Size<TTarget>::Type limit);
template <typename T, typename TValue, typename TPos> inline void assignValue(T & me, TPos pos, TValue const & _value);
template <typename T> SEQAN_HOST_DEVICE inline typename Reference<T const>::Type back(T const & me);
template <typename T> SEQAN_HOST_DEVICE inline typename Reference<T>::Type back(T & me);
template <typename T> SEQAN_HOST_DEVICE inline typename Iterator<T, typename DefaultGetIteratorSpec<T>::Type>::Type begin(T & me);
template <typename T> SEQAN_HOST_DEVICE inline typename Iterator<T const, typename DefaultGetIteratorSpec<T>::Type>::Type begin(T const & me);
template <typename T, typename TSpec> SEQAN_HOST_DEVICE inline typename Iterator<T, Tag<TSpec> const>::Type begin(T & me, Tag<TSpec> const tag_);
template <typename T, typename TSpec> SEQAN_HOST_DEVICE inline typename Iterator<T const, Tag<TSpec> const>::Type begin(T const & me, Tag<TSpec> const tag_);
template <typename T> inline typename Position<T>::Type beginPosition(T &);
template <typename T> inline typename Position<T>::Type beginPosition(T const &);
template <typename T> SEQAN_HOST_DEVICE inline typename Size<T const>::Type capacity(T const & me);
template <typename T, typename TSize> inline TSize computeGenerousCapacity(T const & , TSize capacity);
template <typename T> SEQAN_HOST_DEVICE inline bool empty(T const & me);
template <typename T> SEQAN_HOST_DEVICE inline typename Iterator<T, typename DefaultGetIteratorSpec<T>::Type>::Type end(T & me);
template <typename T> SEQAN_HOST_DEVICE inline typename Iterator<T const, typename DefaultGetIteratorSpec<T>::Type>::Type end(T const & me);
template <typename T, typename TSpec> SEQAN_HOST_DEVICE inline typename Iterator<T, Tag<TSpec> const>::Type end(T & me, Tag<TSpec> const tag_);
template <typename T, typename TSpec> SEQAN_HOST_DEVICE inline typename Iterator<T const, Tag<TSpec> const>::Type end(T const & me, Tag<TSpec> const tag_);
template <typename T> inline typename Position<T>::Type endPosition(T & me);
template <typename T> inline typename Position<T>::Type endPosition(T const & me);
template <typename T, typename TBeginPosition, typename TEndPosition> inline void erase(T & me, TBeginPosition pos, TEndPosition pos_end);
template <typename T, typename TPosition> inline void erase(T & me, TPosition pos);
template <typename T, typename TBeginPosition, typename TEndPosition> inline void erase(T const & me, TBeginPosition pos, TEndPosition pos_end);
template <typename T, typename TPosition> inline void erase(T const & me, TPosition pos);
template <typename T> inline void eraseBack(T & me);
template <typename T> inline typename Reference<T>::Type front(T & me);
template <typename T> inline typename Reference<T const>::Type front(T const & me);
template <typename T> inline void const * getObjectId(T const & me);
template <typename T, typename TPos> inline typename GetValue<T>::Type getValue(T & me, TPos pos);
template <typename T, typename TPos> inline typename GetValue<T const>::Type getValue(T const & me, TPos pos);
template <typename T, typename TPosition, typename TSeq, typename TExpand> inline void insert(T & me, TPosition pos, TSeq const & insertSeq, Tag<TExpand>);
template <typename T, typename TPosition, typename TSeq, typename TExpand> inline void insert(T const & me, TPosition pos, TSeq const & insertSeq, Tag<TExpand>);
template <typename T, typename TPosition, typename TSeq> inline void insert(T & me, TPosition pos, TSeq const & insertSeq);
template <typename T, typename TPosition, typename TSeq> inline void insert(T const & me, TPosition pos, TSeq const & insertSeq);
template <typename T, typename TPosition, typename TValue> inline void insertValue(T & me, TPosition pos, TValue const & _value);
template <typename T, typename TPosition, typename TValue> inline void insertValue(T const & me, TPosition pos, TValue const & _value);
template <typename T, typename TPos> inline typename Iterator<T, typename DefaultGetIteratorSpec<T>::Type>::Type iter(T & me, TPos pos);
template <typename T, typename TPos> inline typename Iterator<T const, typename DefaultGetIteratorSpec<T>::Type>::Type iter(T const & me, TPos pos);
template <typename T, typename TPos, typename TTag> inline typename Iterator<T, Tag<TTag> const>::Type iter(T & me, TPos pos, Tag<TTag> const tag_);
template <typename T, typename TPos, typename TTag> inline typename Iterator<T const, Tag<TTag> const>::Type iter(T const & me, TPos pos, Tag<TTag> const tag_);
template <typename T> inline typename Size<T>::Type length(T const & );
template <typename T, typename TValue, typename TPos> inline void moveValue(T & me, TPos pos, TValue const & _value);
template <typename T, typename TValue, typename TPos> inline void moveValue(T const & me, TPos pos, TValue const & _value);
template <typename TTarget, typename TPositionBegin, typename TPositionEnd, typename TSource> inline void replace(TTarget & target, TPositionBegin pos_begin, TPositionEnd pos_end, TSource & source);
template <typename TTarget, typename TPositionBegin, typename TPositionEnd, typename TSource> inline void replace(TTarget const & target, TPositionBegin pos_begin, TPositionEnd pos_end, TSource & source);
template <typename TTarget, typename TPositionBegin, typename TPositionEnd, typename TSource> inline void replace(TTarget & target, TPositionBegin pos_begin, TPositionEnd pos_end, TSource const & source);
template <typename TTarget, typename TPositionBegin, typename TPositionEnd, typename TSource> inline void replace(TTarget const & target, TPositionBegin pos_begin, TPositionEnd pos_end, TSource const & source);
template <typename TTarget, typename TPositionBegin, typename TPositionEnd, typename TSource> inline void replace(TTarget & target, TPositionBegin pos_begin, TPositionEnd pos_end, TSource & source, typename Size<TTarget>::Type limit);
template <typename TTarget, typename TPositionBegin, typename TPositionEnd, typename TSource> inline void replace(TTarget const & target, TPositionBegin pos_begin, TPositionEnd pos_end, TSource & source, typename Size<TTarget>::Type limit);
template <typename TTarget, typename TPositionBegin, typename TPositionEnd, typename TSource> inline void replace(TTarget & target, TPositionBegin pos_begin, TPositionEnd pos_end, TSource const & source, typename Size<TTarget>::Type limit);
template <typename TTarget, typename TPositionBegin, typename TPositionEnd, typename TSource> inline void replace(TTarget const & target, TPositionBegin pos_begin, TPositionEnd pos_end, TSource const & source, typename Size<TTarget>::Type limit);
template <typename T, typename TSize, typename TExpand> inline typename Size<T>::Type reserve(T & me, TSize const & new_capacity, Tag<TExpand> tag);
template <typename T, typename TSize> inline typename Size<T>::Type reserve(T & me, TSize const & new_capacity);
template <typename T, typename TSize> inline typename Size<T>::Type resize(T & me, TSize new_length);
template <typename T, typename TSize, typename TValue> inline typename Size<T>::Type resize(T & me, TSize new_length, TValue const & val);
template <typename T, typename TSize, typename TBeginPosition, typename TEndPosition> inline TSize resizeSpace(T & me, TSize size, TBeginPosition pos_begin, TEndPosition pos_end);
template <typename T, typename TSize, typename TBeginPosition, typename TEndPosition, typename TLimit> inline TSize resizeSpace(T & me, TSize size, TBeginPosition pos_begin, TEndPosition pos_end, TLimit limit);
template <typename T1, typename T2> inline bool shareResources(T1 const & obj1, T2 const & obj2);
template <typename T> inline void shrinkToFit(T & me);
template <typename T, typename TPos> SEQAN_HOST_DEVICE inline typename Reference<T>::Type value(T & me, TPos );
template <typename T, typename TPos> SEQAN_HOST_DEVICE inline typename Reference<T const>::Type value(T const & me, TPos );

// --------------------------------------------------------------------------
// Forwards For std::vector
// --------------------------------------------------------------------------

template <typename TChar, typename TAlloc, typename TSource> inline void append(std::vector<TChar, TAlloc> & target, TSource const & source, Generous);
template <typename TChar, typename TAlloc, typename TSource> inline void append(std::vector<TChar, TAlloc> & target, TSource const & source, typename Size< std::vector<TChar, TAlloc> >::Type limit, Generous);
template <typename TChar, typename TAlloc, typename TSource> inline void append(std::vector<TChar, TAlloc> & target, TSource const & source, Limit);
template <typename TChar, typename TAlloc, typename TSource> inline void append(std::vector<TChar, TAlloc> & target, TSource const & source, typename Size< std::vector<TChar, TAlloc> >::Type limit, Limit);
template <typename TChar, typename TAlloc, typename TValue, typename TTag> inline void appendValue(std::vector<TChar, TAlloc> & me, TValue const & _value, TTag);
template <typename TChar, typename TAlloc, typename TValue> inline void appendValue(std::vector<TChar, TAlloc> & me, TValue const & _value, Limit);
template <typename TChar, typename TAlloc, typename TSource> inline void assign(std::vector<TChar, TAlloc> & target, TSource & source);
template <typename TChar, typename TAlloc, typename TSource> inline void assign(std::vector<TChar, TAlloc> & target, TSource const & source);
template <typename TChar, typename TAlloc, typename TSource, typename TSize> inline void assign(std::vector<TChar, TAlloc> & target, TSource & source, TSize limit);
template <typename TChar, typename TAlloc, typename TSource, typename TSize> inline void assign(std::vector<TChar, TAlloc> & target, TSource const & source, TSize limit);
template <typename TChar, typename TAlloc, typename TSource> inline void assign(std::vector<TChar, TAlloc> & target, TSource & source, Generous);
template <typename TChar, typename TAlloc, typename TSource> inline void assign(std::vector<TChar, TAlloc> & target, TSource const & source, Generous);
template <typename TChar, typename TAlloc, typename TSource> inline void assign(std::vector<TChar, TAlloc> & target, TSource & source, typename Size< std::vector<TChar, TAlloc> >::Type limit, Generous);
template <typename TChar, typename TAlloc, typename TSource> inline void assign(std::vector<TChar, TAlloc> & target, TSource const & source, typename Size< std::vector<TChar, TAlloc> >::Type limit, Generous);
template <typename TChar, typename TAlloc, typename TSource> inline void assign(std::vector<TChar, TAlloc> & target, TSource & source, Limit);
template <typename TChar, typename TAlloc, typename TSource> inline void assign(std::vector<TChar, TAlloc> & target, TSource const & source, Limit);
template <typename TChar, typename TAlloc, typename TSource> inline void assign(std::vector<TChar, TAlloc> & target, TSource & source, typename Size< std::vector<TChar, TAlloc> >::Type limit, Limit);
template <typename TChar, typename TAlloc, typename TSource> inline void assign(std::vector<TChar, TAlloc> & target, TSource const & source, typename Size< std::vector<TChar, TAlloc> >::Type limit, Limit);
template <typename TChar, typename TAlloc, typename TSource> inline void assign_std_vector_Generous_impl(std::vector<TChar, TAlloc> & target, TSource & source, typename Size< std::vector<TChar, TAlloc> >::Type limit);
template <typename TChar, typename TAlloc> inline typename Iterator< std::vector<TChar, TAlloc>, Standard>::Type begin(std::vector<TChar, TAlloc> & me, Standard);
template <typename TChar, typename TAlloc> inline typename Iterator< std::vector<TChar, TAlloc> const, Standard>::Type begin(std::vector<TChar, TAlloc> const & me, Standard);
template <typename TChar, typename TAlloc> inline typename Size< std::vector<TChar, TAlloc> >::Type capacity(std::vector<TChar, TAlloc> const & me);
template <typename TChar, typename TAlloc> inline void clear(std::vector<TChar, TAlloc> & me);
template <typename TChar, typename TAlloc> inline bool empty(std::vector<TChar, TAlloc> const & me);
template <typename TChar, typename TAlloc> inline typename Iterator< std::vector<TChar, TAlloc>, Standard>::Type end(std::vector<TChar, TAlloc> & me, Standard);
template <typename TChar, typename TAlloc> inline typename Iterator< std::vector<TChar, TAlloc> const, Standard>::Type end(std::vector<TChar, TAlloc> const & me, Standard);
template <typename TChar, typename TAlloc, typename TSize, typename TExpand> inline typename Size< std::vector<TChar, TAlloc> >::Type fill( std::vector<TChar, TAlloc> & me, TSize new_length, TChar const & val, Tag<TExpand>);
template <typename TChar, typename TAlloc> inline void const * getObjectId(std::vector<TChar, TAlloc> const & me);
template <typename TChar, typename TAlloc> inline typename Size< std::vector<TChar, TAlloc> >::Type length(std::vector<TChar, TAlloc> const & me);
template <typename TChar, typename TAlloc, typename TSource> inline void replace(std::vector<TChar, TAlloc> & target, typename Position< std::vector<TChar, TAlloc> >::Type pos_begin, typename Position< std::vector<TChar, TAlloc> >::Type pos_end, TSource const & source, Generous);
template <typename TChar, typename TAlloc, typename TSource> inline void replace(std::vector<TChar, TAlloc> & target, typename Position< std::vector<TChar, TAlloc> >::Type pos_begin, typename Position< std::vector<TChar, TAlloc> >::Type pos_end, TSource const & source, typename Size< std::vector<TChar, TAlloc> >::Type limit, Generous);
template <typename TChar, typename TAlloc, typename TSource> inline void replace(std::vector<TChar, TAlloc> & target, typename Position< std::vector<TChar, TAlloc> >::Type pos_begin, typename Position< std::vector<TChar, TAlloc> >::Type pos_end, TSource const & source, Limit);
template <typename TChar, typename TAlloc, typename TSource> inline void replace(std::vector<TChar, TAlloc> & target, typename Position< std::vector<TChar, TAlloc> >::Type pos_begin, typename Position< std::vector<TChar, TAlloc> >::Type pos_end, TSource const & source, typename Size< std::vector<TChar, TAlloc> >::Type limit, Limit);
template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource, typename TExpand> inline void replace(std::vector<TChar, TAlloc> & target, typename Iterator< std::vector<TChar, TAlloc>, Rooted>::Type pos_begin, typename Iterator< std::vector<TChar, TAlloc>, Rooted>::Type pos_end, TSource & source, Tag<TExpand> tag);
template <typename TChar, typename TAlloc, typename TSize, typename TExpand> inline typename Size< std::vector<TChar, TAlloc> >::Type reserve( std::vector<TChar, TAlloc> & seq, TSize new_capacity, Tag<TExpand> tag);
template <typename TChar, typename TAlloc, typename TSize> inline typename Size< std::vector<TChar, TAlloc> >::Type reserve( std::vector<TChar, TAlloc> & seq, TSize new_capacity, Insist const &);
template <typename TChar, typename TAlloc, typename TSize> inline typename Size< std::vector<TChar, TAlloc> >::Type reserve( std::vector<TChar, TAlloc> & seq, TSize new_capacity, Limit const &);
template <typename TChar, typename TAlloc, typename TSize, typename TExpand> inline typename Size< std::vector<TChar, TAlloc> >::Type resize( std::vector<TChar, TAlloc> & me, TSize new_length, Tag<TExpand>);
template <typename TChar, typename TAlloc, typename TPos> inline typename GetValue< std::vector<TChar, TAlloc> >::Type value(std::vector<TChar, TAlloc> & me, TPos pos);
template <typename TChar, typename TAlloc, typename TPos> inline typename GetValue< std::vector<TChar, TAlloc> const>::Type value(std::vector<TChar, TAlloc> const & me, TPos pos);

// --------------------------------------------------------------------------
// Forwards For std::string
// --------------------------------------------------------------------------

template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource> inline void append(std::basic_string<TChar, TCharTraits, TAlloc> & target, TSource const & source, Generous);
template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource> inline void append(std::basic_string<TChar, TCharTraits, TAlloc> & target, TSource const & source, typename Size< std::basic_string<TChar, TCharTraits, TAlloc> >::Type limit, Generous);
template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource> inline void append(std::basic_string<TChar, TCharTraits, TAlloc> & target, TSource const & source, Limit);
template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource> inline void append(std::basic_string<TChar, TCharTraits, TAlloc> & target, TSource const & source, typename Size< std::basic_string<TChar, TCharTraits, TAlloc> >::Type limit, Limit);
template <typename TChar, typename TCharTraits, typename TAlloc, typename TValue, typename TTag> inline void appendValue(std::basic_string<TChar, TCharTraits, TAlloc> & me, TValue const & _value, TTag);
template <typename TChar, typename TCharTraits, typename TAlloc, typename TValue> inline void appendValue(std::basic_string<TChar, TCharTraits, TAlloc> & me, TValue const & _value, Limit);
template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource> inline void assign(std::basic_string<TChar, TCharTraits, TAlloc> & target, TSource & source);
template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource> inline void assign(std::basic_string<TChar, TCharTraits, TAlloc> & target, TSource const & source);
template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource, typename TSize> inline void assign(std::basic_string<TChar, TCharTraits, TAlloc> & target, TSource & source, TSize limit);
template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource, typename TSize> inline void assign(std::basic_string<TChar, TCharTraits, TAlloc> & target, TSource const & source, TSize limit);
template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource> inline void assign(std::basic_string<TChar, TCharTraits, TAlloc> & target, TSource & source, Generous);
template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource> inline void assign(std::basic_string<TChar, TCharTraits, TAlloc> & target, TSource const & source, Generous);
template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource> inline void assign(std::basic_string<TChar, TCharTraits, TAlloc> & target, TSource & source, typename Size< std::basic_string<TChar, TCharTraits, TAlloc> >::Type limit, Generous);
template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource> inline void assign(std::basic_string<TChar, TCharTraits, TAlloc> & target, TSource const & source, typename Size< std::basic_string<TChar, TCharTraits, TAlloc> >::Type limit, Generous);
template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource> inline void assign(std::basic_string<TChar, TCharTraits, TAlloc> & target, TSource & source, Limit);
template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource> inline void assign(std::basic_string<TChar, TCharTraits, TAlloc> & target, TSource const & source, Limit);
template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource> inline void assign(std::basic_string<TChar, TCharTraits, TAlloc> & target, TSource & source, typename Size< std::basic_string<TChar, TCharTraits, TAlloc> >::Type limit, Limit);
template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource> inline void assign(std::basic_string<TChar, TCharTraits, TAlloc> & target, TSource const & source, typename Size< std::basic_string<TChar, TCharTraits, TAlloc> >::Type limit, Limit);
template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource> inline void assign_std_string_Generous_impl(std::basic_string<TChar, TCharTraits, TAlloc> & target, TSource & source, typename Size< std::basic_string<TChar, TCharTraits, TAlloc> >::Type limit);
template <typename TChar, typename TCharTraits, typename TAlloc> inline typename Iterator< std::basic_string<TChar, TCharTraits, TAlloc>, Standard>::Type begin(std::basic_string<TChar, TCharTraits, TAlloc> & me, Standard);
template <typename TChar, typename TCharTraits, typename TAlloc> inline typename Iterator< std::basic_string<TChar, TCharTraits, TAlloc> const, Standard>::Type begin(std::basic_string<TChar, TCharTraits, TAlloc> const & me, Standard);
template <typename TChar, typename TCharTraits, typename TAlloc> inline typename Size< std::basic_string<TChar, TCharTraits, TAlloc> >::Type capacity(std::basic_string<TChar, TCharTraits, TAlloc> const & me);
template <typename TChar, typename TCharTraits, typename TAlloc> inline void clear(std::basic_string<TChar, TCharTraits, TAlloc> & me);
template <typename TChar, typename TCharTraits, typename TAlloc> inline bool empty(std::basic_string<TChar, TCharTraits, TAlloc> const & me);
template <typename TChar, typename TCharTraits, typename TAlloc> inline typename Iterator< std::basic_string<TChar, TCharTraits, TAlloc>, Standard>::Type end(std::basic_string<TChar, TCharTraits, TAlloc> & me, Standard);
template <typename TChar, typename TCharTraits, typename TAlloc> inline typename Iterator< std::basic_string<TChar, TCharTraits, TAlloc> const, Standard>::Type end(std::basic_string<TChar, TCharTraits, TAlloc> const & me, Standard);
template <typename TChar, typename TCharTraits, typename TAlloc> inline void const * getObjectId(std::basic_string<TChar, TCharTraits, TAlloc> const & me);
template <typename TChar, typename TCharTraits, typename TAlloc> inline typename std::basic_string<TChar, TCharTraits, TAlloc>::size_type length(std::basic_string<TChar, TCharTraits, TAlloc> const & me);
template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource> inline void replace(std::basic_string<TChar, TCharTraits, TAlloc> & target, typename Position< std::basic_string<TChar, TCharTraits, TAlloc> >::Type pos_begin, typename Position< std::basic_string<TChar, TCharTraits, TAlloc> >::Type pos_end, TSource const & source, Generous);
template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource> inline void replace(std::basic_string<TChar, TCharTraits, TAlloc> & target, typename Position< std::basic_string<TChar, TCharTraits, TAlloc> >::Type pos_begin, typename Position< std::basic_string<TChar, TCharTraits, TAlloc> >::Type pos_end, TSource const & source, typename Size< std::basic_string<TChar, TCharTraits, TAlloc> >::Type limit, Generous);
template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource> inline void replace(std::basic_string<TChar, TCharTraits, TAlloc> & target, typename Position< std::basic_string<TChar, TCharTraits, TAlloc> >::Type pos_begin, typename Position< std::basic_string<TChar, TCharTraits, TAlloc> >::Type pos_end, TSource const & source, Limit);
template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource> inline void replace(std::basic_string<TChar, TCharTraits, TAlloc> & target, typename Position< std::basic_string<TChar, TCharTraits, TAlloc> >::Type pos_begin, typename Position< std::basic_string<TChar, TCharTraits, TAlloc> >::Type pos_end, TSource const & source, typename Size< std::basic_string<TChar, TCharTraits, TAlloc> >::Type limit, Limit);
template <typename TChar, typename TCharTraits, typename TAlloc, typename TSize, typename TExpand> inline typename Size< std::basic_string<TChar, TCharTraits, TAlloc> >::Type reserve( std::basic_string<TChar, TCharTraits, TAlloc> & seq, TSize new_capacity, Tag<TExpand> tag);
template <typename TChar, typename TCharTraits, typename TAlloc, typename TSize> inline typename Size< std::basic_string<TChar, TCharTraits, TAlloc> >::Type reserve( std::basic_string<TChar, TCharTraits, TAlloc> & seq, TSize new_capacity, Insist const &);
template <typename TChar, typename TCharTraits, typename TAlloc, typename TSize> inline typename Size< std::basic_string<TChar, TCharTraits, TAlloc> >::Type reserve( std::basic_string<TChar, TCharTraits, TAlloc> & seq, TSize new_capacity, Limit const &);
template <typename TChar, typename TCharTraits, typename TAlloc, typename TSize, typename TExpand> inline typename Size< std::basic_string<TChar, TCharTraits, TAlloc> >::Type resize( std::basic_string<TChar, TCharTraits, TAlloc> & me, TSize new_length, Tag<TExpand>);
template <typename TChar, typename TCharTraits, typename TAlloc, typename TSize, typename TExpand> inline typename Size< std::basic_string<TChar, TCharTraits, TAlloc> >::Type resize( std::basic_string<TChar, TCharTraits, TAlloc> & me, TSize new_length, TChar const & val, Tag<TExpand>);
template <typename TChar, typename TCharTraits, typename TAlloc, typename TPos> inline typename GetValue< std::basic_string<TChar, TCharTraits, TAlloc> >::Type value(std::basic_string<TChar, TCharTraits, TAlloc> & me, TPos pos);
template <typename TChar, typename TCharTraits, typename TAlloc, typename TPos> inline typename GetValue< std::basic_string<TChar, TCharTraits, TAlloc> const>::Type value(std::basic_string<TChar, TCharTraits, TAlloc> const & me, TPos pos);

// --------------------------------------------------------------------------
// Forwards For arrays and pointers.
// --------------------------------------------------------------------------

template <typename TValue> struct DefaultOverflowExplicit;
template <typename TValue> struct DefaultOverflowImplicit;
template <typename TValue> struct IsContiguous;
template <typename TValue, typename TExpand> inline size_t _clearSpace(TValue * me, size_t size, Tag<TExpand>);
template <typename TValue, typename TExpand> inline size_t _clearSpace(TValue * me, size_t size, size_t limit, Tag<TExpand>);
template <typename TValue, typename TPosition, typename TExpand> inline size_t _clearSpace(TValue * me, size_t size, TPosition pos_begin, TPosition pos_end, Tag<TExpand>);
template <typename TValue, typename TPosition, typename TExpand> inline size_t _clearSpace(TValue * me, size_t size, TPosition pos_begin, TPosition pos_end, size_t limit, Tag<TExpand>);
template <typename TValue> inline void _setLength(TValue * me, size_t new_length);
template <typename TTargetValue, typename TSource, typename TExpand> inline void append(TTargetValue * target, TSource const & source, Tag<TExpand>);
template <typename TTargetValue, typename TSource, typename TExpand> inline void append(TTargetValue * target, TSource const & source, size_t limit, Tag<TExpand>);
template <typename TTargetValue, typename TSourceValue, typename TExpand> inline void append(TTargetValue * target, TSourceValue const * source, Tag<TExpand>);
template <typename TTargetValue, typename TSourceValue, typename TExpand> inline void append(TTargetValue * target, TSourceValue const * source, size_t limit, Tag<TExpand>);
template <typename TTargetValue, typename TSource> inline typename EnableIf<IsCharType<TTargetValue> >::Type assign(TTargetValue * target, TSource & source);
template <typename TTargetValue, typename TSource> inline typename EnableIf<IsCharType<TTargetValue> >::Type assign(TTargetValue * target, TSource const & source);
template <typename TTargetValue, typename TSource, typename TExpand> inline void assign(TTargetValue * target, TSource const & source, Tag<TExpand>);
template <typename TTargetValue, typename TSource, typename TExpand> inline void assign(TTargetValue * target, TSource const & source, size_t limit, Tag<TExpand>);
template <typename TTargetValue, typename TSourceValue, typename TExpand> inline void assign(TTargetValue * target, TSourceValue const * source, Tag<TExpand>);
template <typename TTargetValue, typename TSourceValue, typename TExpand> inline void assign(TTargetValue * target, TSourceValue const * source, size_t limit, Tag<TExpand>);
template <typename TValue, typename TPos> inline void assignValue(TValue * me, TPos pos, TValue const & _value);
template <typename TValue> inline bool atEnd(TValue * pos);
template <typename TValue> inline bool atEnd(TValue const * pos, TValue const * );
template <typename T> inline typename Iterator<T *, typename DefaultGetIteratorSpec<T>::Type>::Type begin(T * me);
template <typename TValue> inline typename Iterator<TValue *, Standard>::Type begin(TValue * me, Standard);
template <typename TValue> SEQAN_HOST_DEVICE inline typename Iterator<TValue const *, Standard>::Type begin(TValue const * me, Standard);
template <typename TValue, typename TSpec> inline typename Iterator<TValue *, Tag<TSpec> const>::Type begin(TValue * me, Tag<TSpec> const);
template <typename TValue, typename TSpec> inline typename Iterator<TValue const *, Tag<TSpec> const>::Type begin(TValue const * me, Tag<TSpec> const);
template <typename TValue> inline void clear(TValue * me);
template <typename TValue> inline bool empty(TValue * me);
template <typename TValue> inline typename Iterator<TValue *, Standard>::Type end(TValue * me, Standard);
template <typename TValue> inline typename Iterator<TValue const *, Standard>::Type end(TValue const * me, Standard);
template <typename TValue, typename TSpec> inline typename Iterator<TValue *, Tag<TSpec> const>::Type end(TValue * me, Tag<TSpec> const tag_);
template <typename TValue, typename TSpec> inline typename Iterator<TValue const *, Tag<TSpec> const>::Type end(TValue const * me, Tag<TSpec> const tag_);
template <typename TLeftValue, typename TRight > inline bool isEqual(TLeftValue * left, TRight const & right);
template <typename TLeftValue, typename TRight> inline bool isGreater(TLeftValue * left, TRight const & right);
template <typename TLeftValue, typename TRight> inline bool isGreaterOrEqual(TLeftValue * left, TRight const & right);
template <typename TLeftValue, typename TRight> inline bool isLess(TLeftValue * left, TRight const & right);
template <typename TLeftValue, typename TRight> inline bool isLessOrEqual(TLeftValue * left, TRight const & right);
template <typename TLeftValue, typename TRight > inline bool isNotEqual(TLeftValue * left, TRight const & right);
template <typename TValue> inline size_t length(TValue * me);
template <typename TValue> inline size_t length(TValue const * me);
inline size_t length(char * me);
inline size_t length(char const * me);
template <typename TTargetValue, typename TSource> inline void move(TTargetValue * & target, TSource & source);
template <typename TTargetValue, typename TSource> inline void move(TTargetValue * & target, TSource const & source);
template <typename TValue, typename TPos> inline void moveValue(TValue * me, TPos pos, TValue const & _value);
template <typename TTargetValue, typename TSource, typename TExpand> inline void replace(TTargetValue * target, size_t pos_begin, size_t pos_end, TSource const & source, Tag<TExpand>);
template <typename TTargetValue, typename TSource, typename TExpand> inline void replace(TTargetValue * target, size_t pos_begin, size_t pos_end, TSource const & source, size_t limit, Tag<TExpand>);
template <typename TTargetValue, typename TSourceValue, typename TExpand> inline void replace(TTargetValue * target, size_t pos_begin, size_t pos_end, TSourceValue const * source, Tag<TExpand>);
template <typename TTargetValue, typename TSourceValue, typename TExpand> inline void replace(TTargetValue * target, size_t pos_begin, size_t pos_end, TSourceValue const * source, size_t limit, Tag<TExpand>);
template <typename TValue, typename TSize, typename TExpand> inline size_t resize( TValue * me, TSize new_length, Tag<TExpand>);
template <typename TValue, typename TSize, typename TExpand> inline size_t resize( TValue * me, TSize new_length, TValue const & val, Tag<TExpand>);
template <typename TValue, typename TPos> inline TValue & value(TValue * me, TPos pos);
template <typename TValue, typename TPos> inline TValue const & value(TValue const * me, TPos pos);

// --------------------------------------------------------------------------
// Forwards For thrust::device_vector
// --------------------------------------------------------------------------

#ifdef PLATFORM_CUDA
template <typename TChar,  typename TAlloc> inline void const * getObjectId(thrust::device_vector<TChar, TAlloc> const & me);
template <typename TChar,  typename TAlloc> inline typename Iterator<thrust::device_vector<TChar, TAlloc>, Standard>::Type begin(thrust::device_vector<TChar, TAlloc> & me, Standard);
template <typename TChar,  typename TAlloc> inline typename Iterator<thrust::device_vector<TChar, TAlloc> const, Standard>::Type begin(thrust::device_vector<TChar, TAlloc> const & me, Standard);
template <typename TChar, typename TAlloc> inline typename Iterator<thrust::device_vector<TChar, TAlloc>, Standard>::Type end(thrust::device_vector<TChar, TAlloc> & me, Standard);
template <typename TChar,  typename TAlloc> inline typename Iterator<thrust::device_vector<TChar, TAlloc> const, Standard>::Type end(thrust::device_vector<TChar, TAlloc> const & me, Standard);
template <typename TChar,  typename TAlloc, typename TPos> inline typename GetValue<thrust::device_vector<TChar, TAlloc> >::Type value(thrust::device_vector<TChar, TAlloc> & me, TPos pos);
template <typename TChar,  typename TAlloc, typename TPos> inline typename GetValue<thrust::device_vector<TChar, TAlloc> const>::Type value(thrust::device_vector<TChar, TAlloc> const & me, TPos pos);
template <typename TChar, typename TAlloc> inline typename Size<thrust::device_vector<TChar, TAlloc> >::Type length(thrust::device_vector<TChar, TAlloc> const & me);
template <typename TChar, typename TAlloc> inline typename Size<thrust::device_vector<TChar, TAlloc> >::Type capacity(thrust::device_vector<TChar, TAlloc> const & me);
template <typename TChar, typename TAlloc> inline bool empty(thrust::device_vector<TChar, TAlloc> const & me);
template <typename TChar,  typename TAlloc> inline void clear(thrust::device_vector<TChar, TAlloc> & me);
template <typename TChar, typename TAlloc> inline typename Reference<thrust::device_vector<TChar, TAlloc> >::Type back(thrust::device_vector<TChar, TAlloc> & list);
template <typename TChar, typename TAlloc> inline typename Reference<thrust::device_vector<TChar, TAlloc> const>::Type back(thrust::device_vector<TChar, TAlloc> const & list);
template <typename TChar,  typename TAlloc, typename TSource> inline void assign(thrust::device_vector<TChar, TAlloc> & target, TSource & source);
template <typename TChar,  typename TAlloc, typename TSource> inline void assign(thrust::device_vector<TChar, TAlloc> & target, TSource const & source);
template <typename TChar,  typename TAlloc, typename TSource, typename TSize> inline void assign(thrust::device_vector<TChar, TAlloc> & target, TSource & source, TSize limit);
template <typename TChar,  typename TAlloc, typename TSource, typename TSize>
inline void assign(thrust::device_vector<TChar, TAlloc> & target, TSource const & source, TSize limit);
template <typename TChar,  typename TAlloc, typename TSource> inline void assign(thrust::device_vector<TChar, TAlloc> & target, TSource & source, Generous);
template <typename TChar, typename TAlloc, typename TSource> inline void assign(thrust::device_vector<TChar, TAlloc> & target, TSource const & source, Generous);
template <typename TChar,  typename TAlloc, typename TSource> inline void assign_std_vector_Generous_impl(thrust::device_vector<TChar, TAlloc> & target, TSource & source, typename Size<thrust::device_vector<TChar, TAlloc> >::Type limit);
template <typename TChar,  typename TAlloc, typename TSource> inline void assign(thrust::device_vector<TChar, TAlloc> & target, TSource & source, typename Size<thrust::device_vector<TChar, TAlloc> >::Type limit, Generous);
template <typename TChar,  typename TAlloc, typename TSource> inline void assign(thrust::device_vector<TChar, TAlloc> & target, TSource const & source, typename Size<thrust::device_vector<TChar, TAlloc> >::Type limit, Generous);
template <typename TChar, typename TAlloc, typename TSource> inline void assign(thrust::device_vector<TChar, TAlloc> & target, TSource & source, Limit);
template <typename TChar, typename TAlloc, typename TSource> inline void assign(thrust::device_vector<TChar, TAlloc> & target, TSource const & source, Limit);
template <typename TChar, typename TAlloc, typename TSource> inline void assign(thrust::device_vector<TChar, TAlloc> & target, TSource & source, typename Size<thrust::device_vector<TChar, TAlloc> >::Type limit, Limit);
template <typename TChar, typename TAlloc, typename TSource> inline void assign(thrust::device_vector<TChar, TAlloc> & target, TSource const & source, typename Size<thrust::device_vector<TChar, TAlloc> >::Type limit, Limit);
template <typename TChar, typename TAlloc, typename TSource> inline void append(thrust::device_vector<TChar, TAlloc> & target, TSource const & source, Generous);
template <typename TChar,  typename TAlloc, typename TSource> inline void append(thrust::device_vector<TChar, TAlloc> & target, TSource const & source, typename Size<thrust::device_vector<TChar, TAlloc> >::Type limit, Generous);
template <typename TChar, typename TAlloc, typename TSource> inline void append(thrust::device_vector<TChar, TAlloc> & target, TSource const & source, Limit);
template <typename TChar, typename TAlloc, typename TSource> inline void append(thrust::device_vector<TChar, TAlloc> & target, TSource const & source, typename Size<thrust::device_vector<TChar, TAlloc> >::Type limit, Limit);
template <typename TChar, typename TAlloc, typename TValue, typename TTag> inline void appendValue(thrust::device_vector<TChar, TAlloc> & me, TValue const & _value, TTag);
template <typename TChar, typename TAlloc, typename TValue> inline void appendValue(thrust::device_vector<TChar, TAlloc> & me, TValue const & _value, Limit);
template <typename TChar,  typename TAlloc, typename TSource> inline void replace(thrust::device_vector<TChar, TAlloc> & target, typename Position<thrust::device_vector<TChar, TAlloc> >::Type pos_begin, typename Position<thrust::device_vector<TChar, TAlloc> >::Type pos_end, TSource const & source, Generous);
template <typename TChar, typename TAlloc, typename TSource> inline void replace(thrust::device_vector<TChar, TAlloc> & target, typename Position<thrust::device_vector<TChar, TAlloc> >::Type pos_begin, typename Position<thrust::device_vector<TChar, TAlloc> >::Type pos_end, TSource const & source, typename Size<thrust::device_vector<TChar, TAlloc> >::Type limit, Generous);
template <typename TChar,  typename TAlloc, typename TSource> inline void replace(thrust::device_vector<TChar, TAlloc> & target, typename Position<thrust::device_vector<TChar, TAlloc> >::Type pos_begin, typename Position<thrust::device_vector<TChar, TAlloc> >::Type pos_end, TSource const & source, Limit);
template <typename TChar, typename TAlloc, typename TSource> inline void replace(thrust::device_vector<TChar, TAlloc> & target, typename Position<thrust::device_vector<TChar, TAlloc> >::Type pos_begin, typename Position<thrust::device_vector<TChar, TAlloc> >::Type pos_end, TSource const & source, typename Size<thrust::device_vector<TChar, TAlloc> >::Type limit, Limit);
template<typename TChar, typename TCharTraits, typename TAlloc, typename TSource, typename TExpand> inline void replace(thrust::device_vector<TChar, TAlloc> & target, typename Iterator<thrust::device_vector<TChar, TAlloc>, Rooted>::Type pos_begin, typename Iterator<thrust::device_vector<TChar, TAlloc>, Rooted>::Type pos_end, TSource & source, Tag<TExpand> const tag);
template <typename TChar,  typename TAlloc, typename TSize, typename TExpand> inline typename Size<thrust::device_vector<TChar, TAlloc> >::Type reserve(thrust::device_vector<TChar, TAlloc> & seq, TSize new_capacity, Tag<TExpand> const & tag);
template <typename TChar, typename TAlloc, typename TSize> inline typename Size<thrust::device_vector<TChar, TAlloc> >::Type reserve(thrust::device_vector<TChar, TAlloc> & seq, TSize new_capacity, Insist const &);
template <typename TChar,  typename TAlloc, typename TSize> inline typename Size<thrust::device_vector<TChar, TAlloc> >::Type reserve(thrust::device_vector<TChar, TAlloc> & seq, TSize new_capacity, Limit const &);
template <typename TChar,  typename TAlloc, typename TSize, typename TExpand> inline typename Size<thrust::device_vector<TChar, TAlloc> >::Type resize(thrust::device_vector<TChar, TAlloc> & me, TSize new_length, Tag<TExpand> const &);
template <typename TChar, typename TAlloc, typename TSize, typename TExpand> inline typename Size<thrust::device_vector<TChar, TAlloc> >::Type fill(thrust::device_vector<TChar, TAlloc> & me, TSize new_length, TChar const & val, Tag<TExpand> const &);
#endif

}  // namespace SEQAN_NAMESPACE_MAIN

#endif  // #if !defined(_MSC_VER) || _MSC_VER <= 1600

#endif

