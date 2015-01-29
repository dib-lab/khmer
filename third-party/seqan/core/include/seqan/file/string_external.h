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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// Access a file like a string by implementing a custom paging/lru mechanism.
// ==========================================================================

#ifndef SEQAN_HEADER_STRING_EXTERNAL_H
#define SEQAN_HEADER_STRING_EXTERNAL_H

/* IOREV
 * _nottested_
 * _doc_
 *
 *
 * mostly documented (doc for some functions missing)
 * not tested in any test case nor used in any app right now
 * -> needs testing, especially different iterators
 */



//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{

/**
.Spec.External String:
..cat:Strings
..general:Class.String
..summary:String that is stored in external memory.
..signature:String<TValue, External<> >
..signature:String<TValue, External<TConfig> >
..param.TValue:The value type, that is the type of the items/characters stored in the string.
...remarks:Use @Metafunction.Value@ to get the value type for a given class.
..param.TConfig:A structure to configure the external string.
...type:Tag.ExternalConfig
...type:Tag.ExternalConfigLarge
...type:Tag.ExternalConfigSize
...default:@Tag.ExternalConfigLarge@
..remarks:The External String enables to access sequences larger than the available internal memory (RAM)
by using external memory (e.g. Hard disk, Network storage, ...) via a @Class.File@ object.
Sequences of nearly arbitrary size can be accessed even larger than the logically addressable memory,
i.e. they can in particular contain more than 2^32 elements on a 32bit system (see Tag.ExternalConfigLarge).
See the @Memfunc.External String#String@ constructor for more details.
..remarks:This String also supports fast appending and removing of values at the end (see @Spec.Block String@, @Function.appendValue@)
..remarks:The External String implements a LRU mechanism to swap out pages.
The External String's @Metafunction.Iterator@ detects a forward or backward iteration and asynchronously prefetches pages that
certainly will be accessed and automatically swaps out pages that certainly won't be accessed any more in the iteration
process.
..remarks:The String is implemented like a virtual memory manager.
It divides its character sequence into pages of a fixed length (e.g. 1MB) and maintains a
page table with information for each page (e.g. resides in memory or was swapped out, is dirty and needs to be saved, ...).
Besides the page table the String also contains a size-limited list of page frames. A page frame is reserved internal
memory for a page. When accessing values of a page that is stored in external memory, the page is loaded to a page frame
first. In case that there is no page frame free, another page is swapped out before to free a page frame.
..include:seqan/file.h
*/

/**
.Tag.ExternalConfig:
..cat:Strings
..summary:Standard configuration for the @Spec.External String@.
..signature:String<TValue, External< ExternalConfig<> > >
..signature:String<TValue, External< ExternalConfig<TFile[, pageSize[, frames]]> > >
..param.TFile:Type of file the External String will be based on.
...type:Class.File
..param.pageSize:A positive integer that specifies the number of values in one page.
...remarks:This should be a power of 2, to speed up transfer and calculations.
...default:2^20
..param.frames:A positive integer that specifies the maximum number of pages that should reside in internal memory.
...remarks:To enable prefetching and automatic swap-out, $frames$ should be greater than 1.
...default:2
..remarks:When using this configuration, the @Metafunction.Size@ type of the @Spec.External String@ is $unsigned int$.
Thus, with this configuration at most 4.294.967.296 values can be stored in an @Spec.External String@ on a 32bit system.
For a larger size type, use @Tag.ExternalConfigLarge@.
..include:seqan/file.h
*/
    // standard external string
    // size is uint32
    template < typename TFile_ = File<>,				// default file type
               unsigned PAGESIZE_ = 4 * 1024 * 1024,	// 1MTypes per default
			   unsigned FRAMES_ = 2 >					// simultanous frames
    struct ExternalConfig {
//IOREV _bug_ doc says default page size is 2^20, but it is 2^22
        typedef TFile_ TFile;
        typedef unsigned TSize;
        enum { PAGESIZE = PAGESIZE_ };
        enum { FRAMES = FRAMES_ };
    };

/**
.Tag.ExternalConfigLarge:
..cat:Strings
..summary:Large size type configuration for the @Spec.External String@.
..signature:String<TValue, External< ExternalConfigLarge<> > >
..signature:String<TValue, External< ExternalConfigLarge<TFile[, pageSize[, frames]]> > >
..param.TFile:Type of file the External String will be based on.
...type:Class.File
..param.pageSize:A positive integer that specifies the number of values in one page.
...remarks:This should be a power of 2, to speed up transfer and calculations.
...default:2^20
..param.frames:A positive integer that specifies the maximum number of pages that should reside in internal memory.
...remarks:To enable prefetching and automatic swap-out, $frames$ should be greater than 1.
...default:2
..remarks:When using this configuration, the @Metafunction.Size@ type of the @Spec.External String@ is the @Metafunction.Size@
type of $TFile$. Normally this is a 64bit integer. For a smaller size type, use @Tag.ExternalConfig@.
..remarks:Some data structures store size type values (e.g. suffix arrays in @Class.Index@). To save memory,
you should think of using @Tag.ExternalConfig@.
..include:seqan/file.h
*/
    // the same as ExternalConfig
    // but size type is size type of TFile_ (i.e. uint64)
    //
    // ATTENTION:
    // pipes use the size type 
    // uint64 blows up your suffix arrays, lcp-tables, ...
    template < typename TFile_ = File<>,				// default file type
               unsigned PAGESIZE_ = 4 * 1024 * 1024,	// 1MTypes per default
			   unsigned FRAMES_ = 2 >					// simultanous frames
    struct ExternalConfigLarge {
//IOREV contains warning in code comments, need to investigate
        typedef TFile_ TFile;
        typedef typename MakeUnsigned<typename Size<TFile_>::Type>::Type TSize;
        enum { PAGESIZE = PAGESIZE_ };
        enum { FRAMES = FRAMES_ };
    };

/**
.Tag.ExternalConfigSize:
..cat:Strings
..summary:Arbitrary size type configuration for the @Spec.External String@.
..signature:String<TValue, External< ExternalConfigSize< TSize > > >
..signature:String<TValue, External< ExternalConfigSize< TSize, TFile[, pageSize[, frames]]> > >
..param.TSize:Size type the External String will return via @Metafunction.Size@.
..param.TFile:Type of file the External String will be based on.
...type:Class.File
..param.pageSize:A positive integer that specifies the number of values in one page.
...remarks:This should be a power of 2, to speed up transfer and calculations.
...default:2^20
..param.frames:A positive integer that specifies the maximum number of pages that should reside in internal memory.
...remarks:To enable prefetching and automatic swap-out, $frames$ should be greater than 1.
...default:2
..remarks:When using this configuration, the @Metafunction.Size@ type of the @Spec.External String@ is the @Metafunction.Size@
type of $TFile$. Normally this is a 64bit integer. For a smaller size type, use @Tag.ExternalConfig@.
..remarks:Some data structures store size type values (e.g. suffix arrays in @Class.Index@). To save memory,
you should think of using @Tag.ExternalConfig@.
..include:seqan/file.h
*/
    // custom size type
    template < typename TSize_,
		       typename TFile_ = File<>,				// default file type
               unsigned PAGESIZE_ = 1 * 1024 * 1024,	// 1MTypes per default
			   unsigned FRAMES_ = 2 >					// simultanous frames
    struct ExternalConfigSize {
//IOREV
		typedef TSize_ TSize;
        typedef TFile_ TFile;
        enum { PAGESIZE = PAGESIZE_ };
        enum { FRAMES = FRAMES_ };
    };

    template < typename TConfig = ExternalConfigLarge<> >
    struct External {};
//IOREV


    //////////////////////////////////////////////////////////////////////////////
	// random External String iterator
	template < typename TExtString >
	struct ExtStringIterator
	{
//IOREV
		typedef ExtStringIterator						TIterator;
        typedef ExtStringIterator						TStdIterator;

        typedef typename Value<TExtString>::Type		TValue;
		typedef typename Size<TExtString>::Type			TSize;
		typedef typename Difference<TExtString>::Type	TDifference;
		typedef typename TExtString::TVolatilePtr		TVolatilePtr;

        enum { PAGESIZE = TExtString::PAGESIZE };

		TSize		offset;
		TExtString	*extString;


	    ExtStringIterator():
			offset(0),
			extString(NULL) {}

	    explicit ExtStringIterator(TExtString *_extString, TSize _offset):
			extString(_extString),
			offset(_offset) {}

		//////////////////////////////////////////////////////////////////////////////
		// iterator conversion interface

        ExtStringIterator(const TStdIterator &I):
			offset(I.offset),
			extString(I.extString) {}

		//////////////////////////////////////////////////////////////////////////////
		// iterator arithmetic

        inline TDifference operator- (const TIterator &I) const {
			return offset - I.offset;
		}
		
		inline TIterator operator- (TDifference delta) const {
			return TIterator(extString, offset - delta);
		}
		
		inline TIterator& operator-= (TDifference delta) const {
			offset -= delta;
			return *this;
		}
		
		inline TIterator operator+ (TDifference delta) const {
			return TIterator(extString, offset + delta);
		}
		
		inline TIterator& operator+= (TDifference delta) const {
			offset += delta;
			return *this;
		}
		
		inline TValue & operator* () const {
			return (*extString)[offset];
		}
    
/*		inline TValue const & operator* () const {
			return (*extString)[offset];
		}
*/  
		inline TIterator& operator++ () {
			++offset; return *this;
		}

		inline TIterator operator++ (int) {
			TIterator before = *this;
			++offset; return before;
		}

		inline TIterator& operator-- () {
			--offset; return *this;
		}

		inline TIterator operator-- (int) {
			TIterator before = *this;
			--offset; return before;
		}

		inline bool operator== (const TIterator &I) const {
			SEQAN_ASSERT_EQ(extString, I.extString);
			return offset == I.offset;
		}

		inline bool operator!= (const TIterator &I) const {
			SEQAN_ASSERT_EQ(extString, I.extString);
			return offset != I.offset;
		}

		inline bool operator< (const TIterator &I) const {
			SEQAN_ASSERT_EQ(extString, I.extString);
			return offset < I.offset;
		}
	};


	//////////////////////////////////////////////////////////////////////////////
	// const random External String iterator
	template < typename TExtString >
	struct ExtStringConstIterator
	{
//IOREV
		typedef ExtStringConstIterator					TIterator;
        typedef ExtStringIterator<TExtString>			TStdIterator;
        typedef ExtStringConstIterator					TStdConstIterator;

        typedef typename Value<TExtString>::Type		TValue;
		typedef typename Size<TExtString>::Type			TSize;
		typedef typename Difference<TExtString>::Type	TDifference;
		typedef typename TExtString::TVolatilePtr	    TVolatilePtr;

        enum { PAGESIZE = TExtString::PAGESIZE };

		TSize		offset;
		TExtString	*extString;
		
		
	    ExtStringConstIterator():
			offset(0),
			extString(NULL) {}

	    explicit ExtStringConstIterator(TExtString *_extString, TSize _offset):
			extString(_extString),
			offset(_offset) {}

		//////////////////////////////////////////////////////////////////////////////
		// iterator conversion interface

		ExtStringConstIterator(const TStdIterator &I):
			offset(I.offset),
			extString(I.extString) {}

		ExtStringConstIterator(const TStdConstIterator &I):
			offset(I.offset),
			extString(I.extString) {}

		//////////////////////////////////////////////////////////////////////////////
		// iterator arithmetic

		inline TDifference operator- (const TIterator &I) const {
			return offset - I.offset;
		}
		
		inline TIterator operator- (TDifference delta) const {
			return TIterator(extString, offset - delta);
		}
		
		inline TIterator& operator-= (TDifference delta) const {
			offset -= delta;
			return *this;
		}
		
		inline TIterator operator+ (TDifference delta) const {
			return TIterator(extString, offset + delta);
		}
		
		inline TIterator& operator+= (TDifference delta) const {
			offset += delta;
			return *this;
		}
		
		inline TValue const & operator* () const {
			return (*extString)[offset];
		}
    
		inline TIterator& operator++ () {
			++offset; return *this;
		}

		inline TIterator operator++ (int) {
			TIterator before = *this;
			++offset; return before;
		}

		inline TIterator& operator-- () {
			--offset; return *this;
		}

		inline TIterator operator-- (int) {
			TIterator before = *this;
			--offset; return before;
		}

		inline bool operator== (const TIterator &I) const {
			SEQAN_ASSERT_EQ(extString, I.extString);
			return offset == I.offset;
		}

		inline bool operator!= (const TIterator &I) const {
			SEQAN_ASSERT_EQ(extString, I.extString);
			return offset != I.offset;
		}

		inline bool operator< (const TIterator &I) const {
			SEQAN_ASSERT_EQ(extString, I.extString);
			return offset < I.offset;
		}		
	};


	//////////////////////////////////////////////////////////////////////////////
	// forward External String iterator
	template < typename TExtString >
    struct ExtStringFwdIterator
	{
//IOREV
		typedef ExtStringFwdIterator					TIterator;
        typedef ExtStringIterator<TExtString>			TStdIterator;
        typedef ExtStringConstIterator<TExtString>		TStdConstIterator;

		typedef typename Value<TExtString>::Type		TValue;
		typedef typename Size<TExtString>::Type			TSize;
		typedef typename Difference<TExtString>::Type	TDifference;
		typedef typename TExtString::TVolatilePtr	    TVolatilePtr;

        enum { PAGESIZE = TExtString::PAGESIZE };


		TExtString		*extString;

        bool            dirty;
		int     		pageNo;
		unsigned		pageOfs;
        int             prefetch;   // -n .. prefetch n pages downwards, n .. prefetch n pages upwards, 0 .. disabled
		TVolatilePtr	begin;
		
	    ExtStringFwdIterator():
			extString(NULL),
			pageNo(0),
			pageOfs(0),
            prefetch(0),
			begin(NULL) {}

		ExtStringFwdIterator(const TIterator &I):
			extString(I.extString),
			pageNo(I.pageNo),
			pageOfs(I.pageOfs),
            prefetch(I.prefetch),
			begin(NULL) {}

	    explicit ExtStringFwdIterator(TExtString *_extString, TSize _offset):
			extString(_extString),
			pageNo(_offset / PAGESIZE),
			pageOfs(_offset % PAGESIZE),
            prefetch(0),
			begin(NULL) {}

		explicit ExtStringFwdIterator(TExtString *_extString, TSize _pageNo, TSize _pageOfs):
			extString(_extString),
			pageNo(_pageNo),
			pageOfs(_pageOfs),
            prefetch(0),
			begin(NULL) {}

		~ExtStringFwdIterator() {
			invalidate();
		}

		//////////////////////////////////////////////////////////////////////////////
		// iterator conversion interface

		ExtStringFwdIterator(const TStdIterator &I):
			extString(I.extString),
			pageNo(I.offset / PAGESIZE),
            pageOfs(I.offset % PAGESIZE),
            prefetch(0),
			begin(NULL) {}

		inline TIterator& operator=(TStdIterator const & Right_) {
			invalidate();
			pageNo = Right_.offset / PAGESIZE;
			pageOfs = Right_.offset % PAGESIZE;
            extString = Right_.extString;
			return *this;
		}

        inline TSize position() const
        {
            return (TSize)pageNo * (TSize)PAGESIZE + pageOfs;
        }

        inline operator TStdIterator() const {
            return TStdIterator(extString, position());
        }

        inline operator TStdConstIterator() const {
            return TStdConstIterator(extString, position());
        }

		inline TIterator& operator=(TIterator const & Right_) {
			invalidate();
			extString = Right_.extString;
			pageNo = Right_.pageNo;
			pageOfs = Right_.pageOfs;
            prefetch = Right_.prefetch;
			return *this;
		}

		//////////////////////////////////////////////////////////////////////////////
		// iterator arithmetic

		inline TDifference operator- (const TIterator &I) const
        {
            return position() - I.position();
		}
		
		inline TIterator operator- (TDifference delta) const {
			TDifference dPNo  = delta / PAGESIZE;
			TDifference dPOfs = delta % PAGESIZE;
			if (pageOfs >= dPOfs)
				return TIterator(extString, pageNo - dPNo, pageOfs - dPOfs);
			else
				return TIterator(extString, pageNo - dPNo - 1, PAGESIZE + pageOfs - dPOfs);
		}
		
		inline TIterator& operator-= (TDifference delta) {
			TDifference dPNo  = delta / PAGESIZE;
			TDifference dPOfs = delta % PAGESIZE;
			if (pageOfs < dPOfs) {
				++dPNo;
				pageOfs = PAGESIZE + pageOfs - dPOfs;
			} else
				pageOfs -= dPOfs;
			if (dPNo) invalidate(0);
			pageNo -= dPNo;
			return *this;
		}
		
		inline TIterator operator+ (TDifference delta) const {
			TDifference dPNo  = delta / PAGESIZE;
			TDifference nPOfs = pageOfs + delta % PAGESIZE;
			if (nPOfs < PAGESIZE)
				return TIterator(extString, pageNo + dPNo, nPOfs);
			else
				return TIterator(extString, pageNo + dPNo + 1, nPOfs - PAGESIZE);
		}
		
		inline TIterator& operator+= (TDifference delta) {
			TDifference dPNo  = delta / PAGESIZE;
			TDifference nPOfs = pageOfs + delta % PAGESIZE;
			if (nPOfs >= PAGESIZE) {
				++dPNo;
				nPOfs -= PAGESIZE;
			}
			if (dPNo) invalidate(0);
			pageNo += dPNo;
			pageOfs = nPOfs;
			return *this;
		}
		
		inline void validate() const 
        {
			typename TExtString::TPageFrame &pf = extString->getSharedPage(pageNo, prefetch);
            const_cast<TIterator*>(this)->dirty = pf.dirty;
			const_cast<TIterator*>(this)->begin = pf.begin;
		}

        inline void invalidate(int _prefetch = 0) const 
        {
            if (begin) 
            {
                const_cast<TIterator*>(this)->begin = NULL;
				extString->releasePage(pageNo, (prefetch != 0) || (_prefetch != 0));
                const_cast<TIterator*>(this)->prefetch = _prefetch;
            }
		}

		inline TValue & operator* () const 
        {
			if (!begin) validate();
            // synchronize PageFrame dirty flag on dirty false->true change
            if (!dirty) {
                const_cast<TIterator*>(this)->dirty = true;
    			extString->getPage(pageNo).dirty = true;
            }
			return const_cast<TIterator*>(this)->begin[pageOfs];
		}
/*    
		inline TValue const & operator* () const {
			if (!begin) validate();
			return begin[pageOfs];
		}
*/    
		inline TIterator& operator++ () {
			if (++pageOfs == PAGESIZE) {
				invalidate(1);
				pageOfs = 0;
				++pageNo;
			}
			return *this;
		}

		inline TIterator operator++ (int) {
			TIterator before = *this;
			if (++pageOfs == PAGESIZE) {
				invalidate(1);
				pageOfs = 0;
				++pageNo;
			}
			return before;
		}

		inline TIterator& operator-- () {
			if (pageOfs)
				--pageOfs;
			else {
				invalidate(-1);
				pageOfs = PAGESIZE - 1;
				--pageNo;
			}
			return *this;
		}

		inline TIterator operator-- (int) {
			TIterator before = *this;
			if (pageOfs)
				--pageOfs;
			else {
				invalidate(-1);
				pageOfs = PAGESIZE - 1;
				--pageNo;
			}
			return before;
		}

		inline bool operator== (const TIterator &I) const {
			SEQAN_ASSERT_EQ(extString, I.extString);
			return pageNo == I.pageNo && pageOfs == I.pageOfs;
		}

		inline bool operator!= (const TIterator &I) const {
			SEQAN_ASSERT_EQ(extString, I.extString);
			return pageNo != I.pageNo || pageOfs != I.pageOfs;
		}

		inline bool operator< (const TIterator &I) const {
			SEQAN_ASSERT_EQ(extString, I.extString);
			return pageNo < I.pageNo || (pageNo == I.pageNo && pageOfs < I.pageOfs);
		}
    };


	//////////////////////////////////////////////////////////////////////////////
	// const forward External String iterator
	template < typename TExtString >
    struct ExtStringFwdConstIterator
	{
//IOREV _nottested_ 
		typedef ExtStringFwdConstIterator				TIterator;
        typedef ExtStringIterator<TExtString>			TStdIterator;
        typedef ExtStringConstIterator<TExtString>		TStdConstIterator;
        typedef ExtStringFwdIterator<TExtString>		TFwdIterator;

		typedef typename Value<TExtString>::Type		TValue;
		typedef typename Size<TExtString>::Type			TSize;
		typedef typename Difference<TExtString>::Type	TDifference;
		typedef typename TExtString::TVolatilePtr	    TVolatilePtr;

		enum { PAGESIZE = TExtString::PAGESIZE };

		
		TExtString		*extString;

		int     		pageNo;
		unsigned		pageOfs;
        int             prefetch;   // -n .. prefetch n pages downwards, n .. prefetch n pages upwards, 0 .. disabled
		TVolatilePtr	begin;
		

        ExtStringFwdConstIterator():
			extString(NULL),
			pageNo(0),
			pageOfs(0),
            prefetch(0),
			begin(NULL) {}

		ExtStringFwdConstIterator(const TIterator &I):
			extString(I.extString),
			pageNo(I.pageNo),
			pageOfs(I.pageOfs),
            prefetch(I.prefetch),
			begin(NULL) {}

		ExtStringFwdConstIterator(const TFwdIterator &I):
			extString(I.extString),
			pageNo(I.pageNo),
			pageOfs(I.pageOfs),
            prefetch(I.prefetch),
			begin(NULL) {}

		~ExtStringFwdConstIterator() {
			invalidate();
		}

	    ExtStringFwdConstIterator(TExtString *_extString, TSize _offset):
			extString(_extString),
			pageNo(_offset / PAGESIZE),
			pageOfs(_offset % PAGESIZE),
            prefetch(0),
			begin(NULL) {}

		ExtStringFwdConstIterator(TExtString *_extString, TSize _pageNo, TSize _pageOfs):
			extString(_extString),
			pageNo(_pageNo),
			pageOfs(_pageOfs),
            prefetch(0),
			begin(NULL) {}

		//////////////////////////////////////////////////////////////////////////////
		// iterator conversion interface

		ExtStringFwdConstIterator(TStdIterator &I):
			extString(I.extString),
			pageNo(I.offset / PAGESIZE),
            pageOfs(I.offset % PAGESIZE),
            prefetch(0),
			begin(NULL) {}

        ExtStringFwdConstIterator(TStdConstIterator const &I):
			extString(I.extString),
			pageNo(I.offset / PAGESIZE),
            pageOfs(I.offset % PAGESIZE),
            prefetch(0),
			begin(NULL) {}

		inline TIterator& operator=(TStdIterator const & Right_) {
			invalidate();
			pageNo = Right_.offset / PAGESIZE;
			pageOfs = Right_.offset % PAGESIZE;
            extString = Right_.extString;
			return *this;
		}

		inline TIterator& operator=(TStdConstIterator const & Right_) {
			invalidate();
			pageNo = Right_.offset / PAGESIZE;
			pageOfs = Right_.offset % PAGESIZE;
            extString = Right_.extString;
			return *this;
		}

        inline TSize position() const
        {
            return (TSize)pageNo * (TSize)PAGESIZE + pageOfs;
        }

        inline operator TStdConstIterator() const {
            return TStdConstIterator(extString, position());
        }

		inline TIterator& operator=(TIterator const & Right_) {
			invalidate();
			extString = Right_.extString;
			pageNo = Right_.pageNo;
			pageOfs = Right_.pageOfs;
            prefetch = Right_.prefetch;
			return *this;
		}

		inline TIterator& operator=(TFwdIterator const & Right_) {
			invalidate();
			extString = Right_.extString;
			pageNo = Right_.pageNo;
			pageOfs = Right_.pageOfs;
            prefetch = Right_.prefetch;
			return *this;
		}

		//////////////////////////////////////////////////////////////////////////////
		// iterator arithmetic

		inline TDifference operator- (const TIterator &I) const
        {
            return position() - I.position();
		}
		
		inline TIterator operator- (TDifference delta) const {
			TDifference dPNo  = delta / PAGESIZE;
			TDifference dPOfs = delta % PAGESIZE;
			if (pageOfs >= dPOfs)
				return TIterator(extString, pageNo - dPNo, pageOfs - dPOfs);
			else
				return TIterator(extString, pageNo - dPNo - 1, PAGESIZE + pageOfs - dPOfs);
		}
		
		inline TIterator& operator-= (TDifference delta) {
			TDifference dPNo  = delta / PAGESIZE;
			TDifference dPOfs = delta % PAGESIZE;
			if (pageOfs < dPOfs) {
				++dPNo;
				pageOfs = PAGESIZE + pageOfs - dPOfs;
			} else
				pageOfs -= dPOfs;
			if (dPNo) invalidate(0);
			pageNo -= dPNo;
			return *this;
		}

		inline TIterator operator+ (TDifference delta) const {
			TDifference dPNo  = delta / PAGESIZE;
			TDifference nPOfs = pageOfs + delta % PAGESIZE;
			if (nPOfs < PAGESIZE)
				return TIterator(extString, pageNo + dPNo, nPOfs);
			else
				return TIterator(extString, pageNo + dPNo + 1, nPOfs - PAGESIZE);
		}
		
		inline TIterator& operator+= (TDifference delta) {
			TDifference dPNo  = delta / PAGESIZE;
			TDifference nPOfs = pageOfs + delta % PAGESIZE;
			if (nPOfs >= PAGESIZE) {
				++dPNo;
				nPOfs -= PAGESIZE;
			}
			if (dPNo) invalidate(0);
			pageNo += dPNo;
			pageOfs = nPOfs;
			return *this;
		}
		
		inline void validate() const {
			typename TExtString::TPageFrame &pf = extString->getSharedPage(pageNo, prefetch);
			const_cast<TIterator*>(this)->begin = pf.begin;
		}

        inline void invalidate(int _prefetch = 0) const {
            if (begin) {
                const_cast<TIterator*>(this)->begin = NULL;
				extString->releasePage(pageNo, (prefetch != 0) || (_prefetch != 0));
                const_cast<TIterator*>(this)->prefetch = _prefetch;
            }
		}

		inline TValue const & operator* () const {
			if (!begin) validate();
			return begin[pageOfs];
		}
    
		inline TIterator& operator++ () {
			if (++pageOfs == PAGESIZE) {
				invalidate(1);
				pageOfs = 0;
				++pageNo;
			}
			return *this;
		}

		inline TIterator operator++ (int) {
			TIterator before = *this;
			if (++pageOfs == PAGESIZE) {
				invalidate(1);
				pageOfs = 0;
				++pageNo;
			}
			return before;
		}

		inline TIterator& operator-- () {
			if (pageOfs)
				--pageOfs;
			else {
				invalidate(-1);
				pageOfs = PAGESIZE - 1;
				--pageNo;
			}
			return *this;
		}

		inline TIterator operator-- (int) {
			TIterator before = *this;
			if (pageOfs)
				--pageOfs;
			else {
				invalidate(-1);
				pageOfs = PAGESIZE - 1;
				--pageNo;
			}
			return before;
		}

		inline bool operator== (const TIterator &I) const {
			SEQAN_ASSERT_EQ(extString, I.extString);
			return pageNo == I.pageNo && pageOfs == I.pageOfs;
		}

		inline bool operator!= (const TIterator &I) const {
			SEQAN_ASSERT_EQ(extString, I.extString);
			return pageNo != I.pageNo || pageOfs != I.pageOfs;
		}

		inline bool operator< (const TIterator &I) const {
			SEQAN_ASSERT_EQ(extString, I.extString);
			return pageNo < I.pageNo || (pageNo == I.pageNo && pageOfs < I.pageOfs);
		}

	};



	//////////////////////////////////////////////////////////////////////////////
	// iterator metafunctions

	template < typename TString >
    struct Container< ExtStringIterator<TString> >			{ typedef TString Type; };
//IOREV
	template < typename TString >
    struct Container< ExtStringConstIterator<TString> >		{ typedef TString Type; };
//IOREV
	template < typename TString >
    struct Container< ExtStringFwdIterator<TString> >		{ typedef TString Type; };
//IOREV
	template < typename TString >
    struct Container< ExtStringFwdConstIterator<TString> >	{ typedef TString Type; };
//IOREV

	template < typename TString >
    struct Value< ExtStringIterator<TString> >				{ typedef typename Value<TString>::Type Type; };
//IOREV
	template < typename TString >
    struct Value< ExtStringConstIterator<TString> >			{ typedef typename Value<TString>::Type Type; };
//IOREV
	template < typename TString >
    struct Value< ExtStringFwdIterator<TString> >			{ typedef typename Value<TString>::Type Type; };
//IOREV
	template < typename TString >
    struct Value< ExtStringFwdConstIterator<TString> >		{ typedef typename Value<TString>::Type Type; };
//IOREV

	template < typename TString >
	struct Reference< ExtStringConstIterator<TString> >:
		public Reference<TString const> {};
//IOREV

	template < typename TString >
	struct Reference< ExtStringFwdConstIterator<TString> >:
		public Reference<TString const> {};
//IOREV

	template < typename TString >
    struct Size< ExtStringIterator<TString> >				{ typedef typename Size<TString>::Type Type; };
//IOREV
	template < typename TString >
    struct Size< ExtStringConstIterator<TString> >			{ typedef typename Size<TString>::Type Type; };
//IOREV
	template < typename TString >
    struct Size< ExtStringFwdIterator<TString> >			{ typedef typename Size<TString>::Type Type; };
//IOREV
	template < typename TString >
    struct Size< ExtStringFwdConstIterator<TString> >		{ typedef typename Size<TString>::Type Type; };
//IOREV

	template < typename TString >
    struct Position< ExtStringIterator<TString> >			{ typedef typename Position<TString>::Type Type; };
//IOREV
	template < typename TString >
    struct Position< ExtStringConstIterator<TString> >		{ typedef typename Position<TString>::Type Type; };
//IOREV
	template < typename TString >
    struct Position< ExtStringFwdIterator<TString> >		{ typedef typename Position<TString>::Type Type; };
//IOREV
	template < typename TString >
    struct Position< ExtStringFwdConstIterator<TString> >	{ typedef typename Position<TString>::Type Type; };
//IOREV

	template < typename TString >
    struct Difference< ExtStringIterator<TString> >			{ typedef typename Difference<TString>::Type Type; };
//IOREV
	template < typename TString >
    struct Difference< ExtStringConstIterator<TString> >	{ typedef typename Difference<TString>::Type Type; };
//IOREV
	template < typename TString >
    struct Difference< ExtStringFwdIterator<TString> >		{ typedef typename Difference<TString>::Type Type; };
//IOREV
	template < typename TString >
    struct Difference< ExtStringFwdConstIterator<TString> > { typedef typename Difference<TString>::Type Type; };
//IOREV


	//////////////////////////////////////////////////////////////////////////////
    // global interface

	template <typename TExtString>
	inline TExtString &	container(ExtStringIterator<TExtString> &it) { return *(it.extString); }
//IOREV
	template <typename TExtString>
	inline TExtString &	container(ExtStringIterator<TExtString> const &it) { return *(it.extString); }
//IOREV

	template <typename TExtString>
	inline TExtString &	container(ExtStringConstIterator<TExtString> &it) { return *(it.extString); }
//IOREV
	template <typename TExtString>
	inline TExtString &	container(ExtStringConstIterator<TExtString> const &it) { return *(it.extString); }
//IOREV

	template <typename TExtString>
	inline TExtString &	container(ExtStringFwdIterator<TExtString> &it) { return *(it.extString); }
//IOREV
	template <typename TExtString>
	inline TExtString &	container(ExtStringFwdIterator<TExtString> const &it) { return *(it.extString); }
//IOREV

	template <typename TExtString>
	inline TExtString &	container(ExtStringFwdConstIterator<TExtString> &it) { return *(it.extString); }
//IOREV
	template <typename TExtString>
	inline TExtString &	container(ExtStringFwdConstIterator<TExtString> const &it) { return *(it.extString); }
//IOREV
//____________________________________________________________________________

	template <typename TExtString>
	inline bool	atBegin(ExtStringIterator<TExtString> &it) { return it.offset == 0; }
//IOREV
	template <typename TExtString>
	inline bool	atBegin(ExtStringIterator<TExtString> const &it) { return it.offset == 0; }
//IOREV

	template <typename TExtString>
	inline bool	atBegin(ExtStringConstIterator<TExtString> &it) { return it.offset == 0; }
//IOREV
	template <typename TExtString>
	inline bool	atBegin(ExtStringConstIterator<TExtString> const &it) { return it.offset == 0; }
//IOREV

	template <typename TExtString>
	inline bool	atBegin(ExtStringFwdIterator<TExtString> &it) { 
//IOREV
		return it.pageNo == 0 && it.pageOfs == 0;
	}
	template <typename TExtString>
	inline bool	atBegin(ExtStringFwdIterator<TExtString> const &it) {
//IOREV
		return it.pageNo == 0 && it.pageOfs == 0;
	}

	template <typename TExtString>
	inline bool	atBegin(ExtStringFwdConstIterator<TExtString> &it) { 
//IOREV
		return it.pageNo == 0 && it.pageOfs == 0;
	}
	template <typename TExtString>
	inline bool	atBegin(ExtStringFwdConstIterator<TExtString> const &it) {
//IOREV
		return it.pageNo == 0 && it.pageOfs == 0;
	}
//____________________________________________________________________________

	template <typename TExtString>
	inline bool	atEnd(ExtStringIterator<TExtString> &it) { return it.offset == it.extString->data_size; }
//IOREV
	template <typename TExtString>
	inline bool	atEnd(ExtStringIterator<TExtString> const &it) { return it.offset == it.extString->data_size; }
//IOREV

	template <typename TExtString>
	inline bool	atEnd(ExtStringConstIterator<TExtString> &it) { return it.offset == it.extString->data_size; }
//IOREV
	template <typename TExtString>
	inline bool	atEnd(ExtStringConstIterator<TExtString> const &it) { return it.offset == it.extString->data_size; }
//IOREV

	template <typename TExtString>
	inline bool	atEnd(ExtStringFwdIterator<TExtString> &it) { 
//IOREV
		return TExtString::PAGESIZE * it.pageNo + it.pageOfs == it.extString->data_size;
	}
	template <typename TExtString>
	inline bool	atEnd(ExtStringFwdIterator<TExtString> const &it) {
//IOREV
		return TExtString::PAGESIZE * it.pageNo + it.pageOfs == it.extString->data_size;
	}

	template <typename TExtString>
	inline bool	atEnd(ExtStringFwdConstIterator<TExtString> &it) { 
//IOREV
		return TExtString::PAGESIZE * it.pageNo + it.pageOfs == it.extString->data_size;
	}
	template <typename TExtString>
	inline bool	atEnd(ExtStringFwdConstIterator<TExtString> const &it) {
//IOREV
		return TExtString::PAGESIZE * it.pageNo + it.pageOfs == it.extString->data_size;
	}


	
	
	
	
	//////////////////////////////////////////////////////////////////////////////
    // External String
    //////////////////////////////////////////////////////////////////////////////

    template < typename TValue,
               typename TConfig >
	class String<TValue, External<TConfig> >
	{
//IOREV _doc_ contains TODOs by holtgrew
	public:
        enum { FRAMES    = TConfig::FRAMES,
               PAGESIZE = TConfig::PAGESIZE };

        typedef typename TConfig::TFile                                 TFile;
        typedef typename TConfig::TSize                                 TSize;

		typedef String<int>                                             TPageTable;
		typedef Buffer<TValue, PageFrame<TFile, Fixed<PAGESIZE> > >    TPageFrame;
		typedef PageContainer<TPageFrame, FRAMES>                       TCache;
		typedef VolatilePtr<TValue>                                     TVolatilePtr;

		TPageTable			pager;
		TCache				cache;
		TFile				file;
        bool                _temporary, _ownFile;
		TSize				data_size;
        int                 lastDiskPage;       // the last page on disk and in mem 
        unsigned            lastDiskPageSize;   // can be smaller than PAGESIZE

		String(TSize size = 0):
            file(NULL),
			data_size(0)
        {
            _temporary = true;
            _ownFile = false;
            lastDiskPage = 0;       // actually, these values need not to be initialized
            lastDiskPageSize = 0;   // here, because of "write before read"

			resize(*this, size);
        }

/*
	private:	// making these C'tors private clashes with Holder<external String>
		String(String &) {}
		String(String const &) {}

	public:
*/
		String(String &) { SEQAN_ASSERT_FAIL("Aborted attempt to copy a String<..,External<..> >"); }  // TODO(holtgrew): Actually, this should be an ABORT
		String(String const &) { SEQAN_ASSERT_FAIL("Aborted attempt to copy a String<..,External<..> >"); }  // TODO(holtgrew): Actually, this should be an ABORT

/**
.Memfunc.External String#String:
..class:Spec.External String
..summary:Constructor
..signature:String<TValue, External<TConfig> > ()
..signature:String<TValue, External<TConfig> > (file)
..signature:String<TValue, External<TConfig> > (fileName[, openMode])
..param.file:The @Spec.External String@ will use the file associated with $file$.
...remarks:You must ensure that $file$ is open, as the string won't call @Function.open@ and @Function.close@ on it.
...type:Class.File
..param.fileName:The @Spec.External String@ will @Function.open@ the file with the path $fileName$.
..param.openMode:File mode for @Function.open@.
..remarks:When a file or file name is given, this file will be used for the @Spec.External String@.
If the file exists, this file will be used and determines the strings length and content.
If the file doesn't exist, a new and empty file will be created and used for the string.
In both cases, the string won't delete the file in the destructor.
..remarks:
When no file is given (default c'tor) the string will be empty and no file is used until the string 
needs to swap out page frames. Then a temporary file will be used which will be deleted when the string is destroyed.
..remarks:
Instead of giving $file$ or $fileName$ to the constructor, you could also use the default constructor and call @Function.open@
or @Function.openTemp@ afterwards to reach the same behaviour.
*/
		String(TFile &_file)
        {
			open(*this, _file);
        }

		String(const char *fileName, int openMode = DefaultOpenMode<TFile>::VALUE):
			file(NULL)
        {
			open(*this, fileName, openMode);
        }

		~String() 
		{
			close(*this);
		}

		inline TValue & operator[] (TSize offset) {
			TPageFrame &pf = getPage(offset / PAGESIZE);
			pf.dirty = true;
			return pf[offset % PAGESIZE];
		}

		inline TValue const & operator[] (TSize offset) const {
			return const_cast<String*>(this)->getPage(offset / PAGESIZE)[offset % PAGESIZE];
		}

	    template <typename TSource>
	    inline String & operator= (TSource const & source)
	    {
		    assign(*this, source);
		    return *this;
	    }

	    inline String & operator= (String const & source)
	    {
		    assign(*this, source);
		    return *this;
	    }

        inline operator bool() 
		{
            return file;
        }

		//////////////////////////////////////////////////////////////////////////////
		// swapping interface

		// when a page has to be swapped out and file is not open, open a temporary file
		inline void _ensureFileIsOpen() 
		{
			if (!file) 
			{
				_temporary = true;
				if (!(_ownFile = openTemp(file)))
					::std::cerr << "External String couldn't open temporary file" << ::std::endl;
			}
		}

		// for debugging
        void _dumpCache() 
		{
            for(int i = 0; i < length(cache); ++i) 
			{
                TPageFrame &pf = cache[i];
                ::std::cerr << "[" << pf.pageNo << "]";
                if (pf.dirty)
                    ::std::cerr << "*";
                else
                    ::std::cerr << " ";

                if (pf.status == READY)
                    ::std::cerr << "   ";
                else
                    ::std::cerr << ".  ";
            }
            ::std::cerr << ::std::endl;
        }

        // return a priority for a page frame (the higher is more persistent)
        inline typename TPageFrame::Priority getPriority(int /*pageNo*/) const 
		{
/*            if (keepFirst && pageNo < (int)(length(cache)) - 10) // save 1 for random access
                return TPageFrame::PERMANENT_LEVEL;
            else*/
                return TPageFrame::NORMAL_LEVEL;
        }

		// write page to disk if dirty and remove from page table now or after finishing IO
		inline void flush(TPageFrame &pf) 
		{
            if (pf.status == READY && pf.dirty) {    // write if dirty and not i/o transferring
				nukeCopies(pf.begin);				            // proceeding writes should wait and set dirty bit

                if (pf.priority > TPageFrame::NORMAL_LEVEL && pf.priority <= TPageFrame::ITERATOR_LEVEL)
					cache.upgrade(pf, TPageFrame::PREFETCH_LEVEL);

				_ensureFileIsOpen();
				if (pf.pageNo != (int)(data_size / (TSize)PAGESIZE))
    				writePage(pf, pf.pageNo, file);
                else {
                    lastDiskPage = data_size / PAGESIZE;
                    lastDiskPageSize = data_size % PAGESIZE;
				    writeLastPage(pf, pf.pageNo, file, lastDiskPageSize);
                }
                pf.dataStatus = TPageFrame::ON_DISK;
			}
		}

		// write page synchronously to disk if dirty and remove from page table
		inline void swapOutAndWait(TPageFrame &pf) 
		{
			nukeCopies(pf.begin);      				// proceeding writes should wait and set dirty bit

            if (pf.status != READY)
            {
				pager[pf.pageNo] = TPageFrame::ON_DISK;		// page is not dirty and on disk
				bool waitResult = waitFor(pf);              // after finishing I/O transfer

                // TODO(weese): Throw an I/O exception
                if (!waitResult)
                    SEQAN_FAIL("%s operation could not be completed: \"%s\"", _pageFrameStatusString(pf), strerror(errno));

                pf.pageNo = -1;                             // cut back link
                return;
            }

			if (pf.dirty) {                                 // write if dirty
				_ensureFileIsOpen();
                if (pf.pageNo != (int)(data_size / (TSize)PAGESIZE)) {
    				writePage(pf, pf.pageNo, file);
                    if (pf.pageNo >= lastDiskPage)
                        lastDiskPage = -1;       			// make lastDiskPage(Size) invalid because file size is aligned
                } else {
				    writeLastPage(pf, pf.pageNo, file, data_size % PAGESIZE);
                    lastDiskPage = data_size / PAGESIZE;
                    lastDiskPageSize = data_size % PAGESIZE;
                }
				pager[pf.pageNo] = TPageFrame::ON_DISK;		// page is marked to be on disk
				bool waitResult = waitFor(pf);              // after finishing I/O transfer

                // TODO(weese): Throw an I/O exception
                if (!waitResult)
                    SEQAN_FAIL("%s operation could not be completed: \"%s\"", _pageFrameStatusString(pf), strerror(errno));
			} else
				pager[pf.pageNo] = pf.dataStatus;			// restore original data status

            pf.pageNo = -1;                                 // cut back link
		}

		struct testIODone : public ::std::unary_function<TPageFrame&,bool> 
		{
			String &me;
			testIODone(String &_me): me(_me) {}

			inline bool operator() (TPageFrame &pf)
            {
                PageFrameStatus oldStatus = pf.status;
                bool inProgress;
                bool waitResult = waitFor(pf, 0, inProgress);

                // TODO(weese): Throw an I/O exception
                if (!waitResult)
                    SEQAN_FAIL("%s operation could not be completed: \"%s\"", _pageFrameStatusString(pf), strerror(errno));

                if (!inProgress && (oldStatus != READY))
                {
                    if (pf.pageNo >= me.lastDiskPage)
                        me.lastDiskPage = -1;    // make lastDiskPage(Size) invalid because file size is aligned
                }
                return !inProgress;
			}
		};

        inline TPageFrame &getPage(
            int pageNo,
            typename TPageFrame::Priority maxLevel,
            typename TPageFrame::Priority newLevel,
            int prefetchPages)
        {
			int frameNo = pager[pageNo];
			if (frameNo >= 0)					// cache hit
            {
				TPageFrame &pf = cache[frameNo];
				cache.upgrade(
                    pf, 
                    _max(pf.priority, newLevel));    		// update lru order

                PageFrameStatus oldStatus = pf.status;
				bool waitResult = waitFor(pf);              // wait for I/O transfer to complete

                // TODO(weese): Throw an I/O exception
                if (!waitResult)
                    SEQAN_FAIL("%s operation could not be completed: \"%s\"", _pageFrameStatusString(pf), strerror(errno));

				if (oldStatus != READY)
                    if (pf.pageNo >= lastDiskPage)
                        lastDiskPage = -1;       			// make lastDiskPage(Size) invalid because file size is aligned

                if (prefetchPages > 0) prefetch(pageNo + 1, pageNo + 1 + prefetchPages, frameNo);
                else if (prefetchPages < 0) prefetch(pageNo + prefetchPages, pageNo, frameNo);

				return pf;

			} else {							// cache miss

				typename TPageFrame::DataStatus dataStatus = static_cast<typename TPageFrame::DataStatus>(frameNo);
				frameNo = cache.mru(testIODone(*this), maxLevel);   // try to get an undirty and READY pageframe
				if (frameNo < 0)							// if there is none,
					frameNo = cache.mruDirty();				// get the most recently used dirty frame
				TPageFrame &pf = cache[frameNo];

				// *** frame is choosen ***

				if (pf.begin)
					swapOutAndWait(pf);						// write synchronously to disk, if page is dirty
				else
					allocPage(pf, file);                    // allocate memory if page is virgin

				// *** frame is free now ***

				pf.dataStatus = dataStatus;
				if (dataStatus == TPageFrame::ON_DISK)
				{
                    if (pageNo != lastDiskPage)
					    readPage(pageNo, pf, file);
                    else
                        readLastPage(pageNo, pf, file, lastDiskPageSize);
				}
				pager[pageNo] = frameNo;					// assign new page to page table
				pf.pageNo = pageNo;							// set back link
				cache.upgrade(
                    pf,
                    _max(getPriority(pageNo), newLevel));    // update lru order

                if (prefetchPages > 0) prefetch(pageNo + 1, pageNo + 1 + prefetchPages, frameNo);
                else if (prefetchPages < 0) prefetch(pageNo + prefetchPages, pageNo, frameNo);
                
				bool waitResult = waitFor(pf);              // wait for I/O transfer to complete

                // TODO(weese): Throw an I/O exception
                if (!waitResult)
                    SEQAN_FAIL("%s operation could not be completed: \"%s\"", _pageFrameStatusString(pf), strerror(errno));

				return pf;
			}
		}
        
		inline TPageFrame &getPage(int pageNo)
        {
			return getPage(pageNo, TPageFrame::NORMAL_LEVEL, TPageFrame::NORMAL_LEVEL, 0);
		}

        // prefetch is non-blocking and should speed up swapping
		inline void prefetch(int pageBegin, int pageEnd, int except = -1) 
		{
            if (!file) return;
            if (pageBegin < 0)					pageBegin = 0;
            if (pageEnd >= (int)length(pager))	pageEnd = (int)length(pager) - 1;
            for(int pageNo = pageBegin; pageNo < pageEnd; ++pageNo) {
			    int frameNo = pager[pageNo];
				typename TPageFrame::DataStatus dataStatus = static_cast<typename TPageFrame::DataStatus>(frameNo);
                if (dataStatus == TPageFrame::ON_DISK &&             // prefetch only if page is on disk
                    pageNo != lastDiskPage)                         // reading the last page is blocking
                {   
				    frameNo = cache.mru(
                        testIODone(*this),
                        TPageFrame::NORMAL_LEVEL);                   // choose undirty and ready page

                    if (frameNo < 0 || frameNo == except) return;   // no lowlevel-page left for prefetching
				    TPageFrame &pf = cache[frameNo];
                    #ifdef SEQAN_VERBOSE
						::std::cerr << "prefetch: page " << pageNo << ::std::endl;
                    #endif

                    // *** frame is choosen ***

				    if (pf.begin)
					    swapOutAndWait(pf);						    // write synchronously to disk, if page is dirty
				    else
					    allocPage(pf, file);                        // allocate memory if page is virgin

    				// *** frame is free now ***

    				pf.dataStatus = dataStatus;
                    readPage(pageNo, pf, file);
				    pager[pageNo] = frameNo;					    // assign new page to page table
    				pf.pageNo = pageNo;							    // set back link
                    cache.upgrade(pf, TPageFrame::PREFETCH_LEVEL);  // update lru order
                }
            }
		}
		
	    template < typename T >
		inline static int _prefetchIffAsync(int /*prefetchPages*/, T const &) {
			return 0;
		}
		
		template < typename TSpec >
		inline static int _prefetchIffAsync(int prefetchPages, File<Async<TSpec> > const &) {
			return prefetchPages;
		}

		inline TPageFrame &getSharedPage(int pageNo, int prefetchPages = 0) 
		{
			return getPage(
                pageNo, 
                TPageFrame::PREFETCH_LEVEL, 
                TPageFrame::ITERATOR_LEVEL,
                _prefetchIffAsync(prefetchPages, file));
		}

		inline void releasePage(int pageNo, bool writeThrough = false) 
		{
			int frameNo = pager[pageNo];
			if (frameNo >= 0)							// release only cached pages
			{								        
				TPageFrame &pf = cache[frameNo];
				if (pf.begin.isLonely() && pf.priority <= TPageFrame::ITERATOR_LEVEL) 
				{
					cache.upgrade(pf, _max(getPriority(pageNo), TPageFrame::NORMAL_LEVEL));
                    if (writeThrough) 
					{
                        #ifdef SEQAN_VERBOSE
                            if (pf.dirty)
								::std::cerr << "writeThrough: page " << pageNo << ::std::endl;
                        #endif
					    flush(pf);							        // write if dirty
                    }
				}
			}
		}
        
        inline void rename(unsigned frameNo) 
		{
			TPageFrame &pf = cache[frameNo];
            cache.rename(frameNo);                                  // update lru entry
            if (pf.pageNo >= 0)
                pager[pf.pageNo] = frameNo;					        // update back link
        }

        // change the number of in-mem pageframes
        // more pages mean less swapping, 
        // less pages mean more free mem
        inline void resizeCache(unsigned newFrames) 
		{
            unsigned oldFrames = length(cache);
            if (data_size)
                newFrames = _min(newFrames, (unsigned) enclosingBlocks(data_size, (unsigned)PAGESIZE));
            if (newFrames < oldFrames) {
                flush(*this);
                for(unsigned i = newFrames; i < oldFrames; ++i) {
    			    int frameNo = cache.mruDirty();             // get the most recently used frame (can be dirty)
                    if (frameNo < 0) break;
				    TPageFrame &pf = cache[frameNo];

				    // *** frame is choosen ***

                    if (pf.begin) {
					    swapOutAndWait(pf);						// write synchronously to disk, if page is dirty
        				freePage(pf, file);                     // free memory
                    }

                    cache.erase(frameNo);                       // erase page frame from cache

                    for(unsigned j = frameNo; j < length(cache); ++j)
                        rename(j);                              // update remaining pages
                }
            } else if (oldFrames < newFrames) {
                resize(cache, newFrames);
            }
        }

    };


    //////////////////////////////////////////////////////////////////////////////
    // meta-function interface

    template < typename TValue, typename TConfig >
    struct Size< String<TValue, External<TConfig> > >
    {
//IOREV
        typedef typename String<TValue, External<TConfig> >::TSize Type;
    };

    template < typename TValue, typename TConfig >
    struct Difference< String<TValue, External<TConfig> > >
    {
//IOREV
		typedef typename MakeSigned_<typename String<TValue, External<TConfig> >::TSize>::Type Type;
    };

    template < typename TValue, typename TConfig, typename TSpec >
    struct Iterator< String<TValue, External<TConfig> > const, Tag<TSpec> const >
    {
//IOREV
        typedef ExtStringFwdConstIterator< String<TValue, External<TConfig> > > Type;
    };

    template < typename TValue, typename TConfig, typename TSpec >
    struct Iterator< String<TValue, External<TConfig> >, Tag<TSpec> const > 
    {
//IOREV
        typedef ExtStringFwdIterator< String<TValue, External<TConfig> > > Type;
    };
//____________________________________________________________________________

    template < typename TValue, typename TConfig >
	struct DefaultOverflowExplicit<String<TValue, External<TConfig> > >
	{
//IOREV
		typedef Generous Type;
	};

    template < typename TValue, typename TConfig >
	struct DefaultOverflowImplicit<String<TValue, External<TConfig> > >
	{
//IOREV
		typedef Generous Type;
	};
//____________________________________________________________________________

/*
    template < typename TValue, typename TConfig >
	struct DefaultIteratorSpec< String<TValue, External<TConfig> > > {
		typedef Standard Type;
	};
	
    template < typename TValue, typename TConfig >
	struct DefaultIteratorSpec< String<TValue, External<TConfig> > const > {
		typedef Standard Type;
	};
*/

    template < typename TValue, typename TConfig >
	struct AllowsFastRandomAccess< String<TValue, External<TConfig> > >
	{
//IOREV
		typedef False Type;
		enum { VALUE = false };
	};


	//////////////////////////////////////////////////////////////////////////////
    // global interface

//____________________________________________________________________________

    template < typename TValue, typename TConfig >
    inline void 
    clear(String<TValue, External<TConfig> > &me) 
	{
//IOREV
		clear(me.pager);
        resize(me, 0);
    }
//____________________________________________________________________________

	// wait until IO of every page is finished
    template < typename TValue, typename TConfig >
	inline void 
	waitForAll(String<TValue, External<TConfig> > &me)
	{
//IOREV _nodoc_
		typedef typename String<TValue, External<TConfig> >::TCache	TCache;
		typedef typename Iterator<TCache, Standard>::Type			TIter;

		TIter f = begin(me.cache, Standard());
		TIter fEnd = end(me.cache, Standard());

		for(; f != fEnd ; ++f)
        {
            bool waitResult = waitFor(*f);              // wait for I/O transfer to complete

            if (!waitResult)
                SEQAN_FAIL("%s operation could not be completed: \"%s\"", _pageFrameStatusString(*f), strerror(errno));
        }
	}
	
/**
.Function.flush:
..signature:flush(string)
..param.string:An external string. All dirty pages are flushed to disk.
...type:Spec.External String
..include:seqan/file.h
*/
    template < typename TValue, typename TConfig >
    inline void 
	flush(String<TValue, External<TConfig> > &me) 
	{
//IOREV _doc_
		typedef typename String<TValue, External<TConfig> >::TCache	TCache;
		typedef typename Iterator<TCache, Standard>::Type			TIter;

		// write all dirty pages to disk
		if (me.file)
		{
			TIter f = begin(me.cache, Standard());
			TIter fEnd = end(me.cache, Standard());

			for(; f != fEnd ; ++f)
				if ((*f).begin) me.flush(*f);
			waitForAll(me);
		}
    }

	// cancel all transactions
    template < typename TValue, typename TConfig >
	inline void 
	cancel(String<TValue, External<TConfig> > &me)
	{
//IOREV _nodoc_
		typedef typename String<TValue, External<TConfig> >::TCache	TCache;
		typedef typename Iterator<TCache, Standard>::Type			TIter;

		if (me.file) 
		{
			TIter f = begin(me.cache, Standard());
			TIter fEnd = end(me.cache, Standard());

			for(; f != fEnd ; ++f)
				if ((*f).begin) cancel(*f, me.file);
		}
	}

	// cancel all transactions and free allocated pages
    template < typename TValue, typename TConfig >
	inline void 
	cancelAndFree(String<TValue, External<TConfig> > &me)
	{
//IOREV _nodoc_
		typedef String<TValue, External<TConfig> >					TExtString;
		typedef typename TExtString::TPageFrame						TPageFrame;
		typedef typename String<TValue, External<TConfig> >::TCache	TCache;
		typedef typename Iterator<TCache, Standard>::Type			TIter;

        TIter f = begin(me.cache, Standard());
        TIter fEnd = end(me.cache, Standard());

        for(; f != fEnd ; ++f) 
        {
            if ((me.file) && (*f).begin) cancel(*f, me.file);
            if ((*f).pageNo >= 0) 
            {
                me.pager[(*f).pageNo] = (*f).dataStatus;
                (*f).pageNo = TPageFrame::UNINITIALIZED;
            }
//			::std::cerr << *f << ::std::endl;
            if ((*f).begin)
                freePage(*f, me.file);
  		}
	}

	// flush and free all allocated pages
    template < typename TValue, typename TConfig >
	inline void 
	flushAndFree(String<TValue, External<TConfig> > &me)
	{
//IOREV _nodoc_
		typedef String<TValue, External<TConfig> >			TExtString;
		typedef typename TExtString::TPageFrame				TPageFrame;
		typedef typename TExtString::TCache					TCache;
		typedef typename Iterator<TCache, Standard>::Type	TIter;

        flush(me);

        TIter f = begin(me.cache, Standard());
        TIter fEnd = end(me.cache, Standard());

        for(; f != fEnd ; ++f) 
        {
            if ((*f).pageNo >= 0) 
            {
                me.pager[(*f).pageNo] = (*f).dataStatus;
                (*f).pageNo = TPageFrame::UNINITIALIZED;
            }
//			::std::cerr << *f << ::std::endl;
            if ((*f).begin) freePage(*f, me.file);
        }
	}
//____________________________________________________________________________
/**
.Function.open:
..signature:open(string, fileName[, openMode]))
..param.string:A persistent string, e.g. a @Spec.External String@ or @Spec.MMap String@.
...type:Spec.External String
..include:seqan/file.h
*/

	template < typename TValue, typename TConfig >
    inline bool 
    open(String<TValue, External<TConfig> > &me, const char *fileName, int openMode) 
	{
//IOREV _doc_
		typedef String<TValue, External<TConfig> >			TExtString;
		typedef typename TExtString::TPageFrame				TPageFrame;

		me._temporary = false;
		if ((me._ownFile = open(me.file, fileName, openMode)))
            me.data_size = size(me.file) / sizeof(TValue);
        else
            me.data_size = 0;

		resize(me.pager, enclosingBlocks(me.data_size, 
			(unsigned)me.PAGESIZE), (me.data_size)? 
				TPageFrame::ON_DISK: 
				TPageFrame::UNINITIALIZED);

        me.lastDiskPage = me.data_size / me.PAGESIZE;
        me.lastDiskPageSize = me.data_size % me.PAGESIZE;
		return me._ownFile;
    }

	template < typename TValue, typename TConfig >
    inline bool 
    open(String<TValue, External<TConfig> > &me, const char *fileName) 
	{
//IOREV _doc_
		typedef String<TValue, External<TConfig> >	TExtString;
		typedef typename TExtString::TFile			TFile;

		return open(me, fileName, DefaultOpenMode<TFile>::VALUE);
    }

	template < typename TValue, typename TConfig >
    inline bool 
    open(String<TValue, External<TConfig> > &me, typename TConfig::TFile file) 
	{
//IOREV _doc_
		typedef String<TValue, External<TConfig> >	TExtString;
		typedef typename TExtString::TPageFrame		TPageFrame;

		me.file = file;
        me._temporary = false;
        me._ownFile = false;
        if (me.file)
            me.data_size = size(me.file) / sizeof(TValue);
        else
            me.data_size = 0;

		resize(me.pager, enclosingBlocks(me.data_size, 
			(unsigned)me.PAGESIZE), (me.data_size)? 
				TPageFrame::ON_DISK: 
				TPageFrame::UNINITIALIZED);

        me.lastDiskPage = me.data_size / me.PAGESIZE;
        me.lastDiskPageSize = me.data_size % me.PAGESIZE;
		return me._file;
    }

/**
.Function.openTemp:
..signature:openTemp(string)
..param.string:An external string.
...type:Spec.External String
..include:seqan/file.h
*/
	template < typename TValue, typename TConfig >
    inline bool 
    openTemp(String<TValue, External<TConfig> > &me) 
	{
//IOREV _doc_
        me._temporary = true;
        me.lastDiskPage = 0;
        me.lastDiskPageSize = 0;
		clear(me.pager);
		return me._ownFile = openTemp(me.file);
    }
//____________________________________________________________________________

	template < typename TValue, typename TConfig >
    inline bool 
    save(String<TValue, External<TConfig> > const &/*me*/, const char * /*fileName*/, int /*openMode*/) {
//IOREV _nodoc_ shouldn't we flush here? in case of abnormal termination...
		// External Strings are persistent, thus there is no need to save them
		//ExtStringsDontNeedToBeSaved error;
		return true;
	}

	template < typename TValue, typename TConfig >
    inline bool 
    save(String<TValue, External<TConfig> > const &/*me*/, const char * /*fileName*/) {
//IOREV _nodoc_ shouldn't we flush here? in case of abnormal termination...
		// External Strings are persistent, thus there is no need to save them
		//ExtStringsDontNeedToBeSaved error;
		return true;
	}

	template < typename TValue, typename TConfig >
    inline bool 
    save(String<TValue, External<TConfig> > const &/*me*/, typename TConfig::TFile /*file*/) {
//IOREV _nodoc_ shouldn't we flush here? in case of abnormal termination...
		// External Strings are persistent, thus there is no need to save them
		//ExtStringsDontNeedToBeSaved error;
		return true;
	}
//____________________________________________________________________________
/**
.Function.close:
..signature:close(string)
..param.string:An external string.
...type:Spec.External String
..include:seqan/file.h
*/
	template < typename TValue, typename TConfig >
    inline bool 
    close(String<TValue, External<TConfig> > &me) 
	{
//IOREV _doc_
		// close associated file
		if (me._temporary)
			cancelAndFree(me);
		else
			flushAndFree(me);
		clear(me.pager);

		if (me._ownFile) 
		{
			me._ownFile = false;
			return close(me.file);
		} 
		else
			return true;
    }
//____________________________________________________________________________
///.Function.length.param.object.type:Class.Shape

	template < typename TValue, typename TConfig >
    inline typename Size< String<TValue, External<TConfig> > >::Type
    length(String<TValue, External<TConfig> > const &me)
    {
//IOREV
        return me.data_size;
    }

	template < typename TValue, typename TConfig >
    inline typename Size< String<TValue, External<TConfig> > >::Type
    capacity(String<TValue, External<TConfig> > const &me)
    {
//IOREV
		typedef typename Size< String<TValue, External<TConfig> > >::Type TSize;
        return (TSize)capacity(me.pager) * (TSize)me.PAGESIZE;
    }
//____________________________________________________________________________

    template < typename TValue, typename TConfig, typename TNewSize, typename TExpand >
    inline typename Size< String<TValue, External<TConfig> > >::Type
    resize(
	    String<TValue, External<TConfig> > &me,
		TNewSize new_length,
		Tag<TExpand> expand)
	{
//IOREV
		typedef String<TValue, External<TConfig> >	TString;
		typedef typename TString::TPageFrame		TPageFrame;
		typedef typename Size<TString>::Type		TSize;

		resize(me.pager, enclosingBlocks(new_length, (unsigned)me.PAGESIZE), TPageFrame::UNINITIALIZED, expand);
        if ((TSize)new_length < me.data_size && me.file) 
		{
			// wait for all pending transfers
            waitForAll(me);

			// before shrinking the file size
            resize(me.file, (TSize)new_length * (TSize)sizeof(TValue));
            me.lastDiskPage = new_length / me.PAGESIZE;
            me.lastDiskPageSize = new_length % me.PAGESIZE;
        }
		me.data_size = new_length;
		return length(me);
	}
//____________________________________________________________________________

    template < typename TValue, typename TConfig, typename TSize, typename TExpand >
    inline typename Size< String<TValue, External<TConfig> > >::Type
    reserve(
	    String<TValue, External<TConfig> > &me,
		TSize new_capacity,
		Tag<TExpand> expand)
	{
//IOREV
		reserve(me.pager, enclosingBlocks(new_capacity, (unsigned)me.PAGESIZE), expand);
		return capacity(me);
	}
//____________________________________________________________________________

    template < typename TValue, typename TConfig, typename TSpec >
    inline typename Iterator<String<TValue, External<TConfig> >, Tag<TSpec> const>::Type
    begin(String<TValue, External<TConfig> > &me, Tag<TSpec> const) 
	{
//IOREV
		typedef String<TValue, External<TConfig> > TString;
		return typename Iterator<TString, Tag<TSpec> const>::Type (&me, 0);
    }

    template < typename TValue, typename TConfig, typename TSpec >
    inline typename Iterator<String<TValue, External<TConfig> > const, Tag<TSpec> const>::Type
    begin(String<TValue, External<TConfig> > const &me, Tag<TSpec> const) 
	{
//IOREV
		typedef String<TValue, External<TConfig> > TString;
		return typename Iterator<TString const, Tag<TSpec> const>::Type (const_cast<TString*>(&me), 0);
    }

    template < typename TValue, typename TConfig, typename TSpec >
    inline typename Iterator<String<TValue, External<TConfig> >, Tag<TSpec> const>::Type
    end(String<TValue, External<TConfig> > &me, Tag<TSpec> const) 
	{
//IOREV
		typedef String<TValue, External<TConfig> > TString;
		return typename Iterator<TString, Tag<TSpec> const>::Type (&me, length(me));
    }

    template < typename TValue, typename TConfig, typename TSpec >
    inline typename Iterator<String<TValue, External<TConfig> > const, Tag<TSpec> const>::Type
    end(String<TValue, External<TConfig> > const &me, Tag<TSpec> const) 
	{
//IOREV
		typedef String<TValue, External<TConfig> > TString;
		return typename Iterator<TString const, Tag<TSpec> const>::Type (const_cast<TString*>(&me), length(me));
    }
//____________________________________________________________________________

    template < typename TValue, typename TConfig, typename TPos >
    inline typename Reference<String<TValue, External<TConfig> > >::Type 
    value(String<TValue, External<TConfig> > &me, TPos pos)
    {
//IOREV
	    return me[pos];
    }

    template < typename TValue, typename TConfig, typename TPos >
    inline typename Reference<String<TValue, External<TConfig> > const>::Type 
    value(String<TValue, External<TConfig> > const &me, TPos pos)
    {
//IOREV
	    return me[pos];
    }
//____________________________________________________________________________

	template < typename TValue, typename TConfig, typename TExpand >
	inline void
	appendValue(String<TValue, External<TConfig> > &me, 
				TValue const &Val_,
				Tag<TExpand> expand)
	{
//IOREV
		resize(me, me.data_size + 1, expand);
		back(me) = Val_;
	}

//____________________________________________________________________________
// stack interface

    template < typename TValue, typename TConfig >
    inline void
    push(String<TValue, External<TConfig> > &me, TValue const &Val_)
    {
//IOREV
		appendValue(me, Val_);
    }

    template < typename TValue, typename TConfig >
    inline void
    push_back(String<TValue, External<TConfig> > &me, TValue const &Val_)
    {
//IOREV _nodoc_
		appendValue(me, Val_);
    }

    template < typename TValue, typename TConfig >
	inline void pop_back(String<TValue, External<TConfig> > &me)
    {
//IOREV _nodoc_
		resize(me, me.data_size - 1);
	}

    template < typename TValue, typename TConfig >
	inline TValue &
	front(String<TValue, External<TConfig> > &me)
    {
//IOREV
		return me[0];
	}

    template < typename TValue, typename TConfig >
	inline TValue const &
	front(String<TValue, External<TConfig> > const &me)
    {
//IOREV
		return me[0];
	}

    template < typename TValue, typename TConfig >
	inline TValue &
	back(String<TValue, External<TConfig> > &me)
    {
//IOREV
		return me[me.data_size - 1];
	}

    template < typename TValue, typename TConfig >
	inline TValue const &
	back(String<TValue, External<TConfig> > const &me)
    {
//IOREV
		return me[me.data_size - 1];
	}
//____________________________________________________________________________

	template < typename TValue, typename TConfig, typename TSource, typename TExpand >
	inline void
	append(String<TValue, External<TConfig> > &target, 
				TSource const &source,
				Tag<TExpand> expand)
	{
//IOREV doc says, resize() my invalidate iterators, therefore it_target might be bogus
		typedef String<TValue, External<TConfig> >					TTarget;
        typedef typename Iterator<TSource const, Standard>::Type	ISource;
        typedef typename Iterator<TTarget, Standard>::Type			ITarget;

        ITarget it_target       = end(target, Standard());
		
		resize(target, length(target) + length(source), expand);
		
        ISource it_source       = begin(source, Standard());
        ISource it_source_end   = end(source, Standard());
		
		for (; it_source != it_source_end; ++it_source, ++it_target)
			*it_target = *it_source;
	}

	template < typename TValue, typename TConfig, typename TSourceValue, typename TExpand >
	inline void
	append(String<TValue, External<TConfig> > &target, 
				TSourceValue * source,
				Tag<TExpand> expand)
	{
//IOREV doc says, resize() my invalidate iterators, therefore it_target might be bogus
		typedef String<TValue, External<TConfig> >					TTarget;
        typedef typename Iterator<TSourceValue *, Standard>::Type	ISource;
        typedef typename Iterator<TTarget, Standard>::Type			ITarget;

        ITarget it_target       = end(target, Standard());
		
		resize(target, length(target) + length(source), expand);
		
        ISource it_source       = begin(source, Standard());
        ISource it_source_end   = end(source, Standard());
		
		for (; it_source != it_source_end; ++it_source, ++it_target)
			*it_target = *it_source;
	}


/*
    template < typename TSpec >
	std::ostream& operator<<(std::ostream &out, String<char, TSpec > &p) {

        typename Iterator< String<char, TSpec > >::Type _cur = begin(p), _end = end(p);
        while (_cur != _end) {
		    out << *_cur;
            ++_cur;
        }
		return out;
	}

    template < typename TValue, typename TSpec >
	std::ostream& operator<<(std::ostream &out, String<TValue, TSpec > &p) {

        typename Iterator< String<TValue, TSpec > >::Type _cur = begin(p), _end = end(p);
        while (_cur != _end) {
		    out << *_cur << " ";
            ++_cur;
        }
		return out;
	}
*/
//____________________________________________________________________________
// sequence -> external string

    template < typename TValue,
               typename TConfig,
               typename TSource,
			   typename TExpand >
    inline void assign(
		String<TValue, External<TConfig> > &target, 
		TSource const &source, 
		Tag<TExpand>) 
	{
//IOREV
		typedef String<TValue, External<TConfig> >					TTarget;
        typedef typename Iterator<TSource const, Standard>::Type	ISource;
        typedef typename Iterator<TTarget, Standard>::Type			ITarget;

        resize(target, length(source));

        ISource it_source       = begin(source, Standard());
        ISource it_source_end   = end(source, Standard());
        ITarget it_target       = begin(target, Standard());
		
		for (; it_source != it_source_end; ++it_source, ++it_target)
			*it_target = *it_source;
    }

    template < typename TValue,
               typename TConfig,
               typename TSourceValue,
			   typename TExpand >
    inline void assign(
		String<TValue, External<TConfig> > &target, 
		TSourceValue * source,
		Tag<TExpand>) 
	{
//IOREV
		typedef String<TValue, External<TConfig> >					TTarget;
        typedef typename Iterator<TSourceValue *, Standard>::Type	ISource;
        typedef typename Iterator<TTarget, Standard>::Type			ITarget;

        resize(target, length(source));

        ISource it_source       = begin(source, Standard());
        ISource it_source_end   = end(source, Standard());
        ITarget it_target       = begin(target, Standard());
		
		for (; it_source != it_source_end; ++it_source, ++it_target)
			*it_target = *it_source;
    }

//____________________________________________________________________________

	template < typename TValue, typename TConfig >
    inline void const * 
    getObjectId(String<TValue, External<TConfig> > const &me)
    {
//IOREV
        return &(*begin(me.pager));
    }

}

#endif
