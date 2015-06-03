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
// Memory map a whole file and use it as a string.
// ==========================================================================
// NOTE(esiragusa): this file is deprecated.

#ifndef SEQAN_HEADER_STRING_MMAP_H
#define SEQAN_HEADER_STRING_MMAP_H


/* IOREV
 * _tested_
 * _windows_
 *
 *
 * tested in library/demos/howto/efficiently_import_sequences.cpp and stellar
 *
 * relation to file_format_mmap.h unclear
 *
 * relation to string_external unclear, what benifit does string_mmap provide?
 *
 */


//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{

/*!
 * @class MMapConfig
 * @headerfile <seqan/file.h>
 * @brief Configuration for @link MMapString MMapStrings @endlink.
 *
 * @signature template <[typename TFile, [typename TSize]]>
 *            class MMapConfig;
 *
 * @tparam TFile The file type to use for the memory mapped string.  Defaults to <tt>File&lt;&gt;</tt>
 * @tparam TSize The size type to use.  Defaults to <tt>size_t</tt>.
 *
 *
 * @typedef MMapConfig::TFile;
 * @signature typedef (..) MMapConfig::TFile;
 * @brief The type to use for files.
 *
 * @typedef MMapConfig::TSize;
 * @signature typedef (..) MMapConfig::TSize;
 * @brief The type to use for sizes.
 */

    template < typename TFile_ = File<>,                // default file type
               typename TSize_ = size_t >                // size type
    struct MMapConfig {
//IOREV _nodoc_ doc says using MMap<MMapConfig> is correct, whats this for?
        typedef TFile_ TFile;
        typedef TSize_ TSize;
    };

    template < typename TConfig = MMapConfig<> >
    struct MMap {};
//IOREV


    //////////////////////////////////////////////////////////////////////////////
    // Memory Mapped String
    //////////////////////////////////////////////////////////////////////////////

/*!
 * @class MMapString MMap String
 * @headerfile <seqan/file.h>
 * @extends String
 * @brief String that is stored in external memory using direct memory mapping.
 *
 * @signature template <typename TValue[, typename TConfig]>
 *            class String<TValue, MMap<TConfig> >;
 *
 * @tparam TValue  The value type to use for the items/characters.
 * @tparam TConfig The configuration to use for the underlying file.  Default: @link MMapConfig @endlink.
 *
 * The MMap String enables to access sequences larger than the available physical memory (RAM) by using external memory
 * (e.g. Hard disk, Network storage, ...) mapped into memory.  The size of the string is limited by external memory and
 * the logical address space (4GB on a 32bit OS).  See the String constructor for more details.
 *
 * This string also supports fast appending and removing of values at the end.
 */

/*!
 * @fn MMapString::String
 * @brief Constructor.
 *
 * @signature String::String();
 * @signature String::String(file);
 * @signature String::String(fileName[, openMode]);
 *
 * @param[in,out] file     The @link File @endlink to use for reading and writing.  You must ensture that
 *                         <tt>file</tt> is open as the string will not call <tt>open</tt> and <tt>close</tt>
 *                         on the file.
 * @param[in]     fileName The path to open.  Type: <tt>char const *</tt>
 * @param[in]     openMode The open mode.
 */

    template < typename TValue,
               typename TConfig >
    class String<TValue, MMap<TConfig> >
    {
//IOREV
    public:

        typedef typename TConfig::TFile        TFile;
        typedef typename TConfig::TSize        TSize;

        TValue                  *data_begin;
        TValue                  *data_end;

        FileMapping<>           mapping;
        FileMappingAdvise       advise;

        explicit
        String(TSize size = 0):
            data_begin(NULL),
            data_end(NULL),
            advise(MAP_NORMAL)
        {
            resize(*this, size);
        }

        // associate a file
        explicit
        String(TFile &_file):
            data_begin(NULL),
            data_end(NULL),
            advise(MAP_NORMAL)
        {
            open(*this, _file);
        }

        // associate a file given filename [and open mode]
        explicit
        String(const char *fileName, int openMode = DefaultOpenMode<TFile>::VALUE):
            data_begin(NULL),
            data_end(NULL),
            advise(MAP_NORMAL)
        {
            open(*this, fileName, openMode);
        }

        // copy the contents from another string
        String(String const & source):
            data_begin(NULL),
            data_end(NULL),
            advise(MAP_NORMAL)
        {
            assign(*this, source);
        }
        template <typename TSource>
        String(TSource const & source):
            data_begin(NULL),
            data_end(NULL),
            advise(MAP_NORMAL)
        {
            assign(*this, source);
        }

        template <typename TSource>
        String & operator =(TSource const & source)
        {
            assign(*this, source);
            return *this;
        }
        String & operator =(String const & source)
        {
            assign(*this, source);
            return *this;
        }

        ~String()
        {
            close(*this);
        }

//____________________________________________________________________________

        template <typename TPos>
        inline typename Reference<String>::Type
        operator [] (TPos pos)
        {
            return value(*this, pos);
        }

        template <typename TPos>
        inline typename Reference<String const>::Type
        operator [] (TPos pos) const
        {
            return value(*this, pos);
        }

//____________________________________________________________________________

        inline operator bool()
        {
            return mapping;
        }

//____________________________________________________________________________

};

    template < typename TValue, typename TConfig >
    inline typename Iterator<String<TValue, MMap<TConfig> > , Standard>::Type
    begin(String<TValue, MMap<TConfig> > & me,
        Standard)
    {
//IOREV
SEQAN_CHECKPOINT
        return me.data_begin;
    }
    template < typename TValue, typename TConfig >
    inline typename Iterator<String<TValue, MMap<TConfig> >  const, Standard>::Type
    begin(String<TValue, MMap<TConfig> > const & me,
        Standard)
    {
//IOREV
SEQAN_CHECKPOINT
        return me.data_begin;
    }

//____________________________________________________________________________

    template < typename TValue, typename TConfig >
    inline typename Iterator<String<TValue, MMap<TConfig> > , Standard>::Type
    end(String<TValue, MMap<TConfig> > & me,
        Standard)
    {
//IOREV
SEQAN_CHECKPOINT
        return me.data_end;
    }
    template < typename TValue, typename TConfig >
    inline typename Iterator<String<TValue, MMap<TConfig> >  const, Standard>::Type
    end(String<TValue, MMap<TConfig> > const & me,
        Standard)
    {
//IOREV
SEQAN_CHECKPOINT
        return me.data_end;
    }

//____________________________________________________________________________

    template < typename TValue, typename TConfig >
    inline typename Size<String<TValue, MMap<TConfig> > >::Type
    capacity(String<TValue, MMap<TConfig> > const & me)
    {
//IOREV
SEQAN_CHECKPOINT
        return length(me.mapping) / sizeof(TValue);
    }

//____________________________________________________________________________

    template < typename TValue, typename TConfig >
    inline void
    _setLength(
        String<TValue, MMap<TConfig> > & me,
        size_t new_length)
    {
//IOREV
SEQAN_CHECKPOINT
        me.data_end = me.data_begin + new_length;
    }

    //////////////////////////////////////////////////////////////////////////////
    // meta-function interface

    template < typename TValue, typename TConfig >
    struct Size< String<TValue, MMap<TConfig> > >
    {
//IOREV
        typedef typename TConfig::TSize Type;
    };

    template < typename TValue, typename TConfig >
    struct Difference< String<TValue, MMap<TConfig> > >
    {
//IOREV
        typedef typename MakeSigned_<typename TConfig::TSize>::Type Type;
    };
//____________________________________________________________________________

    template < typename TValue, typename TConfig >
    struct DefaultOverflowExplicit<String<TValue, MMap<TConfig> > >
    {
//IOREV
        typedef Generous Type;
    };

    template < typename TValue, typename TConfig >
    struct DefaultOverflowImplicit<String<TValue, MMap<TConfig> > >
    {
//IOREV
        typedef Generous Type;
    };
//____________________________________________________________________________

    template < typename TValue, typename TConfig >
    struct IsContiguous< String<TValue, MMap<TConfig> > >
    {
//IOREV
        typedef True Type;
        enum { VALUE = true };
    };

    template < typename TValue, typename TConfig >
    struct AllowsFastRandomAccess< String<TValue, MMap<TConfig> > >
    {
//IOREV
        typedef False Type;
        enum { VALUE = false };
    };


    //////////////////////////////////////////////////////////////////////////////
    // global interface

//____________________________________________________________________________

/*!
 * @fn MMapString#flush
 * @brief Waits for all open read/write requests to complete.
 *
 * @signature void flush(str);
 *
 * @param[in,out] str The MMapString to flush.
 */

    template < typename TValue, typename TConfig >
    inline bool
    flush(String<TValue, MMap<TConfig> > &me)
    {
        typedef typename Size<typename TConfig::TFile>::Type TFileSize;
        return flushFileSegment(
            me.mapping,
            me.data_begin,
            0,
            (TFileSize)capacity(me) * (TFileSize)sizeof(TValue));
    }

/*!
 * @fn MMapString#mmapAdvise
 * @brief Call advise function for memory mapped files.
 *
 * @signature bool mmapAdvise(str, scheme[, beginPos, size]);
 *
 * @param[in,out] str        The MMapString to call advise in.
 * @param[in]     scheme     The memory access scheme to use.  Type: FileMappingAdvise.
 * @param[in]     beginPos   Begin position in the string for the advise call.
 * @param[in]     size       Size of the range used for the advise call.
 *
 * @return bool <tt>true</tt> if the advise was successful, <tt>false</tt> otherwise.
 */

    template <typename TValue, typename TConfig, typename TPos, typename TSize>
    inline bool
    mmapAdvise(String<TValue, MMap<TConfig> > &me, FileMappingAdvise advise, TPos beginPos, TSize size)
    {
        typedef typename Size<typename TConfig::TFile>::Type TFileSize;
        me.advise = advise;
        return adviseFileSegment(
            me.mapping,
            advise,
            me.data_begin,
            (TFileSize)beginPos * (TFileSize)sizeof(TValue),
            (TFileSize)size * (TFileSize)sizeof(TValue));
    }

    template <typename TValue, typename TConfig, typename TPos>
    inline int
    mmapAdvise(String<TValue, MMap<TConfig> > &me, FileMappingAdvise advise, TPos beginPos)
    {
        return mmapAdvise(me, advise, beginPos, capacity(me) - beginPos);
    }

    template <typename TValue, typename TConfig>
    inline int
    mmapAdvise(String<TValue, MMap<TConfig> > &me, FileMappingAdvise advise)
    {
        return mmapAdvise(me, advise, 0, capacity(me));
    }

//____________________________________________________________________________

    // cancel all transactions
    template < typename TValue, typename TConfig >
    inline void
    cancel(String<TValue, MMap<TConfig> > &me)
    {
        typedef typename Size<typename TConfig::TFile>::Type TFileSize;
        cancelFileSegment(
            me.mapping,
            me.data_begin,
            0,
            (TFileSize)capacity(me) * (TFileSize)sizeof(TValue));
    }

//____________________________________________________________________________

    // flush and free all allocated pages
    template < typename TValue, typename TConfig >
    inline bool
    flushAndFree(String<TValue, MMap<TConfig> > &me)
    {
        return flush(me) && mmapAdvise(me, MAP_DONTNEED);
    }

//____________________________________________________________________________

    template < typename TValue, typename TConfig >
    inline bool
    _map(String<TValue, MMap<TConfig> > &me, size_t new_capacity)
    {
        typedef typename Size<typename TConfig::TFile>::Type TFileSize;
        if (new_capacity > 0)
        {
            _ensureFileIsOpen(me);
            if (capacity(me) < new_capacity)
                resize(me.mapping, (TFileSize)new_capacity * (TFileSize)sizeof(TValue));
            me.data_begin = static_cast<TValue*>(mapFileSegment(me.mapping, 0, length(me.mapping)));
            if (me.data_begin == NULL)
            {
                me.data_end = NULL;
                return false;
            }
            adviseFileSegment(me.mapping, me.advise, me.data_begin, 0, length(me.mapping));
        }
        else
            resize(me.mapping, 0);
        _setLength(me, new_capacity);
        return true;
    }

    template < typename TValue, typename TConfig >
    inline bool
    _unmap(String<TValue, MMap<TConfig> > &me)
    {
        bool result = true;
        if (me.data_begin != NULL)
        {
            result = unmapFileSegment(me.mapping, me.data_begin, length(me.mapping));
            me.data_begin = NULL;
        }
        me.data_end = NULL;
        return result;
    }

    template < typename TValue, typename TConfig, typename TCapSize >
    inline bool
    _remap(String<TValue, MMap<TConfig> > &me, TCapSize new_capacity)
    {
        typedef typename Size<String<TValue, MMap<TConfig> > >::Type    TSize;
        typedef typename Size<typename TConfig::TFile>::Type            TFileSize;

        bool result = true;

#ifndef PLATFORM_WINDOWS
        // Windows doesn't allow to resize the file while having a mapped file segment
        // Thus, the following part is only supported on Linux/BSD/Mac OS
        TSize old_capacity = capacity(me);
        if (me.data_begin && new_capacity > 0)
        {
            // if file gets bigger, resize first
            if (old_capacity < new_capacity)
                resize(me.mapping, (TFileSize)new_capacity * (TFileSize)sizeof(TValue));

            me.data_begin = static_cast<TValue*>(remapFileSegment(
                me.mapping,
                me.data_begin,
                0,
                (TFileSize)old_capacity * (TFileSize)sizeof(TValue),
                (TFileSize)new_capacity * (TFileSize)sizeof(TValue)));

            // if file gets smaller, resize at last
            if (old_capacity > new_capacity)
                resize(me.mapping, (TFileSize)new_capacity * (TFileSize)sizeof(TValue));

            if (me.data_begin == NULL)
            {
                me.data_end = NULL;
                return false;
            }
            return true;
        }
#endif
        result &= _unmap(me);
        result &= _map(me, new_capacity);
        return result;
    }

    template < typename TValue, typename TConfig >
    inline void
    clear(String<TValue, MMap<TConfig> > &me)
    {
//IOREV
        cancel(me);
        _unmap(me);
        resize(me.mapping, 0);
    }
//____________________________________________________________________________

    template < typename TValue, typename TConfig, typename TSize >
    inline typename Value<String<TValue, MMap<TConfig> > >::Type *
    _allocateStorage(String<TValue, MMap<TConfig> > &me, TSize new_capacity)
    {
//IOREV
        _map(me, _computeSizeForCapacity(me, new_capacity));
        return NULL;
    }

    template < typename TValue, typename TConfig, typename TSize >
    inline typename Value<String<TValue, MMap<TConfig> > >::Type *
    _reallocateStorage(
        String<TValue, MMap<TConfig> > &me,
        TSize new_capacity)
    {
//IOREV
        TSize size = _computeSizeForCapacity(me, new_capacity);
        _remap(me, size);
        return NULL;
    }

    template < typename TValue, typename TConfig, typename TSize >
    inline void
    _deallocateStorage(String<TValue, MMap<TConfig> > &/*me*/, TValue * /*ptr*/, TSize /*capacity*/)
    {
//IOREV
    }
//____________________________________________________________________________
/*!
 * @fn MMapString#open
 * @brief Open the MMapString's underlying file from a path.
 *
 * @signature bool open(str, fileName[, openMode]);
 *
 * @param[in,out] str      The MMapString to open.
 * @param[in]     fileName Path to the file to open. Type: <tt>char const *</tt>.
 * @param[in]     openMode The open mode. Type: <tt>int</tt>.
 *
 * @return bool <tt>true</tt> if the operation succeeded and <tt>false</tt> otherwise.
 */

    template < typename TValue, typename TConfig >
    inline bool
    open(String<TValue, MMap<TConfig> > &me, const char *fileName, int openMode)
    {
//IOREV
        close(me);
        if (open(me.mapping, fileName, openMode))
            return _map(me, capacity(me));
        return false;
    }

    template < typename TValue, typename TConfig >
    inline bool
    open(String<TValue, MMap<TConfig> > &me, const char *fileName)
    {
//IOREV
        typedef typename String<TValue, MMap<TConfig> >::TFile    TFile;
        return open(me, fileName, DefaultOpenMode<TFile>::VALUE);
    }

    template < typename TValue, typename TConfig >
    inline bool
    open(String<TValue, MMap<TConfig> > &me, typename TConfig::TFile file)
    {
//IOREV
        close(me);
        if (open(me.mapping, file))
            return _map(me, capacity(me));
        return false;
    }

/*!
 * @fn MMapString#openTemp
 * @brief Open an MMapString using an temporary file.
 *
 * @signature bool openTemp(str);
 *
 * @param[in,out] str The MMapString to open using temporary file.
 *
 * @return bool <tt>true</tt> if opening succeeded, <tt>false</tt> otherwise.
 */

    template < typename TValue, typename TConfig >
    inline bool
    openTemp(String<TValue, MMap<TConfig> > &me)
    {
//IOREV
        close(me);
        return openTemp(me.mapping);
    }
//____________________________________________________________________________

    template < typename TValue, typename TConfig >
    inline void
    _ensureFileIsOpen(String<TValue, MMap<TConfig> > &me)
    {
//IOREV
        if (!me)
        {
            if (!openTemp(me.mapping))
                SEQAN_FAIL("Memory Mapped String couldn't open temporary file");
        }
    }
//____________________________________________________________________________

    template < typename TValue, typename TConfig >
    inline bool
    save(String<TValue, MMap<TConfig> > const &/*me*/, const char * /*fileName*/, int /*openMode*/) {
//IOREV _nodoc_ shouldn't we flush here? in case of abnormal termination...
        // Memory Mapped Strings are persistent, thus there is no need to save them
        //MMapStringsDontNeedToBeSaved error;
        return true;
    }

    template < typename TValue, typename TConfig >
    inline bool
    save(String<TValue, MMap<TConfig> > const &/*me*/, const char * /*fileName*/) {
//IOREV _nodoc_ shouldn't we flush here? in case of abnormal termination...
        // Memory Mapped Strings are persistent, thus there is no need to save them
        //MMapStringsDontNeedToBeSaved error;
        return true;
    }

    template < typename TValue, typename TConfig >
    inline bool
    save(String<TValue, MMap<TConfig> > const &/*me*/, typename TConfig::TFile /*file*/) {
//IOREV _nodoc_ shouldn't we flush here? in case of abnormal termination...
        // Memory Mapped Strings are persistent, thus there is no need to save them
        //MMapStringsDontNeedToBeSaved error;
        return true;
    }
//____________________________________________________________________________
/*!
 * @fn MMapString#close
 * @brief Close the MMapString.
 *
 * @signature bool close(str);
 *
 * @param[in] str The MMapString to close the file of.
 *
 * @return bool <tt>true</tt> if the closing succeeded, <tt>false</tt> otherwise.
 */

    template < typename TValue, typename TConfig >
    inline bool
    close(String<TValue, MMap<TConfig> > &me)
    {
        typedef typename Size<typename TConfig::TFile>::Type TFileSize;

        if (me)
        {
            TFileSize finalLen = (TFileSize)length(me) * (TFileSize)sizeof(TValue);

            // close associated file
            if (me.mapping.temporary)
                cancel(me);

            _unmap(me);
            closeAndResize(me.mapping, finalLen);
        }
        return true;
    }


    template < typename TValue, typename TConfig >
    inline void
    shrinkToFit(String<TValue, MMap<TConfig> > &me)
    {
        _remap(me, length(me));
    }

//////////////////////////////////////////////////////////////////////////////



} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
