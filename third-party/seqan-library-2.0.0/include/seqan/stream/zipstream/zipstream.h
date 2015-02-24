// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/io/zipstream.hpp
/// @brief  Altered zipstream library header
/// @author Jonathan de Halleux (dehalleux@pelikhan.com)
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author David Kim (dekim@u.washington.edu)
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// zipstream Library License:
// --------------------------
//
// The zlib/libpng License Copyright (c) 2003 Jonathan de Halleux.
//
// This software is provided 'as-is', without any express or implied warranty. In no event will the authors be held liable for any damages arising from the use of this software.
//
// Permission is granted to anyone to use this software for any purpose, including commercial applications, and to alter it and redistribute it freely, subject to the following restrictions:
//
// 1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
//
// 2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
//
// 3. This notice may not be removed or altered from any source distribution
//
// Author: Jonathan de Halleux, dehalleux@pelikhan.com, 2003


#ifndef INCLUDED_utility_io_zipstream_HPP
#define INCLUDED_utility_io_zipstream_HPP


// Zlib headers
#include <zlib.h>
#include "zutil.h"



// C++ headers
#include <algorithm>
#include <iostream>
#include <vector>


namespace zlib_stream {


/// @brief Default gzip buffer size, change this to suite your needs
const size_t default_buffer_size = 921600; // Was 102400; Was 4096;



/// Compression strategy, see zlib doc.
enum EStrategy
{
    StrategyFiltered = 1,
    StrategyHuffmanOnly = 2,
    DefaultStrategy = 0
};


/// @brief A stream decorator that takes raw input and zips it to a ostream.
/// @note  The class wraps up the inflate method of the zlib library 1.1.4 http://www.gzip.org/zlib/
template<
    typename Elem,
    typename Tr = std::char_traits< Elem >,
    typename ElemA = std::allocator< Elem >,
    typename ByteT = unsigned char,
    typename ByteAT = std::allocator< ByteT >
>
class basic_zip_streambuf :
    public std::basic_streambuf< Elem, Tr >
{

public:

    typedef  std::basic_streambuf< Elem, Tr >  basic_streambuf_type;
    typedef  std::basic_ostream< Elem, Tr > &  ostream_reference;
    typedef  Elem  char_type;
    typedef  ElemA  char_allocator_type;
    typedef  ByteT  byte_type;
    typedef  ByteAT  byte_allocator_type;
    typedef  byte_type *  byte_buffer_type;
    typedef  std::vector< byte_type, byte_allocator_type >  byte_vector_type;
    typedef  std::vector< char_type, char_allocator_type >  char_vector_type;
    typedef  Tr  traits_type;
    typedef  typename Tr::int_type  int_type;

    using basic_streambuf_type::epptr;
    using basic_streambuf_type::pbase;
    using basic_streambuf_type::pptr;

    /// @brief Construct a zip stream
    /// @note  More info on the following parameters can be found in the zlib documentation
    basic_zip_streambuf(
        ostream_reference ostream_,
        size_t level_,
        EStrategy strategy_,
        size_t window_size_,
        size_t memory_level_,
        size_t buffer_size_
    );

    ~basic_zip_streambuf();

    int sync();
    int_type overflow( int_type c );

    /// @brief flushes the zip buffer and output buffer
    std::streamsize flush();

    /// @brief flushes the zip buffer and output buffer and finalize the zip stream
    /// @details This method should be called at the end of the compression.
    std::streamsize flush_finalize();

    /// @brief resets the zip stream and zeros the crc
    /// @details This method should be called after flush_finalize()
    /// @deatils to allow future writes
    void reset_state();

    /// @brief returns a reference to the output stream
    ostream_reference get_ostream() const { return m_ostream; }

    /// @brief returns the latest zlib error status
    int get_zerr() const { return m_err; }

    /// @brief returns the crc of the input data compressed so far
    uLong get_crc() const { return m_crc; }

    /// @brief returns the size (bytes) of the input data compressed so far
    uLong get_in_size() const { return m_zip_stream.total_in; }

    /// @brief returns the size (bytes) of the compressed data so far
    uLong get_out_size() const { return m_zip_stream.total_out; }

private:

    bool zip_to_stream( char_type*, std::streamsize );
    size_t fill_input_buffer();

    /// @brief flush the zip buffer using a particular mode and flush output buffer
    std::streamsize flush( int flush_mode );

    ostream_reference m_ostream;
    z_stream m_zip_stream;
    int m_err;
    byte_vector_type m_output_buffer;
    char_vector_type m_buffer;
    uLong m_crc;

}; // basic_zip_streambuf


/// @brief A stream decorator that takes compressed input and unzips it to a istream.
/// @note  The class wraps up the deflate method of the zlib library 1.1.4 http://www.gzip.org/zlib/
template<
    typename Elem,
    typename Tr = std::char_traits< Elem >,
    typename ElemA = std::allocator< Elem >,
    typename ByteT = unsigned char,
    typename ByteAT = std::allocator< ByteT >
>
class basic_unzip_streambuf :
    public std::basic_streambuf< Elem, Tr >
{

public:

    typedef  std::basic_streambuf< Elem, Tr >  basic_streambuf_type;
    typedef  std::basic_istream< Elem, Tr > &  istream_reference;
    typedef  Elem  char_type;
    typedef  ElemA  char_allocator_type;
    typedef  ByteT  byte_type;
    typedef  ByteAT  byte_allocator_type;
    typedef  byte_type *  byte_buffer_type;
    typedef  std::vector< byte_type, byte_allocator_type >  byte_vector_type;
    typedef  std::vector< char_type, char_allocator_type >  char_vector_type;
    typedef  typename Tr::int_type  int_type;

    using basic_streambuf_type::eback;
    using basic_streambuf_type::egptr;
    using basic_streambuf_type::gptr;

    /// @brief Construct an unzip stream
    /// @note  More info on the following parameters can be found in the zlib documentation
    basic_unzip_streambuf(
        istream_reference istream_,
        size_t window_size_,
        size_t read_buffer_size_,
        size_t input_buffer_size_
    );

    ~basic_unzip_streambuf();

    int_type underflow();

    /// @brief returns the compressed input istream
    istream_reference get_istream() { return m_istream; }

    /// @brief returns the zlib stream structure
    z_stream & get_zip_stream() { return m_zip_stream; }

    /// @brief returns the latest zlib error state
    int get_zerr() const { return m_err; }

    /// @brief returns the crc of the uncompressed data so far
    uLong get_crc() const { return m_crc; }

    /// @brief returns the number of uncompressed bytes
    uLong get_out_size() const { return m_zip_stream.total_out; }

    /// @brief returns the number of read compressed bytes
    uLong get_in_size() const { return m_zip_stream.total_in; }

private:

    void put_back_from_zip_stream();

    std::streamsize unzip_from_stream( char_type*, std::streamsize );

    size_t fill_input_buffer();

    istream_reference m_istream;
    z_stream m_zip_stream;
    int m_err;
    byte_vector_type m_input_buffer;
    char_vector_type m_buffer;
    uLong m_crc;

}; // basic_unzip_streambuf


/// @brief Base class for zip ostreams
/// @note  Contains a basic_zip_streambuf
template<
    typename Elem,
    typename Tr = std::char_traits< Elem >,
    typename ElemA = std::allocator< Elem >,
    typename ByteT = unsigned char,
    typename ByteAT = std::allocator< ByteT >
>
class basic_zip_ostreambase :
    virtual public std::basic_ios< Elem, Tr >
{

public:

    typedef  std::basic_ostream<Elem, Tr> &  ostream_reference;
    typedef  basic_zip_streambuf<
        Elem,
        Tr,
        ElemA,
        ByteT,
        ByteAT
    >  zip_streambuf_type;

    /// @brief Construct a zip stream
    /// @note  More info on the following parameters can be found in the zlib documentation.
    basic_zip_ostreambase(
        ostream_reference ostream_,
        size_t level_,
        EStrategy strategy_,
        size_t window_size_,
        size_t memory_level_,
        size_t buffer_size_
    ) :
        m_buf( ostream_, level_, strategy_, window_size_, memory_level_, buffer_size_ )
    {
        this->init( &m_buf );
    }

    /// @brief returns the underlying zip ostream object
    zip_streambuf_type * rdbuf() { return &m_buf; }

    /// @brief returns the zlib error state
    int get_zerr() const { return m_buf.get_err(); }

    /// @brief returns the uncompressed data crc
    uLong get_crc() const { return m_buf.get_crc(); }

    /// @brief returns the compressed data size
    uLong get_out_size() const { return m_buf.get_out_size(); }

    /// @brief returns the uncompressed data size
    uLong get_in_size() const { return m_buf.get_in_size(); }

private:

    zip_streambuf_type m_buf;

}; // basic_zip_ostreambase


/// @brief Base class for unzip istreams
/// @note  Contains a basic_unzip_streambuf
template<
    typename Elem,
    typename Tr = std::char_traits< Elem >,
    typename ElemA = std::allocator< Elem >,
    typename ByteT = unsigned char,
    typename ByteAT = std::allocator< ByteT >
>
class basic_zip_istreambase :
    virtual public std::basic_ios< Elem, Tr >
{

public:

    typedef  std::basic_istream< Elem, Tr > &  istream_reference;
    typedef  basic_unzip_streambuf<
        Elem,
        Tr,
        ElemA,
        ByteT,
        ByteAT
    >  unzip_streambuf_type;

    basic_zip_istreambase(
        istream_reference ostream_,
        size_t window_size_,
        size_t read_buffer_size_,
        size_t input_buffer_size_
    ) :
        m_buf( ostream_, window_size_, read_buffer_size_, input_buffer_size_ )
    {
        this->init( &m_buf );
    }

    /// @brief returns the underlying unzip istream object
    unzip_streambuf_type * rdbuf() { return &m_buf; }

    /// @brief returns the zlib error state
    int get_zerr() const { return m_buf.get_zerr(); }

    /// @brief returns the uncompressed data crc
    uLong get_crc() const { return m_buf.get_crc(); }

    /// @brief returns the uncompressed data size
    uLong get_out_size() const { return m_buf.get_out_size(); }

    /// @brief returns the compressed data size
    uLong get_in_size() const { return m_buf.get_in_size(); }

private:

    unzip_streambuf_type m_buf;

}; // basic_zip_istreambase


/// @brief A zipper ostream
///
/// @remarks
///
/// This class is a ostream decorator that behaves 'almost' like any other ostream.
///
/// At construction, it takes any ostream that shall be used to output of the compressed data.
///
/// When finished, you need to call the special method zflush or call the destructor
/// to flush all the intermidiate streams.
///
/// Example:
/// \code
/// // creating the target zip string, could be a fstream
/// ostringstream ostringstream_;
/// // creating the zip layer
/// zip_ostream zipper(ostringstream_);
///
///
/// // writing data
/// zipper<<f<<" "<<d<<" "<<ui<<" "<<ul<<" "<<us<<" "<<c<<" "<<dum;
/// // zip ostream needs special flushing...
/// zipper.zflush();
/// \endcode
template<
    typename Elem,
    typename Tr = std::char_traits< Elem >,
    typename ElemA = std::allocator< Elem >,
    typename ByteT = unsigned char,
    typename ByteAT = std::allocator< ByteT >
>
class basic_zip_ostream :
    public basic_zip_ostreambase< Elem, Tr, ElemA, ByteT, ByteAT >,
    public std::basic_ostream< Elem, Tr >
{

public:

    typedef  basic_zip_ostreambase< Elem, Tr, ElemA, ByteT, ByteAT >  zip_ostreambase_type;
    typedef  std::basic_ostream< Elem, Tr >  ostream_type;
    typedef  std::basic_ostream< Elem, Tr > &  ostream_reference;
    typedef  Elem  char_type;

    using ostream_type::flush;
    using zip_ostreambase_type::rdbuf;

    /// @brief Constructs a zipper ostream decorator
    ///
    /// @param ostream_ ostream where the compressed output is written
    /// @param is_gzip_ true if gzip header and footer have to be added
    /// @param level_ level of compression 0, bad and fast, 9, good and slower,
    /// @param strategy_ compression strategy
    /// @param window_size_ see zlib doc
    /// @param memory_level_ see zlib doc
    /// @param buffer_size_ the buffer size used to zip data
    ///
    /// @note  When is_gzip_ is true, a gzip header and footer is automatically added
    basic_zip_ostream(
        ostream_reference ostream_,
        // int open_mode = std::ios::out,
        bool is_gzip_ = true,
        size_t level_ = Z_DEFAULT_COMPRESSION,
        EStrategy strategy_ = DefaultStrategy,
        size_t window_size_ = 15,
        size_t memory_level_ = 8,
        size_t buffer_size_ = default_buffer_size
    ) :
        zip_ostreambase_type(
            ostream_,
            level_,
            strategy_,
            window_size_,
            memory_level_,
            buffer_size_
        ),
        ostream_type( rdbuf() ),
        m_is_gzip( is_gzip_ ),
        m_zip_stream_finalized( false )
    {
        if ( m_is_gzip ) add_header();
    }

    ~basic_zip_ostream()
    {
        // adding a footer is not necessary here, as it will be
        // taken care of during the last zflush_finalize()
        // called by the higher level close() routines
        zflush_finalize();
    }

    /// @brief returns true if it is a gzip
    bool is_gzip() const { return m_is_gzip; }

    /// @brief flush inner buffer and zipper buffer
    basic_zip_ostream< Elem, Tr > &
    zflush()
    {
        flush(); rdbuf()->flush(); return *this;
    }

    /// @brief flush inner and zipper buffers and finalize zip stream
    basic_zip_ostream< Elem, Tr > &
    zflush_finalize()
    {

        flush(); rdbuf()->flush_finalize();

        if ( m_is_gzip && ( rdbuf()->get_zerr() == Z_STREAM_END ) && ( !m_zip_stream_finalized ) ) {
            add_footer(); // writes crc trailer to end the current zip stream
            flush();
            m_zip_stream_finalized = true;
        }

        return *this;
    }

    /// @brief stream output
    /// @details if zip stream has been finalized, will reset
    /// @details the stream and add header if necessary
    template< typename T >
    inline
    basic_zip_ostream &
    operator <<( T const & t )
    {
        reset_zip_stream();
        static_cast< std::ostream & >( *this ) << t;
        return *this;
    }

    /// @brief write char
    /// @details if zip stream has been finalized, will reset
    /// @details the stream and add header if necessary
    inline
    basic_zip_ostream &
    put( char const c )
    {
        reset_zip_stream();
        static_cast< std::ostream & >( *this ).put( c );
        return *this;
    }

    /// @brief write a string
    /// @details if zip stream has been finalized, will reset
    /// @details the stream and add header if necessary
    inline
    basic_zip_ostream &
    write( char const * str, std::streamsize const count )
    {
        reset_zip_stream();
        static_cast< std::ostream & >( *this ).write( str, count );
        return *this;
    }

private:

    /// @brief if end of stream, reset the zip stream and add header
    inline
    bool
    reset_zip_stream()
    {
        if ( rdbuf()->get_zerr() == Z_STREAM_END ) {
            rdbuf()->reset_state();
            add_header();
            m_zip_stream_finalized = false;
            return true;
        }
        return false;
    }

    static void put_long_as_uint32( ostream_reference out_, unsigned long x_ );

    void add_header();
    void add_footer();

    bool m_is_gzip;

    /// @brief tracks to see if zip stream was finalized
    /// @details set to true during zflush_finalize()
    /// @details set to false during reset_state()
    bool m_zip_stream_finalized;

#ifdef _WIN32
private:
    void _Add_vtordisp1() { } // Required to avoid VC++ warning C4250
    void _Add_vtordisp2() { } // Required to avoid VC++ warning C4250
#endif
}; // basic_zip_ostream


/// @brief A zipper istream
///
/// @remarks
///
/// This class is a istream decorator that behaves 'almost' like any other ostream.
///
/// At construction, it takes any istream that shall be used to input of the compressed data.
///
/// Simlpe example:
/// \code
/// // create a stream on zip string
/// istringstream istringstream_( ostringstream_.str());
/// // create unzipper istream
/// zip_istream unzipper( istringstream_);
///
/// // read and unzip
/// unzipper>>f_r>>d_r>>ui_r>>ul_r>>us_r>>c_r>>dum_r;
/// \endcode
template<
    typename Elem,
    typename Tr = std::char_traits< Elem >,
    typename ElemA = std::allocator< Elem >,
    typename ByteT = unsigned char,
    typename ByteAT = std::allocator< ByteT >
>
class basic_zip_istream :
    public basic_zip_istreambase< Elem, Tr, ElemA, ByteT, ByteAT >,
    public std::basic_istream< Elem, Tr >
{

public:

    typedef  basic_zip_istreambase< Elem, Tr, ElemA, ByteT, ByteAT >  zip_istreambase_type;
    typedef  std::basic_istream< Elem, Tr >  istream_type;
    typedef  std::basic_istream< Elem, Tr > &  istream_reference;
    typedef  Elem  char_type;
    typedef  unsigned char  byte_type;

    using zip_istreambase_type::get_crc;
    using zip_istreambase_type::get_out_size;
    using zip_istreambase_type::rdbuf;

    /// @brief Construct a unzipper stream
    ///
    /// @param istream_ input buffer
    /// @param window_size_
    /// @param read_buffer_size_
    /// @param input_buffer_size_
    basic_zip_istream(
        istream_reference istream_,
        size_t window_size_ = 15,
        size_t read_buffer_size_ = default_buffer_size,
        size_t input_buffer_size_ = default_buffer_size
    ) :
        zip_istreambase_type( istream_, window_size_, read_buffer_size_, input_buffer_size_ ),
        istream_type( rdbuf() ),
        m_is_gzip( false ),
        m_gzip_crc( 0 ),
        m_gzip_data_size( 0 )
    {
        if ( rdbuf()->get_zerr() == Z_OK ) check_header();
    }

    /// @brief returns true if it is a gzip file
    bool is_gzip() const { return m_is_gzip; }

    /// @brief reads the gzip footer
    void read_footer();

    /// @brief return crc check result
    ///
    /// @note  When you have finished reading the compressed data,
    ///        call read_footer to read the uncompressed data crc.
    /// @note  This method compares it to the crc of the uncompressed data.
    ///
    /// @return true if crc check is succesful
    bool check_crc() const { return get_crc() == m_gzip_crc; }

    /// @brief return data size check
    bool check_data_size() const { return get_out_size() == m_gzip_data_size; }

    /// @brief return the crc value in the file
    uLong get_gzip_crc() const { return m_gzip_crc; }

    /// @brief return the data size in the file
    uLong get_gzip_data_size() const { return m_gzip_data_size; }

protected:

    static void read_uint32( istream_reference in_, unsigned long & x_ );

    int check_header();

    bool m_is_gzip;
    uLong m_gzip_crc;
    uLong m_gzip_data_size;

#ifdef _WIN32
private:
    void _Add_vtordisp1() { } // Required to avoid VC++ warning C4250
    void _Add_vtordisp2() { } // Required to avoid VC++ warning C4250
#endif
}; // basic_zip_istream


// Types
typedef  basic_zip_ostream< char >  zip_ostream;
typedef  basic_zip_ostream< wchar_t >  zip_wostream;
typedef  basic_zip_istream< char >  zip_istream;
typedef  basic_zip_istream< wchar_t >  zip_wistream;


} // namespace zlib_stream


// Implementation
#include "zipstream_impl.h"


#endif // INCLUDED_utility_io_zipstream_HPP
