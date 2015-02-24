// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/io/zipstream.ipp
/// @brief  Altered zipstream library implementation
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


#ifndef INCLUDED_utility_io_zipstream_IPP
#define INCLUDED_utility_io_zipstream_IPP


// Unit headers
#include "zipstream.h"

// C++ headers
#include <sstream>
#include <cstdio>
#include <istream>

namespace zlib_stream {


namespace detail {

    const int gz_magic[2] = { 0x1f, 0x8b }; // gzip magic header

    // gzip flag byte
    const int gz_ascii_flag  = 0x01; // bit 0 set: file probably ascii text
    const int gz_head_crc    = 0x02; // bit 1 set: header CRC present
    const int gz_extra_field = 0x04; // bit 2 set: extra field present
    const int gz_orig_name   = 0x08; // bit 3 set: original file name present
    const int gz_comment     = 0x10; // bit 4 set: file comment present
    const int gz_reserved    = 0xE0; // bits 5..7: reserved

} // namespace detail


    template<
        typename Elem,
        typename Tr,
        typename ElemA,
        typename ByteT,
        typename ByteAT
    >
    basic_zip_streambuf< Elem, Tr, ElemA, ByteT, ByteAT >::basic_zip_streambuf(
        ostream_reference ostream_,
        size_t level_,
        EStrategy strategy_,
        size_t window_size_,
        size_t memory_level_,
        size_t buffer_size_
    ) :
        m_ostream( ostream_ ),
        m_output_buffer( buffer_size_, 0 ),
        m_buffer( buffer_size_, 0 ),
        m_crc( 0 )
    {
        m_zip_stream.zalloc = (alloc_func)0;
        m_zip_stream.zfree = (free_func)0;

        m_zip_stream.next_in = NULL;
        m_zip_stream.avail_in = 0;
        m_zip_stream.avail_out = 0;
        m_zip_stream.next_out = NULL;

        m_err = deflateInit2(
            &m_zip_stream,
            std::min( 9, static_cast< int >( level_ ) ),
            Z_DEFLATED,
            - static_cast< int >(window_size_), // <-- changed
            std::min( 9, static_cast< int >( memory_level_ ) ),
            static_cast< int >( strategy_ )
        );

        this->setp( &(m_buffer[0]), &(m_buffer[m_buffer.size()-1]) );
    }

    template<
        typename Elem,
        typename Tr,
        typename ElemA,
        typename ByteT,
        typename ByteAT
    >
    basic_zip_streambuf< Elem, Tr, ElemA, ByteT, ByteAT >::~basic_zip_streambuf()
    {
        flush_finalize();
        m_ostream.flush();
        m_err=deflateEnd( &m_zip_stream );
    }

    template<
        typename Elem,
        typename Tr,
        typename ElemA,
        typename ByteT,
        typename ByteAT
    >
    int
    basic_zip_streambuf< Elem, Tr, ElemA, ByteT, ByteAT >::sync()
    {
        if ( pptr() && pptr() > pbase() ) {
            if ( traits_type::eq_int_type( overflow( traits_type::eof() ), traits_type::eof() ) ) return -1;
        }

        return 0;
    }

    template<
        typename Elem,
        typename Tr,
        typename ElemA,
        typename ByteT,
        typename ByteAT
    >
    typename basic_zip_streambuf< Elem, Tr, ElemA, ByteT, ByteAT >::int_type
    basic_zip_streambuf< Elem, Tr, ElemA, ByteT, ByteAT >::overflow(
        typename basic_zip_streambuf< Elem, Tr, ElemA, ByteT, ByteAT >::int_type c
    )
    {
        bool const test_eof = traits_type::eq_int_type( c, traits_type::eof() );
        int w = static_cast< int >( pptr() - pbase() );
        if ( !test_eof ) {
            *pptr() = c;
            ++w;
        }
        if ( zip_to_stream( pbase(), w ) ) {
            this->setp( pbase(), epptr() );
            return traits_type::not_eof( c );
        } else {
            return traits_type::eof();
        }
    }

    template<
        typename Elem,
        typename Tr,
        typename ElemA,
        typename ByteT,
        typename ByteAT
    >
    bool
    basic_zip_streambuf< Elem, Tr, ElemA, ByteT, ByteAT >::zip_to_stream(
        typename basic_zip_streambuf< Elem, Tr, ElemA, ByteT, ByteAT >::char_type * buffer_,
        std::streamsize buffer_size_
    )
    {
        std::streamsize written_byte_size = 0, total_written_byte_size = 0;

        m_zip_stream.next_in = (byte_buffer_type)buffer_;
        m_zip_stream.avail_in = static_cast< uInt >( buffer_size_ * sizeof(char_type) );
        m_zip_stream.avail_out = static_cast< uInt >( m_output_buffer.size() );
        m_zip_stream.next_out = &(m_output_buffer[0]);
        size_t remainder = 0;

        // updating crc
        m_crc = crc32(
            m_crc,
            m_zip_stream.next_in,
            m_zip_stream.avail_in
        );

        do {
            m_err = deflate( &m_zip_stream, 0 );

            if ( m_err == Z_OK || m_err == Z_STREAM_END ) {
                written_byte_size =
                    static_cast< std::streamsize >( m_output_buffer.size() )
                    - m_zip_stream.avail_out;
                total_written_byte_size += written_byte_size;
                // ouput buffer is full, dumping to ostream
                //                std::cerr << "dump to m_ostream " << std::endl;
                m_ostream.write(
                    (const char_type*) &(m_output_buffer[0]),
                    static_cast< std::streamsize >( written_byte_size / sizeof( char_type ) )
                );
                //                std::cerr << "dump to m_ostream -- done " << std::endl;
                // checking if some bytes were not written.
                if ( ( remainder = written_byte_size%sizeof( char_type ) ) != 0 ) {
                    // copy to the beginning of the stream
                    //                    std::cerr << "before memcopy in zipstream.ipp: 206" << std::endl;
                    std::memmove(
                        &(m_output_buffer[0]),
                        &(m_output_buffer[written_byte_size-remainder]),
                        remainder
                    );
                    //                    std::cerr << "after memcopy in zipstream.ipp: 206" << std::endl;
                }

                m_zip_stream.avail_out =
                    static_cast< uInt >( m_output_buffer.size() - remainder );
                m_zip_stream.next_out = &m_output_buffer[remainder];
                //                std::cerr << "after static_cast: line 219" << std::endl;
            }
        } while ( m_zip_stream.avail_in != 0 && m_err == Z_OK );
        //        std::cerr << "line 222" << std::endl;
        return m_err == Z_OK;
    }

    template<
        typename Elem,
        typename Tr,
        typename ElemA,
        typename ByteT,
        typename ByteAT
    >
    std::streamsize basic_zip_streambuf< Elem, Tr, ElemA, ByteT, ByteAT >::flush(
        int flush_mode
    )
    {
        std::streamsize written_byte_size = 0, total_written_byte_size = 0;

        int const buffer_size = static_cast< int >( pptr() - pbase() ); // amount of data currently in buffer

        m_zip_stream.next_in = (byte_buffer_type)pbase();
        m_zip_stream.avail_in = static_cast< uInt >( buffer_size * sizeof(char_type) );
        m_zip_stream.avail_out = static_cast< uInt >( m_output_buffer.size() );
        m_zip_stream.next_out = &(m_output_buffer[0]);
        size_t remainder = 0;

        //        std::cerr << "flush() " << std::endl;

        // updating crc
        m_crc = crc32(
            m_crc,
            m_zip_stream.next_in,
            m_zip_stream.avail_in
        );

        do {
            m_err = deflate( &m_zip_stream, flush_mode );
            if ( m_err == Z_OK || m_err == Z_STREAM_END ) {
                written_byte_size =
                    static_cast< std::streamsize >( m_output_buffer.size() )
                    - m_zip_stream.avail_out;
                total_written_byte_size += written_byte_size;
                // ouput buffer is full, dumping to ostream
                m_ostream.write(
                    (const char_type*) &(m_output_buffer[0]),
                    static_cast< std::streamsize >( written_byte_size / sizeof( char_type ) * sizeof( byte_type ) )
                );

                // checking if some bytes were not written.
                if ( ( remainder = written_byte_size%sizeof( char_type ) ) != 0 ) {
                    // copy to the beginning of the stream
                    std::memmove(
                        &(m_output_buffer[0]),
                        &(m_output_buffer[written_byte_size-remainder]),
                        remainder
                    );

                }

                m_zip_stream.avail_out = static_cast< uInt >( m_output_buffer.size() - remainder );
                m_zip_stream.next_out = &m_output_buffer[remainder];
            }
        } while (m_err == Z_OK);

        m_ostream.flush();
        //        std::cerr << "flush() -- done" << std::endl;

        return total_written_byte_size;
    }


    template<
        typename Elem,
        typename Tr,
        typename ElemA,
        typename ByteT,
        typename ByteAT
    >
    std::streamsize basic_zip_streambuf< Elem, Tr, ElemA, ByteT, ByteAT >::flush()
    {
        std::streamsize const total_written_byte_size = flush( Z_SYNC_FLUSH );
        return total_written_byte_size;
    }


    template<
        typename Elem,
        typename Tr,
        typename ElemA,
        typename ByteT,
        typename ByteAT
    >
    std::streamsize basic_zip_streambuf< Elem, Tr, ElemA, ByteT, ByteAT >::flush_finalize()
    {
        std::streamsize const total_written_byte_size = flush( Z_FINISH );
        return total_written_byte_size;
    }


    template<
        typename Elem,
        typename Tr,
        typename ElemA,
        typename ByteT,
        typename ByteAT
    >
    void basic_zip_streambuf< Elem, Tr, ElemA, ByteT, ByteAT >::reset_state()
    {
        m_crc = 0;

        m_zip_stream.next_in = NULL;
        m_zip_stream.avail_in = 0;
        m_zip_stream.avail_out = 0;
        m_zip_stream.next_out = NULL;

        m_err = deflateReset( &m_zip_stream );
    }


    template<
        typename Elem,
        typename Tr,
        typename ElemA,
        typename ByteT,
        typename ByteAT
    >
    basic_unzip_streambuf< Elem, Tr, ElemA, ByteT, ByteAT >::basic_unzip_streambuf(
        istream_reference istream_,
        size_t window_size_,
        size_t read_buffer_size_,
        size_t input_buffer_size_
    ) :
        m_istream( istream_ ),
        m_input_buffer( input_buffer_size_ ),
        m_buffer( read_buffer_size_ ),
        m_crc( 0 )
    {
        // setting zalloc, zfree and opaque
        m_zip_stream.zalloc = (alloc_func)0;
        m_zip_stream.zfree = (free_func)0;

        m_zip_stream.next_in = NULL;
        m_zip_stream.avail_in = 0;
        m_zip_stream.avail_out = 0;
        m_zip_stream.next_out = NULL;

        m_err = inflateInit2( &m_zip_stream, -static_cast< int >( window_size_ ) );

        this->setg(
            &( m_buffer[0] ) + 4, // beginning of putback area
            &( m_buffer[0] ) + 4, // read position
            &( m_buffer[0] ) + 4  // end position
        );
    }

    template<
        typename Elem,
        typename Tr,
        typename ElemA,
        typename ByteT,
        typename ByteAT
    >
    size_t basic_unzip_streambuf< Elem, Tr, ElemA, ByteT, ByteAT >::fill_input_buffer()
    {
        m_zip_stream.next_in = &( m_input_buffer[0] );
        m_istream.read(
            (char_type*)( &( m_input_buffer[0] ) ),
            static_cast< std::streamsize >( m_input_buffer.size() / sizeof( char_type ) )
        );
        return m_zip_stream.avail_in = m_istream.gcount() * sizeof( char_type );
    }


    template<
        typename Elem,
        typename Tr,
        typename ElemA,
        typename ByteT,
        typename ByteAT
    >
    basic_unzip_streambuf< Elem, Tr, ElemA, ByteT, ByteAT >::~basic_unzip_streambuf()
    {
        inflateEnd( &m_zip_stream );
    }

    template<
        typename Elem,
        typename Tr,
        typename ElemA,
        typename ByteT,
        typename ByteAT
    >
    typename basic_unzip_streambuf< Elem, Tr, ElemA, ByteT, ByteAT >::int_type
    basic_unzip_streambuf< Elem, Tr, ElemA, ByteT, ByteAT >::underflow()
    {
        if ( gptr() && ( gptr() < egptr() ) ) return *reinterpret_cast< unsigned char * >( gptr() );

        int n_putback = static_cast< int >( gptr() - eback() );
        if ( n_putback > 4 ) n_putback = 4;
        std::memmove(
            &( m_buffer[0] ) + ( 4 - n_putback ),
            gptr() - n_putback,
            n_putback * sizeof( char_type )
        );

        int num = unzip_from_stream(
            &( m_buffer[0] ) + 4,
            static_cast< std::streamsize >( ( m_buffer.size() - 4 ) * sizeof( char_type ) )
        );
        if ( num <= 0 ) return EOF; // ERROR or EOF

        // reset buffer pointers
        this->setg(
            &( m_buffer[0] ) + ( 4 - n_putback ), // beginning of putback area
            &( m_buffer[0] ) + 4,                 // read position
            &( m_buffer[0] ) + 4 + num            // end of buffer
        );

        // return next character
        return *reinterpret_cast< unsigned char * >( gptr() );
    }

    template<
        typename Elem,
        typename Tr,
        typename ElemA,
        typename ByteT,
        typename ByteAT
    >
    std::streamsize basic_unzip_streambuf< Elem, Tr, ElemA, ByteT, ByteAT >::unzip_from_stream(
        char_type * buffer_,
        std::streamsize buffer_size_
    )
    {
        m_zip_stream.next_out = (byte_buffer_type)buffer_;
        m_zip_stream.avail_out = static_cast< uInt >( buffer_size_ * sizeof( char_type ) );
        size_t count = m_zip_stream.avail_in;

        do {
            if ( m_zip_stream.avail_in == 0 ) count = fill_input_buffer();

            if ( m_zip_stream.avail_in ) {
                m_err = inflate( &m_zip_stream, Z_SYNC_FLUSH );
            }
        } while ( m_err == Z_OK && m_zip_stream.avail_out != 0 && count != 0 );

        // updating crc
        m_crc = crc32(
            m_crc,
            (byte_buffer_type)buffer_,
            buffer_size_ - m_zip_stream.avail_out / sizeof( char_type )
        );
        std::streamsize n_read = buffer_size_ - m_zip_stream.avail_out / sizeof( char_type );

        // check if it is the end
        if ( m_err == Z_STREAM_END ) put_back_from_zip_stream();

        return n_read;
    }

    template<
        typename Elem,
        typename Tr,
        typename ElemA,
        typename ByteT,
        typename ByteAT
    >
    void basic_unzip_streambuf< Elem, Tr, ElemA, ByteT, ByteAT >::put_back_from_zip_stream()
    {
        if ( m_zip_stream.avail_in == 0 ) return;

        m_istream.clear( std::ios::goodbit );
        m_istream.seekg(
            -static_cast< int >( m_zip_stream.avail_in ),
            std::ios_base::cur
        );

        m_zip_stream.avail_in = 0;
    }

    template<
        typename Elem,
        typename Tr,
        typename ElemA,
        typename ByteT,
        typename ByteAT
    >
    int
    basic_zip_istream< Elem, Tr, ElemA, ByteT, ByteAT >::check_header()
    {
        int method; // method byte
        int flags;  // flags byte
        uInt len;
        int c;
        int err = 0;
        z_stream & zip_stream = rdbuf()->get_zip_stream();

        // Check the gzip magic header
        for ( len = 0; len < 2; ++len ) {
            c = (int)rdbuf()->get_istream().get();
            if ( c != detail::gz_magic[len] ) {
                if ( len != 0 ) rdbuf()->get_istream().unget();
                if ( c != EOF ) {
                    rdbuf()->get_istream().unget();
                }

                err = ( zip_stream.avail_in != 0 ? Z_OK : Z_STREAM_END );
                m_is_gzip = false;
                return err;
            }
        }

        m_is_gzip = true;
        method = (int)rdbuf()->get_istream().get();
        flags = (int)rdbuf()->get_istream().get();
        if ( method != Z_DEFLATED || ( flags & detail::gz_reserved ) != 0 ) {
            err = Z_DATA_ERROR;
            return err;
        }

        // Discard time, xflags and OS code:
        for ( len = 0; len < 6; ++len ) rdbuf()->get_istream().get();

        if ( ( flags & detail::gz_extra_field ) != 0 ) {
            // skip the extra field
            len  =  (uInt)rdbuf()->get_istream().get();
            len += ((uInt)rdbuf()->get_istream().get()) << 8;
            // len is garbage if EOF but the loop below will quit anyway
            while ( len-- != 0 && rdbuf()->get_istream().get() != EOF ) {}
        }
        if ( ( flags & detail::gz_orig_name ) != 0 ) {
            // skip the original file name
            while ( ( c = rdbuf()->get_istream().get() ) != 0 && c != EOF ) {}
        }
        if ( ( flags & detail::gz_comment ) != 0 ) {
            // skip the .gz file comment
            while ( ( c = rdbuf()->get_istream().get() ) != 0 && c != EOF ) {}
        }
        if ( ( flags & detail::gz_head_crc ) != 0 ) { // skip the header crc
            for ( len = 0; len < 2; ++len ) rdbuf()->get_istream().get();
        }
        err = ( rdbuf()->get_istream().eof() ? Z_DATA_ERROR : Z_OK );

        return err;
    }

    template<
        typename Elem,
        typename Tr,
        typename ElemA,
        typename ByteT,
        typename ByteAT
    >
    void basic_zip_istream< Elem, Tr, ElemA, ByteT, ByteAT >::read_footer()
    {
        if ( m_is_gzip ) {
            read_uint32( rdbuf()->get_istream(), m_gzip_crc );
            read_uint32( rdbuf()->get_istream(), m_gzip_data_size );
        }
    }

    template<
        typename Elem,
        typename Tr,
        typename ElemA,
        typename ByteT,
        typename ByteAT
    >
    void basic_zip_ostream< Elem, Tr, ElemA, ByteT, ByteAT >::put_long_as_uint32(
        typename basic_zip_ostream< Elem, Tr, ElemA, ByteT, ByteAT >::ostream_reference out_,
        unsigned long x_
    )
    {
        // yab: 20090414, modified to conform to gzip standard where
        // trailer crc and length must both be 32-bit, otherwise there
        // is breakage in systems where 'unsigned long' is not 32-bit
        // and external archiving programs end up complaining.
        char b1, b2, b3, b4; // assuming char is 8 bits
        b1 = 0xFF & x_;
        b2 = 0xFF & ( x_ >> 8 );
        b3 = 0xFF & ( x_ >> 16 );
        b4 = 0xFF & ( x_ >> 24 );

        out_.write( &b1, 1 );
        out_.write( &b2, 1 );
        out_.write( &b3, 1 );
        out_.write( &b4, 1 );
    }

    template<
        typename Elem,
        typename Tr,
        typename ElemA,
        typename ByteT,
        typename ByteAT
    >
    void basic_zip_istream< Elem, Tr, ElemA, ByteT, ByteAT >::read_uint32(
        istream_reference in_,
        unsigned long & x_
    )
    {
        // yab: 20090414, modified to conform to gzip standard where
        // trailer crc and length must both be 32-bit, otherwise there
        // is breakage in systems where 'unsigned long' is not 32-bit
        // and external archiving programs end up complaining.
        unsigned char b1, b2, b3, b4; // assuming char is 8 bits

        in_.read( reinterpret_cast< char * >( &b1 ), 1 );
        in_.read( reinterpret_cast< char * >( &b2 ), 1 );
        in_.read( reinterpret_cast< char * >( &b3 ), 1 );
        in_.read( reinterpret_cast< char * >( &b4 ), 1 );

        x_ = b1 + ( b2 << 8 ) + ( b3 << 16 ) + ( b4 << 24 );
    }

    template<
        typename Elem,
        typename Tr,
        typename ElemA,
        typename ByteT,
        typename ByteAT
    >
    void basic_zip_ostream< Elem, Tr, ElemA, ByteT, ByteAT >::add_header()
    {
        char_type zero = 0;

        rdbuf()->get_ostream()
            .put( static_cast< char_type >( detail::gz_magic[0] ) )
            .put( static_cast< char_type >( detail::gz_magic[1] ) )
            .put( static_cast< char_type >( Z_DEFLATED ) )
            .put( zero ) // flags
            .put( zero ).put( zero ).put( zero ).put( zero ) // time
            .put( zero ) // xflags
            .put( static_cast< char_type >( OS_CODE ) );
    }

    template<
        typename Elem,
        typename Tr,
        typename ElemA,
        typename ByteT,
        typename ByteAT
    >
    void basic_zip_ostream< Elem, Tr, ElemA, ByteT, ByteAT >::add_footer()
    {
        put_long_as_uint32( rdbuf()->get_ostream(), rdbuf()->get_crc() );
        put_long_as_uint32( rdbuf()->get_ostream(), rdbuf()->get_in_size() );
    }


} // namespace zlib_stream


#endif // INCLUDED_utility_io_zipstream_IPP
