/*
bzip2stream Library License:
--------------------------

The zlib/libpng License Copyright (c) 2003 Jonathan de Halleux.

This software is provided 'as-is', without any express or implied warranty. In no event will the authors be held liable for any damages arising from the use of this software.

Permission is granted to anyone to use this software for any purpose, including commercial applications, and to alter it and redistribute it freely, subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.

2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.

3. This notice may not be removed or altered from any source distribution

Author: Jonathan de Halleux, dehalleux@pelikhan.com, 2003
*/
#ifndef BZIP2STREAM_IPP
#define BZIP2STREAM_IPP

#include "bzip2stream.h"

namespace bzip2_stream{

    template<
        typename Elem,
        typename Tr,
        typename ElemA,
        typename ByteT,
        typename ByteAT
    >
    basic_bzip2_streambuf<
        Elem,Tr,ElemA,ByteT,ByteAT
        >:: basic_bzip2_streambuf(
        ostream_reference ostream_,
        size_t block_size_100k_,
        size_t verbosity_,
        size_t work_factor_,
        size_t buffer_size_
        )
    :
        m_ostream(ostream_),
        m_output_buffer(buffer_size_,0),
        m_buffer(buffer_size_,0)
    {
        m_bzip2_stream.bzalloc=NULL;
        m_bzip2_stream.bzfree=NULL;

        m_bzip2_stream.next_in=NULL;
        m_bzip2_stream.avail_in=0;
        m_bzip2_stream.avail_out=0;
        m_bzip2_stream.next_out=NULL;

        m_err=BZ2_bzCompressInit(
            &m_bzip2_stream,
            std::min( 9, static_cast<int>(block_size_100k_) ),
            std::min( 4, static_cast<int>(verbosity_) ),
            std::min( 250, static_cast<int>(work_factor_) )
            );

        this->setp( &(m_buffer[0]), &(m_buffer[m_buffer.size()-1]));
    }

    template<
        typename Elem,
        typename Tr,
        typename ElemA,
        typename ByteT,
        typename ByteAT
    >
    basic_bzip2_streambuf<
        Elem,Tr,ElemA,ByteT,ByteAT
        >::~basic_bzip2_streambuf()
    {
        flush(BZ_FINISH);
        m_ostream.flush();
        m_err=BZ2_bzCompressEnd(&m_bzip2_stream);
    }

    template<
        typename Elem,
        typename Tr,
        typename ElemA,
        typename ByteT,
        typename ByteAT
    >
    int basic_bzip2_streambuf<
        Elem,Tr,ElemA,ByteT,ByteAT
        >::sync ()
    {
        if ( this->pptr() && this->pptr() > this->pbase())
        {
            int c = overflow( EOF);

            if ( c == EOF)
                return -1;
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
    typename basic_bzip2_streambuf<
        Elem,Tr,ElemA,ByteT,ByteAT
        >::int_type
            basic_bzip2_streambuf<
                Elem,Tr,ElemA,ByteT,ByteAT
                >::overflow (
                    typename basic_bzip2_streambuf<
                        Elem,Tr,ElemA,ByteT,ByteAT
                        >::int_type c
                    )
    {
        int w = static_cast<int>(this->pptr() - this->pbase());
        if (c != EOF) {
             *this->pptr() = c;
             ++w;
         }
         if ( bzip2_to_stream( this->pbase(), w)) {
             this->setp( this->pbase(), this->epptr());
             return c;
         } else
             return EOF;
    }

    template<
        typename Elem,
        typename Tr,
        typename ElemA,
        typename ByteT,
        typename ByteAT
    >
    bool basic_bzip2_streambuf<
        Elem,Tr,ElemA,ByteT,ByteAT
        >::bzip2_to_stream(
            typename basic_bzip2_streambuf<
                Elem,Tr,ElemA,ByteT,ByteAT
                >::char_type* buffer_,
            std::streamsize buffer_size_
            )
    {
        std::streamsize written_byte_size=0, total_written_byte_size = 0;

        m_bzip2_stream.next_in=(byte_buffer_type)buffer_;
        m_bzip2_stream.avail_in=buffer_size_*sizeof(char_type);
        m_bzip2_stream.avail_out=static_cast<unsigned int>(m_output_buffer.size());
        m_bzip2_stream.next_out=&(m_output_buffer[0]);
        size_t remainder=0;

        do
        {
            m_err = BZ2_bzCompress (&m_bzip2_stream, BZ_RUN );

            if (m_err == BZ_RUN_OK  || m_err == BZ_STREAM_END)
            {
                written_byte_size= static_cast<std::streamsize>(m_output_buffer.size()) - m_bzip2_stream.avail_out;
                total_written_byte_size+=written_byte_size;
                // ouput buffer is full, dumping to ostream
                m_ostream.write(
                    (const char_type*) &(m_output_buffer[0]),
                    static_cast<std::streamsize>( written_byte_size/sizeof(char_type) )
                    );

                // checking if some bytes were not written.
                if ( (remainder = written_byte_size%sizeof(char_type))!=0)
                {
                    // copy to the beginning of the stream
                    std::memmove(
                        &(m_output_buffer[0]),
                        &(m_output_buffer[written_byte_size-remainder]),
                        remainder);

                }

                m_bzip2_stream.avail_out=static_cast<unsigned int>(m_output_buffer.size()-remainder);
                m_bzip2_stream.next_out=&m_output_buffer[remainder];
            }
        }
        while (m_bzip2_stream.avail_in != 0 && m_err == BZ_RUN_OK);

        return m_err == BZ_RUN_OK || m_err == BZ_FLUSH_OK;
    }

    template<
        typename Elem,
        typename Tr,
        typename ElemA,
        typename ByteT,
        typename ByteAT
    >
    std::streamsize basic_bzip2_streambuf<
        Elem,Tr,ElemA,ByteT,ByteAT
        >::flush(int flush_mode)
    {
        std::streamsize written_byte_size=0, total_written_byte_size=0;

        int const buffer_size = static_cast< int >( pptr() - pbase() ); // amount of data currently in buffer

        m_bzip2_stream.next_in=(byte_buffer_type)pbase();
        m_bzip2_stream.avail_in=static_cast< uInt >(buffer_size*sizeof(char_type));
        m_bzip2_stream.avail_out=static_cast< uInt >(m_output_buffer.size());
        m_bzip2_stream.next_out=&(m_output_buffer[0]);
        size_t remainder=0;

        do
        {
            m_err = BZ2_bzCompress (&m_bzip2_stream, flush_mode);
            if (m_err == BZ_FINISH_OK || m_err == BZ_STREAM_END)
            {
                written_byte_size=
                    static_cast<std::streamsize>(m_output_buffer.size())
                    - m_bzip2_stream.avail_out;
                total_written_byte_size+=written_byte_size;
                // ouput buffer is full, dumping to ostream
                m_ostream.write(
                    (const char_type*) &(m_output_buffer[0]),
                    static_cast<std::streamsize>( written_byte_size/sizeof(char_type)*sizeof(char) )
                    );

                // checking if some bytes were not written.
                if ( (remainder = written_byte_size%sizeof(char_type))!=0)
                {
                    // copy to the beginning of the stream
                    std::memmove(
                        &(m_output_buffer[0]),
                        &(m_output_buffer[written_byte_size-remainder]),
                        remainder);

                }

                m_bzip2_stream.avail_out=static_cast<unsigned int>(m_output_buffer.size()-remainder);
                m_bzip2_stream.next_out=&(m_output_buffer[remainder]);
            }
        } while (m_err == BZ_FINISH_OK);

        m_ostream.flush();

        return total_written_byte_size;
    }

    template<
        typename Elem,
        typename Tr,
        typename ElemA,
        typename ByteT,
        typename ByteAT
    >
    basic_unbzip2_streambuf<
        Elem,Tr,ElemA,ByteT,ByteAT
        >::basic_unbzip2_streambuf(
            istream_reference istream_,
            size_t verbosity_,
            bool small_,
            size_t read_buffer_size_,
            size_t input_buffer_size_
    )
    :
        m_istream(istream_),
        m_input_buffer(input_buffer_size_),
        m_buffer(read_buffer_size_)
    {
        // setting zalloc, zfree and opaque
        m_bzip2_stream.bzalloc=NULL;
        m_bzip2_stream.bzfree=NULL;

        m_bzip2_stream.next_in=NULL;
        m_bzip2_stream.avail_in=0;
        m_bzip2_stream.avail_out=0;
        m_bzip2_stream.next_out=NULL;


        m_err=BZ2_bzDecompressInit (
            &m_bzip2_stream,
            std::min(4, static_cast<int>(verbosity_)),
            static_cast<int>(small_)
        );

        this->setg(
            &(m_buffer[0])+4,     // beginning of putback area
            &(m_buffer[0])+4,     // read position
            &(m_buffer[0])+4);    // end position
    }

    template<
        typename Elem,
        typename Tr,
        typename ElemA,
        typename ByteT,
        typename ByteAT
    >
    size_t basic_unbzip2_streambuf<
        Elem,Tr,ElemA,ByteT,ByteAT
        >::fill_input_buffer()
    {
        m_bzip2_stream.next_in=&(m_input_buffer[0]);
        m_istream.read(
            (char_type*)(&(m_input_buffer[0])),
            static_cast<std::streamsize>(m_input_buffer.size()/sizeof(char_type))
            );
        return m_bzip2_stream.avail_in=m_istream.gcount()*sizeof(char_type);
    }

    template<
        typename Elem,
        typename Tr,
        typename ElemA,
        typename ByteT,
        typename ByteAT
    >
    void basic_unbzip2_streambuf<
        Elem,Tr,ElemA,ByteT,ByteAT
        >::put_back_from_bzip2_stream()
    {
        if (m_bzip2_stream.avail_in==0)
            return;

        m_istream.clear( std::ios::goodbit );
        m_istream.seekg(
            -static_cast<int>(m_bzip2_stream.avail_in),
            std::ios_base::cur
            );

        m_bzip2_stream.avail_in=0;
    }


    template<
        typename Elem,
        typename Tr,
        typename ElemA,
        typename ByteT,
        typename ByteAT
    >
    basic_unbzip2_streambuf<
        Elem,Tr,ElemA,ByteT,ByteAT
        >::~basic_unbzip2_streambuf()
    {
        BZ2_bzDecompressEnd(&m_bzip2_stream);
    }

    template<
        typename Elem,
        typename Tr,
        typename ElemA,
        typename ByteT,
        typename ByteAT
    >
    typename basic_unbzip2_streambuf<
        Elem,Tr,ElemA,ByteT,ByteAT
        >::int_type
            basic_unbzip2_streambuf<
                Elem,Tr,ElemA,ByteT,ByteAT
                >::underflow()
    {
        if ( this->gptr() && ( this->gptr() < this->egptr()))
            return * reinterpret_cast<unsigned char *>( this->gptr());

       int n_putback = static_cast<int>(this->gptr() - this->eback());
       if ( n_putback > 4)
          n_putback = 4;
       std::memmove(
            &(m_buffer[0]) + (4 - n_putback),
            this->gptr() - n_putback,
            n_putback*sizeof(char_type)
            );

       int num = unbzip2_from_stream(
           &(m_buffer[0])+4,
           static_cast<std::streamsize>((m_buffer.size()-4)*sizeof(char_type))
           );
        if (num <= 0) // ERROR or EOF
           return EOF;

        // reset buffer pointers
        this->setg(
              &(m_buffer[0]) + (4 - n_putback),   // beginning of putback area
              &(m_buffer[0]) + 4,                 // read position
              &(m_buffer[0]) + 4 + num);          // end of buffer

         // return next character
         return* reinterpret_cast<unsigned char *>( this->gptr());
     }


    template<
        typename Elem,
        typename Tr,
        typename ElemA,
        typename ByteT,
        typename ByteAT
    >
    std::streamsize basic_unbzip2_streambuf<
        Elem,Tr,ElemA,ByteT,ByteAT
        >::unbzip2_from_stream(
            char_type* buffer_,
            std::streamsize buffer_size_
            )
    {
        m_bzip2_stream.next_out=(byte_buffer_type)buffer_;
        m_bzip2_stream.avail_out=buffer_size_*sizeof(char_type);
        size_t count =m_bzip2_stream.avail_in;

        do
        {
            if (m_bzip2_stream.avail_in==0)
                count=fill_input_buffer();

            if (m_bzip2_stream.avail_in)
            {
                m_err = BZ2_bzDecompress( &m_bzip2_stream );
            }
        } while (m_err==BZ_OK && m_bzip2_stream.avail_out != 0 && count != 0);

        if (m_err == BZ_STREAM_END)
            put_back_from_bzip2_stream();

        return buffer_size_ - m_bzip2_stream.avail_out/sizeof(char_type);
    }


} // zlib_sream

#endif
