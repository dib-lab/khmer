/*
zipstream Library License:
--------------------------

The zlib/libpng License Copyright (c) 2003 Jonathan de Halleux.

This software is provided 'as-is', without any express or implied warranty. In no event will the authors be held liable for any damages arising from the use of this software.

Permission is granted to anyone to use this software for any purpose, including commercial applications, and to alter it and redistribute it freely, subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.

2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.

3. This notice may not be removed or altered from any source distribution

Author: Jonathan de Halleux, dehalleux@pelikhan.com, 2003   (original zlib stream)
Author: David Weese, dave.weese@gmail.com, 2014             (extension to parallel block-wise compression in bgzf format)
*/
#ifndef BGZFSTREAM_IPP
#define BGZFSTREAM_IPP

#include "bgzfstream.h"
#include <sstream>

namespace seqan {

namespace detail{
    const int gz_magic[2] = {0x1f, 0x8b}; /* gzip magic header */

    /* gzip flag byte */
    const int gz_ascii_flag =  0x01; /* bit 0 set: file probably ascii text */
    const int gz_head_crc    = 0x02; /* bit 1 set: header CRC present */
    const int gz_extra_field = 0x04; /* bit 2 set: extra field present */
    const int gz_orig_name  =  0x08; /* bit 3 set: original file name present */
    const int gz_comment    =  0x10; /* bit 4 set: file comment present */
    const int gz_reserved   =  0xE0; /* bits 5..7: reserved */

}

    template<
        typename Elem,
        typename Tr,
        typename ElemA,
        typename ByteT,
        typename ByteAT
    >
    int _checkGZHeader(basic_unbgzf_streambuf<Elem, Tr, ElemA, ByteT, ByteAT> *buf)
    {
        int method; /* method byte */
        int flags;  /* flags byte */
        uInt len;
        int c;
        int err=0;
        z_stream& bgzf_stream = buf->get_bgzf_stream();

        std::basic_istream<Elem, Tr> &istream = buf->get_istream();
        bool m_is_gzip;

        /* Check the gzip magic header */
         for (len = 0; len < 2; len++)
         {
            c = (int)istream.get();
            if (c != detail::gz_magic[len])
            {
                if (len != 0)
                    istream.unget();
                if (c!= EOF)
                {
                    istream.unget();
                }

                err = bgzf_stream.avail_in != 0 ? Z_OK : Z_STREAM_END;
                m_is_gzip = false;
                return err;
            }
        }

        m_is_gzip = true;
        method = (int)istream.get();
        flags = (int)istream.get();
        if (method != Z_DEFLATED || (flags & detail::gz_reserved) != 0)
        {
            err = Z_DATA_ERROR;
            return err;
        }

        /* Discard time, xflags and OS code: */
        for (len = 0; len < 6; len++)
            istream.get();

        if ((flags & detail::gz_extra_field) != 0)
        {
            /* skip the extra field */
            len  =  (uInt)istream.get();
            len += ((uInt)istream.get())<<8;
            /* len is garbage if EOF but the loop below will quit anyway */
            while (len-- != 0 && istream.get() != EOF) ;
        }
        if ((flags & detail::gz_orig_name) != 0)
        {
            /* skip the original file name */
            while ((c = istream.get()) != 0 && c != EOF) ;
        }
        if ((flags & detail::gz_comment) != 0)
        {
            /* skip the .gz file comment */
            while ((c = istream.get()) != 0 && c != EOF) ;
        }
        if ((flags & detail::gz_head_crc) != 0)
        {  /* skip the header crc */
            for (len = 0; len < 2; len++)
                istream.get();
        }
        err = istream.eof() ? Z_DATA_ERROR : Z_OK;

        return err;
    }


    template<
        typename Elem,
        typename Tr,
        typename ElemA,
        typename ByteT,
        typename ByteAT
    >
    int basic_bgzf_istream<
        Elem,Tr,ElemA,ByteT,ByteAT
        >::check_header()
    {
        return _checkGZHeader(this->rdbuf());
    }

    template<
        typename Elem,
        typename Tr
    >
    void _putBinaryLong(std::basic_ostream<Elem,Tr> & out_, unsigned long x_)
    {
        static const int size_ul = sizeof(unsigned long);
        static const int size_c = sizeof(typename Tr::char_type);
        static const int n_end = size_ul/size_c;
        out_.write(reinterpret_cast<typename Tr::char_type const*>(&x_), n_end);
    }

    template<
        typename Elem,
        typename Tr,
        typename ElemA,
        typename ByteT,
        typename ByteAT
    >
    void basic_bgzf_istream<
            Elem,Tr,ElemA,ByteT,ByteAT
            >::read_long(
                istream_reference in_,
            unsigned long& x_
            )
    {
        static const int size_ul = sizeof(unsigned long);
        static const int size_c = sizeof(typename Tr::char_type);
        static const int n_end = size_ul/size_c;
        in_.read(reinterpret_cast<char*>(&x_),n_end);
    }

    template<
        typename Elem,
        typename Tr,
        typename ElemA,
        typename ByteT,
        typename ByteAT
    >
    void _addGZHeader(basic_bgzf_streambuf<Elem, Tr, ElemA, ByteT, ByteAT> *buf)
    {
        typename Tr::char_type zero=0;

        buf->get_ostream()
            .put(static_cast<typename Tr::char_type>(detail::gz_magic[0]))
            .put(static_cast<typename Tr::char_type>(detail::gz_magic[1]))
            .put(static_cast<typename Tr::char_type>(Z_DEFLATED))
            .put(zero) //flags
            .put(zero).put(zero).put(zero).put(zero) // time
            .put(zero) //xflags
            .put(static_cast<typename Tr::char_type>(OS_CODE));
    }

    template<
        typename Elem,
        typename Tr,
        typename ElemA,
        typename ByteT,
        typename ByteAT
    >
    void _addGZFooter(basic_bgzf_streambuf<Elem, Tr, ElemA, ByteT, ByteAT> *buf)
    {
        _putBinaryLong( buf->get_ostream(), buf->get_crc() );
        _putBinaryLong( buf->get_ostream(), buf->get_in_size() );
    }

}  // namespace seqan

#endif
