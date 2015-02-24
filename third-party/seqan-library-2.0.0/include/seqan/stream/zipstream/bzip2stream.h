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

#ifndef BZIP2STREAM_HPP
#define BZIP2STREAM_HPP

#define BZ_NO_STDIO

#include <bzlib.h>
#include <vector>
#include <iostream>
#include <algorithm>

namespace bzip2_stream{

const size_t default_buffer_size = 4096;

template<
    typename Elem,
    typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = char,
    typename ByteAT = std::allocator<ByteT>
>
class basic_bzip2_streambuf :
    public std::basic_streambuf<Elem, Tr>
{
public:
    typedef std::basic_streambuf< Elem, Tr > basic_streambuf_type;
    typedef std::basic_ostream<Elem, Tr>& ostream_reference;
    typedef ElemA char_allocator_type;
    typedef ByteT byte_type;
    typedef ByteAT byte_allocator_type;
    typedef byte_type* byte_buffer_type;
    typedef typename Tr::char_type char_type;
    typedef typename Tr::int_type int_type;
    typedef std::vector<byte_type, byte_allocator_type > byte_vector_type;
    typedef std::vector<char_type, char_allocator_type > char_vector_type;

    using basic_streambuf_type::epptr;
    using basic_streambuf_type::pbase;
    using basic_streambuf_type::pptr;

    basic_bzip2_streambuf(
        ostream_reference ostream_,
        size_t block_size_100k_ ,
        size_t verbosity_ ,
        size_t work_factor_,
        size_t buffer_size_
        );

    ~basic_bzip2_streambuf();

    int sync ();
    int_type overflow (int_type c);

    std::streamsize flush(int flush_mode);
    int get_zerr() const
    {    return m_err;};
    __uint64 get_in_size() const
    {
        return ((__uint64)m_bzip2_stream.total_in_hi32 << 32)
                + m_bzip2_stream.total_in_lo32;
    }
    __uint64 get_out_size() const
    {
        return ((__uint64)m_bzip2_stream.total_out_hi32 << 32)
                + m_bzip2_stream.total_out_lo32;
    }
private:
    bool bzip2_to_stream( char_type*, std::streamsize);
    size_t fill_input_buffer();

    ostream_reference m_ostream;
    bz_stream m_bzip2_stream;
    int m_err;
    byte_vector_type m_output_buffer;
    char_vector_type m_buffer;
};

template<
    typename Elem,
    typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = char,
    typename ByteAT = std::allocator<ByteT>
>
class basic_unbzip2_streambuf :
    public std::basic_streambuf<Elem, Tr>
{
public:
    typedef std::basic_istream<Elem, Tr>& istream_reference;
    typedef ElemA char_allocator_type;
    typedef ByteT byte_type;
    typedef ByteAT byte_allocator_type;
    typedef byte_type* byte_buffer_type;
    typedef typename Tr::char_type char_type;
    typedef typename Tr::int_type int_type;
    typedef std::vector<byte_type, byte_allocator_type > byte_vector_type;
    typedef std::vector<char_type, char_allocator_type > char_vector_type;

    basic_unbzip2_streambuf(
        istream_reference istream_,
        size_t verbosity_,
        bool small_,
        size_t read_buffer_size_,
        size_t input_buffer_size_
        );

    ~basic_unbzip2_streambuf();

    int_type underflow();

    istream_reference get_istream()        {    return m_istream;};
    bz_stream& get_bzip2_stream()        {    return m_bzip2_stream;};
    int get_zerr() const                {    return m_err;};
private:
    std::streamsize unbzip2_from_stream( char_type*, std::streamsize);
    void put_back_from_bzip2_stream();
    size_t fill_input_buffer();

    istream_reference m_istream;
    bz_stream m_bzip2_stream;
    int m_err;
    byte_vector_type m_input_buffer;
    char_vector_type m_buffer;
};

template<
    typename Elem,
    typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = char,
    typename ByteAT = std::allocator<ByteT>
>
class basic_bzip2_ostreambase : virtual public std::basic_ios<Elem,Tr>
{
public:
    typedef std::basic_ostream<Elem, Tr>& ostream_reference;
    typedef basic_bzip2_streambuf<
        Elem,Tr,ElemA,ByteT,ByteAT> bzip2_streambuf_type;

    basic_bzip2_ostreambase(
        ostream_reference ostream_,
        size_t block_size_100k_ ,
        size_t verbosity_ ,
        size_t work_factor_,
        size_t buffer_size_
        )
        : m_buf(ostream_,block_size_100k_, verbosity_, work_factor_, buffer_size_)
    {
        this->init(&m_buf );
    };

    bzip2_streambuf_type* rdbuf() { return &m_buf; };

private:
    bzip2_streambuf_type m_buf;
};

template<
    typename Elem,
    typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = char,
    typename ByteAT = std::allocator<ByteT>
>
class basic_bzip2_istreambase : virtual public std::basic_ios<Elem,Tr>
{
public:
    typedef std::basic_istream<Elem, Tr>& istream_reference;
    typedef basic_unbzip2_streambuf<
        Elem,Tr,ElemA,ByteT,ByteAT> unbzip2_streambuf_type;

    basic_bzip2_istreambase(
        istream_reference ostream_,
        size_t verbosity_,
        bool small_,
        size_t read_buffer_size_,
        size_t input_buffer_size_
        )
        : m_buf(
            ostream_,
            verbosity_,
            small_,
            read_buffer_size_,
            input_buffer_size_
            )
    {
        this->init(&m_buf );
    };

    unbzip2_streambuf_type* rdbuf() { return &m_buf; };

private:
    unbzip2_streambuf_type m_buf;
};

template<
    typename Elem,
    typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = char,
    typename ByteAT = std::allocator<ByteT>
>
class basic_bzip2_ostream :
    public basic_bzip2_ostreambase<Elem,Tr,ElemA,ByteT,ByteAT>,
    public std::basic_ostream<Elem,Tr>
{
public:
    typedef basic_bzip2_ostreambase<
        Elem,Tr,ElemA,ByteT,ByteAT> bzip2_ostreambase_type;
    typedef std::basic_ostream<Elem,Tr> ostream_type;
    typedef ostream_type& ostream_reference;

    basic_bzip2_ostream(
        ostream_reference ostream_,
        size_t block_size_100k_ = 9,
        size_t verbosity_ = 0,
        size_t work_factor_ = 30,
        size_t buffer_size_ = default_buffer_size
        )
    :
        bzip2_ostreambase_type(ostream_,block_size_100k_, verbosity_, work_factor_,buffer_size_),
        ostream_type(bzip2_ostreambase_type::rdbuf())
    {

    };

    basic_bzip2_ostream& add_header();
    basic_bzip2_ostream& zflush()
    {
        this->flush(); this->rdbuf()->flush(); return *this;
    };

#ifdef _WIN32
private:
    void _Add_vtordisp1() { } // Required to avoid VC++ warning C4250
    void _Add_vtordisp2() { } // Required to avoid VC++ warning C4250
#endif
};

template<
    typename Elem,
    typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = char,
    typename ByteAT = std::allocator<ByteT>
>
class basic_bzip2_istream :
    public basic_bzip2_istreambase<Elem,Tr,ElemA,ByteT,ByteAT>,
    public std::basic_istream<Elem,Tr>
{
public:
    typedef basic_bzip2_istreambase<
        Elem,Tr,ElemA,ByteT,ByteAT> bzip2_istreambase_type;
    typedef std::basic_istream<Elem,Tr> istream_type;
    typedef istream_type& istream_reference;
    typedef unsigned char byte_type;

    basic_bzip2_istream(
        istream_reference istream_,
        size_t verbosity_ = 0,
        bool small_ = false,
        size_t read_buffer_size_ = default_buffer_size,
        size_t input_buffer_size_ = default_buffer_size
        )
      :
        bzip2_istreambase_type(istream_,verbosity_, small_, read_buffer_size_, input_buffer_size_),
        istream_type(bzip2_istreambase_type::rdbuf())
    {};
#ifdef _WIN32
private:
    void _Add_vtordisp1() { } // Required to avoid VC++ warning C4250
    void _Add_vtordisp2() { } // Required to avoid VC++ warning C4250
#endif
};

typedef basic_bzip2_ostream<char> bzip2_ostream;
typedef basic_bzip2_ostream<wchar_t> bzip2_wostream;
typedef basic_bzip2_istream<char> bzip2_istream;
typedef basic_bzip2_istream<wchar_t> bzip2_wistream;

} // bzip2_stream

#include "bzip2stream_impl.h"

#endif

