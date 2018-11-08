/***************************************************************************
 *  include/stxxl/bits/common/binary_buffer.h
 *
 *  Classes binary_buffer and binary_reader to construct data blocks with
 *  variable length content. Programs construct blocks using
 *  binary_buffer::put<type>() and read them using
 *  binary_reader::get<type>(). The operation sequences should match.
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2013-2014 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_COMMON_BINARY_BUFFER_HEADER
#define STXXL_COMMON_BINARY_BUFFER_HEADER

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

#include <stxxl/bits/namespace.h>
#include <stxxl/bits/common/types.h>

STXXL_BEGIN_NAMESPACE

//! \addtogroup support
//! \{

/*!
 * binary_buffer represents a dynamically growable area of memory, which can be
 * modified by appending integral data types via put() and other basic
 * operations.
 */
class binary_buffer
{
protected:
    //! Allocated buffer pointer.
    char* m_data;

    //! Size of valid data.
    size_t m_size;

    //! Total capacity of buffer.
    size_t m_capacity;

public:
    //! Create a new empty object
    inline binary_buffer()
        : m_data(NULL), m_size(0), m_capacity(0)
    { }

    //! Copy-Constructor, duplicates memory content.
    inline binary_buffer(const binary_buffer& other)
        : m_data(NULL), m_size(0), m_capacity(0)
    {
        assign(other);
    }

    //! Constructor, copy memory area.
    inline binary_buffer(const void* data, size_t n)
        : m_data(NULL), m_size(0), m_capacity(0)
    {
        assign(data, n);
    }

    //! Constructor, create object with n bytes pre-allocated.
    inline binary_buffer(size_t n)
        : m_data(NULL), m_size(0), m_capacity(0)
    {
        alloc(n);
    }

    //! Constructor from std::string, copies string content.
    inline binary_buffer(const std::string& str)
        : m_data(NULL), m_size(0), m_capacity(0)
    {
        assign(str.data(), str.size());
    }

    //! Destroys the memory space.
    inline ~binary_buffer()
    {
        dealloc();
    }

    //! Return a pointer to the currently kept memory area.
    inline const char * data() const
    {
        return m_data;
    }

    //! Return a writeable pointer to the currently kept memory area.
    inline char * data()
    {
        return m_data;
    }

    //! Return the currently used length in bytes.
    inline size_t size() const
    {
        return m_size;
    }

    //! Return the currently allocated buffer capacity.
    inline size_t capacity() const
    {
        return m_capacity;
    }

    //! Explicit conversion to std::string (copies memory of course).
    inline std::string str() const
    {
        return std::string(reinterpret_cast<const char*>(m_data), m_size);
    }

    //! Set the valid bytes in the buffer, use if the buffer is filled
    //! directly.
    inline binary_buffer & set_size(size_t n)
    {
        assert(n <= m_capacity);
        m_size = n;

        return *this;
    }

    //! Make sure that at least n bytes are allocated.
    inline binary_buffer & alloc(size_t n)
    {
        if (m_capacity < n)
        {
            m_capacity = n;
            m_data = static_cast<char*>(realloc(m_data, m_capacity));
        }

        return *this;
    }

    //! Deallocates the kept memory space (we use dealloc() instead of free()
    //! as a name, because sometimes "free" is replaced by the preprocessor)
    inline binary_buffer & dealloc()
    {
        if (m_data) free(m_data);
        m_data = NULL;
        m_size = m_capacity = 0;

        return *this;
    }

    //! Detach the memory from the object, returns the memory pointer.
    inline const char * detach()
    {
        const char* data = m_data;
        m_data = NULL;
        m_size = m_capacity = 0;
        return data;
    }

    //! Clears the memory contents, does not deallocate the memory.
    inline binary_buffer & clear()
    {
        m_size = 0;
        return *this;
    }

    //! Copy a memory range into the buffer, overwrites all current
    //! data. Roughly equivalent to clear() followed by append().
    inline binary_buffer & assign(const void* data, size_t len)
    {
        if (len > m_capacity) alloc(len);

        memcpy(m_data, data, len);
        m_size = len;

        return *this;
    }

    //! Copy the contents of another buffer object into this buffer, overwrites
    //! all current data. Roughly equivalent to clear() followed by append().
    inline binary_buffer & assign(const binary_buffer& other)
    {
        if (&other != this)
            assign(other.data(), other.size());

        return *this;
    }

    //! Assignment operator: copy other's memory range into buffer.
    inline binary_buffer& operator = (const binary_buffer& other)
    {
        if (&other != this)
            assign(other.data(), other.size());

        return *this;
    }

    //! Align the size of the buffer to a multiple of n. Fills up with 0s.
    inline binary_buffer & align(size_t n)
    {
        assert(n > 0);
        size_t rem = m_size % n;
        if (rem != 0)
        {
            size_t add = n - rem;
            if (m_size + add > m_capacity) dynalloc(m_size + add);
            memset(m_data + m_size, 0, add);
            m_size += add;
        }
        assert((m_size % n) == 0);

        return *this;
    }

    //! Dynamically allocate more memory. At least n bytes will be available,
    //! probably more to compensate future growth.
    inline binary_buffer & dynalloc(size_t n)
    {
        if (m_capacity < n)
        {
            // place to adapt the buffer growing algorithm as need.
            size_t newsize = m_capacity;

            while (newsize < n) {
                if (newsize < 256) newsize = 512;
                else if (newsize < 1024 * 1024) newsize = 2 * newsize;
                else newsize += 1024 * 1024;
            }

            alloc(newsize);
        }

        return *this;
    }

    // *** Appending Write Functions ***

    //! Append a memory range to the buffer
    inline binary_buffer & append(const void* data, size_t len)
    {
        if (m_size + len > m_capacity) dynalloc(m_size + len);

        memcpy(m_data + m_size, data, len);
        m_size += len;

        return *this;
    }

    //! Append the contents of a different buffer object to this one.
    inline binary_buffer & append(const class binary_buffer& bb)
    {
        return append(bb.data(), bb.size());
    }

    //! Append to contents of a std::string, excluding the null (which isn't
    //! contained in the string size anyway).
    inline binary_buffer & append(const std::string& s)
    {
        return append(s.data(), s.size());
    }

    //! Put (append) a single item of the template type T to the buffer. Be
    //! careful with implicit type conversions!
    template <typename Type>
    inline binary_buffer & put(const Type item)
    {
        if (m_size + sizeof(Type) > m_capacity) dynalloc(m_size + sizeof(Type));

        *reinterpret_cast<Type*>(m_data + m_size) = item;
        m_size += sizeof(Type);

        return *this;
    }

    //! Append a varint to the buffer.
    inline binary_buffer & put_varint(uint32 v)
    {
        if (v < 128) {
            put<uint8>(uint8(v));
        }
        else if (v < 128 * 128) {
            put<uint8>((uint8)(((v >> 0) & 0x7F) | 0x80));
            put<uint8>((uint8)((v >> 7) & 0x7F));
        }
        else if (v < 128 * 128 * 128) {
            put<uint8>((uint8)(((v >> 0) & 0x7F) | 0x80));
            put<uint8>((uint8)(((v >> 7) & 0x7F) | 0x80));
            put<uint8>((uint8)((v >> 14) & 0x7F));
        }
        else if (v < 128 * 128 * 128 * 128) {
            put<uint8>((uint8)(((v >> 0) & 0x7F) | 0x80));
            put<uint8>((uint8)(((v >> 7) & 0x7F) | 0x80));
            put<uint8>((uint8)(((v >> 14) & 0x7F) | 0x80));
            put<uint8>((uint8)((v >> 21) & 0x7F));
        }
        else {
            put<uint8>((uint8)(((v >> 0) & 0x7F) | 0x80));
            put<uint8>((uint8)(((v >> 7) & 0x7F) | 0x80));
            put<uint8>((uint8)(((v >> 14) & 0x7F) | 0x80));
            put<uint8>((uint8)(((v >> 21) & 0x7F) | 0x80));
            put<uint8>((uint8)((v >> 28) & 0x7F));
        }

        return *this;
    }

    //! Append a varint to the buffer.
    inline binary_buffer & put_varint(int v)
    {
        return put_varint((uint32)v);
    }

    //! Append a varint to the buffer.
    inline binary_buffer & put_varint(uint64 v)
    {
        if (v < 128) {
            put<uint8>(uint8(v));
        }
        else if (v < 128 * 128) {
            put<uint8>((uint8)(((v >> 00) & 0x7F) | 0x80));
            put<uint8>((uint8)((v >> 07) & 0x7F));
        }
        else if (v < 128 * 128 * 128) {
            put<uint8>((uint8)(((v >> 00) & 0x7F) | 0x80));
            put<uint8>((uint8)(((v >> 07) & 0x7F) | 0x80));
            put<uint8>((uint8)((v >> 14) & 0x7F));
        }
        else if (v < 128 * 128 * 128 * 128) {
            put<uint8>((uint8)(((v >> 00) & 0x7F) | 0x80));
            put<uint8>((uint8)(((v >> 07) & 0x7F) | 0x80));
            put<uint8>((uint8)(((v >> 14) & 0x7F) | 0x80));
            put<uint8>((uint8)((v >> 21) & 0x7F));
        }
        else if (v < ((uint64)128) * 128 * 128 * 128 * 128) {
            put<uint8>((uint8)(((v >> 00) & 0x7F) | 0x80));
            put<uint8>((uint8)(((v >> 07) & 0x7F) | 0x80));
            put<uint8>((uint8)(((v >> 14) & 0x7F) | 0x80));
            put<uint8>((uint8)(((v >> 21) & 0x7F) | 0x80));
            put<uint8>((uint8)((v >> 28) & 0x7F));
        }
        else if (v < ((uint64)128) * 128 * 128 * 128 * 128 * 128) {
            put<uint8>((uint8)(((v >> 00) & 0x7F) | 0x80));
            put<uint8>((uint8)(((v >> 07) & 0x7F) | 0x80));
            put<uint8>((uint8)(((v >> 14) & 0x7F) | 0x80));
            put<uint8>((uint8)(((v >> 21) & 0x7F) | 0x80));
            put<uint8>((uint8)(((v >> 28) & 0x7F) | 0x80));
            put<uint8>((uint8)((v >> 35) & 0x7F));
        }
        else if (v < ((uint64)128) * 128 * 128 * 128 * 128 * 128 * 128) {
            put<uint8>((uint8)(((v >> 00) & 0x7F) | 0x80));
            put<uint8>((uint8)(((v >> 07) & 0x7F) | 0x80));
            put<uint8>((uint8)(((v >> 14) & 0x7F) | 0x80));
            put<uint8>((uint8)(((v >> 21) & 0x7F) | 0x80));
            put<uint8>((uint8)(((v >> 28) & 0x7F) | 0x80));
            put<uint8>((uint8)(((v >> 35) & 0x7F) | 0x80));
            put<uint8>((uint8)((v >> 42) & 0x7F));
        }
        else if (v < ((uint64)128) * 128 * 128 * 128 * 128 * 128 * 128 * 128) {
            put<uint8>((uint8)(((v >> 00) & 0x7F) | 0x80));
            put<uint8>((uint8)(((v >> 07) & 0x7F) | 0x80));
            put<uint8>((uint8)(((v >> 14) & 0x7F) | 0x80));
            put<uint8>((uint8)(((v >> 21) & 0x7F) | 0x80));
            put<uint8>((uint8)(((v >> 28) & 0x7F) | 0x80));
            put<uint8>((uint8)(((v >> 35) & 0x7F) | 0x80));
            put<uint8>((uint8)(((v >> 42) & 0x7F) | 0x80));
            put<uint8>((uint8)((v >> 49) & 0x7F));
        }
        else if (v < ((uint64)128) * 128 * 128 * 128 * 128 * 128 * 128 * 128 * 128) {
            put<uint8>((uint8)(((v >> 00) & 0x7F) | 0x80));
            put<uint8>((uint8)(((v >> 07) & 0x7F) | 0x80));
            put<uint8>((uint8)(((v >> 14) & 0x7F) | 0x80));
            put<uint8>((uint8)(((v >> 21) & 0x7F) | 0x80));
            put<uint8>((uint8)(((v >> 28) & 0x7F) | 0x80));
            put<uint8>((uint8)(((v >> 35) & 0x7F) | 0x80));
            put<uint8>((uint8)(((v >> 42) & 0x7F) | 0x80));
            put<uint8>((uint8)(((v >> 49) & 0x7F) | 0x80));
            put<uint8>((uint8)((v >> 56) & 0x7F));
        }
        else {
            put<uint8>((uint8)(((v >> 00) & 0x7F) | 0x80));
            put<uint8>((uint8)(((v >> 07) & 0x7F) | 0x80));
            put<uint8>((uint8)(((v >> 14) & 0x7F) | 0x80));
            put<uint8>((uint8)(((v >> 21) & 0x7F) | 0x80));
            put<uint8>((uint8)(((v >> 28) & 0x7F) | 0x80));
            put<uint8>((uint8)(((v >> 35) & 0x7F) | 0x80));
            put<uint8>((uint8)(((v >> 42) & 0x7F) | 0x80));
            put<uint8>((uint8)(((v >> 49) & 0x7F) | 0x80));
            put<uint8>((uint8)(((v >> 56) & 0x7F) | 0x80));
            put<uint8>((uint8)((v >> 63) & 0x7F));
        }

        return *this;
    }

    //! Put a string by saving it's length followed by the data itself.
    inline binary_buffer & put_string(const char* data, size_t len)
    {
        return put_varint((uint32)len).append(data, len);
    }

    //! Put a string by saving it's length followed by the data itself.
    inline binary_buffer & put_string(const std::string& str)
    {
        return put_string(str.data(), str.size());
    }

    //! Put a binary_buffer by saving it's length followed by the data itself.
    inline binary_buffer & put_string(const binary_buffer& bb)
    {
        return put_string(bb.data(), bb.size());
    }
};

/*!
 * binary_buffer_ref represents a memory area as pointer and valid length. It
 * is not deallocated or otherwise managed. This class can be used to pass
 * around references to binary_buffer objects.
 */
class binary_buffer_ref
{
protected:
    //! Allocated buffer pointer.
    const char* m_data;

    //! Size of valid data.
    size_t m_size;

public:
    //! Constructor, assign memory area from binary_buffer.
    binary_buffer_ref(const binary_buffer& bb)
        : m_data(bb.data()), m_size(bb.size())
    { }

    //! Constructor, assign memory area from pointer and length.
    binary_buffer_ref(const void* data, size_t n)
        : m_data(reinterpret_cast<const char*>(data)), m_size(n)
    { }

    //! Constructor, assign memory area from string, does NOT copy.
    inline binary_buffer_ref(const std::string& str)
        : m_data(str.data()), m_size(str.size())
    { }

    //! Return a pointer to the currently kept memory area.
    const void * data() const
    { return m_data; }

    //! Return the currently valid length in bytes.
    size_t size() const
    { return m_size; }

    //! Explicit conversion to std::string (copies memory of course).
    inline std::string str() const
    { return std::string(reinterpret_cast<const char*>(m_data), m_size); }

    //! Compare contents of two binary_buffer_refs.
    bool operator == (const binary_buffer_ref& br) const
    {
        if (m_size != br.m_size) return false;
        return memcmp(m_data, br.m_data, m_size) == 0;
    }

    //! Compare contents of two binary_buffer_refs.
    bool operator != (const binary_buffer_ref& br) const
    {
        if (m_size != br.m_size) return true;
        return memcmp(m_data, br.m_data, m_size) != 0;
    }
};

/*!
 * binary_reader represents a binary_buffer_ref with an additional cursor with which
 * the memory can be read incrementally.
 */
class binary_reader : public binary_buffer_ref
{
protected:
    //! Current read cursor
    size_t m_curr;

public:
    //! Constructor, assign memory area from binary_buffer.
    inline binary_reader(const binary_buffer_ref& br)
        : binary_buffer_ref(br), m_curr(0)
    { }

    //! Constructor, assign memory area from pointer and length.
    inline binary_reader(const void* data, size_t n)
        : binary_buffer_ref(data, n), m_curr(0)
    { }

    //! Constructor, assign memory area from string, does NOT copy.
    inline binary_reader(const std::string& str)
        : binary_buffer_ref(str), m_curr(0)
    { }

    //! Return the current read cursor.
    inline size_t curr() const
    {
        return m_curr;
    }

    //! Reset the read cursor.
    inline binary_reader & rewind()
    {
        m_curr = 0;
        return *this;
    }

    //! Check that n bytes are available at the cursor.
    inline bool cursor_available(size_t n) const
    {
        return (m_curr + n <= m_size);
    }

    //! Throws a std::underflow_error unless n bytes are available at the
    //! cursor.
    inline void check_available(size_t n) const
    {
        if (!cursor_available(n))
            throw (std::underflow_error("binary_reader underrun"));
    }

    //! Return true if the cursor is at the end of the buffer.
    inline bool empty() const
    {
        return (m_curr == m_size);
    }

    //! Advance the cursor given number of bytes without reading them.
    inline binary_reader & skip(size_t n)
    {
        check_available(n);
        m_curr += n;

        return *this;
    }

    //! Fetch a number of unstructured bytes from the buffer, advancing the
    //! cursor.
    inline binary_reader & read(void* outdata, size_t datalen)
    {
        check_available(datalen);
        memcpy(outdata, m_data + m_curr, datalen);
        m_curr += datalen;

        return *this;
    }

    //! Fetch a number of unstructured bytes from the buffer as std::string,
    //! advancing the cursor.
    inline std::string read(size_t datalen)
    {
        check_available(datalen);
        std::string out(m_data + m_curr, datalen);
        m_curr += datalen;
        return out;
    }

    //! Fetch a single item of the template type Type from the buffer,
    //! advancing the cursor. Be careful with implicit type conversions!
    template <typename Type>
    inline Type get()
    {
        check_available(sizeof(Type));

        Type ret = *reinterpret_cast<const Type*>(m_data + m_curr);
        m_curr += sizeof(Type);

        return ret;
    }

    //! Fetch a varint with up to 32-bit from the buffer at the cursor.
    inline uint32 get_varint()
    {
        uint32 u, v = get<uint8>();
        if (!(v & 0x80)) return v;
        v &= 0x7F;
        u = get<uint8>(), v |= (u & 0x7F) << 7;
        if (!(u & 0x80)) return v;
        u = get<uint8>(), v |= (u & 0x7F) << 14;
        if (!(u & 0x80)) return v;
        u = get<uint8>(), v |= (u & 0x7F) << 21;
        if (!(u & 0x80)) return v;
        u = get<uint8>();
        if (u & 0xF0)
            throw (std::overflow_error("Overflow during varint decoding."));
        v |= (u & 0x7F) << 28;
        return v;
    }

    //! Fetch a 64-bit varint from the buffer at the cursor.
    inline uint64 get_varint64()
    {
        uint64 u, v = get<uint8>();
        if (!(v & 0x80)) return v;
        v &= 0x7F;
        u = get<uint8>(), v |= (u & 0x7F) << 7;
        if (!(u & 0x80)) return v;
        u = get<uint8>(), v |= (u & 0x7F) << 14;
        if (!(u & 0x80)) return v;
        u = get<uint8>(), v |= (u & 0x7F) << 21;
        if (!(u & 0x80)) return v;
        u = get<uint8>(), v |= (u & 0x7F) << 28;
        if (!(u & 0x80)) return v;
        u = get<uint8>(), v |= (u & 0x7F) << 35;
        if (!(u & 0x80)) return v;
        u = get<uint8>(), v |= (u & 0x7F) << 42;
        if (!(u & 0x80)) return v;
        u = get<uint8>(), v |= (u & 0x7F) << 49;
        if (!(u & 0x80)) return v;
        u = get<uint8>(), v |= (u & 0x7F) << 56;
        if (!(u & 0x80)) return v;
        u = get<uint8>();
        if (u & 0xFE)
            throw (std::overflow_error("Overflow during varint64 decoding."));
        v |= (u & 0x7F) << 63;
        return v;
    }

    //! Fetch a string which was put via put_string().
    inline std::string get_string()
    {
        uint32 len = get_varint();
        return read(len);
    }

    //! Fetch a binary_buffer_ref to a binary string or blob which was put via
    //! put_string(). Does NOT copy the data.
    inline binary_buffer_ref get_binary_buffer_ref()
    {
        uint32 len = get_varint();
        // save object
        binary_buffer_ref br(m_data + m_curr, len);
        // skip over sub block data
        skip(len);
        return br;
    }
};

//! \}

STXXL_END_NAMESPACE

#endif // !STXXL_COMMON_BINARY_BUFFER_HEADER
