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

//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

#ifndef SEQAN_HEADER_FILE_DIRECTORY_H
#define SEQAN_HEADER_FILE_DIRECTORY_H

#ifdef PLATFORM_WINDOWS
# include <io.h>
#else
# include <dirent.h>
#endif


/* IOREV
 *
 * _nodoc_
 * _tested_
 * _windows_
 *
 *
 * used by file/file_format_mmap.h and some apps
 * no documentation at all
 *
 */



//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{

#ifdef PLATFORM_WINDOWS

    class Directory
    {
//IOREV _nodoc_ _windows_
    protected:

        intptr_t            handle;
        struct _finddata_t    entry;
        bool                _atEnd;

        friend inline bool open(Directory &dir, char const *dirName);
        friend inline bool close(Directory &dir);
        friend inline char const * value(Directory &dir);
        friend inline char const * value(Directory const &dir);
        friend inline Directory & goNext(Directory &dir);
        friend inline bool atEnd(Directory &dir);
        friend inline bool atEnd(Directory const &dir);

    public:

        Directory()
        {
            handle = 0;
            _atEnd = true;
        }

        Directory(char const *dirName)
        {
            open(*this, dirName);
        }

        ~Directory()
        {
            close(*this);
        }

        inline char const * operator* () const
        {
            return value(*this);
        }

        inline Directory & operator++ ()
        {
            return goNext(*this);
        }

        inline operator bool () const
        {
            return !_atEnd;
        }
    };

//////////////////////////////////////////////////////////////////////////////

    inline bool
    open(Directory &dir, char const *dirName)
    {
//IOREV _nodoc_
        CharString selection = dirName;
        append(selection, "\\*");
        dir._atEnd = ((dir.handle = _findfirst(toCString(selection), &dir.entry)) == -1L);
        return !dir._atEnd;
    }

    inline bool
    close(Directory &dir)
    {
//IOREV _nodoc_
        int result = 0;
        if (dir.handle)
            result = _findclose(dir.handle);

        dir._atEnd = true;
        dir.handle = 0;
        return result == 0;
    }

    inline char const *
    value(Directory &dir)
    {
//IOREV _nodoc_
        return dir.entry.name;
    }

    inline char const *
    value(Directory const &dir)
    {
//IOREV _nodoc_
        return dir.entry.name;
    }

    inline Directory &
    goNext(Directory &dir)
    {
//IOREV _nodoc_
        dir._atEnd = (_findnext(dir.handle, &dir.entry) != 0);
        return dir;
    }

    inline bool
    atEnd(Directory &dir)
    {
//IOREV _nodoc_
        return dir._atEnd;
    }

    inline bool
    atEnd(Directory const &dir)
    {
//IOREV _nodoc_
        return dir._atEnd;
    }

//////////////////////////////////////////////////////////////////////////////


#else


//////////////////////////////////////////////////////////////////////////////

    class Directory
    {
//IOREV _nodoc_
    protected:

        DIR        *handle;
        dirent    *it;

        friend inline bool open(Directory &dir, char const *dirName);
        friend inline bool close(Directory &dir);
        friend inline char const * value(Directory &dir);
        friend inline char const * value(Directory const &dir);
        friend inline Directory & goBegin(Directory &dir);
        friend inline Directory & goNext(Directory &dir);
        friend inline bool atEnd(Directory &dir);
        friend inline bool atEnd(Directory const &dir);

    public:

        Directory()
        {
            handle = NULL;
            it = NULL;
        }

        Directory(char const *dirName)
        {
            open(*this, dirName);
        }

        ~Directory()
        {
            close(*this);
        }

        inline char const * operator* () const
        {
            return value(*this);
        }

        inline Directory & operator++ ()
        {
            return goNext(*this);
        }

        inline operator bool () const
        {
            return !atEnd(*this);
        }
    };

//////////////////////////////////////////////////////////////////////////////

    inline bool
    open(Directory &dir, char const *dirName)
    {
//IOREV _nodoc_
        if ((dir.handle = opendir(dirName)) != NULL)
        {
            goNext(dir);
            return true;
        }
        dir.it = NULL;
        return false;
    }

    inline bool
    close(Directory &dir)
    {
//IOREV _nodoc_
        int result = 0;
        if (dir.handle != NULL)
            result = closedir(dir.handle);

        dir.handle = NULL;
        dir.it = NULL;
        return result == 0;
    }

    inline char const *
    value(Directory &dir)
    {
//IOREV _nodoc_
        return dir.it->d_name;
    }

    inline char const *
    value(Directory const &dir)
    {
//IOREV _nodoc_
        return dir.it->d_name;
    }

    inline Directory &
    goBegin(Directory &dir)
    {
//IOREV _nodoc_
        rewinddir(dir.handle);
        return goNext(dir);
    }

    inline Directory &
    goNext(Directory &dir)
    {
//IOREV _nodoc_
        dir.it = readdir(dir.handle);
        return dir;
    }

    inline bool
    atEnd(Directory &dir)
    {
//IOREV _nodoc_
        return dir.it == NULL;
    }

    inline bool
    atEnd(Directory const &dir)
    {
//IOREV _nodoc_
        return dir.it == NULL;
    }

#endif

}

#endif
