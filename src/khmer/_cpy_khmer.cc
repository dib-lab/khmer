/*
This file is part of khmer, https://github.com/dib-lab/khmer/, and is
Copyright (C) 2010-2015, Michigan State University.
Copyright (C) 2015-2016, The Regents of the University of California.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.

    * Neither the name of the Michigan State University nor the names
      of its contributors may be used to endorse or promote products
      derived from this software without specific prior written
      permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
LICENSE (END)

Contact: khmer-project@idyll.org
*/

//
// A module for Python that exports khmer C++ library functions.
//

// Must be first.
#include <Python.h>
#include <string>

#include "khmer/_cpy_khmer.hh"
#include "khmer/_cpy_utils.hh"

using namespace oxli;


//
// Function necessary for Python loading:
//

extern "C" {
    MOD_INIT(_khmer);
}

namespace khmer {

PyObject * forward_hash(PyObject * self, PyObject * args)
{
    const char * kmer;
    WordLength ksize;

    if (!PyArg_ParseTuple(args, "sb", &kmer, &ksize)) {
        return NULL;
    }

    if (ksize > KSIZE_MAX) {
        PyErr_Format(PyExc_ValueError, "k-mer size must be <= %u", KSIZE_MAX);
        return NULL;
    }

    if (strlen(kmer) != ksize) {
        PyErr_Format(PyExc_ValueError, "k-mer size different from ksize");
        return NULL;
    }

    try {
        PyObject * hash = nullptr;
        const HashIntoType h(_hash(kmer, ksize));
        convert_HashIntoType_to_PyObject(h, &hash);
        return hash;
    } catch (oxli_exception &e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        return NULL;
    }
}

PyObject * forward_hash_no_rc(PyObject * self, PyObject * args)
{
    const char * kmer;
    WordLength ksize;

    if (!PyArg_ParseTuple(args, "sb", &kmer, &ksize)) {
        return NULL;
    }

    if (ksize > KSIZE_MAX) {
        PyErr_Format(PyExc_ValueError, "k-mer size must be <= %u", KSIZE_MAX);
        return NULL;
    }

    if (strlen(kmer) != ksize) {
        PyErr_SetString(PyExc_ValueError,
                        "k-mer length must equal the k-size");
        return NULL;
    }

    PyObject * hash = nullptr;
    const HashIntoType h(_hash_forward(kmer, ksize));
    convert_HashIntoType_to_PyObject(h, &hash);
    return hash;
}

PyObject * reverse_hash(PyObject * self, PyObject * args)
{
    PyObject * val;
    HashIntoType hash;
    WordLength ksize;

    if (!PyArg_ParseTuple(args, "Ob", &val, &ksize)) {
        return NULL;
    }

    if (PyLong_Check(val) || PyInt_Check(val)) {
        if (!convert_PyLong_to_HashIntoType(val, hash)) {
            return NULL;
        }
    } else {
        PyErr_SetString(PyExc_TypeError,
                        "Hash value must be an integer.");
        return NULL;
    }

    if (ksize > KSIZE_MAX) {
        PyErr_Format(PyExc_ValueError, "k-mer size must be <= %u", KSIZE_MAX);
        return NULL;
    }

    return PyUnicode_FromString(_revhash(hash, ksize).c_str());
}

PyObject * murmur3_forward_hash(PyObject * self, PyObject * args)
{
    const char * kmer;

    if (!PyArg_ParseTuple(args, "s", &kmer)) {
        return NULL;
    }

    PyObject * hash = nullptr;
    const HashIntoType h(_hash_murmur(kmer, strlen(kmer)));
    convert_HashIntoType_to_PyObject(h, &hash);
    return hash;
}

PyObject * murmur3_forward_hash_no_rc(PyObject * self, PyObject * args)
{
    const char * kmer;

    if (!PyArg_ParseTuple(args, "s", &kmer)) {
        return NULL;
    }

    PyObject * hash = nullptr;
    const HashIntoType h(_hash_murmur_forward(kmer, strlen(kmer)));
    convert_HashIntoType_to_PyObject(h, &hash);
    return hash;
}

PyObject * reverse_complement(PyObject * self, PyObject * args)
{
    const char * sequence;
    if (!PyArg_ParseTuple(args, "s", &sequence)) {
        return NULL;
    }

    std::string s(sequence);
    try {
        s = _revcomp(s);
    } catch (oxli_exception &e) {
        PyErr_SetString(PyExc_RuntimeError, e.what());
        return NULL;
    }
    return PyUnicode_FromString(s.c_str());
}

//
// technique for resolving literal below found here:
// https://gcc.gnu.org/onlinedocs/gcc-4.9.1/cpp/Stringification.html
//

PyObject *
get_version_cpp( PyObject * self, PyObject * args )
{
#define xstr(s) str(s)
#define str(s) #s
    std::string dVersion = xstr(VERSION);
    return PyUnicode_FromString(dVersion.c_str());
}

PyMethodDef KhmerMethods[] = {
    {
        "forward_hash",     forward_hash,
        METH_VARARGS,       "",
    },
    {
        "forward_hash_no_rc",   forward_hash_no_rc,
        METH_VARARGS,       "",
    },
    {
        "reverse_hash",     reverse_hash,
        METH_VARARGS,       "",
    },
    {
        "hash_murmur3",
        murmur3_forward_hash,
        METH_VARARGS,
        "Calculate the hash value of a k-mer using MurmurHash3 "
        "(with reverse complement)",
    },
    {
        "hash_no_rc_murmur3",
        murmur3_forward_hash_no_rc,
        METH_VARARGS,
        "Calculate the hash value of a k-mer using MurmurHash3 "
        "(no reverse complement)",
    },
    {
        "reverse_complement",
        reverse_complement,
        METH_VARARGS,
        "Calculate the reverse-complement of the DNA sequence "
        "with alphabet ACGT",
    },
    {
        "get_version_cpp", get_version_cpp,
        METH_VARARGS, "return the VERSION c++ compiler option"
    },
    { NULL, NULL, 0, NULL } // sentinel
};

} // namespace khmer

//
// Module machinery.
//

MOD_INIT(_khmer)
{
    using namespace khmer;
    using namespace oxli;
    using namespace oxli::read_parsers;



    _init_ReadParser_Type_constants();
    if (PyType_Ready( &khmer_ReadParser_Type ) < 0) {
        return MOD_ERROR_VAL;
    }

    if (PyType_Ready(&khmer_Read_Type ) < 0) {
        return MOD_ERROR_VAL;
    }

    if (PyType_Ready(&khmer_ReadPairIterator_Type ) < 0) {
        return MOD_ERROR_VAL;
    }

    PyObject * m;

    MOD_DEF(m, "_khmer", "interface for the khmer module low-level extensions",
            KhmerMethods);

    if (m == NULL) {
        return MOD_ERROR_VAL;
    }

    PyObject * filetype_dict = Py_BuildValue("{s,i,s,i,s,i,s,i,s,i,s,i,s,i}",
                               "COUNTING_HT", SAVED_COUNTING_HT,
                               "HASHBITS", SAVED_HASHBITS,
                               "TAGS", SAVED_TAGS,
                               "STOPTAGS", SAVED_STOPTAGS,
                               "SUBSET", SAVED_SUBSET,
                               "LABELSET", SAVED_LABELSET,
                               "SMALLCOUNT", SAVED_SMALLCOUNT);
    if (PyModule_AddObject( m, "FILETYPES", filetype_dict ) < 0) {
        return MOD_ERROR_VAL;
    }

    Py_INCREF(&khmer_Read_Type);
    if (PyModule_AddObject( m, "Read",
                            (PyObject *)&khmer_Read_Type ) < 0) {
        return MOD_ERROR_VAL;
    }

    Py_INCREF(&khmer_ReadParser_Type);
    if (PyModule_AddObject( m, "ReadParser",
                            (PyObject *)&khmer_ReadParser_Type ) < 0) {
        return MOD_ERROR_VAL;
    }


    return MOD_SUCCESS_VAL(m);
}

// vim: set ft=cpp sts=4 sw=4 tw=79:
