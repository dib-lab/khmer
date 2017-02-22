#include "khmer/_cpy_utils.hh"
#include "oxli/oxli.hh"

using namespace oxli;
using namespace oxli::read_parsers;

namespace khmer {

// Convert a hash to a python long object.
bool convert_HashIntoType_to_PyObject(const HashIntoType &hashval,
        PyObject **value)
{
    *value = PyLong_FromUnsignedLongLong(hashval);
    return true;
}


// Convert a python long to a hash
bool convert_PyLong_to_HashIntoType(PyObject * value, HashIntoType &hashval)
{
    if (PyLong_Check(value)) {
        //(PyLongObject *)
        hashval = PyLong_AsUnsignedLongLong(value);
        return true;
    } else if (PyInt_Check(value)) {
        hashval = PyInt_AsLong(value);
        return true;
    } else {
        PyErr_SetString(PyExc_ValueError, "could not convert to hash");
        return false;
    }
}


// Take a Python object and (try to) convert it to a HashIntoType.
// Note: will set error condition and return false if cannot do.

bool convert_PyObject_to_HashIntoType(PyObject * value,
        HashIntoType& hashval,
        WordLength ksize)
{
    if (PyInt_Check(value) || PyLong_Check(value)) {
        return convert_PyLong_to_HashIntoType(value, hashval);
    } else {
        PyErr_SetString(PyExc_ValueError,
                        "must use a hash");
        return false;
    }
}

// Take a Python object and (try to) convert it to a HashIntoType.
// Note: will set error condition and return false if cannot do.
// Further note: the main difference between this and
// ht_convert_PyObject_to_Kmer is that this will not pass HashIntoType
// numbers through the Kmer class, which means reverse complements
// will not be calculated.  There is a test in test_nodegraph.py
// that checks this.

bool ht_convert_PyObject_to_HashIntoType(PyObject * value,
        HashIntoType& hashval,
        const Hashtable * ht)
{
    if (PyInt_Check(value) || PyLong_Check(value)) {
        return convert_PyLong_to_HashIntoType(value, hashval);
    } else if (PyUnicode_Check(value))  {
        PyObject* val_as_str = PyUnicode_AsEncodedString(value,
                               "utf-8", "strict");
        std::string s = PyBytes_AsString(val_as_str);
        if (strlen(s.c_str()) != ht->ksize()) {
            Py_DECREF(val_as_str);
            PyErr_SetString(PyExc_ValueError,
                            "k-mer length must equal the k-mer size");
            return false;
        }

        try {
            hashval = ht->hash_dna(s.c_str());
        } catch (oxli_exception &e) {
            PyErr_SetString(PyExc_ValueError, e.what());
            Py_DECREF(val_as_str);
            return false;
        }

        Py_DECREF(val_as_str);
        return true;

    } else if (PyBytes_Check(value)) {
        std::string s = PyBytes_AsString(value);
        if (strlen(s.c_str()) != ht->ksize()) {
            PyErr_SetString(PyExc_ValueError,
                            "k-mer length must equal the k-mer size");
            return false;
        }
        try {
            hashval = ht->hash_dna(s.c_str());
        } catch (oxli_exception &e) {
            PyErr_SetString(PyExc_ValueError, e.what());
            return false;
        }
        return true;
    } else {
        PyErr_SetString(PyExc_ValueError,
                        "k-mers must be either a hash or a string");
        return false;
    }
}

// Take a Python object and (try to) convert it to a oxli::Kmer.
// Note: will set error condition and return false if cannot do.

bool ht_convert_PyObject_to_Kmer(PyObject * value,
                                 Kmer& kmer, const Hashtable * ht)
{
    if (PyInt_Check(value) || PyLong_Check(value)) {
        HashIntoType h;
        if (!convert_PyLong_to_HashIntoType(value, h)) {
            return false;
        }
        kmer.set_from_unique_hash(h, ht->ksize());
        return true;
    } else if (PyUnicode_Check(value))  {
        std::string s = PyBytes_AsString(PyUnicode_AsEncodedString(
                                             value, "utf-8", "strict"));
        if (strlen(s.c_str()) != ht->ksize()) {
            PyErr_SetString(PyExc_ValueError,
                            "k-mer length must equal the k-mer size");
            return false;
        }
        kmer = Kmer(s, ht->ksize());
        return true;

    } else if (PyBytes_Check(value)) {
        std::string s = PyBytes_AsString(value);
        if (strlen(s.c_str()) != ht->ksize()) {
            PyErr_SetString(PyExc_ValueError,
                            "k-mer length must equal the k-mer size");
            return false;
        }
        kmer = Kmer(s, ht->ksize());
        return true;
    } else {
        PyErr_SetString(PyExc_ValueError,
                        "k-mers must be either a hash or a string");
        return false;
    }
}


bool convert_Pytablesizes_to_vector(PyListObject * sizes_list_o,
                                           std::vector<uint64_t>& sizes)
{
    Py_ssize_t sizes_list_o_length = PyList_GET_SIZE(sizes_list_o);
    if (sizes_list_o_length < 1) {
        PyErr_SetString(PyExc_ValueError,
                        "tablesizes needs to be one or more numbers");
        return false;
    }
    for (Py_ssize_t i = 0; i < sizes_list_o_length; i++) {
        PyObject * size_o = PyList_GET_ITEM(sizes_list_o, i);
        if (PyLong_Check(size_o)) {
            sizes.push_back(PyLong_AsUnsignedLongLong(size_o));
        } else if (PyInt_Check(size_o)) {
            sizes.push_back(PyInt_AsLong(size_o));
        } else if (PyFloat_Check(size_o)) {
            sizes.push_back(PyFloat_AS_DOUBLE(size_o));
        } else {
            PyErr_SetString(PyExc_TypeError,
                            "2nd argument must be a list of ints, longs, or floats");
            return false;
        }
    }
    return true;
}

FastxParserPtr& _PyObject_to_khmer_ReadParser(PyObject * py_object);

}
