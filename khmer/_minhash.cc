#include "_minhash.hh"

//
// Function necessary for Python loading:
//

extern "C" {
    MOD_INIT(_minhash);
}


#include <string>
#include <map>
#include <exception>
#include <iostream>
#include <fstream>
#include <sstream> // IWYU pragma: keep

using namespace khmer;

PyTypeObject MinHash_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)        /* init & ob_size */
    "_minhash.MinHash",                   /* tp_name */
    sizeof(MinHash_Object),               /* tp_basicsize */
    0,                                    /* tp_itemsize */
    0,                                    /* tp_dealloc */
    0,                                    /* tp_print */
    0,                                    /* tp_getattr */
    0,                                    /* tp_setattr */
    0,                                    /* tp_compare */
    0,                                    /* tp_repr */
    0,                                    /* tp_as_number */
    0,                                    /* tp_as_sequence */
    0,                                    /* tp_as_mapping */
    0,                                    /* tp_hash */
    0,                                    /* tp_call */
    0,                                    /* tp_str */
    0,                                    /* tp_getattro */
    0,                                    /* tp_setattro */
    0,                                    /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,                   /* tp_flags */
    "A MinHash sketch.",                  /* tp_doc */
};

PyObject * build_MinHash_Object(khmer::KmerMinHash * mh)
{
    MinHash_Object * obj = (MinHash_Object *) \
                           PyObject_New(MinHash_Object, &MinHash_Type);
    obj->mh = mh;

    return (PyObject *) obj;
}

PyTypeObject NeighborhoodMinHash_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)            /* init & ob_size */
    "_minhash.NeighborhoodMinHash",           /* tp_name */
    sizeof(NeighborhoodMinHash_Object),       /* tp_basicsize */
    0,                                        /* tp_itemsize */
    0,                                        /* tp_dealloc */
    0,                                        /* tp_print */
    0,                                        /* tp_getattr */
    0,                                        /* tp_setattr */
    0,                                        /* tp_compare */
    0,                                        /* tp_repr */
    0,                                        /* tp_as_number */
    0,                                        /* tp_as_sequence */
    0,                                        /* tp_as_mapping */
    0,                                        /* tp_hash */
    0,                                        /* tp_call */
    0,                                        /* tp_str */
    0,                                        /* tp_getattro */
    0,                                        /* tp_setattro */
    0,                                        /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,                       /* tp_flags */
    "A MinHash sketch on tag neighborhoods.", /* tp_doc */
};

PyObject * build_NeighborhoodMinHash_Object(khmer::NeighborhoodMinHash * nbhd_mh)
{
    NeighborhoodMinHash_Object * obj = (NeighborhoodMinHash_Object *) \
                                       PyObject_New(NeighborhoodMinHash_Object, &NeighborhoodMinHash_Type);
    obj->nbhd_mh = nbhd_mh;

    return (PyObject *) obj;
}

PyTypeObject CombinedMinHash_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)            /* init & ob_size */
    "_minhash.CombinedMinHash",               /* tp_name */
    sizeof(CombinedMinHash_Object),           /* tp_basicsize */
    0,                                        /* tp_itemsize */
    0,                                        /* tp_dealloc */
    0,                                        /* tp_print */
    0,                                        /* tp_getattr */
    0,                                        /* tp_setattr */
    0,                                        /* tp_compare */
    0,                                        /* tp_repr */
    0,                                        /* tp_as_number */
    0,                                        /* tp_as_sequence */
    0,                                        /* tp_as_mapping */
    0,                                        /* tp_hash */
    0,                                        /* tp_call */
    0,                                        /* tp_str */
    0,                                        /* tp_getattro */
    0,                                        /* tp_setattro */
    0,                                        /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,                       /* tp_flags */
    "A MinHash sketch on tag neighborhoods.", /* tp_doc */
};

PyObject * build_CombinedMinHash_Object(khmer::CombinedMinHash * combined_mh)
{
    CombinedMinHash_Object * obj = (CombinedMinHash_Object *) \
                                   PyObject_New(CombinedMinHash_Object, &CombinedMinHash_Type);
    obj->combined_mh = combined_mh;

    return (PyObject *) obj;
}

bool check_IsMinHash(PyObject * mh);
bool check_IsCombinedMinHash(PyObject * combined_mh);

////

static
void
MinHash_dealloc(MinHash_Object * obj)
{
    delete obj->mh;
    obj->mh = NULL;
    Py_TYPE(obj)->tp_free((PyObject*)obj);
}

static
PyObject *
minhash_add_sequence(MinHash_Object * me, PyObject * args)
{
    const char * sequence = NULL;
    if (!PyArg_ParseTuple(args, "s", &sequence)) {
        return NULL;
    }
    KmerMinHash * mh = me->mh;

    mh->add_sequence(sequence);

    Py_INCREF(Py_None);
    return Py_None;
}

static
PyObject *
minhash_add_protein(MinHash_Object * me, PyObject * args)
{
    const char * sequence = NULL;
    if (!PyArg_ParseTuple(args, "s", &sequence)) {
        return NULL;
    }
    KmerMinHash * mh = me->mh;

    unsigned int ksize = mh->ksize / 3;

    if(strlen(sequence) < ksize) {
        Py_INCREF(Py_None);
        return Py_None;
    }

    if (!mh->is_protein) {
        assert(0);
    } else {                      // protein
        std::string seq = sequence;
        for (unsigned int i = 0; i < seq.length() - ksize + 1; i ++) {
            std::string aa = seq.substr(i, ksize);

            mh->add_kmer(aa);
        }
    }

    Py_INCREF(Py_None);
    return Py_None;
}

static
PyObject *
minhash_add_hash(MinHash_Object * me, PyObject * args)
{
    long int hh;
    if (!PyArg_ParseTuple(args, "l", &hh)) {
        return NULL;
    }

    me->mh->add_hash(hh);

    Py_INCREF(Py_None);
    return Py_None;
}

static
PyObject *
minhash_get_mins(MinHash_Object * me, PyObject * args)
{
    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    KmerMinHash * mh = me->mh;
    PyObject * mins_o = PyList_New(mh->mins.size());

    unsigned int j = 0;
    for (CMinHashType::iterator i = mh->mins.begin(); i != mh->mins.end(); ++i) {
        PyList_SET_ITEM(mins_o, j, PyLong_FromUnsignedLongLong(*i));
        j++;
    }
    return(mins_o);
}

static PyObject * minhash___copy__(MinHash_Object * me, PyObject * args)
{
    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    KmerMinHash * mh = me->mh;
    KmerMinHash * new_mh = new KmerMinHash(mh->num, mh->ksize, mh->prime,
                                           mh->is_protein);
    new_mh->merge(*mh);

    return build_MinHash_Object(new_mh);
}

static PyObject * minhash___len__(MinHash_Object * me, PyObject * args)
{
    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    KmerMinHash * mh = me->mh;

    return PyInt_FromLong(mh->mins.size());
}

static PyObject * minhash_count_common(MinHash_Object * me, PyObject * args)
{
    PyObject * other_mh;

    if (!PyArg_ParseTuple(args, "O", &other_mh)) {
        return NULL;
    }

    if (!check_IsMinHash(other_mh)) {
        return NULL;
    }

    unsigned int n = me->mh->count_common(*((MinHash_Object*)other_mh)->mh);
    return PyInt_FromLong(n);
}

static PyObject * minhash_compare(MinHash_Object * me, PyObject * args)
{
    PyObject * other_mh;

    if (!PyArg_ParseTuple(args, "O", &other_mh)) {
        return NULL;
    }

    if (!check_IsMinHash(other_mh)) {
        return NULL;
    }

    unsigned int n = me->mh->count_common(*((MinHash_Object*)other_mh)->mh);
    unsigned int size = me->mh->mins.size();

    return PyFloat_FromDouble(float(n) / float(size));
}

static PyMethodDef MinHash_methods [] = {
    {
        "add_sequence",
        (PyCFunction)minhash_add_sequence, METH_VARARGS,
        "Add kmer into MinHash"
    },
    {
        "add_protein",
        (PyCFunction)minhash_add_protein, METH_VARARGS,
        "Add AA kmer into protein MinHash"
    },
    {
        "add_hash",
        (PyCFunction)minhash_add_hash, METH_VARARGS,
        "Add kmer into MinHash"
    },
    {
        "get_mins",
        (PyCFunction)minhash_get_mins, METH_VARARGS,
        "Get MinHash signature"
    },
    {
        "__copy__",
        (PyCFunction)minhash___copy__, METH_VARARGS,
        "Copy this MinHash object",
    },
    {
        "__len__",                  // @CTB this doesn't work...?
        (PyCFunction)minhash___len__, METH_VARARGS,
        "Number of hashes in this MinHash object",
    },
    {
        "count_common",
        (PyCFunction)minhash_count_common, METH_VARARGS,
        "Get number of hashes in common with other."
    },
    {
        "compare",
        (PyCFunction)minhash_compare, METH_VARARGS,
        "Get the Jaccard similarity between this and other."
    },
    { NULL, NULL, 0, NULL } // sentinel
};

static
PyObject *
MinHash_new(PyTypeObject * subtype, PyObject * args, PyObject * kwds)
{
    PyObject * self     = subtype->tp_alloc( subtype, 1 );
    if (self == NULL) {
        return NULL;
    }

    unsigned int _n, _ksize;
    long int _p = DEFAULT_MINHASH_PRIME;
    PyObject * is_protein_o = NULL;
    if (!PyArg_ParseTuple(args, "II|lO", &_n, &_ksize, &_p, &is_protein_o)) {
        return NULL;
    }

    MinHash_Object * myself = (MinHash_Object *)self;
    bool is_protein = false;
    if (is_protein_o && PyObject_IsTrue(is_protein_o)) {
        is_protein = true;
    }

    myself->mh = new KmerMinHash(_n, _ksize, _p, is_protein);

    return self;
}

bool check_IsMinHash(PyObject * mh)
{
    if (!PyObject_TypeCheck(mh, &MinHash_Type)) {
        return false;
    }
    return true;
}


///////////////////////////////////////////////////////////////////////////////

static
void
NeighborhoodMinHash_dealloc(NeighborhoodMinHash_Object * obj)
{
    delete obj->nbhd_mh;
    obj->nbhd_mh = NULL;
    Py_TYPE(obj)->tp_free((PyObject*)obj);
}

// build_combined_minhashes: find quasi-linear paths in the neighborhoods,
//   turn them into CombinedMinHashes.

void build_combined_minhashes(NeighborhoodMinHash& nbhd_mh,
                              std::vector<CombinedMinHash *>& combined_mhs,
                              unsigned int combined_minhash_size=500,
                              unsigned int nbhd_size_limit=3)
{
    unsigned int k = nbhd_mh.tag_to_mh.begin()->second->ksize;
    unsigned int p = nbhd_mh.tag_to_mh.begin()->second->prime;
    bool prot = nbhd_mh.tag_to_mh.begin()->second->is_protein;

    CombinedMinHash * combined_mh = NULL;

    // iterate over all tags:
    TagToMinHash::const_iterator mhi;
    mhi = nbhd_mh.tag_to_mh.begin();

    while(mhi != nbhd_mh.tag_to_mh.end()) {
        // already merged this 'un? move to next.
        TagToTagSet::iterator posn;
        posn = nbhd_mh.tag_connections.find(mhi->first);
        if (posn == nbhd_mh.tag_connections.end() ||
                posn->second.size() > nbhd_size_limit) {
            mhi++;
            continue;
        }

        // build minhash to merge into:
        if (combined_mh == NULL) {
            combined_mh = new CombinedMinHash;
            combined_mh->mh = new KmerMinHash(combined_minhash_size, k, p, prot);
        }

        // keep track of tags that could be merged into this:
        TagSet to_be_merged = posn->second;

        // merge nbhd minhashes in, tag by tag, until no more to combine.
        bool did_combine = true;
        while (did_combine) {
            did_combine = false;

            // walk through the list of tags to be merged:
            TagSet::iterator ti;
            while(to_be_merged.size()) {
                // grab the first tag & its list of connected tags:
                TagToTagSet::iterator posn2;
                ti = to_be_merged.begin();
                posn2 = nbhd_mh.tag_connections.find(*ti);

                // already merged & removed from nbhd_mh.tag_connections? ok, ignore.
                if (posn2 == nbhd_mh.tag_connections.end()) {
                    to_be_merged.erase(ti);
                    continue;
                }

                // more than three tags in its connections? ignore.
                if (posn2->second.size() > 5) {
                    to_be_merged.erase(ti);
                    continue;
                }

                // finally! merge.
                KmerMinHash * mh = nbhd_mh.tag_to_mh[*ti];
                combined_mh->mh->merge(*mh);
                combined_mh->tags.insert(*ti);

                // track info.
                did_combine = true;

                // add the just-merged one's connected tags to to_be_merged
                to_be_merged.insert(posn2->second.begin(), posn2->second.end());

                // remove:
                nbhd_mh.tag_connections.erase(posn2);
                to_be_merged.erase(ti);
            }
        }

        // store combined_mhs in a list; recycle empty ones.
        if (combined_mh->tags.size()) {
            combined_mhs.push_back(combined_mh);
            combined_mh = NULL;
        }

        mhi++;
    }

    // last bit of cleanup.
    if (combined_mh) {
        if (combined_mh->tags.size()) {
            combined_mhs.push_back(combined_mh);
            combined_mh = NULL;
        } else {
            delete combined_mh;
        }

    }
}

void combine_from_tags(NeighborhoodMinHash& nbhd_mh,
                       TagSet& tagset,
                       CombinedMinHash& combined_mh,
                       unsigned int combined_minhash_size=500)
{
    unsigned int k = nbhd_mh.tag_to_mh.begin()->second->ksize;
    unsigned int p = nbhd_mh.tag_to_mh.begin()->second->prime;
    bool prot = nbhd_mh.tag_to_mh.begin()->second->is_protein;

    combined_mh.mh = new KmerMinHash(combined_minhash_size, k, p, prot);

    TagSet::iterator ti;
    for (ti = tagset.begin(); ti != tagset.end(); ti++) {
        KmerMinHash * mh = nbhd_mh.tag_to_mh[*ti];
        combined_mh.mh->merge(*mh);
        combined_mh.tags.insert(*ti);
    }
}

void save_neighborhood(const char * filename,
                       NeighborhoodMinHash * nbhd_mh)
{
    std::ofstream outfile(filename, std::ios::binary);

    outfile.write(SAVED_SIGNATURE, 4);
    unsigned char version = SAVED_FORMAT_VERSION;
    outfile.write((const char *) &version, 1);

    unsigned char ht_type = SAVED_NEIGHBORHOOD_HASH;
    outfile.write((const char *) &ht_type, 1);

    ///

    unsigned long num = nbhd_mh->tag_connections.size();
    outfile.write((const char *) &num, sizeof(num));

    for (TagToTagSet::iterator ni = nbhd_mh->tag_connections.begin();
            ni != nbhd_mh->tag_connections.end(); ni++) {
        TagSet * ts = &(ni->second);

        HashIntoType tag = ni->first;
        outfile.write((const char *) &tag, sizeof(tag));
        unsigned int ts_size = ts->size();
        outfile.write((const char *) &ts_size, sizeof(ts_size));
        for (TagSet::iterator tsi = ts->begin(); tsi != ts->end(); tsi++) {
            HashIntoType otag = *tsi;
            outfile.write((const char *) &otag, sizeof(otag));
        }
    }

    KmerMinHash * mh = nbhd_mh->tag_to_mh.begin()->second;
    unsigned int mh_max_size = mh->num;
    outfile.write((const char *) &num, sizeof(mh_max_size));

    unsigned int ksize = mh->ksize;
    outfile.write((const char *) &ksize, sizeof(ksize));

    long int prime = mh->prime;
    outfile.write((const char *) &prime, sizeof(prime));

    char is_protein = mh->is_protein ? 1 : 0;
    outfile.write((const char *) &is_protein, sizeof(is_protein));

    for (TagToMinHash::iterator mhi = nbhd_mh->tag_to_mh.begin();
            mhi != nbhd_mh->tag_to_mh.end(); mhi++) {
        HashIntoType tag = mhi->first;
        outfile.write((const char *) &tag, sizeof(tag));

        mh = mhi->second;
        CMinHashType * sketch = &(mh->mins);

        unsigned int mh_size = sketch->size();
        outfile.write((const char*) &mh_size, sizeof(mh_size));

        for (CMinHashType::iterator si = sketch->begin();
                si != sketch->end(); si++) {
            HashIntoType h = *si;
            outfile.write((const char*) &h, sizeof(h));
        }
    }

    if (outfile.fail()) {
        throw khmer_file_exception(strerror(errno));
    }
    outfile.close();
}

void load_neighborhood(const char * infilename_c,
                       NeighborhoodMinHash * nbhd_mh)
{
    std::string infilename(infilename_c);
    std::ifstream infile;

    // configure ifstream to raise exceptions for everything.
    infile.exceptions(std::ifstream::failbit | std::ifstream::badbit |
                      std::ifstream::eofbit);

    try {
        infile.open(infilename, std::ios::binary);
    } catch (std::ifstream::failure &e) {
        std::string err;
        if (!(infile.is_open())) {
            err = "Cannot open file: " + infilename;
        } else {
            err = "Unknown error in opening file: " + infilename;
        }
        throw khmer_file_exception(err);
    }

    unsigned char version, ht_type;

    try {
        char signature[4];
        infile.read(signature, 4);
        infile.read((char *) &version, 1);
        infile.read((char *) &ht_type, 1);
        if (!(std::string(signature, 4) == SAVED_SIGNATURE)) {
            std::ostringstream err;
            err << "Incorrect file signature 0x";
            for(size_t i=0; i < 4; ++i) {
                err << std::hex << (int) signature[i];
            }
            err << " while reading from " << infilename
                << "; should be " << SAVED_SIGNATURE;
            throw khmer_file_exception(err.str());
        } else if (!(version == SAVED_FORMAT_VERSION)) {
            std::ostringstream err;
            err << "Incorrect file format version " << (int) version
                << " while reading from " << infilename
                << "; should be " << (int) SAVED_FORMAT_VERSION;
            throw khmer_file_exception(err.str());
        } else if (!(ht_type == SAVED_NEIGHBORHOOD_HASH)) {
            std::ostringstream err;
            err << "Incorrect file format type " << (int) ht_type
                << " while reading tagset from " << infilename;
            throw khmer_file_exception(err.str());
        }


        unsigned long tag_connections_to_load = 0;
        infile.read((char *) &tag_connections_to_load,
                    sizeof(tag_connections_to_load));

        for (unsigned int i = 0; i < tag_connections_to_load; i++) {
            HashIntoType tag = 0;
            infile.read((char *) &tag, sizeof(tag));
            unsigned int ts_size = 0;
            infile.read((char *) &ts_size, sizeof(ts_size));

            HashIntoType tag2;
            TagSet tagset;
            for (unsigned int j = 0; j < ts_size; j++) {
                infile.read((char *) &tag2, sizeof(tag2));
                tagset.insert(tag2);
            }
            nbhd_mh->tag_connections[tag] = tagset;
        }

        unsigned int mh_max_size = 0;
        infile.read((char *) &mh_max_size, sizeof(mh_max_size));

        unsigned int ksize = 0;
        infile.read((char *) &ksize, sizeof(ksize));

        long int prime = 0;
        infile.read((char *) &prime, sizeof(prime));

        char is_protein_ch;
        infile.read((char *) &is_protein_ch, sizeof(is_protein_ch));
        bool is_protein = is_protein_ch;

        for (unsigned int i = 0; i < tag_connections_to_load; i++) {
            KmerMinHash * mh = new KmerMinHash(mh_max_size, ksize, prime,
                                               is_protein);

            HashIntoType tag = 0;
            infile.read((char *) &tag, sizeof(tag));

            CMinHashType * sketch = &(mh->mins);
            unsigned int mh_size = 0;
            infile.read((char *) &mh_size, sizeof(mh_size));

            HashIntoType h;
            for (unsigned int i = 0; i < mh_size; i++) {

                infile.read((char *) &h, sizeof(h));
                sketch->insert(h);
            }

            nbhd_mh->tag_to_mh[tag] = mh;
        }
    } catch (std::ifstream::failure &e) {
        std::string err = "Error reading data from: " + infilename;
        throw khmer_file_exception(err);
    }

}


static
PyObject *
nbhd_build_combined_minhashes(NeighborhoodMinHash_Object * me,
                              PyObject * args)
{
    unsigned int new_minhash_size;
    if (!PyArg_ParseTuple(args, "I", &new_minhash_size)) {
        return NULL;
    }

    std::vector<CombinedMinHash *> combined_mhs;
    NeighborhoodMinHash * nbhd_mh = me->nbhd_mh;

    Py_BEGIN_ALLOW_THREADS

    build_combined_minhashes(*nbhd_mh, combined_mhs, new_minhash_size);

    std::cout << "went from " << nbhd_mh->tag_to_mh.size() << " to "
              << combined_mhs.size() << " merged mhs.\n";

    Py_END_ALLOW_THREADS

    PyObject * list_of_mhs = PyList_New(combined_mhs.size());
    for (unsigned int i = 0; i < combined_mhs.size(); i++) {
        PyList_SET_ITEM(list_of_mhs, i,
                        build_CombinedMinHash_Object(combined_mhs[i]));
    }


    return list_of_mhs;
}

static
PyObject *
nbhd_combine_from_tags(NeighborhoodMinHash_Object * me, PyObject * args)
{
    unsigned int new_minhash_size;
    PyObject * list_of_tags_o;

    if (!PyArg_ParseTuple(args, "IO",
                          &new_minhash_size, &list_of_tags_o)) {
        return NULL;
    }

    if (!PyList_Check(list_of_tags_o)) {
        return NULL;
    }

    Py_ssize_t len = PyList_Size(list_of_tags_o);
    TagSet tags;
    for (Py_ssize_t i = 0; i < len; i++) {
        PyObject * tag_o = PyList_GetItem(list_of_tags_o, i);
        HashIntoType tag = PyLong_AsUnsignedLongLong(tag_o);
        tags.insert(tag);
    }

    CombinedMinHash * combined_mh = new CombinedMinHash;
    NeighborhoodMinHash * nbhd_mh = me->nbhd_mh;

    Py_BEGIN_ALLOW_THREADS

    combine_from_tags(*nbhd_mh, tags, *combined_mh, new_minhash_size);

    Py_END_ALLOW_THREADS

    return build_CombinedMinHash_Object(combined_mh);
}

static PyObject * nbhd_load(PyObject * self, PyObject * args)
{
    const char * filename = NULL;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    NeighborhoodMinHash * nbhd_mh = new NeighborhoodMinHash;
    load_neighborhood(filename, nbhd_mh);

    return build_NeighborhoodMinHash_Object(nbhd_mh);
}

static PyObject * nbhd_save(NeighborhoodMinHash_Object * me, PyObject * args)
{
    const char * filename = NULL;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    save_neighborhood(filename, me->nbhd_mh);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject * nbhd_search(NeighborhoodMinHash_Object *me,
                              PyObject * args)
{
    double threshold;
    PyObject * search_mh_obj;

    if (!PyArg_ParseTuple(args, "Od", &search_mh_obj, &threshold)) {
        return NULL;
    }

    if (!check_IsMinHash(search_mh_obj)) {
        return NULL;
    }

    NeighborhoodMinHash * nbhd_mh = me->nbhd_mh;
    const KmerMinHash * search_mh = extract_KmerMinHash(search_mh_obj);

    TagSet matching_tags;
    for (TagToMinHash::const_iterator tsi = nbhd_mh->tag_to_mh.begin();
            tsi != nbhd_mh->tag_to_mh.end(); tsi++) {
        unsigned int common = tsi->second->count_common(*search_mh);
        float similarity = float(common) / float(tsi->second->mins.size());
        if (similarity >= threshold) {
            matching_tags.insert(tsi->first);
        }
    }

    PyObject * tags_o = PyList_New(matching_tags.size());
    unsigned int j = 0;
    for (TagSet::iterator ti = matching_tags.begin();
            ti != matching_tags.end(); ++ti) {
        PyList_SET_ITEM(tags_o, j, PyLong_FromUnsignedLongLong(*ti));
        j++;
    }

    return tags_o;
}

static PyObject * nbhd_get_minhash(NeighborhoodMinHash_Object * me,
                                   PyObject * args)
{
    HashIntoType tag;
    if (!PyArg_ParseTuple(args, "K", &tag)) {
        return NULL;
    }

    TagToMinHash * tag_to_mh = &(me->nbhd_mh->tag_to_mh);
    TagToMinHash::iterator tm = tag_to_mh->find(tag);
    if (tm != tag_to_mh->end()) {
        KmerMinHash * mh = tm->second;
        KmerMinHash * new_mh = new KmerMinHash(mh->num, mh->ksize, mh->prime,
                                               mh->is_protein);
        new_mh->merge(*mh);

        return build_MinHash_Object(new_mh);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject * nbhd_get_connected_tags(NeighborhoodMinHash_Object * me,
        PyObject * args)
{
    HashIntoType tag;
    if (!PyArg_ParseTuple(args, "K", &tag)) {
        return NULL;
    }

    TagSet this_conn = me->nbhd_mh->tag_connections[tag];

    PyObject * tags_o = PyList_New(this_conn.size());
    unsigned int i = 0;
    for (TagSet::const_iterator tsi = this_conn.begin();
            tsi != this_conn.end(); tsi++) {
        PyList_SET_ITEM(tags_o, i, PyLong_FromUnsignedLongLong(*tsi));
        i++;
    }

    return tags_o;
}

static PyObject * nbhd_get_all_tags(NeighborhoodMinHash_Object * me,
                                    PyObject * args)
{
    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    TagToTagSet * tag_conn = &(me->nbhd_mh->tag_connections);

    PyObject * tags_o = PyList_New(tag_conn->size());
    unsigned int i = 0;

    TagToTagSet::iterator tsi = tag_conn->begin();
    for (; tsi != tag_conn->end(); tsi++) {
        PyList_SET_ITEM(tags_o, i, PyLong_FromUnsignedLongLong(tsi->first));
        i++;
    }

    return tags_o;
}

static PyMethodDef NeighborhoodMinHash_methods [] = {
    {
        "build_combined_minhashes",
        (PyCFunction)nbhd_build_combined_minhashes,
        METH_VARARGS, "Combine neighborhood minhashes"
    },
    {
        "combine_from_tags",
        (PyCFunction)nbhd_combine_from_tags,
        METH_VARARGS, "Combine nbhd minhashes from tags into one MinHash obj"
    },
    {
        "save",
        (PyCFunction)nbhd_save,
        METH_VARARGS, "Save nbhd hash to disk."
    },
    {
        "get_minhash",
        (PyCFunction)nbhd_get_minhash,
        METH_VARARGS, "Get the MinHash for this tag."
    },
    {
        "get_connected_tags",
        (PyCFunction)nbhd_get_connected_tags,
        METH_VARARGS, "Get the connections for this tag."
    },
    {
        "get_all_tags",
        (PyCFunction)nbhd_get_all_tags,
        METH_VARARGS, "Get a list of all tags."
    },
    {
        "search",
        (PyCFunction)nbhd_search,
        METH_VARARGS, "Find all tags with MinHash comparisons above threshold"
    },
    { NULL, NULL, 0, NULL } // sentinel
};

static
PyObject *
NeighborhoodMinHash_new(PyTypeObject * subtype, PyObject * args,
                        PyObject * kwds)
{
    PyObject * self     = subtype->tp_alloc( subtype, 1 );
    if (self == NULL) {
        return NULL;
    }

    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    NeighborhoodMinHash_Object * myself = (NeighborhoodMinHash_Object *)self;

    myself->nbhd_mh = new NeighborhoodMinHash;

    return self;
}

bool check_IsNeighborhoodMinHash(PyObject * nbhd_mh)
{
    if (!PyObject_TypeCheck(nbhd_mh, &NeighborhoodMinHash_Type)) {
        return false;
    }
    return true;
}

///////////////////////////////////////////////////////////////////////////////

static
void
CombinedMinHash_dealloc(CombinedMinHash_Object * obj)
{
    delete obj->combined_mh;
    obj->combined_mh = NULL;
    Py_TYPE(obj)->tp_free((PyObject*)obj);
}

static
PyObject *
combined_get_tags(CombinedMinHash_Object * me, PyObject * args)
{
    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    CombinedMinHash * combined_mh = me->combined_mh;
    TagSet * tags = &(combined_mh->tags);
    PyObject * tags_o = PyList_New(tags->size());

    unsigned int j = 0;
    for (TagSet::iterator tsi = tags->begin(); tsi != tags->end(); ++tsi) {
        PyList_SET_ITEM(tags_o, j, PyLong_FromUnsignedLongLong(*tsi));
        j++;
    }

    return(tags_o);
}

static
PyObject *
combined_get_tag_neighbors(CombinedMinHash_Object * me, PyObject * args)
{
    PyObject * nbhd_mh_o;
    if (!PyArg_ParseTuple(args, "O", &nbhd_mh_o)) {
        return NULL;
    }

    if (!check_IsNeighborhoodMinHash(nbhd_mh_o)) {
        return NULL;
    }

    CombinedMinHash * combined_mh = me->combined_mh;
    NeighborhoodMinHash * nbhd_mh = extract_NeighborhoodMinHash(nbhd_mh_o);

    TagSet * my_tags = &(combined_mh->tags);
    TagSet neighbor_tags;

    // look at all tags in this combined tag set
    for (TagSet::const_iterator mti = my_tags->begin();
            mti != my_tags->end(); mti++) {

        // ...and grab all their neighbors.
        TagSet * otags = &(nbhd_mh->tag_connections[*mti]);
        for (TagSet::const_iterator oti = otags->begin();
                oti != otags->end(); oti++) {
            if (my_tags->find(*oti) == my_tags->end()) { // not in set already.
                neighbor_tags.insert(*mti);
            }
        }
    }

    PyObject * tags_o = PyList_New(neighbor_tags.size());
    unsigned int j = 0;
    for (TagSet::iterator tsi = neighbor_tags.begin();
            tsi != neighbor_tags.end(); ++tsi) {
        PyList_SET_ITEM(tags_o, j, PyLong_FromUnsignedLongLong(*tsi));
        j++;
    }

    return(tags_o);
}

static
PyObject *
combined_extend(CombinedMinHash_Object * me, PyObject * args)
{
    PyObject * nbhd_mh_o;
    PyObject * list_o;
    if (!PyArg_ParseTuple(args, "OO", &nbhd_mh_o, &list_o)) {
        return NULL;
    }

    if (!check_IsNeighborhoodMinHash(nbhd_mh_o)) {
        return NULL;
    }

    if (!PyList_Check(list_o)) {
        return NULL;
    }

    CombinedMinHash * combined_mh = me->combined_mh;
    NeighborhoodMinHash * nbhd_mh = extract_NeighborhoodMinHash(nbhd_mh_o);

    bool did_extend = false;
    for (int i = 0; i < PyList_Size(list_o); i++) {
        PyObject * o = PyList_GET_ITEM(list_o, i);
        HashIntoType tag = PyLong_AsUnsignedLongLong(o);

        if (nbhd_mh->tag_connections.find(tag) == nbhd_mh->tag_connections.end()) {
            continue;
        }

        combined_mh->tags.insert(tag);
        combined_mh->mh->merge(*nbhd_mh->tag_to_mh[tag]);
        nbhd_mh->tag_connections.erase(tag);
        did_extend = true;
    }

    if (did_extend) {
        Py_INCREF(Py_True);
        return Py_True;
    }
    Py_INCREF(Py_False);
    return Py_False;
}

static PyObject * combined_get_minhash(CombinedMinHash_Object * me,
                                       PyObject * args)
{
    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    KmerMinHash * mh = me->combined_mh->mh;
    KmerMinHash * new_mh = new KmerMinHash(mh->num, mh->ksize, mh->prime,
                                           mh->is_protein);
    new_mh->merge(*mh);

    return build_MinHash_Object(new_mh);
}

static PyObject * combined_combine(CombinedMinHash_Object * me,
                                   PyObject * args)
{
    PyObject * second_o;
    unsigned int new_minhash_size = 0;
    if (!PyArg_ParseTuple(args, "O|I", &second_o, &new_minhash_size)) {
        return NULL;
    }
    if (!check_IsCombinedMinHash(second_o)) {
        return NULL;
    }

    CombinedMinHash * first = me->combined_mh;
    if (new_minhash_size == 0) {
        new_minhash_size = first->mh->num;
    }
    CombinedMinHash * second = ((CombinedMinHash_Object *)second_o)->combined_mh;

    CombinedMinHash * combined = new CombinedMinHash;
    for (TagSet::const_iterator tsi = first->tags.begin();
            tsi != first->tags.end(); tsi++) {
        combined->tags.insert(*tsi);
    }
    for (TagSet::const_iterator tsi = second->tags.begin();
            tsi != second->tags.end(); tsi++) {
        combined->tags.insert(*tsi);
    }

    KmerMinHash * new_mh = new KmerMinHash(new_minhash_size,
                                           first->mh->ksize,
                                           first->mh->prime,
                                           first->mh->is_protein);
    new_mh->merge(*first->mh);
    new_mh->merge(*second->mh);
    combined->mh = new_mh;

    return build_CombinedMinHash_Object(combined);
}

static PyMethodDef CombinedMinHash_methods [] = {
    {
        "get_tags",
        (PyCFunction)combined_get_tags,
        METH_VARARGS, "Get the list of tags in this combined set."
    },
    {
        "get_tag_neighbors",
        (PyCFunction)combined_get_tag_neighbors,
        METH_VARARGS, "Get the list of neighboring tags in this combined set."
    },
    {
        "get_minhash",
        (PyCFunction)combined_get_minhash,
        METH_VARARGS, "Get the MinHash for this combined set of tags."
    },
    {
        "combine",
        (PyCFunction)combined_combine,
        METH_VARARGS, "Combine this and another Combined MinHash into a new one."
    },
    {
        "extend",
        (PyCFunction)combined_extend,
        METH_VARARGS, "Add the given tags into this combined MinHash."
    },
    { NULL, NULL, 0, NULL } // sentinel
};

static
PyObject *
CombinedMinHash_new(PyTypeObject * subtype, PyObject * args,
                    PyObject * kwds)
{
    PyObject * self     = subtype->tp_alloc( subtype, 1 );
    if (self == NULL) {
        return NULL;
    }

    unsigned int _n, _ksize;
    long int _p;
    PyObject * is_protein_o;
    if (!PyArg_ParseTuple(args, "IIlO", &_n, &_ksize, &_p, &is_protein_o)) {
        return NULL;
    }

    CombinedMinHash_Object * myself = (CombinedMinHash_Object *)self;

    myself->combined_mh = new CombinedMinHash;

    return self;
}

bool check_IsCombinedMinHash(PyObject * combined_mh)
{
    if (!PyObject_TypeCheck(combined_mh, &CombinedMinHash_Type)) {
        return false;
    }
    return true;
}

///

static PyMethodDef MinHashModuleMethods[] = {
    {
        "load_neighborhood_minhash",
        (PyCFunction)nbhd_load,
        METH_VARARGS, "load nbhd hash from disk."
    },
    { NULL, NULL, 0, NULL } // sentinel
};

MOD_INIT(_minhash)
{
    MinHash_Type.tp_methods = MinHash_methods;
    MinHash_Type.tp_dealloc = (destructor)MinHash_dealloc;
    MinHash_Type.tp_new = MinHash_new;

    if (PyType_Ready( &MinHash_Type ) < 0) {
        return MOD_ERROR_VAL;
    }

    NeighborhoodMinHash_Type.tp_methods = NeighborhoodMinHash_methods;
    NeighborhoodMinHash_Type.tp_dealloc = (destructor)NeighborhoodMinHash_dealloc;
    NeighborhoodMinHash_Type.tp_new = NeighborhoodMinHash_new;

    if (PyType_Ready( &NeighborhoodMinHash_Type ) < 0) {
        return MOD_ERROR_VAL;
    }

    CombinedMinHash_Type.tp_methods = CombinedMinHash_methods;
    CombinedMinHash_Type.tp_dealloc = (destructor)CombinedMinHash_dealloc;
    CombinedMinHash_Type.tp_new = CombinedMinHash_new;
    if (PyType_Ready( &CombinedMinHash_Type ) < 0) {
        return MOD_ERROR_VAL;
    }

    PyObject * m;

    MOD_DEF(m, "_minhash",
            "interface for the sourmash module low-level extensions",
            MinHashModuleMethods);

    if (m == NULL) {
        return MOD_ERROR_VAL;
    }

    Py_INCREF(&MinHash_Type);
    if (PyModule_AddObject( m, "MinHash",
                            (PyObject *)&MinHash_Type ) < 0) {
        return MOD_ERROR_VAL;
    }

    Py_INCREF(&NeighborhoodMinHash_Type);
    if (PyModule_AddObject( m, "NeighborhoodMinHash",
                            (PyObject *)&NeighborhoodMinHash_Type ) < 0) {
        return MOD_ERROR_VAL;
    }

    Py_INCREF(&CombinedMinHash_Type);
    if (PyModule_AddObject( m, "CombinedMinHash",
                            (PyObject *)&CombinedMinHash_Type ) < 0) {
        return MOD_ERROR_VAL;
    }
    return MOD_SUCCESS_VAL(m);
}
