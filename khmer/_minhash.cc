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

bool check_IsMinHash(PyObject * mh);

////

static
void
MinHash_dealloc(MinHash_Object * obj)
{
  delete obj->mh;
  obj->mh = NULL;
  Py_TYPE(obj)->tp_free((PyObject*)obj);
}

static std::map<std::string, std::string> * _codon_table = NULL;

static std::string _dna_to_aa(const std::string& dna)
{
  if (_codon_table == NULL) {
    _codon_table = new std::map<std::string, std::string>;
    (*_codon_table)["TTT"] = "F";
    (*_codon_table)["TTC"] = "F";
    (*_codon_table)["TTA"] = "L";
    (*_codon_table)["TTG"] = "L";
    
    (*_codon_table)["TCT"] = "S";
    (*_codon_table)["TCC"] = "S";
    (*_codon_table)["TCA"] = "S";
    (*_codon_table)["TCG"] = "S";
    
    (*_codon_table)["TAT"] = "Y";
    (*_codon_table)["TAC"] = "Y";
    (*_codon_table)["TAA"] = "*";
    (*_codon_table)["TAG"] = "*";
    
    (*_codon_table)["TGT"] = "C";
    (*_codon_table)["TGC"] = "C";
    (*_codon_table)["TGA"] = "*";
    (*_codon_table)["TGG"] = "W";
    
    (*_codon_table)["CTT"] = "L";
    (*_codon_table)["CTC"] = "L";
    (*_codon_table)["CTA"] = "L";
    (*_codon_table)["CTG"] = "L";
    
    (*_codon_table)["CCT"] = "P";
    (*_codon_table)["CCC"] = "P";
    (*_codon_table)["CCA"] = "P";
    (*_codon_table)["CCG"] = "P";
    
    (*_codon_table)["CAT"] = "H";
    (*_codon_table)["CAC"] = "H";
    (*_codon_table)["CAA"] = "Q";
    (*_codon_table)["CAG"] = "Q";
    
    (*_codon_table)["CGT"] = "R";
    (*_codon_table)["CGC"] = "R";
    (*_codon_table)["CGA"] = "R";
    (*_codon_table)["CGG"] = "R";
    
    (*_codon_table)["ATT"] = "I";
    (*_codon_table)["ATC"] = "I";
    (*_codon_table)["ATA"] = "I";
    (*_codon_table)["ATG"] = "M";
    
    (*_codon_table)["ACT"] = "T";
    (*_codon_table)["ACC"] = "T";
    (*_codon_table)["ACA"] = "T";
    (*_codon_table)["ACG"] = "T";
    
    (*_codon_table)["AAT"] = "N";
    (*_codon_table)["AAC"] = "N";
    (*_codon_table)["AAA"] = "K";
    (*_codon_table)["AAG"] = "K";
    
    (*_codon_table)["AGT"] = "S";
    (*_codon_table)["AGC"] = "S";
    (*_codon_table)["AGA"] = "R";
    (*_codon_table)["AGG"] = "R";
    
    (*_codon_table)["GTT"] = "V";
    (*_codon_table)["GTC"] = "V";
    (*_codon_table)["GTA"] = "V";
    (*_codon_table)["GTG"] = "V";
    
    (*_codon_table)["GCT"] = "A";
    (*_codon_table)["GCC"] = "A";
    (*_codon_table)["GCA"] = "A";
    (*_codon_table)["GCG"] = "A";
    
    (*_codon_table)["GAT"] = "D";
    (*_codon_table)["GAC"] = "D";
    (*_codon_table)["GAA"] = "E";
    (*_codon_table)["GAG"] = "E";
    
    (*_codon_table)["GGT"] = "G";
    (*_codon_table)["GGC"] = "G";
    (*_codon_table)["GGA"] = "G";
    (*_codon_table)["GGG"] = "G";
  }

  std::string aa;
  for (unsigned int j = 0; j < dna.size(); j += 3) {
    std::string codon = dna.substr(j, 3);
    aa += (*_codon_table)[codon];
  }
  return aa;
}

std::string _revcomp(const std::string& kmer);

static
PyObject *
minhash_add_sequence(MinHash_Object * me, PyObject * args)
{
  const char * sequence = NULL;
  if (!PyArg_ParseTuple(args, "s", &sequence)) {
    return NULL;
  }
  KmerMinHash * mh = me->mh;
  CMinHashType::iterator mins_end;
  
  long int h = 0;
  unsigned int ksize = mh->ksize;

  if (!mh->is_protein) {
    std::string seq = sequence;
    for (unsigned int i = 0; i < seq.length() - ksize + 1; i++) {
      std::string kmer = seq.substr(i, ksize);
      mh->add_kmer(kmer);
    }
    std::string rc = _revcomp(seq);
    for (unsigned int i = 0; i < rc.length() - ksize + 1; i++) {
      std::string kmer = rc.substr(i, ksize);
      mh->add_kmer(kmer);
    }
  } else {                      // protein
    std::string seq = sequence;
    for (unsigned int i = 0; i < seq.length() - ksize + 1; i ++) {
      std::string kmer = seq.substr(i, ksize);
      std::string aa = _dna_to_aa(kmer);

      mh->add_kmer(aa);

      std::string rc = _revcomp(kmer);
      aa = _dna_to_aa(rc);

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
  { "add_sequence",
    (PyCFunction)minhash_add_sequence, METH_VARARGS,
    "Add kmer into MinHash"
  },
  { "add_hash",
    (PyCFunction)minhash_add_hash, METH_VARARGS,
    "Add kmer into MinHash"
  },
  { "get_mins",
    (PyCFunction)minhash_get_mins, METH_VARARGS,
    "Get MinHash signature"
  },
  { "__copy__",
    (PyCFunction)minhash___copy__, METH_VARARGS,
    "Copy this MinHash object",
  },
  { "__len__",                  // @CTB this doesn't work...?
    (PyCFunction)minhash___len__, METH_VARARGS,
    "Number of hashes in this MinHash object",
  },
  { "count_common",
    (PyCFunction)minhash_count_common, METH_VARARGS,
    "Get number of hashes in common with other."
  },
  { "compare",
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
    long int _p;
    PyObject * is_protein_o;
    if (!PyArg_ParseTuple(args, "IIlO", &_n, &_ksize, &_p, &is_protein_o)){
      return NULL;
    }
    
    MinHash_Object * myself = (MinHash_Object *)self;

    myself->mh = new KmerMinHash(_n, _ksize, _p,\
                                 PyObject_IsTrue(is_protein_o));

    return self;
}

static PyTypeObject MinHash_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)        /* init & ob_size */
    "_minhash.MinHash",                    /* tp_name */
    sizeof(MinHash_Object),         /* tp_basicsize */
    0,                                    /* tp_itemsize */
    (destructor)MinHash_dealloc,    /* tp_dealloc */
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
    0,                                    /* tp_traverse */
    0,                                    /* tp_clear */
    0,                                    /* tp_richcompare */
    0,                                    /* tp_weaklistoffset */
    0,                                    /* tp_iter */
    0,                                    /* tp_iternext */
    MinHash_methods,               /* tp_methods */
    0,                                    /* tp_members */
    0,                                    /* tp_getset */
    0,                                         /* tp_base */
    0,                                         /* tp_dict */
    0,                                         /* tp_descr_get */
    0,                                         /* tp_descr_set */
    0,                                         /* tp_dictoffset */
    0,                                         /* tp_init */
    0,                                         /* tp_alloc */
    MinHash_new,                        /* tp_new */
};

PyObject * build_MinHash_Object(KmerMinHash * mh)
{
  MinHash_Object * obj = (MinHash_Object *) \
    PyObject_New(MinHash_Object, &MinHash_Type);
  obj->mh = mh;

  return (PyObject *) obj;
}

bool check_IsMinHash(PyObject * mh)
{
  if (!PyObject_TypeCheck(mh, &MinHash_Type)) {
    return false;
  }
  return true;
}

khmer::KmerMinHash * extract_KmerMinHash(PyObject * mh_obj)
{
  if (!PyObject_TypeCheck(mh_obj, &MinHash_Type)) {
    return NULL;
  }
  MinHash_Object * obj = (MinHash_Object *) mh_obj;
  return obj->mh;
}

std::string _revcomp(const std::string& kmer)
{
    std::string out = kmer;
    size_t ksize = out.size();

    for (size_t i=0; i < ksize; ++i) {
        char complement;

        switch(kmer[i]) {
        case 'A':
            complement = 'T';
            break;
        case 'C':
            complement = 'G';
            break;
        case 'G':
            complement = 'C';
            break;
        case 'T':
            complement = 'A';
            break;
        default:
            throw std::exception();
            break;
        }
        out[ksize - i - 1] = complement;
    }
    return out;
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

void build_combined_minhashes(NeighborhoodMinHash& nbhd_mh,
                              std::vector<CombinedMinHash *>& combined_mhs,
                              unsigned int combine_this_many_tags=10000,
                              unsigned int combined_minhash_size=500,
                              HashIntoType start_tag=0)
{
  unsigned int k = nbhd_mh.tag_to_mh.begin()->second->ksize;
  unsigned int p = nbhd_mh.tag_to_mh.begin()->second->prime;
  bool prot = nbhd_mh.tag_to_mh.begin()->second->is_protein;

  CombinedMinHash * combined_mh = NULL;
  unsigned int total_combined = 0;
  unsigned int combined_tags = 0;

  // iterate over all tags:
  TagToMinHash::const_iterator mhi;
  if (start_tag == 0) {
    mhi = nbhd_mh.tag_to_mh.begin();
  } else {
    mhi = nbhd_mh.tag_to_mh.find(start_tag);
  }
  while(mhi != nbhd_mh.tag_to_mh.end()) {
    HashIntoType start_tag = mhi->first;

    // already merged this 'un? move to next.
    TagToTagSet::iterator posn;
    posn = nbhd_mh.tag_connections.find(mhi->first);
    if (posn == nbhd_mh.tag_connections.end()) {
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
    
    // & clear to-be-merged tag from list
    nbhd_mh.tag_connections.erase(posn);

    // merge nbhd minhashes in, tag by tag, until stop.
    bool did_combine = true;
    while (combined_tags < combine_this_many_tags && did_combine) {
      did_combine = false;

      // walk through the list of tags to be merged:
      TagSet::iterator ti;
      while(combined_tags < combine_this_many_tags && to_be_merged.size()) {
        // grab the first tag & its list of connected tags:
        TagToTagSet::iterator posn2;
        ti = to_be_merged.begin();
        posn2 = nbhd_mh.tag_connections.find(*ti);

        // already merged & removed from nbhd_mh.tag_connections? ok, ignore.
        if (posn2 == nbhd_mh.tag_connections.end()) {
          to_be_merged.erase(ti);
          continue;
        }

        // finally! merge.
        KmerMinHash * mh = nbhd_mh.tag_to_mh[*ti];
        combined_mh->mh->merge(*mh);
        combined_mh->tags.insert(*ti);

        // track info.
        did_combine = true;
        combined_tags++;
        total_combined++;

        // add the just-merged one's connected tags to to_be_merged
        to_be_merged.insert(posn2->second.begin(), posn2->second.end());
        
        // remove:
        nbhd_mh.tag_connections.erase(posn2);
        to_be_merged.erase(ti);
      }
    }
    if (combined_tags >= combine_this_many_tags) {
      combined_mhs.push_back(combined_mh);
      combined_mh = NULL;
      combined_tags = 0;
    }

    mhi++;

    if (start_tag && mhi == nbhd_mh.tag_to_mh.end()) {
      mhi = nbhd_mh.tag_to_mh.begin();
    }
  }
  if (combined_mh) {
    combined_mhs.push_back(combined_mh);
    combined_mh = NULL;
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
  
#if 0

        infile.read((char *) &save_ksize, sizeof(save_ksize));
        if (!(save_ksize == _ksize)) {
            std::ostringstream err;
            err << "Incorrect k-mer size " << save_ksize
                << " while reading tagset from " << infilename;
            throw khmer_file_exception(err.str());
        }

        infile.read((char *) &tagset_size, sizeof(tagset_size));
        infile.read((char *) &_tag_density, sizeof(_tag_density));

        buf = new HashIntoType[tagset_size];

        infile.read((char *) buf, sizeof(HashIntoType) * tagset_size);

        for (unsigned int i = 0; i < tagset_size; i++) {
            all_tags.insert(buf[i]);
        }

        delete[] buf;
    } catch (std::ifstream::failure &e) {
        std::string err = "Error reading data from: " + infilename;
        if (buf != NULL) {
            delete[] buf;
        }
        throw khmer_file_exception(err);
    }
#endif
}


static
PyObject *
nbhd_build_combined_minhashes(NeighborhoodMinHash_Object * me, PyObject * args)
{
    unsigned int num_tags_to_combine, new_minhash_size;
    unsigned int start_tag = 0;
    if (!PyArg_ParseTuple(args, "II|K",
                          &num_tags_to_combine, &new_minhash_size,
                          &start_tag)) {
        return NULL;
    }

    std::vector<CombinedMinHash *> combined_mhs;
    NeighborhoodMinHash * nbhd_mh = me->nbhd_mh;

    Py_BEGIN_ALLOW_THREADS

    build_combined_minhashes(*nbhd_mh, combined_mhs,
                             num_tags_to_combine, new_minhash_size,
                             start_tag=start_tag);

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

static PyMethodDef NeighborhoodMinHash_methods [] = {
  { "build_combined_minhashes",
    (PyCFunction)nbhd_build_combined_minhashes,
    METH_VARARGS, "Combine neighborhood minhashes" },
  { "combine_from_tags",
    (PyCFunction)nbhd_combine_from_tags,
    METH_VARARGS, "Combine nbhd minhashes from tags into one MinHash obj" },
  { "save",
    (PyCFunction)nbhd_save,
    METH_VARARGS, "Save nbhd hash to disk." },
  { "get_minhash",
    (PyCFunction)nbhd_get_minhash,
    METH_VARARGS, "Get the MinHash for this tag." },
  { "search",
    (PyCFunction)nbhd_search,
    METH_VARARGS, "Find all tags with MinHash comparisons above threshold" },
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

    unsigned int _n, _ksize;
    long int _p;
    PyObject * is_protein_o;
    if (!PyArg_ParseTuple(args, "IIlO", &_n, &_ksize, &_p, &is_protein_o)){
      return NULL;
    }
    
    NeighborhoodMinHash_Object * myself = (NeighborhoodMinHash_Object *)self;

    myself->nbhd_mh = new NeighborhoodMinHash;

    return self;
}

static PyTypeObject NeighborhoodMinHash_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)        /* init & ob_size */
    "_minhash.NeighborhoodMinHash",       /* tp_name */
    sizeof(NeighborhoodMinHash_Object),   /* tp_basicsize */
    0,                                    /* tp_itemsize */
    (destructor)NeighborhoodMinHash_dealloc, /* tp_dealloc */
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
    "A MinHash sketch on tag neighborhoods.", /* tp_doc */
    0,                                    /* tp_traverse */
    0,                                    /* tp_clear */
    0,                                    /* tp_richcompare */
    0,                                    /* tp_weaklistoffset */
    0,                                    /* tp_iter */
    0,                                    /* tp_iternext */
    NeighborhoodMinHash_methods,          /* tp_methods */
    0,                                    /* tp_members */
    0,                                    /* tp_getset */
    0,                                    /* tp_base */
    0,                                    /* tp_dict */
    0,                                    /* tp_descr_get */
    0,                                    /* tp_descr_set */
    0,                                    /* tp_dictoffset */
    0,                                    /* tp_init */
    0,                                    /* tp_alloc */
    NeighborhoodMinHash_new,              /* tp_new */
};

PyObject * build_NeighborhoodMinHash_Object(NeighborhoodMinHash * nbhd_mh)
{
  NeighborhoodMinHash_Object * obj = (NeighborhoodMinHash_Object *) \
    PyObject_New(NeighborhoodMinHash_Object, &NeighborhoodMinHash_Type);
  obj->nbhd_mh = nbhd_mh;

  return (PyObject *) obj;
}

khmer::NeighborhoodMinHash * extract_NeighborhoodMinHash(PyObject * nbhd_obj)
{
  if (!PyObject_TypeCheck(nbhd_obj, &NeighborhoodMinHash_Type)) {
    return NULL;
  }
  NeighborhoodMinHash_Object * obj = (NeighborhoodMinHash_Object *) nbhd_obj;
  return obj->nbhd_mh;
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

static PyMethodDef CombinedMinHash_methods [] = {
  { "get_tags",
    (PyCFunction)combined_get_tags,
    METH_VARARGS, "Get the list of tags in this combined set." },
  { "get_minhash",
    (PyCFunction)combined_get_minhash,
    METH_VARARGS, "Get the MinHash for this combined set of tags." },
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
    if (!PyArg_ParseTuple(args, "IIlO", &_n, &_ksize, &_p, &is_protein_o)){
      return NULL;
    }
    
    CombinedMinHash_Object * myself = (CombinedMinHash_Object *)self;

    myself->combined_mh = new CombinedMinHash;

    return self;
}

static PyTypeObject CombinedMinHash_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)        /* init & ob_size */
    "_minhash.CombinedMinHash",       /* tp_name */
    sizeof(CombinedMinHash_Object),   /* tp_basicsize */
    0,                                    /* tp_itemsize */
    (destructor)CombinedMinHash_dealloc, /* tp_dealloc */
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
    "A MinHash sketch on tag neighborhoods.", /* tp_doc */
    0,                                    /* tp_traverse */
    0,                                    /* tp_clear */
    0,                                    /* tp_richcompare */
    0,                                    /* tp_weaklistoffset */
    0,                                    /* tp_iter */
    0,                                    /* tp_iternext */
    CombinedMinHash_methods,          /* tp_methods */
    0,                                    /* tp_members */
    0,                                    /* tp_getset */
    0,                                         /* tp_base */
    0,                                         /* tp_dict */
    0,                                         /* tp_descr_get */
    0,                                         /* tp_descr_set */
    0,                                         /* tp_dictoffset */
    0,                                         /* tp_init */
    0,                                         /* tp_alloc */
    CombinedMinHash_new,                        /* tp_new */
};

PyObject * build_CombinedMinHash_Object(CombinedMinHash * combined_mh)
{
  CombinedMinHash_Object * obj = (CombinedMinHash_Object *) \
    PyObject_New(CombinedMinHash_Object, &CombinedMinHash_Type);
  obj->combined_mh = combined_mh;

  return (PyObject *) obj;
}

khmer::CombinedMinHash * extract_CombinedMinHash(PyObject * combined_obj)
{
  if (!PyObject_TypeCheck(combined_obj, &CombinedMinHash_Type)) {
    return NULL;
  }
  CombinedMinHash_Object * obj = (CombinedMinHash_Object *) combined_obj;
  return obj->combined_mh;
}

///

static PyMethodDef MinHashModuleMethods[] = {
  { "load_neighborhood_minhash",
    (PyCFunction)nbhd_load,
    METH_VARARGS, "load nbhd hash from disk." },
  { NULL, NULL, 0, NULL } // sentinel
};

MOD_INIT(_minhash)
{
    if (PyType_Ready( &MinHash_Type ) < 0) {
        return MOD_ERROR_VAL;
    }
    
    if (PyType_Ready( &NeighborhoodMinHash_Type ) < 0) {
        return MOD_ERROR_VAL;
    }
    
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
    Py_INCREF(&MinHash_Type);
    
    if (PyModule_AddObject( m, "NeighborhoodMinHash",
                            (PyObject *)&NeighborhoodMinHash_Type ) < 0) {
        return MOD_ERROR_VAL;
    }
    Py_INCREF(&MinHash_Type);
    if (PyModule_AddObject( m, "CombinedMinHash",
                            (PyObject *)&CombinedMinHash_Type ) < 0) {
        return MOD_ERROR_VAL;
    }
    return MOD_SUCCESS_VAL(m);
}
