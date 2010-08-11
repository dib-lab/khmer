//
// A module for Python that exports khmer C++ library functions.
//

#include <iostream>

#include "Python.h"
#include "khmer.hh"
#include "ktable.hh"
#include "hashtable.hh"
#include "storage.hh"

//
// Function necessary for Python loading:
//

extern "C" {
  void init_khmer();
}

class _khmer_exception {
private:
  std::string _message;
public:
  _khmer_exception(std::string message) : _message(message) { };
  inline const std::string get_message() const { return _message; };
};

class _khmer_signal : public _khmer_exception {
public:
  _khmer_signal(std::string message) : _khmer_exception(message) { };
};

// Python exception to raise
static PyObject *KhmerError;

// default callback obj;
static PyObject *callback_obj = NULL;

// callback function to pass into C++ functions

void _report_fn(const char * info, void * data, unsigned int n_reads,
		unsigned long long other)
{
  // handle signals etc. (like CTRL-C)
  if (PyErr_CheckSignals() != 0) {
    throw _khmer_signal("PyErr_CheckSignals received a signal");
  }

  // set data to default?
  if (!data && callback_obj) {
    data = callback_obj;
  }

  // if 'data' is set, it is None, or a Python callable
  if (data) {
    PyObject * obj = (PyObject *) data;
    if (obj != Py_None) {
      PyObject * args = Py_BuildValue("siL", info, n_reads, other);

      PyObject * r = PyObject_Call(obj, args, NULL);
      Py_XDECREF(r);
      Py_DECREF(args);
    }
  }

  if (PyErr_Occurred()) {
    throw _khmer_signal("PyErr_Occurred is set");
  }

  // ...allow other Python threads to do stuff...
  Py_BEGIN_ALLOW_THREADS;
  Py_END_ALLOW_THREADS;
}



/***********************************************************************/

//
// KTable object -- exact counting of k-mers.
//

typedef struct {
  PyObject_HEAD
  khmer::KTable * ktable;
} khmer_KTableObject;

static void khmer_ktable_dealloc(PyObject *);

//
// ktable_forward_hash - hash k-mers into numbers
//

static PyObject * ktable_forward_hash(PyObject * self, PyObject * args)
{
  khmer_KTableObject * me = (khmer_KTableObject *) self;
  khmer::KTable * ktable = me->ktable;

  char * kmer;

  if (!PyArg_ParseTuple(args, "s", &kmer)) {
    return NULL;
  }

  if (strlen(kmer) != ktable->ksize()) {
    PyErr_SetString(PyExc_ValueError,
		    "k-mer length must be the same as the hashtable k-size");
    return NULL;
  }

  return PyInt_FromLong(khmer::_hash(kmer, ktable->ksize()));
}

//
// ktable_forward_hash_no_rc -- hash k-mers into numbers, with no RC handling
//

static PyObject * ktable_forward_hash_no_rc(PyObject * self, PyObject * args)
{
  khmer_KTableObject * me = (khmer_KTableObject *) self;
  khmer::KTable * ktable = me->ktable;

  char * kmer;

  if (!PyArg_ParseTuple(args, "s", &kmer)) {
    return NULL;
  }

  if (strlen(kmer) != ktable->ksize()) {
    PyErr_SetString(PyExc_ValueError,
		    "k-mer length must be the same as the hashtable k-size");
    return NULL;
  }

  return PyInt_FromLong(khmer::_hash_forward(kmer, ktable->ksize()));
}

//
// ktable_reverse_hash -- reverse-hash numbers into DNA strings
//

static PyObject * ktable_reverse_hash(PyObject * self, PyObject * args)
{
  khmer_KTableObject * me = (khmer_KTableObject *) self;
  khmer::KTable * ktable = me->ktable;

  unsigned int val;
  
  if (!PyArg_ParseTuple(args, "I", &val)) {
    return NULL;
  }

  return PyString_FromString(khmer::_revhash(val, ktable->ksize()).c_str());
}

//
// ktable_count -- count the given k-mer in the ktable
//

static PyObject * ktable_count(PyObject * self, PyObject * args)
{
  khmer_KTableObject * me = (khmer_KTableObject *) self;
  khmer::KTable * ktable = me->ktable;

  char * kmer;

  if (!PyArg_ParseTuple(args, "s", &kmer)) {
    return NULL;
  }

  if (strlen(kmer) != ktable->ksize()) {
    PyErr_SetString(PyExc_ValueError,
		    "k-mer length must be the same as the hashtable k-size");
    return NULL;
  }

  ktable->count(kmer);

  return PyInt_FromLong(1);
}

//
// ktable_consume -- count all of the k-mers in the given string
//

static PyObject * ktable_consume(PyObject * self, PyObject * args)
{
  khmer_KTableObject * me = (khmer_KTableObject *) self;
  khmer::KTable * ktable = me->ktable;

  char * long_str;

  if (!PyArg_ParseTuple(args, "s", &long_str)) {
    return NULL;
  }

  if (strlen(long_str) < ktable->ksize()) {
    PyErr_SetString(PyExc_ValueError,
		    "string length must >= the hashtable k-mer size");
    return NULL;
  }

  ktable->consume_string(long_str);

  unsigned int n_consumed = strlen(long_str) - ktable->ksize() + 1;

  return PyInt_FromLong(n_consumed);
}

static PyObject * ktable_get(PyObject * self, PyObject * args)
{
  khmer_KTableObject * me = (khmer_KTableObject *) self;
  khmer::KTable * ktable = me->ktable;

  PyObject * arg;

  if (!PyArg_ParseTuple(args, "O", &arg)) {
    return NULL;
  }

  unsigned long count = 0;

  if (PyInt_Check(arg)) {
    long pos = PyInt_AsLong(arg);
    count = ktable->get_count((unsigned int) pos);
  } else if (PyString_Check(arg)) {
    std::string s = PyString_AsString(arg);
    count = ktable->get_count(s.c_str());
  }

  return PyInt_FromLong(count);
}

static PyObject * ktable__getitem__(PyObject * self, Py_ssize_t index)
{
  khmer_KTableObject * me = (khmer_KTableObject *) self;
  khmer::KTable * ktable = me->ktable;

  khmer::ExactCounterType count = ktable->get_count(index);

  return PyInt_FromLong(count);
}

static PyObject * ktable_set(PyObject * self, PyObject * args)
{
  khmer_KTableObject * me = (khmer_KTableObject *) self;
  khmer::KTable * ktable = me->ktable;

  PyObject * arg;
  unsigned long count;

  if (!PyArg_ParseTuple(args, "OI", &arg, &count)) {
    return NULL;
  }

  if (PyInt_Check(arg)) {
    long pos = PyInt_AsLong(arg);
    ktable->set_count((unsigned int) pos, count);
  } else if (PyString_Check(arg)) {
    std::string s = PyString_AsString(arg);
    ktable->set_count(s.c_str(), count);
  }
  
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject * ktable_max_hash(PyObject * self, PyObject * args)
{
  khmer_KTableObject * me = (khmer_KTableObject *) self;
  khmer::KTable * ktable = me->ktable;

  if (!PyArg_ParseTuple(args, "")) {
    return NULL;
  }

  return PyInt_FromLong(ktable->max_hash());
}

static PyObject * ktable_n_entries(PyObject * self, PyObject * args)
{
  khmer_KTableObject * me = (khmer_KTableObject *) self;
  khmer::KTable * ktable = me->ktable;

  if (!PyArg_ParseTuple(args, "")) {
    return NULL;
  }

  return PyInt_FromLong(ktable->n_entries());
}

Py_ssize_t ktable__len__(PyObject * self)
{
  khmer_KTableObject * me = (khmer_KTableObject *) self;
  khmer::KTable * ktable = me->ktable;

  return ktable->n_entries();
}

static PyObject * ktable_ksize(PyObject * self, PyObject * args)
{
  khmer_KTableObject * me = (khmer_KTableObject *) self;
  khmer::KTable * ktable = me->ktable;

  if (!PyArg_ParseTuple(args, "")) {
    return NULL;
  }

  return PyInt_FromLong(ktable->ksize());
}

static PyObject * ktable_clear(PyObject * self, PyObject * args)
{
  khmer_KTableObject * me = (khmer_KTableObject *) self;
  khmer::KTable * ktable = me->ktable;

  if (!PyArg_ParseTuple(args, "")) {
    return NULL;
  }

  ktable->clear();

  Py_INCREF(Py_None);
  return Py_None;
}


// fwd decl --> defined below
static PyObject * ktable_update(PyObject * self, PyObject * args);
static PyObject * ktable_intersect(PyObject * self, PyObject * args);

static PyMethodDef khmer_ktable_methods[] = {
  { "forward_hash", ktable_forward_hash, METH_VARARGS, "Convert string to int" },
  { "forward_hash_no_rc", ktable_forward_hash_no_rc, METH_VARARGS, "Convert string to int, with no reverse complement handling" },
  { "reverse_hash", ktable_reverse_hash, METH_VARARGS, "Convert int to string" },
  { "count", ktable_count, METH_VARARGS, "Count the given kmer" },
  { "consume", ktable_consume, METH_VARARGS, "Count all k-mers in the given string" },
  { "get", ktable_get, METH_VARARGS, "Get the count for the given k-mer" },
  { "max_hash", ktable_max_hash, METH_VARARGS, "Get the maximum hash value"},
  { "n_entries", ktable_n_entries, METH_VARARGS, "Get the number of possible entries"},
  { "ksize", ktable_ksize, METH_VARARGS, "Get k"},
  { "set", ktable_set, METH_VARARGS, "Set counter to a value"},
  { "update", ktable_update, METH_VARARGS, "Combine another ktable with this one"},
  { "intersect", ktable_intersect, METH_VARARGS,
    "Create another ktable containing the intersection of two ktables:"
    "where both ktables have an entry, the counts will be summed."},
  { "clear", ktable_clear, METH_VARARGS, "Set all entries to 0." },

  {NULL, NULL, 0, NULL}           /* sentinel */
};

static PyObject *
khmer_ktable_getattr(PyObject * obj, char * name)
{
  return Py_FindMethod(khmer_ktable_methods, obj, name);
}

#define is_ktable_obj(v)  ((v)->ob_type == &khmer_KTableType)

static PySequenceMethods khmer_KTable_SequenceMethods = {
  ktable__len__,
  0,
  0,
  ktable__getitem__,
  0,
  0,
  0,
  0
};

static PyTypeObject khmer_KTableType = {
    PyObject_HEAD_INIT(NULL)
    0,
    "KTable", sizeof(khmer_KTableObject),
    0,
    khmer_ktable_dealloc,	/*tp_dealloc*/
    0,				/*tp_print*/
    khmer_ktable_getattr,	/*tp_getattr*/
    0,				/*tp_setattr*/
    0,				/*tp_compare*/
    0,				/*tp_repr*/
    0,				/*tp_as_number*/
    &khmer_KTable_SequenceMethods, /*tp_as_sequence*/
    0,				/*tp_as_mapping*/
    0,				/*tp_hash */
    0,				/*tp_call*/
    0,				/*tp_str*/
    0,				/*tp_getattro*/
    0,				/*tp_setattro*/
    0,				/*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT,		/*tp_flags*/
    "ktable object",           /* tp_doc */
};

//
// new_ktable
//

static PyObject* new_ktable(PyObject * self, PyObject * args)
{
  unsigned int size = 0;

  if (!PyArg_ParseTuple(args, "I", &size)) {
    return NULL;
  }

  khmer_KTableObject * ktable_obj = (khmer_KTableObject *) \
    PyObject_New(khmer_KTableObject, &khmer_KTableType);

  ktable_obj->ktable = new khmer::KTable(size);

  return (PyObject *) ktable_obj;
}

//
// khmer_ktable_dealloc -- clean up a table object.
//

static void khmer_ktable_dealloc(PyObject* self)
{
  khmer_KTableObject * obj = (khmer_KTableObject *) self;
  delete obj->ktable;
  obj->ktable = NULL;
  
  PyObject_Del((PyObject *) obj);
}

static PyObject * ktable_update(PyObject * self, PyObject * args)
{
  khmer_KTableObject * me = (khmer_KTableObject *) self;
  khmer::KTable * ktable = me->ktable;

  PyObject * other_o;

  PyArg_ParseTuple(args, "O", &other_o);

  assert(is_ktable_obj(other_o));

  khmer::KTable * other = ((khmer_KTableObject*) other_o)->ktable;

  ktable->update(*other);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject * ktable_intersect(PyObject * self, PyObject * args)
{
  khmer_KTableObject * me = (khmer_KTableObject *) self;
  khmer::KTable * ktable = me->ktable;

  PyObject * other_o;

  PyArg_ParseTuple(args, "O", &other_o);

  assert(is_ktable_obj(other_o));

  khmer::KTable * other = ((khmer_KTableObject*) other_o)->ktable;

  khmer::KTable * intersection = ktable->intersect(*other);

  khmer_KTableObject * ktable_obj = (khmer_KTableObject *) \
    PyObject_New(khmer_KTableObject, &khmer_KTableType);

  ktable_obj->ktable = intersection;

  return (PyObject *) ktable_obj;
}

PyObject * consume_genome(PyObject * self, PyObject * args)
{
  unsigned int size;
  char * genome;

  if (!PyArg_ParseTuple(args, "is", &size, &genome)) {
    return NULL;
  }

  khmer_KTableObject * ktable_obj = (khmer_KTableObject *) \
    PyObject_New(khmer_KTableObject, &khmer_KTableType);

  //  Py_BEGIN_ALLOW_THREADS
    {
      ktable_obj->ktable = new khmer::KTable(size);
      ktable_obj->ktable->consume_string(genome);
    }

    //  Py_END_ALLOW_THREADS

  return (PyObject *) ktable_obj;
}

/***********************************************************************/

//
// KHashtable object
//

typedef struct {
  PyObject_HEAD
  khmer::Hashtable * hashtable;
} khmer_KHashtableObject;

typedef struct {
  PyObject_HEAD
  khmer::ReadMaskTable * mask;
} khmer_ReadMaskObject;

#define is_readmask_obj(v)  ((v)->ob_type == &khmer_ReadMaskType)

static void khmer_readmask_dealloc(PyObject *);
static PyObject * khmer_readmask_getattr(PyObject *, char *);

static PyTypeObject khmer_ReadMaskType = {
    PyObject_HEAD_INIT(NULL)
    0,
    "ReadMask", sizeof(khmer_ReadMaskObject),
    0,
    khmer_readmask_dealloc,	/*tp_dealloc*/
    0,				/*tp_print*/
    khmer_readmask_getattr,	/*tp_getattr*/
    0,				/*tp_setattr*/
    0,				/*tp_compare*/
    0,				/*tp_repr*/
    0,				/*tp_as_number*/
    0,				/*tp_as_sequence*/
    0,				/*tp_as_mapping*/
    0,				/*tp_hash */
    0,				/*tp_call*/
    0,				/*tp_str*/
    0,				/*tp_getattro*/
    0,				/*tp_setattro*/
    0,				/*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT,		/*tp_flags*/
    "readmask object",           /* tp_doc */
};

typedef struct {
  PyObject_HEAD
  khmer::MinMaxTable * mmt;
} khmer_MinMaxObject;


#define is_minmax_obj(v)  ((v)->ob_type == &khmer_MinMaxType)

static void khmer_minmax_dealloc(PyObject* self);
static PyObject * khmer_minmax_getattr(PyObject *, char *);

static PyTypeObject khmer_MinMaxType = {
    PyObject_HEAD_INIT(NULL)
    0,
    "MinMax", sizeof(khmer_MinMaxObject),
    0,
    khmer_minmax_dealloc,	/*tp_dealloc*/
    0,				/*tp_print*/
    khmer_minmax_getattr,	/*tp_getattr*/
    0,				/*tp_setattr*/
    0,				/*tp_compare*/
    0,				/*tp_repr*/
    0,				/*tp_as_number*/
    0,				/*tp_as_sequence*/
    0,				/*tp_as_mapping*/
    0,				/*tp_hash */
    0,				/*tp_call*/
    0,				/*tp_str*/
    0,				/*tp_getattro*/
    0,				/*tp_setattro*/
    0,				/*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT,		/*tp_flags*/
    "minmax object",           /* tp_doc */
};



static void khmer_hashtable_dealloc(PyObject *);

static PyObject * hash_n_occupied(PyObject * self, PyObject * args)
{
  khmer_KHashtableObject * me = (khmer_KHashtableObject *) self;
  khmer::Hashtable * hashtable = me->hashtable;

  khmer::HashIntoType start = 0, stop = 0;

  if (!PyArg_ParseTuple(args, "|LL", &start, &stop)) {
    return NULL;
  }

  khmer::HashIntoType n = hashtable->n_occupied(start, stop);

  return PyInt_FromLong(n);
}

static PyObject * hash_n_entries(PyObject * self, PyObject * args)
{
  khmer_KHashtableObject * me = (khmer_KHashtableObject *) self;
  khmer::Hashtable * hashtable = me->hashtable;

  if (!PyArg_ParseTuple(args, "")) {
    return NULL;
  }

  return PyInt_FromLong(hashtable->n_entries());
}


static PyObject * hash_count(PyObject * self, PyObject * args)
{
  khmer_KHashtableObject * me = (khmer_KHashtableObject *) self;
  khmer::Hashtable * hashtable = me->hashtable;

  char * kmer;

  if (!PyArg_ParseTuple(args, "s", &kmer)) {
    return NULL;
  }

  if (strlen(kmer) != hashtable->ksize()) {
    PyErr_SetString(PyExc_ValueError,
		    "k-mer length must be the same as the hashtable k-size");
    return NULL;
  }

  hashtable->count(kmer);

  return PyInt_FromLong(1);
}

static PyObject * hash_filter_file_connected(PyObject * self, PyObject *args)
{
  khmer_KHashtableObject * me = (khmer_KHashtableObject *) self;
  khmer::Hashtable * hashtable = me->hashtable;

  char * est;
  char * readsfile;
  unsigned int total_reads;

  khmer::ReadMaskTable * readmask;

  if (!PyArg_ParseTuple(args, "ssi", &est, &readsfile, &total_reads)) {
    return NULL;
  }

  try {
    readmask = hashtable->filter_file_connected(est, readsfile, total_reads);
  } catch (_khmer_signal &e) {
    return NULL;
  }

  khmer_ReadMaskObject * readmask_obj = (khmer_ReadMaskObject *) \
    PyObject_New(khmer_ReadMaskObject, &khmer_ReadMaskType);

  readmask_obj->mask = readmask;

  return (PyObject *) readmask_obj;
}

static PyObject * hash_output_fasta_kmer_pos_freq(PyObject * self, PyObject *args)
{
  khmer_KHashtableObject * me = (khmer_KHashtableObject *) self;
  khmer::Hashtable * hashtable = me->hashtable;

  char * infile;
  char * outfile;

  if (!PyArg_ParseTuple(args, "ss", &infile, &outfile)) {
    return NULL;
  }

  hashtable->output_fasta_kmer_pos_freq(infile, outfile);

  return PyInt_FromLong(0);
}

static PyObject * hash_fasta_file_to_minmax(PyObject * self, PyObject *args)
{
  khmer_KHashtableObject * me = (khmer_KHashtableObject *) self;
  khmer::Hashtable * hashtable = me->hashtable;

  char * filename;
  unsigned int total_reads;
  PyObject * readmask_obj = NULL;
  PyObject * callback_obj = NULL;

  if (!PyArg_ParseTuple(args, "si|OO", &filename, &total_reads,
			&readmask_obj, &callback_obj)) {
    return NULL;
  }

  khmer::ReadMaskTable * readmask = NULL;
  if (readmask_obj && readmask_obj != Py_None) {
    if (!is_readmask_obj(readmask_obj)) {
      PyErr_SetString(PyExc_TypeError,
		      "third argument must be None or a readmask object");
      return NULL;
    }
    readmask = ((khmer_ReadMaskObject *) readmask_obj)->mask;
  }

  khmer::MinMaxTable * mmt;
  try {
    mmt = hashtable->fasta_file_to_minmax(filename, total_reads, readmask,
					  _report_fn, callback_obj);
  } catch (_khmer_signal &e) {
    return NULL;
  }
  
  khmer_MinMaxObject * minmax_obj = (khmer_MinMaxObject *) \
    PyObject_New(khmer_MinMaxObject, &khmer_MinMaxType);

  minmax_obj->mmt = mmt;

  return (PyObject *) minmax_obj;
}

static PyObject * hash_filter_fasta_file_limit_n(PyObject * self, PyObject *args)
{
  khmer_KHashtableObject * me = (khmer_KHashtableObject *) self;
  khmer::Hashtable * hashtable = me->hashtable;

  unsigned int threshold, n;
  char * filename;
  PyObject * o1 = NULL, * o2 = NULL;
  PyObject * callback_obj = NULL;

  if (!PyArg_ParseTuple(args, "sOii|OO", &filename, &o1, &threshold, &n, &o2, &callback_obj)) {
    return NULL;
  }

  if (!is_minmax_obj(o1)) {
    PyErr_SetString(PyExc_TypeError,
                    "third argument must be a minmax object");
    return NULL;
  }
  khmer::MinMaxTable * mmt = ((khmer_MinMaxObject *) o1)->mmt;

  khmer::ReadMaskTable * old_readmask = NULL;
  if (o2 && o2 != Py_None) {
    if (!is_readmask_obj(o2)) {
      PyErr_SetString(PyExc_TypeError,
                      "sixth must be None or a readmask object");
      return NULL;
    }
    old_readmask = ((khmer_ReadMaskObject *) o2)->mask;
  }

  khmer::ReadMaskTable * readmask;
  try {
    readmask = hashtable->filter_fasta_file_limit_n(filename, *mmt, threshold,
                                                n, old_readmask,
                                                _report_fn, callback_obj);
  } catch (_khmer_signal &e) {
    return NULL;
  }

  khmer_ReadMaskObject * readmask_obj = (khmer_ReadMaskObject *) \
    PyObject_New(khmer_ReadMaskObject, &khmer_ReadMaskType);

  readmask_obj->mask = readmask;

  return (PyObject *) readmask_obj;

}

static PyObject * hash_filter_fasta_file_any(PyObject * self, PyObject *args)
{
  khmer_KHashtableObject * me = (khmer_KHashtableObject *) self;
  khmer::Hashtable * hashtable = me->hashtable;

  unsigned int threshold;

  PyObject * o1 = NULL, * o2 = NULL;
  PyObject * callback_obj = NULL;

  if (!PyArg_ParseTuple(args, "Oi|OO", &o1, &threshold, &o2, &callback_obj)) {
    return NULL;
  }

  if (!is_minmax_obj(o1)) {
    PyErr_SetString(PyExc_TypeError,
		    "second argument must be a minmax object");
    return NULL;
  }
  khmer::MinMaxTable * mmt = ((khmer_MinMaxObject *) o1)->mmt;

  khmer::ReadMaskTable * old_readmask = NULL;
  if (o2 && o2 != Py_None) {
    if (!is_readmask_obj(o2)) {
      PyErr_SetString(PyExc_TypeError,
		      "fourth argument must be None or a readmask object");
      return NULL;
    }
    old_readmask = ((khmer_ReadMaskObject *) o2)->mask;
  }

  khmer::ReadMaskTable * readmask;
  try {
    readmask = hashtable->filter_fasta_file_any(*mmt, threshold,
						old_readmask,
						_report_fn, callback_obj);
  } catch (_khmer_signal &e) {
    return NULL;
  }

  khmer_ReadMaskObject * readmask_obj = (khmer_ReadMaskObject *) \
    PyObject_New(khmer_ReadMaskObject, &khmer_ReadMaskType);

  readmask_obj->mask = readmask;

  return (PyObject *) readmask_obj;
}

static PyObject * hash_filter_fasta_file_all(PyObject * self, PyObject *args)
{
  khmer_KHashtableObject * me = (khmer_KHashtableObject *) self;
  khmer::Hashtable * hashtable = me->hashtable;

  unsigned int threshold;

  PyObject * o1 = NULL, * o2 = NULL;
  PyObject * callback_obj = NULL;

  if (!PyArg_ParseTuple(args, "Oi|OO", &o1, &threshold, &o2, &callback_obj)) {
    return NULL;
  }

  if (!is_minmax_obj(o1)) {
    PyErr_SetString(PyExc_TypeError,
		    "second argument must be a minmax object");
    return NULL;
  }
  khmer::MinMaxTable * mmt = ((khmer_MinMaxObject *) o1)->mmt;

  khmer::ReadMaskTable * old_readmask = NULL;
  if (o2 && o2 != Py_None) {
    if (!is_readmask_obj(o2)) {
      PyErr_SetString(PyExc_TypeError,
		      "fourth argument must be None or a readmask object");
      return NULL;
    }
    old_readmask = ((khmer_ReadMaskObject *) o2)->mask;
  }

  khmer::ReadMaskTable * readmask;
  try {
    readmask = hashtable->filter_fasta_file_all(*mmt, threshold,
						old_readmask,
						_report_fn, callback_obj);
  } catch (_khmer_signal &e) {
    return NULL;
  }

  khmer_ReadMaskObject * readmask_obj = (khmer_ReadMaskObject *) \
    PyObject_New(khmer_ReadMaskObject, &khmer_ReadMaskType);

  readmask_obj->mask = readmask;

  return (PyObject *) readmask_obj;
}

static PyObject * hash_filter_fasta_file_run(PyObject * self, PyObject *args)
{
  khmer_KHashtableObject * me = (khmer_KHashtableObject *) self;
  khmer::Hashtable * hashtable = me->hashtable;

  char * filename;
  unsigned int threshold;
  unsigned int total_reads;
  unsigned int runlength;

  PyObject * o1 = NULL;
  PyObject * callback_obj = NULL;

  if (!PyArg_ParseTuple(args, "siii|OO", &filename, &total_reads, &threshold,
			&runlength, &o1, &callback_obj)) {
    return NULL;
  }

  khmer::ReadMaskTable * old_readmask = NULL;
  if (o1 && o1 != Py_None) {
    if (!is_readmask_obj(o1)) {
      PyErr_SetString(PyExc_TypeError,
		      "fifth argument must be None or a readmask object");
      return NULL;
    }
    old_readmask = ((khmer_ReadMaskObject *) o1)->mask;
  }

  khmer::ReadMaskTable * readmask;
  try {
    readmask = hashtable->filter_fasta_file_run(filename, total_reads,
						threshold, runlength,
						old_readmask,
						_report_fn, callback_obj);
  } catch (_khmer_signal &e) {
    return NULL;
  }

  khmer_ReadMaskObject * readmask_obj = (khmer_ReadMaskObject *) \
    PyObject_New(khmer_ReadMaskObject, &khmer_ReadMaskType);

  readmask_obj->mask = readmask;

  return (PyObject *) readmask_obj;
}

static PyObject * hash_consume_fasta(PyObject * self, PyObject * args)
{
  khmer_KHashtableObject * me = (khmer_KHashtableObject *) self;
  khmer::Hashtable * hashtable = me->hashtable;

  char * filename;
  PyObject * readmask_obj = NULL;
  PyObject * update_readmask_bool = NULL;
  khmer::HashIntoType lower_bound = 0, upper_bound = 0;
  PyObject * callback_obj = NULL;

  if (!PyArg_ParseTuple(args, "s|iiOOO", &filename, &lower_bound, &upper_bound,
			&readmask_obj, &update_readmask_bool,
			&callback_obj)) {
    return NULL;
  }

  // make sure update_readmask_bool is the right type of object
  if (update_readmask_bool && !PyBool_Check(update_readmask_bool)) {
    PyErr_SetString(PyExc_TypeError, "fifth argument must be True/False");
    return NULL;
  }

  // set C++ parameters accordingly
  bool update_readmask = false;
  khmer::ReadMaskTable * readmask = NULL;

  if (readmask_obj && readmask_obj != Py_None) {
    if (update_readmask_bool == Py_True) {
      update_readmask = true;
    }

    if (!is_readmask_obj(readmask_obj)) {
      PyErr_SetString(PyExc_TypeError,
		      "fourth argument must be None or a readmask object");
      return NULL;
    }
    
    readmask = ((khmer_ReadMaskObject *) readmask_obj)->mask;
  }

  // call the C++ function, and trap signals => Python

  unsigned long long n_consumed;
  unsigned int total_reads;

  try {
    hashtable->consume_fasta(filename, total_reads, n_consumed,
			     lower_bound, upper_bound, &readmask,
			     update_readmask, _report_fn, callback_obj);
  } catch (_khmer_signal &e) {
    return NULL;
  }

  // error checking -- this should still be null!
  if (!update_readmask && !readmask_obj) {
    assert(readmask == NULL);
  }

  return Py_BuildValue("iL", total_reads, n_consumed);
}

static PyObject * hash_consume_fasta_build_readmask(PyObject * self, PyObject * args)
{
  khmer_KHashtableObject * me = (khmer_KHashtableObject *) self;
  khmer::Hashtable * hashtable = me->hashtable;

  char * filename;
  khmer::HashIntoType lower_bound = 0, upper_bound = 0;
  PyObject * callback_obj = NULL;

  if (!PyArg_ParseTuple(args, "s|iiO", &filename, &lower_bound, &upper_bound,
			&callback_obj)) {
    return NULL;
  }

  khmer::ReadMaskTable * readmask = NULL;
  unsigned int total_reads;
  unsigned long long n_consumed;

  // this will allocate 'readmask' and fill it in.
  try {
    hashtable->consume_fasta(filename, total_reads, n_consumed,
			     lower_bound, upper_bound, &readmask, true,
			     _report_fn, callback_obj);
  } catch  (_khmer_signal &e) {
    return NULL;
  }

  if (!readmask) {
    PyErr_SetString(PyExc_RuntimeError,
		    "unexpected error in C++/consume_fasta; die die die.");
    return NULL;
  }

  khmer_ReadMaskObject * readmask_obj = (khmer_ReadMaskObject *) \
    PyObject_New(khmer_ReadMaskObject, &khmer_ReadMaskType);
  readmask_obj->mask = readmask;

  return Py_BuildValue("iLO", total_reads, n_consumed, readmask_obj);
}

static PyObject * hash_consume(PyObject * self, PyObject * args)
{
  khmer_KHashtableObject * me = (khmer_KHashtableObject *) self;
  khmer::Hashtable * hashtable = me->hashtable;

  char * long_str;
  khmer::HashIntoType lower_bound = 0, upper_bound = 0;

  if (!PyArg_ParseTuple(args, "s|ll", &long_str, &lower_bound, &upper_bound)) {
    return NULL;
  }
  
  if (strlen(long_str) < hashtable->ksize()) {
    PyErr_SetString(PyExc_ValueError,
		    "string length must >= the hashtable k-mer size");
    return NULL;
  }

  unsigned int n_consumed;
  n_consumed = hashtable->consume_string(long_str, lower_bound, upper_bound);

  return PyInt_FromLong(n_consumed);
}

static PyObject * hash_get_min_count(PyObject * self, PyObject * args)
{
  khmer_KHashtableObject * me = (khmer_KHashtableObject *) self;
  khmer::Hashtable * hashtable = me->hashtable;
  khmer::HashIntoType lower_bound = 0, upper_bound = 0;

  char * long_str;

  if (!PyArg_ParseTuple(args, "s|ll", &long_str, &lower_bound, &upper_bound)) {
    return NULL;
  }

  if (strlen(long_str) < hashtable->ksize()) {
    PyErr_SetString(PyExc_ValueError,
		    "string length must >= the hashtable k-mer size");
    return NULL;
  }

  khmer::BoundedCounterType c = hashtable->get_min_count(long_str,
							 lower_bound,
							 upper_bound);
  unsigned int N = c;

  return PyInt_FromLong(N);
}

static PyObject * hash_get_max_count(PyObject * self, PyObject * args)
{
  khmer_KHashtableObject * me = (khmer_KHashtableObject *) self;
  khmer::Hashtable * hashtable = me->hashtable;
  khmer::HashIntoType lower_bound = 0, upper_bound = 0;

  char * long_str;

  if (!PyArg_ParseTuple(args, "s|ll", &long_str, &lower_bound, &upper_bound)) {
    return NULL;
  }

  if (strlen(long_str) < hashtable->ksize()) {
    PyErr_SetString(PyExc_ValueError,
		    "string length must >= the hashtable k-mer size");
    return NULL;
  }

  khmer::BoundedCounterType c = hashtable->get_max_count(long_str,
							 lower_bound,
							 upper_bound);
  unsigned int N = c;

  return PyInt_FromLong(N);
}

static PyObject * hash_get(PyObject * self, PyObject * args)
{
  khmer_KHashtableObject * me = (khmer_KHashtableObject *) self;
  khmer::Hashtable * hashtable = me->hashtable;

  PyObject * arg;

  if (!PyArg_ParseTuple(args, "O", &arg)) {
    return NULL;
  }

  unsigned long count = 0;

  if (PyInt_Check(arg)) {
    long pos = PyInt_AsLong(arg);
    count = hashtable->get_count((unsigned int) pos);
  } else if (PyString_Check(arg)) {
    std::string s = PyString_AsString(arg);
    count = hashtable->get_count(s.c_str());
  }

  return PyInt_FromLong(count);
}

static PyObject * hash_abundance_distribution(PyObject * self, PyObject * args)
{
  khmer_KHashtableObject * me = (khmer_KHashtableObject *) self;
  khmer::Hashtable * hashtable = me->hashtable;

  if (!PyArg_ParseTuple(args, "")) {
    return NULL;
  }

  khmer::HashIntoType * dist;
  dist = hashtable->abundance_distribution();
  
  PyObject * x = PyList_New(256);
  for (int i = 0; i < 256; i++) {
    PyList_SET_ITEM(x, i, PyInt_FromLong(dist[i]));
  }

  delete dist;

  return x;
}

static PyObject * hash_fasta_count_kmers_by_position(PyObject * self, PyObject * args)
{
  khmer_KHashtableObject * me = (khmer_KHashtableObject *) self;
  khmer::Hashtable * hashtable = me->hashtable;

  char * inputfile;
  int max_read_len;
  int limit_by = 0;
  PyObject * readmask_obj = NULL;
  PyObject * callback_obj = NULL;

  if (!PyArg_ParseTuple(args, "sii|OO", &inputfile, &max_read_len, &limit_by,
			&readmask_obj, &callback_obj)) {
    return NULL;
  }

  khmer::ReadMaskTable * readmask = NULL;
  if (readmask_obj && readmask_obj != Py_None){
    if (!is_readmask_obj(readmask_obj)) {
      PyErr_SetString(PyExc_TypeError,
		      "fourth argument must be None or a readmask object");
      return NULL;
    }
    readmask = ((khmer_ReadMaskObject *) readmask_obj)->mask;
  }
    

  unsigned long long * counts;
  counts = hashtable->fasta_count_kmers_by_position(inputfile, max_read_len,
						    readmask, limit_by,
						    _report_fn, callback_obj);
					 
  PyObject * x = PyList_New(max_read_len);
  for (int i = 0; i < max_read_len; i++) {
    PyList_SET_ITEM(x, i, PyInt_FromLong(counts[i]));
  }

  delete counts;

  return x;
}

static PyObject * hash_fasta_dump_kmers_by_abundance(PyObject * self, PyObject * args)
{
  khmer_KHashtableObject * me = (khmer_KHashtableObject *) self;
  khmer::Hashtable * hashtable = me->hashtable;

  char * inputfile;
  int limit_by = 0;
  PyObject * readmask_obj = NULL;
  PyObject * callback_obj = NULL;

  if (!PyArg_ParseTuple(args, "si|OO", &inputfile, &limit_by,
			&readmask_obj, &callback_obj)) {
    return NULL;
  }

  khmer::ReadMaskTable * readmask = NULL;
  if (readmask_obj && readmask_obj != Py_None){
    if (!is_readmask_obj(readmask_obj)) {
      PyErr_SetString(PyExc_TypeError,
		      "third argument must be None or a readmask object");
      return NULL;
    }
    readmask = ((khmer_ReadMaskObject *) readmask_obj)->mask;
  }
    

  hashtable->fasta_dump_kmers_by_abundance(inputfile,
					   readmask, limit_by,
					   _report_fn, callback_obj);
					 

  Py_INCREF(Py_None);
  return Py_None;
}

// callback function to pass into dump function

void _dump_report_fn(const char * info, unsigned int count, void * data)
{
  // handle signals etc. (like CTRL-C)
  if (PyErr_CheckSignals() != 0) {
    throw _khmer_signal("PyErr_CheckSignals received a signal");
  }

  // if 'data' is set, it is a Python callable
  if (data) {
    PyObject * obj = (PyObject *) data;
    if (obj != Py_None) {
      PyObject * args = Py_BuildValue("si", info, count);

      PyObject * r = PyObject_Call(obj, args, NULL);
      Py_XDECREF(r);
      Py_DECREF(args);
    }
  }

  if (PyErr_Occurred()) {
    throw _khmer_signal("PyErr_Occurred is set");
  }

  // ...allow other Python threads to do stuff...
  Py_BEGIN_ALLOW_THREADS;
  Py_END_ALLOW_THREADS;
}


static PyObject * hash_dump_kmers_and_counts(PyObject * self, PyObject * args)
{
  khmer_KHashtableObject * me = (khmer_KHashtableObject *) self;
  khmer::Hashtable * hashtable = me->hashtable;

  PyObject * cb = NULL;
  if (!PyArg_ParseTuple(args, "|O", &cb)) {
    return NULL;
  }
  hashtable->dump_kmers_and_counts(_dump_report_fn, cb);
  
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject * hash_calc_connected_graph_size(PyObject * self, PyObject * args)
{
  khmer_KHashtableObject * me = (khmer_KHashtableObject *) self;
  khmer::Hashtable * hashtable = me->hashtable;

  char * _kmer;
  unsigned int max_size = 0;
  if (!PyArg_ParseTuple(args, "s|i", &_kmer, &max_size)) {
    return NULL;
  }

  unsigned long long size = 0;

  Py_BEGIN_ALLOW_THREADS
  khmer::SeenSet keeper;
  hashtable->calc_connected_graph_size(_kmer, size, keeper, max_size);
  Py_END_ALLOW_THREADS

  return PyInt_FromLong(size);
}

static PyObject * hash_trim_graphs(PyObject * self, PyObject * args)
{
  khmer_KHashtableObject * me = (khmer_KHashtableObject *) self;
  khmer::Hashtable * hashtable = me->hashtable;

  unsigned int threshold = 0;
  char * filename = NULL;
  char * outfile = NULL;
  PyObject * callback_obj = NULL;

  if (!PyArg_ParseTuple(args, "sis|O", &filename, &threshold, &outfile,
			&callback_obj)) {
    return NULL;
  }
  
  try {
    hashtable->trim_graphs(filename, outfile, threshold, _report_fn,
			   callback_obj);
  } catch (_khmer_signal &e) {
    return NULL;
  }
  
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject * hash_graphsize_distribution(PyObject * self, PyObject * args)
{
  khmer_KHashtableObject * me = (khmer_KHashtableObject *) self;
  khmer::Hashtable * hashtable = me->hashtable;

  unsigned int threshold = 0;
  if (!PyArg_ParseTuple(args, "i", &threshold)) {
    return NULL;
  }

  khmer::HashIntoType * p = hashtable->graphsize_distribution(threshold);
  
  PyObject * x = PyList_New(threshold);
  for (unsigned int i = 0; i < threshold; i++) {
    if (i > 0) { p[i] /= i; }
    PyList_SET_ITEM(x, i, PyInt_FromLong(p[i]));
  }

  delete p;

  return x;
}

static PyObject * hash_do_exact_partition(PyObject * self, PyObject * args)
{
  khmer_KHashtableObject * me = (khmer_KHashtableObject *) self;
  khmer::Hashtable * hashtable = me->hashtable;

  char * filename = NULL;
  PyObject * callback_obj = NULL;

  if (!PyArg_ParseTuple(args, "s|O", &filename, &callback_obj)) {
    return NULL;
  }

  unsigned int n_partitions = 0;
  try {
    n_partitions = hashtable->do_exact_partition(filename, _report_fn, callback_obj);
  } catch (_khmer_signal &e) {
    return NULL;
  }

  return PyInt_FromLong(n_partitions);
}

static PyObject * hash_do_truncated_partition(PyObject * self, PyObject * args)
{
  khmer_KHashtableObject * me = (khmer_KHashtableObject *) self;
  khmer::Hashtable * hashtable = me->hashtable;

  char * filename = NULL;
  char * output = NULL;
  unsigned int threshold = 0;
  PyObject * callback_obj = NULL;

  if (!PyArg_ParseTuple(args, "ssi|O", &filename, &output, &threshold,
			&callback_obj)) {
    return NULL;
  }

  unsigned int n_partitions = 0;
  try {
    n_partitions = hashtable->do_truncated_partition(filename, output,
						     threshold,
						     _report_fn, callback_obj);
  } catch (_khmer_signal &e) {
    return NULL;
  }

  return PyInt_FromLong(n_partitions);
}

static PyMethodDef khmer_hashtable_methods[] = {
  { "n_occupied", hash_n_occupied, METH_VARARGS, "Count the number of occupied bins" },
  { "n_entries", hash_n_entries, METH_VARARGS, "" },
  { "count", hash_count, METH_VARARGS, "Count the given kmer" },
  { "consume", hash_consume, METH_VARARGS, "Count all k-mers in the given string" },
  { "consume_fasta", hash_consume_fasta, METH_VARARGS, "Count all k-mers in a given file" },
  { "consume_fasta_build_readmask", hash_consume_fasta_build_readmask, METH_VARARGS, "Count all k-mers in a given file, creating a readmask object to mask off bad reads" },
  { "fasta_file_to_minmax", hash_fasta_file_to_minmax, METH_VARARGS, "" },
  { "filter_fasta_file_limit_n", hash_filter_fasta_file_limit_n, METH_VARARGS, "" },
  { "filter_fasta_file_any", hash_filter_fasta_file_any, METH_VARARGS, "" },
  { "filter_fasta_file_all", hash_filter_fasta_file_all, METH_VARARGS, "" },
  { "filter_fasta_file_run", hash_filter_fasta_file_run, METH_VARARGS, "" },
  { "output_fasta_kmer_pos_freq", hash_output_fasta_kmer_pos_freq, METH_VARARGS, "" },
  { "get", hash_get, METH_VARARGS, "Get the count for the given k-mer" },
  { "get_min_count", hash_get_min_count, METH_VARARGS, "Get the smallest count of all the k-mers in the string" },
  { "get_max_count", hash_get_max_count, METH_VARARGS, "Get the largest count of all the k-mers in the string" },
  { "abundance_distribution", hash_abundance_distribution, METH_VARARGS, "" },
  { "fasta_count_kmers_by_position", hash_fasta_count_kmers_by_position, METH_VARARGS, "" },
  { "fasta_dump_kmers_by_abundance", hash_fasta_dump_kmers_by_abundance, METH_VARARGS, "" },
  { "dump_kmers_and_counts", hash_dump_kmers_and_counts, METH_VARARGS, "" },
  { "calc_connected_graph_size", hash_calc_connected_graph_size, METH_VARARGS, "" },
  { "trim_graphs", hash_trim_graphs, METH_VARARGS, "" },
  { "graphsize_distribution", hash_graphsize_distribution, METH_VARARGS, "" },
  { "do_exact_partition", hash_do_exact_partition, METH_VARARGS, "" },
  { "do_truncated_partition", hash_do_truncated_partition, METH_VARARGS, "" },
  { "filter_file_connected", hash_filter_file_connected, METH_VARARGS, "" },
  {NULL, NULL, 0, NULL}           /* sentinel */
};

static PyObject *
khmer_hashtable_getattr(PyObject * obj, char * name)
{
  return Py_FindMethod(khmer_hashtable_methods, obj, name);
}

#define is_hashtable_obj(v)  ((v)->ob_type == &khmer_KHashtableType)

static PyTypeObject khmer_KHashtableType = {
    PyObject_HEAD_INIT(NULL)
    0,
    "KHashtable", sizeof(khmer_KHashtableObject),
    0,
    khmer_hashtable_dealloc,	/*tp_dealloc*/
    0,				/*tp_print*/
    khmer_hashtable_getattr,	/*tp_getattr*/
    0,				/*tp_setattr*/
    0,				/*tp_compare*/
    0,				/*tp_repr*/
    0,				/*tp_as_number*/
    0,				/*tp_as_sequence*/
    0,				/*tp_as_mapping*/
    0,				/*tp_hash */
    0,				/*tp_call*/
    0,				/*tp_str*/
    0,				/*tp_getattro*/
    0,				/*tp_setattro*/
    0,				/*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT,		/*tp_flags*/
    "hashtable object",           /* tp_doc */
};

//
// new_hashtable
//

static PyObject* new_hashtable(PyObject * self, PyObject * args)
{
  unsigned int k = 0;
  unsigned long long size = 0;

  if (!PyArg_ParseTuple(args, "IL", &k, &size)) {
    return NULL;
  }

  khmer_KHashtableObject * khashtable_obj = (khmer_KHashtableObject *) \
    PyObject_New(khmer_KHashtableObject, &khmer_KHashtableType);

  khashtable_obj->hashtable = new khmer::Hashtable(k, size);

  return (PyObject *) khashtable_obj;
}

//
// khmer_hashtable_dealloc -- clean up a table object.
//

static void khmer_hashtable_dealloc(PyObject* self)
{
  khmer_KHashtableObject * obj = (khmer_KHashtableObject *) self;
  delete obj->hashtable;
  obj->hashtable = NULL;
  
  PyObject_Del((PyObject *) obj);
}

//
// ReadMask object
//

static PyObject * readmask_get(PyObject * self, PyObject * args)
{
  khmer_ReadMaskObject * me = (khmer_ReadMaskObject *) self;
  khmer::ReadMaskTable * mask = me->mask;

  khmer::BoundedCounterType val;
  unsigned int index;

  if (!PyArg_ParseTuple(args, "I", &index)) {
    return NULL;
  }

  val = mask->get(index);

  return PyBool_FromLong(val);
}

static PyObject * readmask_n_kept(PyObject * self, PyObject * args)
{
  khmer_ReadMaskObject * me = (khmer_ReadMaskObject *) self;
  khmer::ReadMaskTable * mask = me->mask;

  unsigned int n_kept;

  if (!PyArg_ParseTuple(args, "")) {
    return NULL;
  }

  n_kept = mask->n_kept();

  return PyInt_FromLong(n_kept);
}

static PyObject * readmask_set(PyObject * self, PyObject * args)
{
  khmer_ReadMaskObject * me = (khmer_ReadMaskObject *) self;
  khmer::ReadMaskTable * mask = me->mask;

  unsigned int index;
  unsigned int setval;

  if (!PyArg_ParseTuple(args, "II", &index, &setval)) {
    return NULL;
  }

  if (setval) { mask->set(index, 1); }
  else { mask->set(index, 0); }
  
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject * readmask_and(PyObject * self, PyObject * args)
{
  khmer_ReadMaskObject * me = (khmer_ReadMaskObject *) self;
  khmer::ReadMaskTable * mask = me->mask;

  khmer::BoundedCounterType val;
  unsigned int index;
  unsigned int setval;

  if (!PyArg_ParseTuple(args, "II", &index, &setval)) {
    return NULL;
  }

  val = mask->get(index);

  if (setval && val) { val = 1; } 
  else { val = 0; }

  mask->set(index, val);
  
  return PyBool_FromLong(val);
}


static PyObject * readmask_merge(PyObject * self, PyObject * args)
{
  khmer_ReadMaskObject * me = (khmer_ReadMaskObject *) self;
  PyObject * other_py_obj;

  khmer::ReadMaskTable * mask = me->mask;

  if (!PyArg_ParseTuple(args, "O", &other_py_obj)) {
    return NULL;
  }

  khmer_ReadMaskObject * other = (khmer_ReadMaskObject *) other_py_obj;

  mask->merge(*(other->mask));

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject * readmask_save(PyObject * self, PyObject * args)
{
  khmer_ReadMaskObject * me = (khmer_ReadMaskObject *) self;

  khmer::ReadMaskTable * mask = me->mask;
  char * filename;

  if (!PyArg_ParseTuple(args, "s", &filename)) {
    return NULL;
  }

  mask->save(filename);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject * readmask_load(PyObject * self, PyObject * args)
{
  khmer_ReadMaskObject * me = (khmer_ReadMaskObject *) self;

  khmer::ReadMaskTable * mask = me->mask;
  char * filename;

  if (!PyArg_ParseTuple(args, "s", &filename)) {
    return NULL;
  }

  mask->load(filename);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject * readmask_invert(PyObject * self, PyObject * args)
{
  khmer_ReadMaskObject * me = (khmer_ReadMaskObject *) self;
  khmer::ReadMaskTable * mask = me->mask;

  if (!PyArg_ParseTuple(args, "")) {
    return NULL;
  }

  mask->invert();

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject * readmask_filter_fasta_file(PyObject * self, PyObject *args)
{
  char * inputfilename;
  char * outputfilename;
  khmer::ReadMaskTable * readmask = ((khmer_ReadMaskObject *) self)->mask;
  PyObject * callback_obj = NULL;

  if (!PyArg_ParseTuple(args, "ss|O", &inputfilename, &outputfilename,
			&callback_obj)) {
    return NULL;
  }

  unsigned int n_kept;
  try {
    n_kept = readmask->filter_fasta_file(inputfilename, outputfilename,
					 _report_fn, callback_obj);
  } catch (_khmer_signal &e) {
    return NULL;
  }

  return PyInt_FromLong(n_kept);
}


static PyObject * readmask_tablesize(PyObject * self, PyObject * args)
{
  khmer::ReadMaskTable * readmask = ((khmer_ReadMaskObject *) self)->mask;

  if (!PyArg_ParseTuple(args, "")) {
    return NULL;
  }

  return PyInt_FromLong(readmask->get_tablesize());
}


static PyMethodDef khmer_readmask_methods[] = {
  { "n_kept", readmask_n_kept, METH_VARARGS, "" },
  { "tablesize", readmask_tablesize, METH_VARARGS, "" },
  { "get", readmask_get, METH_VARARGS, "" },
  { "set", readmask_set, METH_VARARGS, "" },
  { "do_and", readmask_and, METH_VARARGS, "" },
  { "merge", readmask_merge, METH_VARARGS, "" },
  { "invert", readmask_invert, METH_VARARGS, "" },
  { "save", readmask_save, METH_VARARGS, "" },
  { "load", readmask_load, METH_VARARGS, "" },
  { "filter_fasta_file", readmask_filter_fasta_file, METH_VARARGS, "" },
  {NULL, NULL, 0, NULL}           /* sentinel */
};

static PyObject *
khmer_readmask_getattr(PyObject * obj, char * name)
{
  return Py_FindMethod(khmer_readmask_methods, obj, name);
}

//
// new_readmask
//

static PyObject* new_readmask(PyObject * self, PyObject * args)
{
  unsigned int size = 0;

  if (!PyArg_ParseTuple(args, "I", &size)) {
    return NULL;
  }

  khmer_ReadMaskObject * readmask_obj = (khmer_ReadMaskObject *) \
    PyObject_New(khmer_ReadMaskObject, &khmer_ReadMaskType);

  readmask_obj->mask = new khmer::ReadMaskTable(size);

  return (PyObject *) readmask_obj;
}

//
// khmer_readmask_dealloc -- clean up a readmask object.
//

static void khmer_readmask_dealloc(PyObject* self)
{
  khmer_ReadMaskObject * obj = (khmer_ReadMaskObject *) self;
  delete obj->mask;
  obj->mask = NULL;
  
  PyObject_Del((PyObject *) obj);
}

//
// MinMaxTable object
//

static PyObject * minmax_clear(PyObject * self, PyObject * args)
{
  khmer_MinMaxObject * me = (khmer_MinMaxObject *) self;
  khmer::MinMaxTable * mmt = me->mmt;

  unsigned int index;

  if (!PyArg_ParseTuple(args, "I", &index)) {
    return NULL;
  }

  mmt->clear(index);
  
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject * minmax_get_min(PyObject * self, PyObject * args)
{
  khmer_MinMaxObject * me = (khmer_MinMaxObject *) self;
  khmer::MinMaxTable * mmt = me->mmt;

  unsigned int index;

  if (!PyArg_ParseTuple(args, "I", &index)) {
    return NULL;
  }

  unsigned int val = mmt->get_min(index);

  return PyInt_FromLong(val);
}

static PyObject * minmax_get_max(PyObject * self, PyObject * args)
{
  khmer_MinMaxObject * me = (khmer_MinMaxObject *) self;
  khmer::MinMaxTable * mmt = me->mmt;

  unsigned int index;

  if (!PyArg_ParseTuple(args, "I", &index)) {
    return NULL;
  }

  unsigned int val = mmt->get_max(index);

  return PyInt_FromLong(val);
}

static PyObject * minmax_add_min(PyObject * self, PyObject * args)
{
  khmer_MinMaxObject * me = (khmer_MinMaxObject *) self;
  khmer::MinMaxTable * mmt = me->mmt;

  unsigned int index;
  unsigned int val;

  if (!PyArg_ParseTuple(args, "II", &index, &val)) {
    return NULL;
  }

  val = mmt->add_min(index, val);
  
  return PyInt_FromLong(val);
}

static PyObject * minmax_add_max(PyObject * self, PyObject * args)
{
  khmer_MinMaxObject * me = (khmer_MinMaxObject *) self;
  khmer::MinMaxTable * mmt = me->mmt;

  unsigned int index;
  unsigned int val;

  if (!PyArg_ParseTuple(args, "II", &index, &val)) {
    return NULL;
  }

  val = mmt->add_max(index, val);
  
  return PyInt_FromLong(val);
}

static PyObject * minmax_merge(PyObject * self, PyObject * args)
{
  khmer_MinMaxObject * me = (khmer_MinMaxObject *) self;
  PyObject * other_py_obj;

  khmer::MinMaxTable * mmt = me->mmt;

  if (!PyArg_ParseTuple(args, "O", &other_py_obj)) {
    return NULL;
  }

  khmer_MinMaxObject * other = (khmer_MinMaxObject *) other_py_obj;

  mmt->merge(*(other->mmt));

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject * minmax_save(PyObject * self, PyObject * args)
{
  khmer_MinMaxObject * me = (khmer_MinMaxObject *) self;

  khmer::MinMaxTable * mmt = me->mmt;
  char * filename;

  if (!PyArg_ParseTuple(args, "s", &filename)) {
    return NULL;
  }

  mmt->save(filename);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject * minmax_load(PyObject * self, PyObject * args)
{
  khmer_MinMaxObject * me = (khmer_MinMaxObject *) self;

  khmer::MinMaxTable * mmt = me->mmt;
  char * filename;

  if (!PyArg_ParseTuple(args, "s", &filename)) {
    return NULL;
  }

  mmt->load(filename);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject * minmax_tablesize(PyObject * self, PyObject * args)
{
  khmer_MinMaxObject * me = (khmer_MinMaxObject *) self;
  khmer::MinMaxTable * mmt = me->mmt;

  if (!PyArg_ParseTuple(args, "")) {
    return NULL;
  }

  return PyInt_FromLong(mmt->get_tablesize());
}

static PyMethodDef khmer_minmax_methods[] = {
  { "tablesize", minmax_tablesize, METH_VARARGS, "" },
  { "get_min", minmax_get_min, METH_VARARGS, "" },
  { "get_max", minmax_get_max, METH_VARARGS, "" },
  { "add_min", minmax_add_min, METH_VARARGS, "" },
  { "add_max", minmax_add_max, METH_VARARGS, "" },
  { "clear", minmax_clear, METH_VARARGS, "" },
  { "merge", minmax_merge, METH_VARARGS, "" },
  { "load", minmax_load, METH_VARARGS, "" },
  { "save", minmax_save, METH_VARARGS, "" },
  {NULL, NULL, 0, NULL}           /* sentinel */
};

static PyObject *
khmer_minmax_getattr(PyObject * obj, char * name)
{
  return Py_FindMethod(khmer_minmax_methods, obj, name);
}

//
// new_minmax
//

static PyObject* new_minmax(PyObject * self, PyObject * args)
{
  unsigned int size = 0;

  if (!PyArg_ParseTuple(args, "I", &size)) {
    return NULL;
  }

  khmer_MinMaxObject * minmax_obj = (khmer_MinMaxObject *) \
    PyObject_New(khmer_MinMaxObject, &khmer_MinMaxType);

  minmax_obj->mmt = new khmer::MinMaxTable(size);

  return (PyObject *) minmax_obj;
}

//
// khmer_minmax_dealloc -- clean up a minmax object.
//

static void khmer_minmax_dealloc(PyObject* self)
{
  khmer_MinMaxObject * obj = (khmer_MinMaxObject *) self;
  delete obj->mmt;
  obj->mmt = NULL;
  
  PyObject_Del((PyObject *) obj);
}


//////////////////////////////
// standalone functions

static PyObject * forward_hash(PyObject * self, PyObject * args)
{
  char * kmer;
  int ksize;

  if (!PyArg_ParseTuple(args, "si", &kmer, &ksize)) {
    return NULL;
  }

  if ((char)ksize != ksize) {
    PyErr_SetString(PyExc_ValueError, "k-mer size must be <= 255");
    return NULL;
  }

  return PyInt_FromLong(khmer::_hash(kmer, ksize));
}

static PyObject * forward_hash_no_rc(PyObject * self, PyObject * args)
{
  char * kmer;
  unsigned int ksize;

  if (!PyArg_ParseTuple(args, "si", &kmer, &ksize)) {
    return NULL;
  }

  if ((unsigned char)ksize != ksize) {
    PyErr_SetString(PyExc_ValueError, "k-mer size must be <= 255");
    return NULL;
  }

  if (strlen(kmer) != ksize) {
    PyErr_SetString(PyExc_ValueError,
		    "k-mer length must be the same as the hashtable k-size");
    return NULL;
  }

  return PyInt_FromLong(khmer::_hash_forward(kmer, ksize));
}

static PyObject * reverse_hash(PyObject * self, PyObject * args)
{
  khmer::HashIntoType val;
  int ksize;
  
  if (!PyArg_ParseTuple(args, "lI", &val, &ksize)) {
    return NULL;
  }

  if ((char)ksize != ksize) {
    PyErr_SetString(PyExc_ValueError, "k-mer size must be <= 255");
    return NULL;
  }

  return PyString_FromString(khmer::_revhash(val, ksize).c_str());
}

static PyObject * set_reporting_callback(PyObject * self, PyObject * args)
{
  PyObject * o;
  
  if (!PyArg_ParseTuple(args, "O", &o)) {
    return NULL;
  }

  Py_XDECREF(callback_obj);
  Py_INCREF(o);
  callback_obj = o;

  Py_INCREF(Py_None);
  return Py_None;
}

//
// Module machinery.
//

static PyMethodDef KhmerMethods[] = {
  { "new_ktable", new_ktable, METH_VARARGS, "Create an empty ktable" },
  { "new_hashtable", new_hashtable, METH_VARARGS, "Create an empty hashtable" },
  { "new_readmask", new_readmask, METH_VARARGS, "Create a new read mask table" },
  { "new_minmax", new_minmax, METH_VARARGS, "Create a new min/max value table" },
  { "consume_genome", consume_genome, METH_VARARGS, "Create a new ktable from a genome" },
  { "forward_hash", forward_hash, METH_VARARGS, "", },
  { "forward_hash_no_rc", forward_hash_no_rc, METH_VARARGS, "", },
  { "reverse_hash", reverse_hash, METH_VARARGS, "", },
  { "set_reporting_callback", set_reporting_callback, METH_VARARGS, "" },
  { NULL, NULL, 0, NULL }
};

DL_EXPORT(void) init_khmer(void)
{
  khmer_KTableType.ob_type = &PyType_Type;
  khmer_KHashtableType.ob_type = &PyType_Type;

  PyObject * m;
  m = Py_InitModule("_khmer", KhmerMethods);

  KhmerError = PyErr_NewException((char *)"_khmer.error", NULL, NULL);
  Py_INCREF(KhmerError);

  PyModule_AddObject(m, "error", KhmerError);
}
