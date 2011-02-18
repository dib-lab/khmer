//
// A module for Python that exports khmer C++ library functions.
//

#include <iostream>

#include "Python.h"
#include "khmer.hh"
#include "ktable.hh"
#include "hashtable.hh"
#include "hashbits.hh"
#include "counting.hh"
#include "storage.hh"
#include "intertable.hh"

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

class _pre_partition_info {
public:
  khmer::HashIntoType kmer;
  khmer::SeenSet tagged_kmers;

  _pre_partition_info(khmer::HashIntoType _kmer) : kmer(_kmer) {};
};

// Python exception to raise
static PyObject *KhmerError;

// default callback obj;
static PyObject *_callback_obj = NULL;

// callback function to pass into C++ functions

void _report_fn(const char * info, void * data, unsigned int n_reads,
		unsigned long long other)
{
  // handle signals etc. (like CTRL-C)
  if (PyErr_CheckSignals() != 0) {
    throw _khmer_signal("PyErr_CheckSignals received a signal");
  }

  // set data to default?
  if (!data && _callback_obj) {
    data = _callback_obj;
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

  return PyLong_FromUnsignedLongLong(khmer::_hash(kmer, ktable->ksize()));
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

  return PyLong_FromUnsignedLongLong(khmer::_hash_forward(kmer, ktable->ksize()));
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
// KCountingHash object
//

typedef struct {
  PyObject_HEAD
  khmer::CountingHash * counting;
} khmer_KCountingHashObject;

typedef struct {
  PyObject_HEAD
  khmer::Hashbits * hashbits;
} khmer_KHashbitsObject;

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



static void khmer_counting_dealloc(PyObject *);
static void khmer_hashbits_dealloc(PyObject *);

static PyObject * hash_n_occupied(PyObject * self, PyObject * args)
{
  khmer_KCountingHashObject * me = (khmer_KCountingHashObject *) self;
  khmer::CountingHash * counting = me->counting;

  khmer::HashIntoType start = 0, stop = 0;

  if (!PyArg_ParseTuple(args, "|LL", &start, &stop)) {
    return NULL;
  }

  khmer::HashIntoType n = counting->n_occupied(start, stop);

  return PyInt_FromLong(n);
}

static PyObject * hash_n_entries(PyObject * self, PyObject * args)
{
  khmer_KCountingHashObject * me = (khmer_KCountingHashObject *) self;
  khmer::CountingHash * counting = me->counting;

  if (!PyArg_ParseTuple(args, "")) {
    return NULL;
  }

  return PyInt_FromLong(counting->n_entries());
}

static PyObject * hash_count(PyObject * self, PyObject * args)
{
  khmer_KCountingHashObject * me = (khmer_KCountingHashObject *) self;
  khmer::CountingHash * counting = me->counting;

  char * kmer;

  if (!PyArg_ParseTuple(args, "s", &kmer)) {
    return NULL;
  }

  if (strlen(kmer) != counting->ksize()) {
    PyErr_SetString(PyExc_ValueError,
		    "k-mer length must be the same as the hashtable k-size");
    return NULL;
  }

  counting->count(kmer);

  return PyInt_FromLong(1);
}

static PyObject * hash_output_fasta_kmer_pos_freq(PyObject * self, PyObject *args)
{
  khmer_KCountingHashObject * me = (khmer_KCountingHashObject *) self;
  khmer::CountingHash * counting = me->counting;

  char * infile;
  char * outfile;

  if (!PyArg_ParseTuple(args, "ss", &infile, &outfile)) {
    return NULL;
  }

  counting->output_fasta_kmer_pos_freq(infile, outfile);

  return PyInt_FromLong(0);
}

static PyObject * hash_fasta_file_to_minmax(PyObject * self, PyObject *args)
{
  khmer_KCountingHashObject * me = (khmer_KCountingHashObject *) self;
  khmer::CountingHash * counting = me->counting;

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
    mmt = counting->fasta_file_to_minmax(filename, total_reads, readmask,
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
  khmer_KCountingHashObject * me = (khmer_KCountingHashObject *) self;
  khmer::CountingHash * counting = me->counting;

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
    readmask = counting->filter_fasta_file_limit_n(filename, *mmt, threshold,
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
  khmer_KCountingHashObject * me = (khmer_KCountingHashObject *) self;
  khmer::CountingHash * counting = me->counting;

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
    readmask = counting->filter_fasta_file_any(*mmt, threshold,
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
  khmer_KCountingHashObject * me = (khmer_KCountingHashObject *) self;
  khmer::CountingHash * counting = me->counting;

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
    readmask = counting->filter_fasta_file_all(*mmt, threshold,
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
  khmer_KCountingHashObject * me = (khmer_KCountingHashObject *) self;
  khmer::CountingHash * counting = me->counting;

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
    readmask = counting->filter_fasta_file_run(filename, total_reads,
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
  khmer_KCountingHashObject * me = (khmer_KCountingHashObject *) self;
  khmer::CountingHash * counting = me->counting;

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
    counting->consume_fasta(filename, total_reads, n_consumed,
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
  khmer_KCountingHashObject * me = (khmer_KCountingHashObject *) self;
  khmer::CountingHash * counting = me->counting;

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
    counting->consume_fasta(filename, total_reads, n_consumed,
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
  khmer_KCountingHashObject * me = (khmer_KCountingHashObject *) self;
  khmer::CountingHash * counting = me->counting;

  char * long_str;
  khmer::HashIntoType lower_bound = 0, upper_bound = 0;

  if (!PyArg_ParseTuple(args, "s|ll", &long_str, &lower_bound, &upper_bound)) {
    return NULL;
  }
  
  if (strlen(long_str) < counting->ksize()) {
    PyErr_SetString(PyExc_ValueError,
		    "string length must >= the hashtable k-mer size");
    return NULL;
  }

  unsigned int n_consumed;
  n_consumed = counting->consume_string(long_str, lower_bound, upper_bound);

  return PyInt_FromLong(n_consumed);
}

static PyObject * hash_get_min_count(PyObject * self, PyObject * args)
{
  khmer_KCountingHashObject * me = (khmer_KCountingHashObject *) self;
  khmer::CountingHash * counting = me->counting;
  khmer::HashIntoType lower_bound = 0, upper_bound = 0;

  char * long_str;

  if (!PyArg_ParseTuple(args, "s|ll", &long_str, &lower_bound, &upper_bound)) {
    return NULL;
  }

  if (strlen(long_str) < counting->ksize()) {
    PyErr_SetString(PyExc_ValueError,
		    "string length must >= the hashtable k-mer size");
    return NULL;
  }

  khmer::BoundedCounterType c = counting->get_min_count(long_str,
							 lower_bound,
							 upper_bound);
  unsigned int N = c;

  return PyInt_FromLong(N);
}

static PyObject * hash_get_max_count(PyObject * self, PyObject * args)
{
  khmer_KCountingHashObject * me = (khmer_KCountingHashObject *) self;
  khmer::CountingHash * counting = me->counting;
  khmer::HashIntoType lower_bound = 0, upper_bound = 0;

  char * long_str;

  if (!PyArg_ParseTuple(args, "s|ll", &long_str, &lower_bound, &upper_bound)) {
    return NULL;
  }

  if (strlen(long_str) < counting->ksize()) {
    PyErr_SetString(PyExc_ValueError,
		    "string length must >= the hashtable k-mer size");
    return NULL;
  }

  khmer::BoundedCounterType c = counting->get_max_count(long_str,
							 lower_bound,
							 upper_bound);
  unsigned int N = c;

  return PyInt_FromLong(N);
}

static PyObject * hash_get_median_count(PyObject * self, PyObject * args)
{
  khmer_KCountingHashObject * me = (khmer_KCountingHashObject *) self;
  khmer::CountingHash * counting = me->counting;

  char * long_str;

  if (!PyArg_ParseTuple(args, "s", &long_str)) {
    return NULL;
  }

  if (strlen(long_str) < counting->ksize()) {
    PyErr_SetString(PyExc_ValueError,
		    "string length must >= the hashtable k-mer size");
    return NULL;
  }

  khmer::BoundedCounterType med = 0;
  float average = 0, stddev = 0;

  counting->get_median_count(long_str, med, average, stddev);

  return Py_BuildValue("iff", med, average, stddev);
}

static PyObject * hash_get_kmer_abund_mean(PyObject * self, PyObject * args)
{
  khmer_KCountingHashObject * me = (khmer_KCountingHashObject *) self;
  khmer::CountingHash * counting = me->counting;

  char * filename = NULL;

  if (!PyArg_ParseTuple(args, "s", &filename)) {
    return NULL;
  }

  unsigned long long total = 0;
  unsigned long long count = 0;
  float mean = 0.0;
  counting->get_kmer_abund_mean(filename, total, count, mean);

  return Py_BuildValue("LLf", total, count, mean);
}

static PyObject * hash_get_kmer_abund_abs_deviation(PyObject * self, PyObject * args)
{
  khmer_KCountingHashObject * me = (khmer_KCountingHashObject *) self;
  khmer::CountingHash * counting = me->counting;

  char * filename = NULL;
  float mean = 0.0;

  if (!PyArg_ParseTuple(args, "sf", &filename, &mean)) {
    return NULL;
  }

  float abs_dev = 0.0;
  counting->get_kmer_abund_abs_deviation(filename, mean, abs_dev);

  return Py_BuildValue("f", abs_dev);
}

static PyObject * hash_get(PyObject * self, PyObject * args)
{
  khmer_KCountingHashObject * me = (khmer_KCountingHashObject *) self;
  khmer::CountingHash * counting = me->counting;

  PyObject * arg;

  if (!PyArg_ParseTuple(args, "O", &arg)) {
    return NULL;
  }

  unsigned long count = 0;

  if (PyInt_Check(arg)) {
    long pos = PyInt_AsLong(arg);
    count = counting->get_count((unsigned int) pos);
  } else if (PyString_Check(arg)) {
    std::string s = PyString_AsString(arg);
    count = counting->get_count(s.c_str());
  }

  return PyInt_FromLong(count);
}

static PyObject * hash_abundance_distribution(PyObject * self, PyObject * args)
{
  khmer_KCountingHashObject * me = (khmer_KCountingHashObject *) self;
  khmer::CountingHash * counting = me->counting;

  char * filename = NULL;
  PyObject * callback_obj = NULL;
  if (!PyArg_ParseTuple(args, "s|O", &filename, &callback_obj)) {
    return NULL;
  }

  khmer::HashIntoType * dist;
  dist = counting->abundance_distribution(filename, _report_fn, callback_obj);
  
  PyObject * x = PyList_New(256);
  for (int i = 0; i < 256; i++) {
    PyList_SET_ITEM(x, i, PyInt_FromLong(dist[i]));
  }

  delete dist;

  return x;
}

static PyObject * hash_fasta_count_kmers_by_position(PyObject * self, PyObject * args)
{
  khmer_KCountingHashObject * me = (khmer_KCountingHashObject *) self;
  khmer::CountingHash * counting = me->counting;

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
  counts = counting->fasta_count_kmers_by_position(inputfile, max_read_len,
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
  khmer_KCountingHashObject * me = (khmer_KCountingHashObject *) self;
  khmer::CountingHash * counting = me->counting;

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
    

  counting->fasta_dump_kmers_by_abundance(inputfile,
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


static PyObject * hash_load(PyObject * self, PyObject * args)
{
  khmer_KCountingHashObject * me = (khmer_KCountingHashObject *) self;
  khmer::CountingHash * counting = me->counting;

  char * filename = NULL;

  if (!PyArg_ParseTuple(args, "s", &filename)) {
    return NULL;
  }

  counting->load(filename);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject * hash_save(PyObject * self, PyObject * args)
{
  khmer_KCountingHashObject * me = (khmer_KCountingHashObject *) self;
  khmer::CountingHash * counting = me->counting;

  char * filename = NULL;

  if (!PyArg_ParseTuple(args, "s", &filename)) {
    return NULL;
  }

  counting->save(filename);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyMethodDef khmer_counting_methods[] = {
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
  { "get_median_count", hash_get_median_count, METH_VARARGS, "Get the median, average, and stddev of the k-mer counts in the string" },
  { "abundance_distribution", hash_abundance_distribution, METH_VARARGS, "" },
  { "fasta_count_kmers_by_position", hash_fasta_count_kmers_by_position, METH_VARARGS, "" },
  { "fasta_dump_kmers_by_abundance", hash_fasta_dump_kmers_by_abundance, METH_VARARGS, "" },
  { "load", hash_load, METH_VARARGS, "" },
  { "save", hash_save, METH_VARARGS, "" },
  { "get_kmer_abund_abs_deviation", hash_get_kmer_abund_abs_deviation, METH_VARARGS, "" },
  { "get_kmer_abund_mean", hash_get_kmer_abund_mean, METH_VARARGS, "" },

  {NULL, NULL, 0, NULL}           /* sentinel */
};

static PyObject *
khmer_counting_getattr(PyObject * obj, char * name)
{
  return Py_FindMethod(khmer_counting_methods, obj, name);
}

#define is_counting_obj(v)  ((v)->ob_type == &khmer_KCountingHashType)

static PyTypeObject khmer_KCountingHashType = {
    PyObject_HEAD_INIT(NULL)
    0,
    "KCountingHash", sizeof(khmer_KCountingHashObject),
    0,
    khmer_counting_dealloc,	/*tp_dealloc*/
    0,				/*tp_print*/
    khmer_counting_getattr,	/*tp_getattr*/
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
    "counting hash object",           /* tp_doc */
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

  khmer_KCountingHashObject * kcounting_obj = (khmer_KCountingHashObject *) \
    PyObject_New(khmer_KCountingHashObject, &khmer_KCountingHashType);

  kcounting_obj->counting = new khmer::CountingHash(k, size);

  return (PyObject *) kcounting_obj;
}

//
// new_counting_hash
//

static PyObject* _new_counting_hash(PyObject * self, PyObject * args)
{
  unsigned int k = 0;
  PyObject* sizes_list_o = NULL;

  if (!PyArg_ParseTuple(args, "IO", &k, &sizes_list_o)) {
    return NULL;
  }

  std::vector<khmer::HashIntoType> sizes;
  for (unsigned int i = 0; i < PyObject_Length(sizes_list_o); i++) {
    PyObject * size_o = PyList_GET_ITEM(sizes_list_o, i);
    sizes.push_back(PyLong_AsLongLong(size_o));
  }

  khmer_KCountingHashObject * kcounting_obj = (khmer_KCountingHashObject *) \
    PyObject_New(khmer_KCountingHashObject, &khmer_KCountingHashType);

  kcounting_obj->counting = new khmer::CountingHash(k, sizes);

  return (PyObject *) kcounting_obj;
}

//
// hashbits stuff
//

static PyObject * hashbits_n_unique_kmers(PyObject * self, PyObject * args)
{
    khmer_KHashbitsObject * me = (khmer_KHashbitsObject *) self;
    khmer::Hashbits * hashbits = me->hashbits;
    
    khmer::HashIntoType start = 0, stop = 0;
    
    if (!PyArg_ParseTuple(args, "|LL", &start, &stop)) {
        return NULL;
    }
    
    khmer::HashIntoType n = hashbits->n_kmers(start, stop);
    
    return PyInt_FromLong(n);
}

static PyObject * hashbits_n_occupied(PyObject * self, PyObject * args)
{
  khmer_KHashbitsObject * me = (khmer_KHashbitsObject *) self;
  khmer::Hashbits * hashbits = me->hashbits;

  khmer::HashIntoType start = 0, stop = 0;

  if (!PyArg_ParseTuple(args, "|LL", &start, &stop)) {
    return NULL;
  }

  khmer::HashIntoType n = hashbits->n_occupied(start, stop);

  return PyInt_FromLong(n);
}

static PyObject * hashbits_n_tags(PyObject * self, PyObject * args)
{
  khmer_KHashbitsObject * me = (khmer_KHashbitsObject *) self;
  khmer::Hashbits * hashbits = me->hashbits;

  if (!PyArg_ParseTuple(args, "")) {
    return NULL;
  }

  return PyInt_FromLong(hashbits->n_tags());
}

static PyObject * hashbits_count(PyObject * self, PyObject * args)
{
  khmer_KHashbitsObject * me = (khmer_KHashbitsObject *) self;
  khmer::Hashbits * hashbits = me->hashbits;

  char * kmer;

  if (!PyArg_ParseTuple(args, "s", &kmer)) {
    return NULL;
  }

  if (strlen(kmer) != hashbits->ksize()) {
    PyErr_SetString(PyExc_ValueError,
		    "k-mer length must be the same as the hashbits k-size");
    return NULL;
  }

  hashbits->count(kmer);

  return PyInt_FromLong(1);
}

static PyObject * hashbits_filter_file_connected(PyObject * self, PyObject *args)
{
  khmer_KHashbitsObject * me = (khmer_KHashbitsObject *) self;
  khmer::Hashbits * hashbits = me->hashbits;

  char * est;
  char * readsfile;
  unsigned int total_reads;

  khmer::ReadMaskTable * readmask;

  if (!PyArg_ParseTuple(args, "ssi", &est, &readsfile, &total_reads)) {
    return NULL;
  }

  try {
    readmask = hashbits->filter_file_connected(est, readsfile, total_reads);
  } catch (_khmer_signal &e) {
    return NULL;
  }

  khmer_ReadMaskObject * readmask_obj = (khmer_ReadMaskObject *) \
    PyObject_New(khmer_ReadMaskObject, &khmer_ReadMaskType);

  readmask_obj->mask = readmask;

  return (PyObject *) readmask_obj;
}


static PyObject * hashbits_consume(PyObject * self, PyObject * args)
{
  khmer_KHashbitsObject * me = (khmer_KHashbitsObject *) self;
  khmer::Hashbits * hashbits = me->hashbits;

  char * long_str;
  khmer::HashIntoType lower_bound = 0, upper_bound = 0;

  if (!PyArg_ParseTuple(args, "s|ll", &long_str, &lower_bound, &upper_bound)) {
    return NULL;
  }
  
  if (strlen(long_str) < hashbits->ksize()) {
    PyErr_SetString(PyExc_ValueError,
		    "string length must >= the hashbits k-mer size");
    return NULL;
  }

  unsigned int n_consumed;
  n_consumed = hashbits->consume_string(long_str, lower_bound, upper_bound);

  return PyInt_FromLong(n_consumed);
}

static PyObject * hashbits_get(PyObject * self, PyObject * args)
{
  khmer_KHashbitsObject * me = (khmer_KHashbitsObject *) self;
  khmer::Hashbits * hashbits = me->hashbits;

  PyObject * arg;

  if (!PyArg_ParseTuple(args, "O", &arg)) {
    return NULL;
  }

  unsigned long count = 0;

  if (PyInt_Check(arg)) {
    long pos = PyInt_AsLong(arg);
    count = hashbits->get_count((unsigned int) pos);
  } else if (PyString_Check(arg)) {
    std::string s = PyString_AsString(arg);
    count = hashbits->get_count(s.c_str());
  }

  return PyInt_FromLong(count);
}

static PyObject * hashbits_calc_connected_graph_size(PyObject * self, PyObject * args)
{
  khmer_KHashbitsObject * me = (khmer_KHashbitsObject *) self;
  khmer::Hashbits * hashbits = me->hashbits;

  char * _kmer;
  unsigned int max_size = 0;
  PyObject * break_on_circum_o = NULL;
  if (!PyArg_ParseTuple(args, "s|iO", &_kmer, &max_size, &break_on_circum_o)) {
    return NULL;
  }

  bool break_on_circum = false;
  if (break_on_circum_o && !PyBool_Check(break_on_circum_o)) {
    PyErr_SetString(PyExc_TypeError, "third argument must be True/False");
    return NULL;
  }

  if (break_on_circum_o == Py_True) {
    break_on_circum = true;
  }

  unsigned long long size = 0;

  Py_BEGIN_ALLOW_THREADS
  khmer::SeenSet keeper;
  hashbits->calc_connected_graph_size(_kmer, size, keeper, max_size,
				      break_on_circum);
  Py_END_ALLOW_THREADS

  return PyInt_FromLong(size);
}

static PyObject * hashbits_trim_graphs(PyObject * self, PyObject * args)
{
  khmer_KHashbitsObject * me = (khmer_KHashbitsObject *) self;
  khmer::Hashbits * hashbits = me->hashbits;

  unsigned int threshold = 0;
  char * filename = NULL;
  char * outfile = NULL;
  PyObject * callback_obj = NULL;

  if (!PyArg_ParseTuple(args, "sis|O", &filename, &threshold, &outfile,
			&callback_obj)) {
    return NULL;
  }
  
  try {
    hashbits->trim_graphs(filename, outfile, threshold, _report_fn,
			   callback_obj);
  } catch (_khmer_signal &e) {
    return NULL;
  }
  
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject * hashbits_graphsize_distribution(PyObject * self, PyObject * args)
{
  khmer_KHashbitsObject * me = (khmer_KHashbitsObject *) self;
  khmer::Hashbits * hashbits = me->hashbits;

  unsigned int threshold = 0;
  if (!PyArg_ParseTuple(args, "i", &threshold)) {
    return NULL;
  }

  khmer::HashIntoType * p = hashbits->graphsize_distribution(threshold);
  
  PyObject * x = PyList_New(threshold);
  for (unsigned int i = 0; i < threshold; i++) {
    if (i > 0) { p[i] /= i; }
    PyList_SET_ITEM(x, i, PyInt_FromLong(p[i]));
  }

  delete p;

  return x;
}

static PyObject * hashbits_connectivity_distribution(PyObject * self, PyObject * args)
{
  khmer_KHashbitsObject * me = (khmer_KHashbitsObject *) self;
  khmer::Hashbits * hashbits = me->hashbits;

  char * filename = NULL;
  PyObject * callback_obj = NULL;

  if (!PyArg_ParseTuple(args, "s|O", &filename, &callback_obj)) {
    return NULL;
  }

  khmer::HashIntoType dist[9];
  try {
    hashbits->connectivity_distribution(filename, dist, _report_fn, callback_obj);
  } catch (_khmer_signal &e) {
    return NULL;
  }

  PyObject * x = PyList_New(9);
  for (unsigned int i = 0; i < 9; i++) {
    PyList_SET_ITEM(x, i, PyLong_FromUnsignedLongLong(dist[i]));
  }

  return x;
}

static PyObject * hashbits_kmer_degree(PyObject * self, PyObject * args)
{
  khmer_KHashbitsObject * me = (khmer_KHashbitsObject *) self;
  khmer::Hashbits * hashbits = me->hashbits;

  char * kmer_s = NULL;
  PyObject * callback_obj = NULL;

  if (!PyArg_ParseTuple(args, "s|O", &kmer_s, &callback_obj)) {
    return NULL;
  }

  return PyInt_FromLong(hashbits->kmer_degree(kmer_s));
}

void free_subset_partition_info(void * p)
{
  khmer::SubsetPartition * subset_p = (khmer::SubsetPartition *) p;
  delete subset_p;
}

static PyObject * hashbits_do_subset_partition(PyObject * self, PyObject * args)
{
  khmer_KHashbitsObject * me = (khmer_KHashbitsObject *) self;
  khmer::Hashbits * hashbits = me->hashbits;

  PyObject * callback_obj = NULL;
  khmer::HashIntoType start_kmer = 0, end_kmer = 0;

  if (!PyArg_ParseTuple(args, "K|KO", &start_kmer, &end_kmer,
			&callback_obj)) {
    return NULL;
  }

  khmer::SubsetPartition * subset_p = NULL;
  try {
    Py_BEGIN_ALLOW_THREADS
    subset_p = new khmer::SubsetPartition(hashbits);
    subset_p->do_partition(start_kmer, end_kmer, _report_fn, callback_obj);
    Py_END_ALLOW_THREADS
  } catch (_khmer_signal &e) {
    return NULL;
  }

  return PyCObject_FromVoidPtr(subset_p, free_subset_partition_info);
}

static PyObject * hashbits_merge_subset(PyObject * self, PyObject *args)
{
  khmer_KHashbitsObject * me = (khmer_KHashbitsObject *) self;
  khmer::Hashbits * hashbits = me->hashbits;

  PyObject * subset_obj;
  if (!PyArg_ParseTuple(args, "O", &subset_obj)) {
    return NULL;
  }

  if (!PyCObject_Check(subset_obj)) {
    return NULL;
  }

  khmer::SubsetPartition * subset_p;
  subset_p = (khmer::SubsetPartition *) PyCObject_AsVoidPtr(subset_obj);

  hashbits->partition->merge(subset_p);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject * hashbits_merge_from_disk(PyObject * self, PyObject *args)
{
  khmer_KHashbitsObject * me = (khmer_KHashbitsObject *) self;
  khmer::Hashbits * hashbits = me->hashbits;

  char * filename = NULL;
  if (!PyArg_ParseTuple(args, "s", &filename)) {
    return NULL;
  }

  hashbits->partition->merge_from_disk(filename);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject * hashbits_consume_fasta(PyObject * self, PyObject * args)
{
  khmer_KHashbitsObject * me = (khmer_KHashbitsObject *) self;
  khmer::Hashbits * hashbits = me->hashbits;

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
    hashbits->consume_fasta(filename, total_reads, n_consumed,
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

void sig(unsigned int total_reads, unsigned int n_consumed)
{
   std::cout << total_reads << " " << n_consumed << std::endl;
}

static PyObject * hashbits_consume_fasta_and_tag(PyObject * self, PyObject * args)
{
  khmer_KHashbitsObject * me = (khmer_KHashbitsObject *) self;
  khmer::Hashbits * hashbits = me->hashbits;

  char * filename;
  PyObject * callback_obj = NULL;

  if (!PyArg_ParseTuple(args, "s|O", &filename, &callback_obj)) {
    return NULL;
  }

  // call the C++ function, and trap signals => Python

  unsigned long long n_consumed;
  unsigned int total_reads;

  try {
    hashbits->consume_fasta_and_tag(filename, total_reads, n_consumed,
				     _report_fn, callback_obj);
  } catch (_khmer_signal &e) {
    return NULL;
  }

  return Py_BuildValue("iL", total_reads, n_consumed);
}

static PyObject * hashbits_thread_fasta(PyObject * self, PyObject * args)
{
  khmer_KHashbitsObject * me = (khmer_KHashbitsObject *) self;
  khmer::Hashbits * hashbits = me->hashbits;

  char * filename;
  PyObject * callback_obj = NULL;

  if (!PyArg_ParseTuple(args, "s|O", &filename, &callback_obj)) {
    return NULL;
  }

  // call the C++ function, and trap signals => Python

  unsigned long long n_consumed;
  unsigned int total_reads;

  try {
    hashbits->thread_fasta(filename, total_reads, n_consumed,
			   _report_fn, callback_obj);
  } catch (_khmer_signal &e) {
    return NULL;
  }

  return Py_BuildValue("iL", total_reads, n_consumed);
}

static PyObject * hashbits_consume_partitioned_fasta(PyObject * self, PyObject * args)
{
  khmer_KHashbitsObject * me = (khmer_KHashbitsObject *) self;
  khmer::Hashbits * hashbits = me->hashbits;

  char * filename;
  PyObject * callback_obj = NULL;

  if (!PyArg_ParseTuple(args, "s|O", &filename, &callback_obj)) {
    return NULL;
  }

  // call the C++ function, and trap signals => Python

  unsigned long long n_consumed;
  unsigned int total_reads;

  try {
    hashbits->consume_partitioned_fasta(filename, total_reads, n_consumed,
					 _report_fn, callback_obj);
  } catch (_khmer_signal &e) {
    return NULL;
  }

  return Py_BuildValue("iL", total_reads, n_consumed);
}

void free_pre_partition_info(void * p)
{
  _pre_partition_info * ppi = (_pre_partition_info *) p;
  delete ppi;
}

static PyObject * hashbits_find_all_tags(PyObject * self, PyObject *args)
{
  khmer_KHashbitsObject * me = (khmer_KHashbitsObject *) self;
  khmer::Hashbits * hashbits = me->hashbits;

  char * kmer_s = NULL;

  if (!PyArg_ParseTuple(args, "s", &kmer_s)) {
    return NULL;
  }

  if (strlen(kmer_s) < hashbits->ksize()) { // @@
    return NULL;
  }

  _pre_partition_info * ppi = NULL;

  Py_BEGIN_ALLOW_THREADS

  khmer::HashIntoType kmer_f, kmer_r;
  khmer::_hash(kmer_s, hashbits->ksize(), kmer_f, kmer_r);

  ppi = new _pre_partition_info(kmer_f);
  hashbits->partition->find_all_tags(kmer_f, kmer_r,
				     ppi->tagged_kmers,
				     false);
  hashbits->add_kmer_to_tags(kmer_f);

  Py_END_ALLOW_THREADS

  return PyCObject_FromVoidPtr(ppi, free_pre_partition_info);
}

static PyObject * hashbits_assign_partition_id(PyObject * self, PyObject *args)
{
  khmer_KHashbitsObject * me = (khmer_KHashbitsObject *) self;
  khmer::Hashbits * hashbits = me->hashbits;

  PyObject * ppi_obj;
  if (!PyArg_ParseTuple(args, "O", &ppi_obj)) {
    return NULL;
  }

  if (!PyCObject_Check(ppi_obj)) {
    return NULL;
  }

  _pre_partition_info * ppi;
  ppi = (_pre_partition_info *) PyCObject_AsVoidPtr(ppi_obj);
  
  khmer::PartitionID p;
  p = hashbits->partition->assign_partition_id(ppi->kmer,
					       ppi->tagged_kmers);

  return PyInt_FromLong(p);
}

static PyObject * hashbits_output_partitions(PyObject * self, PyObject * args)
{
  khmer_KHashbitsObject * me = (khmer_KHashbitsObject *) self;
  khmer::Hashbits * hashbits = me->hashbits;

  char * filename = NULL;
  char * output = NULL;
  PyObject * callback_obj = NULL;

  if (!PyArg_ParseTuple(args, "ss|O", &filename, &output, &callback_obj)) {
    return NULL;
  }

  unsigned int n_partitions = 0;

  try {
    khmer::SubsetPartition * subset_p = hashbits->partition;
    n_partitions = subset_p->output_partitioned_file(filename,
						     output,
						     false,
						     _report_fn,
						     callback_obj);
  } catch (_khmer_signal &e) {
    return NULL;
  }

  return PyInt_FromLong(n_partitions);
}

static PyObject * hashbits_filter_if_present(PyObject * self, PyObject * args)
{
  khmer_KHashbitsObject * me = (khmer_KHashbitsObject *) self;
  khmer::Hashbits * hashbits = me->hashbits;

  char * filename = NULL;
  char * output = NULL;
  PyObject * callback_obj = NULL;

  if (!PyArg_ParseTuple(args, "ss|O", &filename, &output, &callback_obj)) {
    return NULL;
  }

  try {
    hashbits->filter_if_present(filename, output, _report_fn, callback_obj);
  } catch (_khmer_signal &e) {
    return NULL;
  }
  
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject * hashbits_save_partitionmap(PyObject * self, PyObject * args)
{
  khmer_KHashbitsObject * me = (khmer_KHashbitsObject *) self;
  khmer::Hashbits * hashbits = me->hashbits;

  char * filename = NULL;

  if (!PyArg_ParseTuple(args, "s", &filename)) {
    return NULL;
  }

  hashbits->partition->save_partitionmap(filename);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject * hashbits_load_partitionmap(PyObject * self, PyObject * args)
{
  khmer_KHashbitsObject * me = (khmer_KHashbitsObject *) self;
  khmer::Hashbits * hashbits = me->hashbits;

  char * filename = NULL;

  if (!PyArg_ParseTuple(args, "s", &filename)) {
    return NULL;
  }

  hashbits->partition->load_partitionmap(filename);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject * hashbits__validate_partitionmap(PyObject * self, PyObject * args)
{
  khmer_KHashbitsObject * me = (khmer_KHashbitsObject *) self;
  khmer::Hashbits * hashbits = me->hashbits;

  if (!PyArg_ParseTuple(args, "")) {
    return NULL;
  }

  hashbits->partition->_validate_pmap();

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject * hashbits_count_partitions(PyObject * self, PyObject * args)
{
  khmer_KHashbitsObject * me = (khmer_KHashbitsObject *) self;
  khmer::Hashbits * hashbits = me->hashbits;

  if (!PyArg_ParseTuple(args, "")) {
    return NULL;
  }
  
  unsigned int n_partitions = 0, n_unassigned = 0;
  hashbits->partition->count_partitions(n_partitions, n_unassigned);

  return Py_BuildValue("ii", n_partitions, n_unassigned);
}

static PyObject * hashbits_subset_count_partitions(PyObject * self,
					       PyObject * args)
{
  PyObject * subset_obj = NULL;

  if (!PyArg_ParseTuple(args, "O", &subset_obj)) {
    return NULL;
  }
  
  khmer::SubsetPartition * subset_p;
  subset_p = (khmer::SubsetPartition *) PyCObject_AsVoidPtr(subset_obj);

  unsigned int n_partitions = 0, n_unassigned = 0;
  subset_p->count_partitions(n_partitions, n_unassigned);

  return Py_BuildValue("ii", n_partitions, n_unassigned);
}

static PyObject * hashbits_load(PyObject * self, PyObject * args)
{
  khmer_KHashbitsObject * me = (khmer_KHashbitsObject *) self;
  khmer::Hashbits * hashbits = me->hashbits;

  char * filename = NULL;

  if (!PyArg_ParseTuple(args, "s", &filename)) {
    return NULL;
  }

  hashbits->load(filename);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject * hashbits_save(PyObject * self, PyObject * args)
{
  khmer_KHashbitsObject * me = (khmer_KHashbitsObject *) self;
  khmer::Hashbits * hashbits = me->hashbits;

  char * filename = NULL;

  if (!PyArg_ParseTuple(args, "s", &filename)) {
    return NULL;
  }

  hashbits->save(filename);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject * hashbits_load_tagset(PyObject * self, PyObject * args)
{
  khmer_KHashbitsObject * me = (khmer_KHashbitsObject *) self;
  khmer::Hashbits * hashbits = me->hashbits;

  char * filename = NULL;

  if (!PyArg_ParseTuple(args, "s", &filename)) {
    return NULL;
  }

  hashbits->load_tagset(filename);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject * hashbits_save_tagset(PyObject * self, PyObject * args)
{
  khmer_KHashbitsObject * me = (khmer_KHashbitsObject *) self;
  khmer::Hashbits * hashbits = me->hashbits;

  char * filename = NULL;

  if (!PyArg_ParseTuple(args, "s", &filename)) {
    return NULL;
  }

  hashbits->save_tagset(filename);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject * hashbits_save_subset_partitionmap(PyObject * self, PyObject * args)
{
  char * filename = NULL;
  PyObject * subset_obj = NULL;

  if (!PyArg_ParseTuple(args, "Os", &subset_obj, &filename)) {
    return NULL;
  }

  khmer::SubsetPartition * subset_p;
  subset_p = (khmer::SubsetPartition *) PyCObject_AsVoidPtr(subset_obj);

  Py_BEGIN_ALLOW_THREADS

  subset_p->save_partitionmap(filename);

  Py_END_ALLOW_THREADS

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject * hashbits_load_subset_partitionmap(PyObject * self, PyObject * args)
{
  khmer_KHashbitsObject * me = (khmer_KHashbitsObject *) self;
  khmer::Hashbits * hashbits = me->hashbits;

  char * filename = NULL;

  if (!PyArg_ParseTuple(args, "s", &filename)) {
    return NULL;
  }

  khmer::SubsetPartition * subset_p;
  subset_p = new khmer::SubsetPartition(hashbits);

  Py_BEGIN_ALLOW_THREADS

  subset_p->load_partitionmap(filename);

  Py_END_ALLOW_THREADS

  return PyCObject_FromVoidPtr(subset_p, free_subset_partition_info);
}

static PyObject * hashbits_merge2_subset(PyObject * self, PyObject * args)
{
  // khmer_KHashbitsObject * me = (khmer_KHashbitsObject *) self;
  // khmer::Hashbits * hashbits = me->hashbits;

  PyObject * subset1_obj, * subset2_obj;

  if (!PyArg_ParseTuple(args, "OO", &subset1_obj, &subset2_obj)) {
    return NULL;
  }

  khmer::SubsetPartition * subset1_p;
  khmer::SubsetPartition * subset2_p;
  subset1_p = (khmer::SubsetPartition *) PyCObject_AsVoidPtr(subset1_obj);
  subset2_p = (khmer::SubsetPartition *) PyCObject_AsVoidPtr(subset2_obj);

  Py_BEGIN_ALLOW_THREADS

    subset1_p->merge(subset2_p);

  Py_END_ALLOW_THREADS

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject * hashbits_merge2_from_disk(PyObject * self, PyObject * args)
{
  // khmer_KHashbitsObject * me = (khmer_KHashbitsObject *) self;
  // khmer::Hashbits * hashbits = me->hashbits;

  PyObject * subset1_obj;
  char * filename = NULL;

  if (!PyArg_ParseTuple(args, "Os", &subset1_obj, &filename)) {
    return NULL;
  }

  khmer::SubsetPartition * subset1_p;
  subset1_p = (khmer::SubsetPartition *) PyCObject_AsVoidPtr(subset1_obj);

  Py_BEGIN_ALLOW_THREADS

    subset1_p->merge_from_disk(filename);

  Py_END_ALLOW_THREADS

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject * hashbits__validate_subset_partitionmap(PyObject * self, PyObject * args)
{
  PyObject * subset_obj = NULL;

  if (!PyArg_ParseTuple(args, "O", &subset_obj)) {
    return NULL;
  }

  khmer::SubsetPartition * subset_p;
  subset_p = (khmer::SubsetPartition *) PyCObject_AsVoidPtr(subset_obj);
  subset_p->_validate_pmap();

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject * hashbits_divide_tags_into_subsets(PyObject * self, PyObject * args)
{
  khmer_KHashbitsObject * me = (khmer_KHashbitsObject *) self;
  khmer::Hashbits * hashbits = me->hashbits;

  unsigned int subset_size = 0;

  if (!PyArg_ParseTuple(args, "i", &subset_size)) {
    return NULL;
  }

  khmer::SeenSet divvy;
  hashbits->divide_tags_into_subsets(subset_size, divvy);

  PyObject * x = PyList_New(divvy.size());
  unsigned int i = 0;
  for (khmer::SeenSet::const_iterator si = divvy.begin(); si != divvy.end();
       si++, i++) {
    PyList_SET_ITEM(x, i, PyLong_FromUnsignedLongLong(*si));
  }

  return x;
}

void free_tag_map(void * p)
{
  khmer::TagCountMap * m = (khmer::TagCountMap *) p;
  delete m;
}

static PyObject * hashbits_new_tagmap(PyObject * self, PyObject * args)
{
  // khmer_KHashbitsObject * me = (khmer_KHashbitsObject *) self;
  // khmer::Hashbits * hashbits = me->hashbits;

  if (!PyArg_ParseTuple(args, "")) {
    return NULL;
  }

  khmer::TagCountMap * tag_map = new khmer::TagCountMap();

  return PyCObject_FromVoidPtr(tag_map, free_tag_map);
}

static PyObject * hashbits_subset_maxify_partition_size(PyObject * self, PyObject * args)
{
  PyObject * subset_obj = NULL;
  PyObject * tagmap_obj = NULL;

  if (!PyArg_ParseTuple(args, "OO", &subset_obj, &tagmap_obj)) {
    return NULL;
  }

  khmer::SubsetPartition * subset_p;
  subset_p = (khmer::SubsetPartition *) PyCObject_AsVoidPtr(subset_obj);

  khmer::TagCountMap * tag_count;
  tag_count = (khmer::TagCountMap *) PyCObject_AsVoidPtr(tagmap_obj);

  subset_p->maxify_partition_size(*tag_count);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject * hashbits_discard_tags(PyObject * self, PyObject * args)
{
  khmer_KHashbitsObject * me = (khmer_KHashbitsObject *) self;
  khmer::Hashbits * hashbits = me->hashbits;

  PyObject * tagmap_obj = NULL;
  unsigned int threshold = 0;

  if (!PyArg_ParseTuple(args, "Oi", &tagmap_obj, &threshold)) {
    return NULL;
  }

  khmer::TagCountMap * tag_count;
  tag_count = (khmer::TagCountMap *) PyCObject_AsVoidPtr(tagmap_obj);

  hashbits->discard_tags(*tag_count, threshold);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject * hashbits_subset_filter_against_tags(PyObject * self, PyObject * args)
{
  PyObject * subset_obj = NULL;
  PyObject * tagmap_obj = NULL;

  if (!PyArg_ParseTuple(args, "OO", &subset_obj, &tagmap_obj)) {
    return NULL;
  }

  khmer::SubsetPartition * subset_p;
  subset_p = (khmer::SubsetPartition *) PyCObject_AsVoidPtr(subset_obj);

  khmer::TagCountMap * tag_count;
  tag_count = (khmer::TagCountMap *) PyCObject_AsVoidPtr(tagmap_obj);

  subset_p->filter_against_tags(*tag_count);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject * hashbits_set_partition_id(PyObject * self, PyObject * args)
{
  khmer_KHashbitsObject * me = (khmer_KHashbitsObject *) self;
  khmer::Hashbits * hashbits = me->hashbits;

  char * kmer = NULL;
  khmer::PartitionID p = 0;

  if (!PyArg_ParseTuple(args, "si", &kmer, &p)) {
    return NULL;
  }

  hashbits->partition->set_partition_id(kmer, p);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject * hashbits_join_partitions(PyObject * self, PyObject * args)
{
  khmer_KHashbitsObject * me = (khmer_KHashbitsObject *) self;
  khmer::Hashbits * hashbits = me->hashbits;

  khmer::PartitionID p1 = 0, p2 = 0;

  if (!PyArg_ParseTuple(args, "ii", &p1, &p2)) {
    return NULL;
  }

  p1 = hashbits->partition->join_partitions(p1, p2);

  return PyInt_FromLong(p1);
}

static PyObject * hashbits_get_partition_id(PyObject * self, PyObject * args)
{
  khmer_KHashbitsObject * me = (khmer_KHashbitsObject *) self;
  khmer::Hashbits * hashbits = me->hashbits;

  char * kmer = NULL;

  if (!PyArg_ParseTuple(args, "s", &kmer)) {
    return NULL;
  }

  khmer::PartitionID partition_id;
  partition_id = hashbits->partition->get_partition_id(kmer);

  return PyInt_FromLong(partition_id);
}

static PyObject * hashbits_count_kmers_within_radius(PyObject * self, PyObject * args)
{
  khmer_KHashbitsObject * me = (khmer_KHashbitsObject *) self;
  khmer::Hashbits * hashbits = me->hashbits;

  char * kmer = NULL;
  unsigned long radius = 0;
  unsigned long max_count = 0;

  if (!PyArg_ParseTuple(args, "sL|L", &kmer, &radius, &max_count)) {
    return NULL;
  }

  unsigned int n;

  Py_BEGIN_ALLOW_THREADS

  khmer::HashIntoType kmer_f, kmer_r;
  khmer::_hash(kmer, hashbits->ksize(), kmer_f, kmer_r);
  n = hashbits->count_kmers_within_radius(kmer_f, kmer_r, radius,
						       max_count);

  Py_END_ALLOW_THREADS

  return PyLong_FromUnsignedLong(n);
}

static PyObject * hashbits_count_kmers_on_radius(PyObject * self, PyObject * args)
{
  khmer_KHashbitsObject * me = (khmer_KHashbitsObject *) self;
  khmer::Hashbits * hashbits = me->hashbits;

  char * kmer = NULL;
  unsigned long radius = 0;
  unsigned long max_volume = 0;

  if (!PyArg_ParseTuple(args, "sL|L", &kmer, &radius, &max_volume)) {
    return NULL;
  }

  unsigned int n;

  Py_BEGIN_ALLOW_THREADS

  khmer::HashIntoType kmer_f, kmer_r;
  khmer::_hash(kmer, hashbits->ksize(), kmer_f, kmer_r);
  n = hashbits->count_kmers_on_radius(kmer_f, kmer_r, radius, max_volume);

  Py_END_ALLOW_THREADS

  return PyLong_FromUnsignedLong(n);
}

static PyObject * hashbits_find_radius_for_volume(PyObject * self, PyObject * args)
{
  khmer_KHashbitsObject * me = (khmer_KHashbitsObject *) self;
  khmer::Hashbits * hashbits = me->hashbits;

  char * kmer = NULL;
  unsigned long max_count = 0;
  unsigned long max_radius = 0;

  if (!PyArg_ParseTuple(args, "sLL", &kmer, &max_count, &max_radius)) {
    return NULL;
  }

  unsigned int n;

  Py_BEGIN_ALLOW_THREADS

  khmer::HashIntoType kmer_f, kmer_r;
  khmer::_hash(kmer, hashbits->ksize(), kmer_f, kmer_r);
  n = hashbits->find_radius_for_volume(kmer_f, kmer_r, max_count,
						    max_radius);

  Py_END_ALLOW_THREADS

  return PyLong_FromUnsignedLong(n);
}

static PyMethodDef khmer_hashbits_methods[] = {
  { "n_occupied", hashbits_n_occupied, METH_VARARGS, "Count the number of occupied bins" },
  { "n_unique_kmers", hashbits_n_unique_kmers,  METH_VARARGS, "Count the number of unique kmers" },
  { "count", hashbits_count, METH_VARARGS, "Count the given kmer" },
  { "consume", hashbits_consume, METH_VARARGS, "Count all k-mers in the given string" },
  { "get", hashbits_get, METH_VARARGS, "Get the count for the given k-mer" },
  { "calc_connected_graph_size", hashbits_calc_connected_graph_size, METH_VARARGS, "" },
  { "trim_graphs", hashbits_trim_graphs, METH_VARARGS, "" },
  { "graphsize_distribution", hashbits_graphsize_distribution, METH_VARARGS, "" },
  { "kmer_degree", hashbits_kmer_degree, METH_VARARGS, "" },
  { "connectivity_distribution", hashbits_connectivity_distribution, METH_VARARGS, "" },
  { "do_subset_partition", hashbits_do_subset_partition, METH_VARARGS, "" },
  { "filter_file_connected", hashbits_filter_file_connected, METH_VARARGS, "" },
  { "find_all_tags", hashbits_find_all_tags, METH_VARARGS, "" },
  { "assign_partition_id", hashbits_assign_partition_id, METH_VARARGS, "" },
  { "output_partitions", hashbits_output_partitions, METH_VARARGS, "" },
  { "filter_if_present", hashbits_filter_if_present, METH_VARARGS, "" },
  { "load", hashbits_load, METH_VARARGS, "" },
  { "save", hashbits_save, METH_VARARGS, "" },
  { "load_tagset", hashbits_load_tagset, METH_VARARGS, "" },
  { "save_tagset", hashbits_save_tagset, METH_VARARGS, "" },
  { "n_tags", hashbits_n_tags, METH_VARARGS, "" },
  { "divide_tags_into_subsets", hashbits_divide_tags_into_subsets, METH_VARARGS, "" },
  { "load_partitionmap", hashbits_load_partitionmap, METH_VARARGS, "" },
  { "save_partitionmap", hashbits_save_partitionmap, METH_VARARGS, "" },
  { "_validate_partitionmap", hashbits__validate_partitionmap, METH_VARARGS, "" },
  { "consume_fasta", hashbits_consume_fasta, METH_VARARGS, "Count all k-mers in a given file" },
  { "consume_fasta_and_tag", hashbits_consume_fasta_and_tag, METH_VARARGS, "Count all k-mers in a given file" },
  { "thread_fasta", hashbits_thread_fasta, METH_VARARGS, "Count all k-mers in a given file" },
  { "consume_partitioned_fasta", hashbits_consume_partitioned_fasta, METH_VARARGS, "Count all k-mers in a given file" },
  { "merge_subset", hashbits_merge_subset, METH_VARARGS, "" },
  { "merge_subset_from_disk", hashbits_merge_from_disk, METH_VARARGS, "" },
  { "count_partitions", hashbits_count_partitions, METH_VARARGS, "" },
  { "subset_count_partitions", hashbits_subset_count_partitions, METH_VARARGS, "" },
  { "save_subset_partitionmap", hashbits_save_subset_partitionmap, METH_VARARGS },
  { "load_subset_partitionmap", hashbits_load_subset_partitionmap, METH_VARARGS },
  { "merge2_subset", hashbits_merge2_subset, METH_VARARGS },
  { "merge2_subset_from_disk", hashbits_merge2_from_disk, METH_VARARGS },
  { "_validate_subset_partitionmap", hashbits__validate_subset_partitionmap, METH_VARARGS, "" },
  { "new_tagmap", hashbits_new_tagmap, METH_VARARGS, "" },
  { "subset_filter_against_tags", hashbits_subset_filter_against_tags, METH_VARARGS, "" },
  { "discard_tags", hashbits_discard_tags, METH_VARARGS, "" },
  { "subset_maxify_partition_size", hashbits_subset_maxify_partition_size, METH_VARARGS, "" },
  { "set_partition_id", hashbits_set_partition_id, METH_VARARGS, "" },
  { "join_partitions", hashbits_join_partitions, METH_VARARGS, "" },
  { "get_partition_id", hashbits_get_partition_id, METH_VARARGS, "" },
  { "count_kmers_within_radius", hashbits_count_kmers_within_radius, METH_VARARGS, "" },
  { "count_kmers_on_radius", hashbits_count_kmers_on_radius, METH_VARARGS, "" },
  { "find_radius_for_volume", hashbits_find_radius_for_volume, METH_VARARGS, "" },
  {NULL, NULL, 0, NULL}           /* sentinel */
};

static PyObject *
khmer_hashbits_getattr(PyObject * obj, char * name)
{
  return Py_FindMethod(khmer_hashbits_methods, obj, name);
}

#define is_hashbits_obj(v)  ((v)->ob_type == &khmer_KHashbitsType)

static PyTypeObject khmer_KHashbitsType = {
    PyObject_HEAD_INIT(NULL)
    0,
    "KHashbits", sizeof(khmer_KHashbitsObject),
    0,
    khmer_hashbits_dealloc,	/*tp_dealloc*/
    0,				/*tp_print*/
    khmer_hashbits_getattr,	/*tp_getattr*/
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
    "hashbits object",           /* tp_doc */
};

//
// new_hashbits
//

static PyObject* _new_hashbits(PyObject * self, PyObject * args)
{
  unsigned int k = 0;
  PyObject* sizes_list_o = NULL;

  if (!PyArg_ParseTuple(args, "IO", &k, &sizes_list_o)) {
    return NULL;
  }

  std::vector<khmer::HashIntoType> sizes;
  for (unsigned int i = 0; i < PyObject_Length(sizes_list_o); i++) {
    PyObject * size_o = PyList_GET_ITEM(sizes_list_o, i);
    sizes.push_back(PyLong_AsLongLong(size_o));
  }

  khmer_KHashbitsObject * khashbits_obj = (khmer_KHashbitsObject *) \
    PyObject_New(khmer_KHashbitsObject, &khmer_KHashbitsType);

  khashbits_obj->hashbits = new khmer::Hashbits(k, sizes);

  return (PyObject *) khashbits_obj;
}

//
// khmer_counting_dealloc -- clean up a counting hash object.
//

static void khmer_counting_dealloc(PyObject* self)
{
  khmer_KCountingHashObject * obj = (khmer_KCountingHashObject *) self;
  delete obj->counting;
  obj->counting = NULL;
  
  PyObject_Del((PyObject *) obj);
}

//
// khmer_hashbits_dealloc -- clean up a hashbits object.
//

static void khmer_hashbits_dealloc(PyObject* self)
{
  khmer_KHashbitsObject * obj = (khmer_KHashbitsObject *) self;
  delete obj->hashbits;
  obj->hashbits = NULL;
  
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

  return PyLong_FromUnsignedLongLong(khmer::_hash(kmer, ksize));
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

  return PyLong_FromUnsignedLongLong(khmer::_hash_forward(kmer, ksize));
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

  Py_XDECREF(_callback_obj);
  Py_INCREF(o);
  _callback_obj = o;

  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject * do_intersection_partition(PyObject * self, PyObject * args)
{
  int ksize;
  khmer::HashIntoType tablesize;
  char * infile = NULL;
  char * outfile = NULL;
  
  if (!PyArg_ParseTuple(args, "iLss", &ksize, &tablesize, &infile,
			&outfile)) {
    return NULL;
  }

  khmer::IntersectTable * iitable;
  iitable = new khmer::IntersectTable(ksize, tablesize);
  iitable->do_partition(infile, outfile);

  delete iitable;

  Py_INCREF(Py_None);
  return Py_None;
}


//
// Module machinery.
//

static PyMethodDef KhmerMethods[] = {
  { "new_ktable", new_ktable, METH_VARARGS, "Create an empty ktable" },
  { "new_hashtable", new_hashtable, METH_VARARGS, "Create an empty single-table counting hash" },
  { "_new_counting_hash", _new_counting_hash, METH_VARARGS, "Create an empty counting hash" },
  { "_new_hashbits", _new_hashbits, METH_VARARGS, "Create an empty hashbits table" },
  { "new_readmask", new_readmask, METH_VARARGS, "Create a new read mask table" },
  { "new_minmax", new_minmax, METH_VARARGS, "Create a new min/max value table" },
  { "consume_genome", consume_genome, METH_VARARGS, "Create a new ktable from a genome" },
  { "forward_hash", forward_hash, METH_VARARGS, "", },
  { "forward_hash_no_rc", forward_hash_no_rc, METH_VARARGS, "", },
  { "reverse_hash", reverse_hash, METH_VARARGS, "", },
  { "set_reporting_callback", set_reporting_callback, METH_VARARGS, "" },
  { "do_intersection_partition", do_intersection_partition, METH_VARARGS, "" },
  { NULL, NULL, 0, NULL }
};

DL_EXPORT(void) init_khmer(void)
{
  khmer_KTableType.ob_type = &PyType_Type;
  khmer_KCountingHashType.ob_type = &PyType_Type;

  PyObject * m;
  m = Py_InitModule("_khmer", KhmerMethods);

  KhmerError = PyErr_NewException((char *)"_khmer.error", NULL, NULL);
  Py_INCREF(KhmerError);

  PyModule_AddObject(m, "error", KhmerError);
}
