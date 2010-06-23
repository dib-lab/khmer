//
// A module for Python that exports two functions, 'find_iupac' and 'find_pwm'.
//

#include "Python.h"
#include "seqfuncs.hh"
#include "khmer.hh"

//
// Function necessary for Python loading:
//

extern "C" {
  void init_khmer();
}

// exception to raise
static PyObject *KhmerError;

/***********************************************************************/

//
// KTable object
//

typedef struct {
  PyObject_HEAD
  khmer::KTable * ktable;
} khmer_KTableObject;

static void khmer_ktable_dealloc(PyObject *);

static PyObject * forward_hash(PyObject * self, PyObject * args)
{
  khmer_KTableObject * me = (khmer_KTableObject *) self;
  khmer::KTable * ktable = me->ktable;

  char * kmer;

  if (!PyArg_ParseTuple(args, "s", &kmer)) {
    return NULL;
  }

  if (strlen(kmer) != ktable->ksize()) {
    // @CTB
    return NULL;
  }

  return PyInt_FromLong(khmer::_hash(kmer, ktable->ksize()));
}

static PyObject * reverse_hash(PyObject * self, PyObject * args)
{
  khmer_KTableObject * me = (khmer_KTableObject *) self;
  khmer::KTable * ktable = me->ktable;

  unsigned int val;
  
  if (!PyArg_ParseTuple(args, "I", &val)) {
    return NULL;
  }

  return PyString_FromString(khmer::_revhash(val, ktable->ksize()).c_str());
}

static PyObject * count(PyObject * self, PyObject * args)
{
  khmer_KTableObject * me = (khmer_KTableObject *) self;
  khmer::KTable * ktable = me->ktable;

  char * kmer;

  if (!PyArg_ParseTuple(args, "s", &kmer)) {
    return NULL;
  }

  if (strlen(kmer) != ktable->ksize()) {
    // @CTB
    return NULL;
  }

  ktable->count(kmer);

  return PyInt_FromLong(1);
}

static PyObject * consume(PyObject * self, PyObject * args)
{
  khmer_KTableObject * me = (khmer_KTableObject *) self;
  khmer::KTable * ktable = me->ktable;

  char * long_str;

  if (!PyArg_ParseTuple(args, "s", &long_str)) {
    return NULL;
  }

  if (strlen(long_str) <= ktable->ksize()) {
    // @CTB
    return NULL;
  }

  ktable->consume_string(long_str);

  unsigned int n_consumed = strlen(long_str) - ktable->ksize() + 1;

  return PyInt_FromLong(n_consumed);
}

static PyObject * get(PyObject * self, PyObject * args)
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

static PyObject * set(PyObject * self, PyObject * args)
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

static PyObject * max_hash(PyObject * self, PyObject * args)
{
  khmer_KTableObject * me = (khmer_KTableObject *) self;
  khmer::KTable * ktable = me->ktable;

  if (!PyArg_ParseTuple(args, "")) {
    return NULL;
  }

  return PyInt_FromLong(ktable->max_hash());
}

static PyObject * n_entries(PyObject * self, PyObject * args)
{
  khmer_KTableObject * me = (khmer_KTableObject *) self;
  khmer::KTable * ktable = me->ktable;

  if (!PyArg_ParseTuple(args, "")) {
    return NULL;
  }

  return PyInt_FromLong(ktable->n_entries());
}

static PyObject * ksize(PyObject * self, PyObject * args)
{
  khmer_KTableObject * me = (khmer_KTableObject *) self;
  khmer::KTable * ktable = me->ktable;

  if (!PyArg_ParseTuple(args, "")) {
    return NULL;
  }

  return PyInt_FromLong(ktable->ksize());
}

static PyObject * clear(PyObject * self, PyObject * args)
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
static PyObject * update(PyObject * self, PyObject * args);
static PyObject * intersect(PyObject * self, PyObject * args);

static PyMethodDef khmer_ktable_methods[] = {
  { "forward_hash", forward_hash, METH_VARARGS, "Convert string to int" },
  { "reverse_hash", reverse_hash, METH_VARARGS, "Convert int to string" },
  { "count", count, METH_VARARGS, "Count the given kmer" },
  { "consume", consume, METH_VARARGS, "Count all k-mers in the given string" },
  { "get", get, METH_VARARGS, "Get the count for the given k-mer" },
  { "max_hash", max_hash, METH_VARARGS, "Get the maximum hash value"},
  { "n_entries", n_entries, METH_VARARGS, "Get the number of possible entries"},
  { "ksize", ksize, METH_VARARGS, "Get k"},
  { "set", set, METH_VARARGS, "Set counter to a value"},
  { "update", update, METH_VARARGS, "Combine another ktable with this one"},
  { "intersect", intersect, METH_VARARGS,
    "Create another ktable containing the intersection of two ktables:"
    "where both ktables have an entry, the counts will be summed."},
  { "clear", clear, METH_VARARGS, "Set all entries to 0." },

  {NULL, NULL, 0, NULL}           /* sentinel */
};

static PyObject *
khmer_ktable_getattr(PyObject * obj, char * name)
{
  return Py_FindMethod(khmer_ktable_methods, obj, name);
}

#define is_ktable_obj(v)  ((v)->ob_type == &khmer_KTableType)

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
    0,				/*tp_as_sequence*/
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

static PyObject * update(PyObject * self, PyObject * args)
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

static PyObject * intersect(PyObject * self, PyObject * args)
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

static void khmer_hashtable_dealloc(PyObject *);

static PyObject * hash_count(PyObject * self, PyObject * args)
{
  khmer_KHashtableObject * me = (khmer_KHashtableObject *) self;
  khmer::Hashtable * hashtable = me->hashtable;

  char * kmer;

  if (!PyArg_ParseTuple(args, "s", &kmer)) {
    return NULL;
  }

  if (strlen(kmer) != hashtable->ksize()) {
    // @CTB
    return NULL;
  }

  hashtable->count(kmer);

  return PyInt_FromLong(1);
}

static PyObject * hash_filter_fasta_file(PyObject * self, PyObject *args)
{
  khmer_KHashtableObject * me = (khmer_KHashtableObject *) self;
  khmer::Hashtable * hashtable = me->hashtable;

  char * long_str1;
  char * long_str2;
  int i, j;

  if (!PyArg_ParseTuple(args, "ssii", &long_str1, &long_str2, &i, &j)) {
    return NULL;
  }

  hashtable->filter_fasta_file(long_str1, long_str2, i, j);

  return PyInt_FromLong(0);
}

static PyObject * hash_consume_fasta(PyObject * self, PyObject * args)
{
  khmer_KHashtableObject * me = (khmer_KHashtableObject *) self;
  khmer::Hashtable * hashtable = me->hashtable;

  char * long_str;

  if (!PyArg_ParseTuple(args, "s", &long_str)) {
    return NULL;
  }

  hashtable->consume_fasta(long_str);

  unsigned int n_consumed = strlen(long_str) - hashtable->ksize() + 1;

  return PyInt_FromLong(n_consumed);
}

static PyObject * hash_consume(PyObject * self, PyObject * args)
{
  khmer_KHashtableObject * me = (khmer_KHashtableObject *) self;
  khmer::Hashtable * hashtable = me->hashtable;

  char * long_str;

  if (!PyArg_ParseTuple(args, "s", &long_str)) {
    return NULL;
  }

  if (strlen(long_str) < hashtable->ksize()) {
    // @CTB
    return NULL;
  }

  hashtable->consume_string(long_str);

  unsigned int n_consumed = strlen(long_str) - hashtable->ksize() + 1;

  return PyInt_FromLong(n_consumed);
}

static PyObject * hash_get_min_count(PyObject * self, PyObject * args)
{
  khmer_KHashtableObject * me = (khmer_KHashtableObject *) self;
  khmer::Hashtable * hashtable = me->hashtable;

  char * long_str;

  if (!PyArg_ParseTuple(args, "s", &long_str)) {
    return NULL;
  }

  if (strlen(long_str) < hashtable->ksize()) {
    // @CTB
    return NULL;
  }

  khmer::BoundedCounterType c = hashtable->get_min_count(long_str);
  unsigned int N = c;

  return PyInt_FromLong(N);
}

static PyObject * hash_get_max_count(PyObject * self, PyObject * args)
{
  khmer_KHashtableObject * me = (khmer_KHashtableObject *) self;
  khmer::Hashtable * hashtable = me->hashtable;

  char * long_str;

  if (!PyArg_ParseTuple(args, "s", &long_str)) {
    return NULL;
  }

  if (strlen(long_str) < hashtable->ksize()) {
    // @CTB
    return NULL;
  }

  khmer::BoundedCounterType c = hashtable->get_max_count(long_str);
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

static PyMethodDef khmer_hashtable_methods[] = {
  { "count", hash_count, METH_VARARGS, "Count the given kmer" },
  { "consume", hash_consume, METH_VARARGS, "Count all k-mers in the given string" },
  { "consume_fasta", hash_consume_fasta, METH_VARARGS, "Count all k-mers in a given file" },
  { "filter_fasta_file", hash_filter_fasta_file, METH_VARARGS, "Filter and trim reads file"},
  { "get", hash_get, METH_VARARGS, "Get the count for the given k-mer" },
  { "get_min_count", hash_get_min_count, METH_VARARGS, "Get the smallest count of all the k-mers in the string" },
  { "get_max_count", hash_get_max_count, METH_VARARGS, "Get the largest count of all the k-mers in the string" },

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
// Module machinery.
//

static PyMethodDef KhmerMethods[] = {
  { "new_ktable", new_ktable, METH_VARARGS, "Create an empty ktable" },
  { "new_hashtable", new_hashtable, METH_VARARGS, "Create an empty hashtable" },
  { "consume_genome", consume_genome, METH_VARARGS, "Create a new ktable from a genome" },
  { NULL, NULL, 0, NULL }
};

DL_EXPORT(void) init_khmer(void)
{
  khmer_KTableType.ob_type = &PyType_Type;
  khmer_KHashtableType.ob_type = &PyType_Type;

  PyObject * m;
  m = Py_InitModule("_khmer", KhmerMethods);

  KhmerError = PyErr_NewException("_khmer.error", NULL, NULL);
  Py_INCREF(KhmerError);

  PyModule_AddObject(m, "error", KhmerError);
}
