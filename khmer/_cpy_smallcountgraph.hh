/*
This file is part of khmer, https://github.com/dib-lab/khmer/, and is
Copyright (C) 2016, The Regents of the University of California.

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
// KSmallCountgraph object
//

typedef struct {
    khmer_KHashgraph_Object khashgraph;
    SmallCountgraph * countgraph;
} khmer_KSmallCountgraph_Object;

static void khmer_smallcountgraph_dealloc(khmer_KSmallCountgraph_Object * obj);

static
PyObject *
smallcount_get_raw_tables(khmer_KSmallCountgraph_Object * self, PyObject * args)
{
    SmallCountgraph * countgraph = self->countgraph;

    khmer::Byte ** table_ptrs = countgraph->get_raw_tables();
    std::vector<uint64_t> sizes = countgraph->get_tablesizes();

    PyObject * raw_tables = PyList_New(sizes.size());
    for (unsigned int i=0; i<sizes.size(); ++i) {
        Py_buffer buffer;
        int res = PyBuffer_FillInfo(&buffer, NULL, table_ptrs[i],
                                    sizes[i] / 2 +1, 0,
                                    PyBUF_FULL_RO);
        if (res == -1) {
            return NULL;
        }
        PyObject * buf = PyMemoryView_FromBuffer(&buffer);
        if(!PyMemoryView_Check(buf)) {
            return NULL;
        }
        PyList_SET_ITEM(raw_tables, i, buf);
    }

    return raw_tables;
}

static PyMethodDef khmer_smallcountgraph_methods[] = {
    {
        "get_raw_tables",
        (PyCFunction)smallcount_get_raw_tables, METH_VARARGS,
        "Get a list of the raw storage tables as memoryview objects."
    },
    {NULL, NULL, 0, NULL}           /* sentinel */
};

static PyObject* khmer_smallcountgraph_new(PyTypeObject * type, PyObject * args,
        PyObject * kwds);

static PyTypeObject khmer_KSmallCountgraph_Type
CPYCHECKER_TYPE_OBJECT_FOR_TYPEDEF("khmer_KSmallCountgraph_Object")
= {
    PyVarObject_HEAD_INIT(NULL, 0)       /* init & ob_size */
    "_khmer.SmallCountgraph",                 /*tp_name*/
    sizeof(khmer_KSmallCountgraph_Object),  /*tp_basicsize*/
    0,                                   /*tp_itemsize*/
    (destructor)khmer_smallcountgraph_dealloc,  /*tp_dealloc*/
    0,                                   /*tp_print*/
    0,                                   /*tp_getattr*/
    0,                                   /*tp_setattr*/
    0,                                   /*tp_compare*/
    0,                                   /*tp_repr*/
    0,                                   /*tp_as_number*/
    0,                                   /*tp_as_sequence*/
    0,                                   /*tp_as_mapping*/
    0,                                   /*tp_hash */
    0,                                   /*tp_call*/
    0,                                   /*tp_str*/
    0,                                   /*tp_getattro*/
    0,                                   /*tp_setattro*/
    0,                                   /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,                  /*tp_flags*/
    "SmallCountgraph hash object",              /* tp_doc */
    0,                                   /* tp_traverse */
    0,                                   /* tp_clear */
    0,                                   /* tp_richcompare */
    0,                                   /* tp_weaklistoffset */
    0,                                   /* tp_iter */
    0,                                   /* tp_iternext */
    khmer_smallcountgraph_methods,              /* tp_methods */
    0,                                   /* tp_members */
    0,                                   /* tp_getset */
    0,                                   /* tp_base */
    0,                                   /* tp_dict */
    0,                                   /* tp_descr_get */
    0,                                   /* tp_descr_set */
    0,                                   /* tp_dictoffset */
    0,                                   /* tp_init */
    0,                                   /* tp_alloc */
    khmer_smallcountgraph_new,                /* tp_new */
};

//
// khmer_smallcountgraph_dealloc -- clean up a countgraph hash object.
//

static void khmer_smallcountgraph_dealloc(khmer_KSmallCountgraph_Object * obj)
{
    delete obj->countgraph;
    obj->countgraph = NULL;
    Py_TYPE(obj)->tp_free((PyObject*)obj);
}

//
// khmer_smallcountgraph_new
//

static PyObject* khmer_smallcountgraph_new(PyTypeObject * type, PyObject * args,
        PyObject * kwds)
{
    khmer_KSmallCountgraph_Object * self;

    self = (khmer_KSmallCountgraph_Object *)type->tp_alloc(type, 0);

    if (self != NULL) {
        WordLength k = 0;
        PyListObject * sizes_list_o = NULL;

        if (!PyArg_ParseTuple(args, "bO!", &k, &PyList_Type, &sizes_list_o)) {
            Py_DECREF(self);
            return NULL;
        }

        std::vector<uint64_t> sizes;
        if (!convert_Pytablesizes_to_vector(sizes_list_o, sizes)) {
            Py_DECREF(self);
            return NULL;
        }

        try {
            self->countgraph = new SmallCountgraph(k, sizes);
        } catch (std::bad_alloc &e) {
            Py_DECREF(self);
            return PyErr_NoMemory();
        }
        self->khashgraph.khashtable.hashtable =
            dynamic_cast<Hashtable*>(self->countgraph);
        self->khashgraph.hashgraph = dynamic_cast<Hashgraph*>(self->countgraph);
    }

    return (PyObject *) self;
}
