#include "khmer/_cpy_hashset.hh"
#include "oxli/hashgraph.hh"

using namespace oxli;

namespace khmer {


PyTypeObject _HashSet_iter_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)        /* init & ob_size */
    "_khmer.HashSet_iter",                /* tp_name */
    sizeof(_HashSet_iterobj),             /* tp_basicsize */
    0,                                    /* tp_itemsize */
    (destructor)_HashSet_iter_dealloc,    /* tp_dealloc */
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
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_ITER, /* tp_flags */
    "iterator object for HashSet objects.", /* tp_doc */
    0,                                    /* tp_traverse */
    0,                                    /* tp_clear */
    0,                                    /* tp_richcompare */
    0,                                    /* tp_weaklistoffset */
    _HashSet_iter,                        /* tp_iter */
    _HashSet_iternext,                    /* tp_iternext */
    0,                                    /* tp_methods */
    0,                                    /* tp_members */
    0,                                    /* tp_getset */
    0,                                    /* tp_base */
    0,                                    /* tp_dict */
    0,                                    /* tp_descr_get */
    0,                                    /* tp_descr_set */
    0,                                    /* tp_dictoffset */
    0,                                    /* tp_init */
    0,                                    /* tp_alloc */
    0,                                    /* tp_new */
};

PyTypeObject khmer_HashSet_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)        /* init & ob_size */
    "_khmer.HashSet",                     /* tp_name */
    sizeof(khmer_HashSet_Object),         /* tp_basicsize */
    0,                                    /* tp_itemsize */
    (destructor)khmer_HashSet_dealloc,    /* tp_dealloc */
    0,                                    /* tp_print */
    0,                                    /* tp_getattr */
    0,                                    /* tp_setattr */
    0,                                    /* tp_compare */
    0,                                    /* tp_repr */
    0,                                    /* tp_as_number */
    khmer_HashSet_seqmethods,             /* tp_as_sequence */
    0,                                    /* tp_as_mapping */
    0,                                    /* tp_hash */
    0,                                    /* tp_call */
    0,                                    /* tp_str */
    0,                                    /* tp_getattro */
    0,                                    /* tp_setattro */
    0,                                    /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_ITER, /* tp_flags */
    "Stores a set of hashed k-mers.",     /* tp_doc */
    0,                                    /* tp_traverse */
    0,                                    /* tp_clear */
    0,                                    /* tp_richcompare */
    0,                                    /* tp_weaklistoffset */
    khmer_HashSet_iter,                   /* tp_iter */
    0,                                    /* tp_iternext */
    khmer_HashSet_methods,                /* tp_methods */
    0,                                    /* tp_members */
    0,                                    /* tp_getset */
    0,                                    /* tp_base */
    0,                                    /* tp_dict */
    0,                                    /* tp_descr_get */
    0,                                    /* tp_descr_set */
    0,                                    /* tp_dictoffset */
    0,                                    /* tp_init */
    0,                                    /* tp_alloc */
    khmer_HashSet_new,                    /* tp_new */
};

PyMethodDef khmer_HashSet_methods[] = {
    {
        "add",
        (PyCFunction)hashset_add, METH_VARARGS,
        "Add element to the HashSet."
    },
    {
        "remove",
        (PyCFunction)hashset_remove, METH_VARARGS,
        "Remove an element from the HashSet."
    },
    {
        "update",
        (PyCFunction)hashset_update, METH_VARARGS,
        "Add a list of elements to the HashSet."
    },
    {NULL, NULL, 0, NULL}           /* sentinel */
};

PySequenceMethods khmer_HashSet_seqmethods[] = {
    (lenfunc)khmer_HashSet_len, /* sq_length */
    (binaryfunc)khmer_HashSet_concat,      /* sq_concat */
    0,                          /* sq_repeat */
    0,                          /* sq_item */
    0,                          /* sq_slice */
    0,                          /* sq_ass_item */
    0,                          /* sq_ass_slice */
    (objobjproc)khmer_HashSet_contains, /* sq_contains */
    (binaryfunc)khmer_HashSet_concat_inplace,      /* sq_inplace_concat */
    0                           /* sq_inplace_repeat */
};


void khmer_HashSet_dealloc(khmer_HashSet_Object * obj)
{
    delete obj->hashes;
    obj->hashes = NULL;
    obj->ksize = 0;
    Py_TYPE(obj)->tp_free((PyObject*)obj);
}

PyObject* khmer_HashSet_new(PyTypeObject * type, PyObject * args,
                                   PyObject * kwds)
{
    khmer_HashSet_Object * self;

    self = (khmer_HashSet_Object *)type->tp_alloc(type, 0);

    if (self != NULL) {
        PyObject * list_o = NULL;
        WordLength k;
        if (!PyArg_ParseTuple(args, "b|O!", &k, &PyList_Type, &list_o)) {
            Py_DECREF(self);
            return NULL;
        }

        try {
            self->hashes = new SeenSet;
            self->ksize = k;
        } catch (std::bad_alloc &e) {
            Py_DECREF(self);
            return PyErr_NoMemory();
        }

        if (list_o) {
            Py_ssize_t size = PyList_Size(list_o);
            for (Py_ssize_t i = 0; i < size; i++) {
                PyObject * item = PyList_GET_ITEM(list_o, i);
                HashIntoType h;

                if (!convert_PyObject_to_HashIntoType(item, h, self->ksize)) {
                    return NULL;
                }
                self->hashes->insert(h);
            }
        }
    }
    return (PyObject *) self;
}


void
_HashSet_iter_dealloc(_HashSet_iterobj * obj)
{
    delete obj->it;
    obj->it = NULL;
    Py_DECREF(obj->parent);
    Py_TYPE(obj)->tp_free((PyObject*)obj);
}

PyObject * _HashSet_iter(PyObject * self)
{
    Py_INCREF(self);
    return self;
}

PyObject * _HashSet_iternext(PyObject * self)
{
    _HashSet_iterobj * iter_obj = (_HashSet_iterobj *) self;
    SeenSet * hashes = iter_obj->parent->hashes;
    if (*iter_obj->it != hashes->end()) {
        PyObject * ret = nullptr;
        convert_HashIntoType_to_PyObject((**iter_obj->it), &ret);
        (*(iter_obj->it))++;
        return ret;
    }

    PyErr_SetString(PyExc_StopIteration, "end of HashSet");
    return NULL;
}


PyObject * khmer_HashSet_iter(PyObject * self)
{
    khmer_HashSet_Object * me = (khmer_HashSet_Object *) self;
    _HashSet_iterobj * iter_obj = (_HashSet_iterobj *)
                                  _HashSet_iter_Type.tp_alloc(&_HashSet_iter_Type, 0);
    if (iter_obj != NULL) {
        Py_INCREF(me);
        iter_obj->parent = me;

        iter_obj->it = new SeenSet::iterator;
        *iter_obj->it = me->hashes->begin();
    }
    return (PyObject *) iter_obj;
}

int khmer_HashSet_len(khmer_HashSet_Object * o)
{
    return (Py_ssize_t) o->hashes->size();
}

PyObject * khmer_HashSet_concat(khmer_HashSet_Object * o,
                                       khmer_HashSet_Object * o2)
{
    if (o->ksize != o2->ksize) {
        PyErr_SetString(PyExc_ValueError,
                        "cannot add HashSets with different ksizes");
        return NULL;
    }
    khmer_HashSet_Object * no = create_HashSet_Object(new SeenSet,
                                o->ksize);
    no->hashes->insert(o->hashes->begin(), o->hashes->end());
    no->hashes->insert(o2->hashes->begin(), o2->hashes->end());

    return (PyObject *) no;
}

PyObject * khmer_HashSet_concat_inplace(khmer_HashSet_Object * o,
        khmer_HashSet_Object * o2)
{
    if (o->ksize != o2->ksize) {
        PyErr_SetString(PyExc_ValueError,
                        "cannot add HashSets with different ksizes");
        return NULL;
    }
    o->hashes->insert(o2->hashes->begin(), o2->hashes->end());

    Py_INCREF(o);
    return (PyObject *) o;
}

int khmer_HashSet_contains(khmer_HashSet_Object * o, PyObject * val)
{
    HashIntoType v;

    if (convert_PyObject_to_HashIntoType(val, v, 0)) {
        if (set_contains(*o->hashes, v)) {
            return 1;
        }
    }
    return 0;
}

PyObject *
hashset_add(khmer_HashSet_Object * me, PyObject * args)
{
    PyObject * hash_obj;
    HashIntoType h;
    if (!PyArg_ParseTuple(args, "O", &hash_obj)) {
        return NULL;
    }

    if (!convert_PyObject_to_HashIntoType(hash_obj, h, 0)) {
        return NULL;
    }
    me->hashes->insert(h);

    Py_INCREF(Py_None);
    return Py_None;
}

PyObject *
hashset_remove(khmer_HashSet_Object * me, PyObject * args)
{
    PyObject * hash_obj;
    HashIntoType h;
    if (!PyArg_ParseTuple(args, "O", &hash_obj)) {
        return NULL;
    }

    if (!convert_PyObject_to_HashIntoType(hash_obj, h, 0)) {
        return NULL;
    }
    SeenSet::iterator it = me->hashes->find(h);
    if (it == me->hashes->end()) {
        PyErr_SetString(PyExc_ValueError, "h not in list");
        return NULL;
    }
    me->hashes->erase(it);

    Py_INCREF(Py_None);
    return Py_None;
}

PyObject *
hashset_update(khmer_HashSet_Object * me, PyObject * args)
{
    PyObject * obj;
    if (!PyArg_ParseTuple(args, "O", &obj)) {
        return NULL;
    }

    PyObject * iterator = PyObject_GetIter(obj);
    if (iterator == NULL) {
        return NULL;
    }
    PyObject * item = PyIter_Next(iterator);
    while(item) {
        HashIntoType h;

        if (!convert_PyObject_to_HashIntoType(item, h, 0)) {
            PyErr_SetString(PyExc_ValueError, "unknown item type for update");
            Py_DECREF(item);
            return NULL;
        }
        me->hashes->insert(h);

        Py_DECREF(item);
        item = PyIter_Next(iterator);
    }
    Py_DECREF(iterator);
    if (PyErr_Occurred()) {
        return NULL;
    }

    Py_INCREF(Py_None);
    return Py_None;
}


khmer_HashSet_Object * create_HashSet_Object(SeenSet * h, WordLength k)
{
    khmer_HashSet_Object * self;

    self = (khmer_HashSet_Object *)
           khmer_HashSet_Type.tp_alloc(&khmer_HashSet_Type, 0);
    if (self != NULL) {
        self->hashes = h;
        self->ksize = k;
    }
    return self;
}

}
