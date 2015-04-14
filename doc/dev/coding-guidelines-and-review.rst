Coding guidelines and code review checklist
===========================================

This document is for anyone who want to contribute code to the khmer
project, and describes our coding standards and code review checklist.

----

Coding standards
----------------

All plain-text files should have line widths of 80 characters or less unless
that is not supported for the particular file format.

For C++, we use `Todd Hoff's coding standard
<http://www.possibility.com/Cpp/CppCodingStandard.html>`__, and
`astyle -A10 / "One True Brace Style"
<http://astyle.sourceforge.net/astyle.html>`__ indentation and
bracing.  Note: @CTB needs Emacs settings that work for this.

Vim users may want to set the ARTISTIC_STYLE_OPTIONS shell variable to "-A10
--max-code-length=80" and run ```:%!astyle``` to reformat. The four space
indentation can be set with::

	set expandtab
	set shiftwidth=4
	set softtabstop=4

For Python, `PEP 8 <http://www.python.org/dev/peps/pep-0008/>`__ is our
standard. The ```pep8``` and ```autopep8``` Makefile targets are helpful. 

Code, scripts, and documentation must have its spelling checked. Vim users can
run::

        :setlocal spell spelllang=en_us

Use `]s` and `[s` to navigate between misspellings and `z=` to suggest a
correctly spelled word. `zg` will add a word as a good word.

GNU `aspell` can also be used to check the spelling in a single file::

        aspell check --mode ccpp $filename

Code Review
-----------

Please read `11 Best Practices for Peer Code Review
<http://smartbear.com/SmartBear/media/pdfs/WP-CC-11-Best-Practices-of-Peer-Code-Review.pdf>`__.

See also `Code reviews: the lab meeting for code
<http://fperez.org/py4science/code_reviews.html>`__ and
`the PyCogent coding guidelines
<http://pycogent.org/coding_guidelines.html>`__.

Checklist
---------

Copy and paste the following into a pull request comment when it is
ready for review::
   
   - [ ] Is it mergeable?
   - [ ] Did it pass the tests?
   - [ ] If it introduces new functionality in scripts/ is it tested?
     Check for code coverage with `make clean diff-cover`
   - [ ] Is it well formatted? Look at `make pep8`, `make diff_pylint_report`,
     `make cppcheck`, and `make doc` output. Use `make format` and manual
     fixing as needed.
   - [ ] Did it change the command-line interface? Only additions are allowed
     without a major version increment. Changing file formats also requires a
     major version number increment.
   - [ ] Is it documented in the ChangeLog?
     http://en.wikipedia.org/wiki/Changelog#Format
   - [ ] Was a spellchecker run on the source code and documentation after
     changes were made?
   - [ ] Is the Copyright year up to date?

**Note** that after you submit the comment you can check and uncheck
the individual boxes on the formatted comment; no need to put x or y
in the middle.

CPython Checklist
-----------------

Here's a checklist for new CPython types with future-proofing for Python 3::

   - [ ] the CPython object name is of the form `khmer_${OBJECTNAME}_Object`
   - [ ] Named struct with `PyObject_HEAD` macro
   - [ ] `static PyTypeObject khmer_${OBJECTNAME}_Type` with the following
     entries
      - [ ] `PyVarObject_HEAD_INIT(NULL, 0)` as the object init (this includes
        the `ob_size` field).
      - [ ] all fields should have their name in a comment for readability
      - [ ] The `tp_name` filed is a dotted name with both the module name and
        the name of the type within the module. Example: `khmer.ReadAligner`
      - [ ] Deallocator defined and cast to `(destructor)` in tp_dealloc
        - [ ] The object's deallocator must be
          `Py_TYPE(obj)->tp_free((PyObject*)obj);`
      - [ ] Do _not_ define a `tp_getattr`
      - [ ] BONUS: write methods to present the state of the object via
        `tp_str` & `tp_repr`
      - [ ] _Do_ pass in the array of methods in `tp_methods`
      - [ ] _Do_ define a new method in `tp_new`
   - [ ] PyMethodDef arrays contain doc strings
      - [ ] Methods are cast to `PyCFunctions`s
   - [ ] Type methods use their type Object in the method signature.
   - [ ] Type creation method decrements the reference to self
     (`Py_DECREF(self);`) before each error-path exit (`return NULL;`)
   - [ ] No factory methods. Example: `khmer_new_readaligner`
   - [ ] Type object is passed to `PyType_Ready` and its return code is checked
     in `init_khmer()`
   - [ ] The reference count for the type object is incremented before adding
     it to the module: `Py_INCREF(&khmer_${OBJECTNAME}_Type);`.
