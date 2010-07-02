from distutils.core import setup, Extension

# the c++ extension module (needs to be linked in with ktable.o ...)
extension_mod = Extension("khmer._khmermodule",
                          ["_khmermodule.cc"],
                          extra_compile_args=['-g'],
                          include_dirs=['../lib',],
                          library_dirs=['../lib',],
                          extra_objects=['../lib/ktable.o',
                                         '../lib/hashtable.o'])


# python modules: only 'khmer'
py_mod = 'khmer'

setup(name = "khmer", version = "0.2",
      description = 'khmer k-mer counting library',
      author = 'C. Titus Brown',
      author_email = 'ctb@msu.edu',
      url = 'http://cartwheel.caltech.edu/',
      packages = [py_mod,],
      ext_modules = [extension_mod,])
