"""
khmer specific YouCompleteMe configuration
"""
# pylint: disable=missing-docstring

import os
import sys
from distutils import core
from distutils import ccompiler

# These are the compilation flags that will be used in case there's no
# compilation database set (by default, one is not set).
# CHANGE THIS LIST OF FLAGS. YES, THIS IS THE DROID YOU HAVE BEEN LOOKING FOR.
FLAGS = [
    # '-Wextra',
    '-Werror',
    # '-Wc++98-compat',
    # '-Wno-long-long',
    # '-Wno-variadic-macros',
    # '-fexceptions',
    # THIS IS IMPORTANT! Without a "-std=<something>" flag, clang won't know
    # which language to use when compiling headers. So it will guess. Badly. So
    # C++ headers will be compiled as C headers. You don't want that so ALWAYS
    # specify a "-std=<something>".
    # For a C project, you would set this to something like 'c99' instead of
    # 'c++11'.
    '-std=c++11',
    # ...and the same thing goes for the magic -x option which specifies the
    # language that the files to be compiled are written in. This is mostly
    # relevant for c++ headers.
    # For a C project, you would set this to 'c' instead of 'c++'.
    '-x',
    'c++',
    '-I/usr/include/c++/4.8',
    '-I/usr/include/x86_64-linux-gnu/c++/4.8',
    '-I/usr/include/c++/4.8/backward',
    '-I/usr/lib/gcc/x86_64-linux-gnu/4.8/include',
    '-I/usr/local/include',
    '-I/usr/lib/gcc/x86_64-linux-gnu/4.8/include-fixed',
    '-I/usr/include/x86_64-linux-gnu',
    '-I/usr/include'
]


def directory_of_this_script():
    return os.path.dirname(os.path.abspath(__file__))


def make_relative_paths_in_flags_absolute(flags, working_directory):
    if not working_directory:
        return list(flags)
    new_flags = []
    make_next_absolute = False
    path_flags = ['-isystem', '-I', '-iquote', '--sysroot=']
    for flag in flags:
        new_flag = flag

        if make_next_absolute:
            make_next_absolute = False
            if not flag.startswith('/'):
                new_flag = os.path.join(working_directory, flag)

        for path_flag in path_flags:
            if flag == path_flag:
                make_next_absolute = True
                break

            if flag.startswith(path_flag):
                path = flag[len(path_flag):]
                new_flag = path_flag + os.path.join(working_directory, path)
                break

        if new_flag:
            new_flags.append(new_flag)
    return new_flags


def is_header_file(filename):
    extension = os.path.splitext(filename)[1]
    return extension in ['.h', '.hxx', '.hpp', '.hh']


def get_prepared_build_ext():
    sys.path.insert(1, directory_of_this_script())
    core._setup_stop_after = "commandline"  # pylint: disable=protected-access
    env = {'__file__': "setup.py"}
    sys.argv = ["setup.py", "build_ext"]
    try:
        script = open(sys.argv[0])
        try:
            exec script in env, env  # pylint: disable=exec-used
        finally:
            script.close()
    except SystemExit:
        pass
    except:
        raise

    core._setup_stop_after = None  # pylint: disable=protected-access
    distribution = core._setup_distribution  # pylint: disable=protected-access
    distribution.dry_run = True
    build_ext = distribution.get_command_obj(command="build_ext")
    build_ext.ensure_finalized()
    build_ext.run()
    return build_ext


def get_setuptools_options():
    build_ext = get_prepared_build_ext()
    compiler = build_ext.compiler
    options = []
    if compiler:
        if compiler.compiler_so:
            options.extend(compiler.compiler_so[1:])
        options.extend(ccompiler.gen_preprocess_options(
            macros=compiler.macros,
            include_dirs=compiler.include_dirs))
        options.extend(ccompiler.gen_lib_options(
            compiler=compiler,
            library_dirs=compiler.library_dirs,
            runtime_library_dirs=compiler.runtime_library_dirs,
            libraries=compiler.libraries))
    if build_ext.extensions:
        for ext in build_ext.extensions:
            options.extend(ext.extra_compile_args)
    return options

FLAGS.extend(get_setuptools_options())

SOURCE_EXTENSIONS = ['.cpp', '.cxx', '.cc', '.c', '.m', '.mm']


def FlagsForFile(  # pylint: disable=unused-argument,invalid-name
        filename, **kwargs):
    relative_to = directory_of_this_script()
    final_flags = make_relative_paths_in_flags_absolute(FLAGS, relative_to)

    return {
        'flags': final_flags,
        'do_cache': True
    }
