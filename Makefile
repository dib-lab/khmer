# Use multiple threads?
# Set this variable to true if you wish the codes to use multiple threads when they can.
WANT_THREADING=true

# Profile?
# Set this variable to true if you wish to profile the codes.
WANT_PROFILING=false

# Which profiling tool to use?
# Assuming you have TAU installed and setup properly, you can instrument codes with it to get detailed multi-threaded profiling.
# Otherwise, gprof is able to give you some information without threading info.
# Choose one of: gprof, TAU
PROFILER_OF_CHOICE=TAU

# Perform extra sanity checking?
# Set this variable to true if you wish the codes to perform extra sanity checking (to the possible detriment of performance).
WANT_EXTRA_SANITY_CHECKING=false

# Compile with debugging symbols?
# Set this variable to true if you wish the codes to be built with debugging symbols (increases code size and does not always produce accurate stepping in a debugger when optimization is turned on).
WANT_DEBUGGING=true

# Compile with tracing logic turned on?
# Set this variable to true if you want to use instrumentation provided in the sources for debugging purposes and are willing to accept the overhead such instrumentation introduces.
WITH_INTERNAL_TRACING=false

# Compile with performance metrics turned on?
# Set this variable to true if you want to use instrumentation provided in the sources for performance measurement purposes and are willing to accept the overhead such instrumentation introduces.
WITH_INTERNAL_METRICS=false

# Use Cython?
# Set this variable to true if you wish to build the Python wrapper with Cython rather than the directly using the Python C API.
USE_CYTHON=false

CXXFLAGS=
CXX_WARNING_FLAGS=-Wall
CXX_OPTIMIZATION_FLAGS=-O3
CXX_SHARED_LIB_FLAGS=-fPIC
CXXFLAGS+= $(CXX_WARNING_FLAGS) $(CXX_OPTIMIZATION_FLAGS) $(CXX_SHARED_LIB_FLAGS)

LIBS=

ifeq ($(WANT_DEBUGGING), true)
CXX_DEBUG_FLAGS=-g
CXXFLAGS+= $(CXX_DEBUG_FLAGS)
else
CXX_DEBUG_FLAGS=
endif

ifeq ($(WANT_EXTRA_SANITY_CHECKING), true)
DEFINE_KHMER_EXTRA_SANITY_CHECKS=-DKHMER_EXTRA_SANITY_CHECKS
CXXFLAGS+= $(DEFINE_KHMER_EXTRA_SANITY_CHECKS)
else
DEFINE_KHMER_EXTRA_SANITY_CHECKS=
endif

ifeq ($(WANT_THREADING), true)
DEFINE_KHMER_THREADED=-DKHMER_THREADED
CXX_THREADING_FLAGS=-fopenmp
THREADING_LIBS=-fopenmp
CXXFLAGS+= $(DEFINE_KHMER_THREADED) $(CXX_THREADING_FLAGS)
LIBS+= $(THREADING_LIBS)
else
DEFINE_KHMER_THREADED=
CXX_THREADING_FLAGS=
THREADING_LIBS=
endif

ifeq ($(WANT_PROFILING), true)
ifeq ($(PROFILER_OF_CHOICE), TAU)
CXX=tau_cxx.sh
endif
ifeq ($(PROFILER_OF_CHOICE), gprof)
CXXFLAGS+= -pg
endif
endif

ifeq ($(WITH_INTERNAL_TRACING), true)
CXXFLAGS+= -DWITH_INTERNAL_TRACING
endif

ifeq ($(WITH_INTERNAL_METRICS), true)
CXXFLAGS+= -DWITH_INTERNAL_METRICS
endif

ifeq ($(USE_CYTHON), true)
CYTHON_ENABLED_BOOL=True
else
CYTHON_ENABLED_BOOL=False
endif

all: lib_files python_files

clean:
	cd lib && make clean
	cd python && make clean
	cd tests && rm -rf khmertest_*

doc: FORCE
	cd doc && make html

lib_files:
	cd lib && \
	make CXX="$(CXX)" CXXFLAGS="$(CXXFLAGS)" LIBS="$(LIBS)"

python_files: lib_files
	cd python && \
	make 	DEFINE_KHMER_THREADED="$(DEFINE_KHMER_THREADED)" \
		DEFINE_KHMER_EXTRA_SANITY_CHECKS="$(DEFINE_KHMER_EXTRA_SANITY_CHECKS)" \
		CXX_DEBUG_FLAGS="$(CXX_DEBUG_FLAGS)" \
		CXX_THREADING_FLAGS="$(CXX_THREADING_FLAGS)" \
		THREADING_LIBS="$(THREADING_LIBS)" \
		CYTHON_ENABLED_BOOL="$(CYTHON_ENABLED_BOOL)"
#	python setup.py build_ext -i

test: all
	nosetests -v -x

FORCE:
