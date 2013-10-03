# Profile?
# Set this variable to true if you wish to profile the codes.
WANT_PROFILING=false

# Which profiling tool to use?
# Assuming you have TAU installed and setup properly, 
# you can instrument codes with it to get detailed multi-threaded profiling.
# Otherwise, gprof is able to give you some information without threading info.
# Choose one of: gprof, TAU
PROFILER_OF_CHOICE=gprof

# Perform extra sanity checking?
# Set this variable to true 
# if you wish the codes to perform extra sanity checking 
# (to the possible detriment of performance).
WANT_EXTRA_SANITY_CHECKING=false

# Compile with debugging symbols?
# Set this variable to true 
# if you wish the codes to be built with debugging symbols 
# (increases code size 
# and does not always produce accurate stepping in a debugger 
# when optimization is turned on).
WANT_DEBUGGING=false

# Compile with tracing logic turned on?
# Set this variable to true if you want to use instrumentation provided 
# in the sources for debugging purposes 
# and are willing to accept the overhead such instrumentation introduces.
WITH_INTERNAL_TRACING=false

# Trace state transitions?
# Set this variable to true if you want to use instrumentation which reports
# on transitions which occur between the states of the various elements of the 
# processing stack.
# 'WITH_INTERNAL_TRACING' must be true for this to have effect.
TRACE_STATE_CHANGES=true

# Trace busywaits?
# Set this variable to true if you want to use instrumentation which reports
# on various busywaits, such as synchronization barriers, spinlock trials, and 
# polling loops.
# Spinlock trials will only be reported if 'TRACE_SPINLOCKS' is also true.
# 'WITH_INTERNAL_TRACING' must be true for this to have effect.
TRACE_BUSYWAITS=true

# Trace spinlocks?
# Set this variable to true if you want to use instrumentation which reports 
# on entries into and exits from spinlocks and spinlock trials.
# Spinlock trials will only be reported if 'TRACE_BUSYWAITS' is also true.
# 'WITH_INTERNAL_TRACING' must be true for this to have effect.
TRACE_SPINLOCKS=true

# Trace memory copies?
# Set this variable to true if you want to use instrumentation which reports
# on the sizes of memory copies between various caches and buffers.
# 'WITH_INTERNAL_TRACING' must be true for this to have effect.
TRACE_MEMCOPIES=true

# Trace data?
# Set this variable to true if you want to use instrumentation which reports 
# on the pieces of data handled by various levels of the processing stack.
# WARNING! This can generate *very* large trace logs - use with caution or 
# lots of free storage.
# 'WITH_INTERNAL_TRACING' must be true for this to have effect.
TRACE_DATA=false

# Compile with performance metrics turned on?
# Set this variable to true if you want to use instrumentation provided 
# in the sources for performance measurement purposes 
# and are willing to accept the overhead such instrumentation introduces.
WITH_INTERNAL_METRICS=false


### NOTE: No user-servicable parts below this line! ###


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

ifeq ($(WANT_PROFILING), true)
ifeq ($(PROFILER_OF_CHOICE), TAU)
CXX=tau_cxx.sh
endif
ifeq ($(PROFILER_OF_CHOICE), gprof)
PROFILING_LIBS=-pg
CXXFLAGS+= -pg
LIBS+= $(PROFILING_LIBS)
endif
endif

ifeq ($(WITH_INTERNAL_TRACING), true)
CXXFLAGS+= -DWITH_INTERNAL_TRACING
ifeq ($(TRACE_STATE_CHANGES), true)
CXXFLAGS+= -DTRACE_STATE_CHANGES
endif
ifeq ($(TRACE_BUSYWAITS), true)
CXXFLAGS+= -DTRACE_BUSYWAITS
endif
ifeq ($(TRACE_SPINLOCKS), true)
CXXFLAGS+= -DTRACE_SPINLOCKS
endif
ifeq ($(TRACE_MEMCOPIES), true)
CXXFLAGS+= -DTRACE_MEMCOPIES
endif
ifeq ($(TRACE_DATA), true)
CXXFLAGS+= -DTRACE_DATA
endif
endif

ifeq ($(WITH_INTERNAL_METRICS), true)
CXXFLAGS+= -DWITH_INTERNAL_METRICS
endif

# Place POSIX threads last in linking order, if needed.
ifneq ($(shell uname), Linux)
LIBS+= -pthread
endif

export CXX
export CXXFLAGS
export LIBS

all:
	python setup.py build

clean:
	python setup.py clean --all
	cd lib && make clean
	cd tests && rm -rf khmertest_*

doc: FORCE
	python setup.py build_sphinx

lib_files:
	cd lib && \
	make

test: all
	python setup.py nosetests

FORCE:
