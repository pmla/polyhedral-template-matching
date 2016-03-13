CXX = gcc

C_SRC_FILES = ptmmodule.c canonical.c graph_data.c convex_hull_incremental.c index_PTM.c\
	alloy_types.c qcprot.c deformation_gradient.c normalize_vertices.c quat.c\
	svdpolar/polar_decomposition.c\
	#selftest.c

HEADER_FILES = $(wildcard *.h)

C_OBJECT_FILES = $(C_SRC_FILES:%.c=%.o)

PYTHONMODULE = ptmmodule.so

PYTHON = python

PYTHONVERSION = $(shell $(PYTHON) -c 'import sys; print "{0[0]}.{0[1]}".format(sys.version_info)')
PYTHONPREFIX = $(shell $(PYTHON) -c 'import sys; print sys.prefix')
PYTHONEXECPREFIX = $(shell $(PYTHON) -c 'import sys; print sys.exec_prefix')
NUMPY_INCLUDE := $(shell $(PYTHON) -c 'import numpy; print numpy.get_include()')/numpy

PYTHONINCLDIR = $(PYTHONPREFIX)/include/python$(PYTHONVERSION)
PYTHONLIBDIR = $(PYTHONEXECPREFIX)/lib/python$(PYTHONVERSION)/config
PYTHONLIB = python$(PYTHONVERSION)

CXXFLAGS = -std=c99 -fPIC -g -O3 -Wall -Wextra

ifeq ($(shell uname),Darwin)
MAKESHARED = -bundle -undefined dynamic_lookup
else
MAKESHARED = -shared
endif

# Rule for linking module
$(PYTHONMODULE): $(C_OBJECT_FILES)
	$(CXX) $(MAKESHARED) -fPIC -g -O2 -o $(PYTHONMODULE) $(C_OBJECT_FILES) -z defs -L$(PYTHONLIBDIR) -l$(PYTHONLIB) -lm

# Lazy again: all object files depend on all header files
%.o: $(HEADER_FILES)

# Rule for compiling C source
%.o: %.c
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) -o $@ -I$(PYTHONINCLDIR) -I$(NUMPY_INCLUDE) $<
