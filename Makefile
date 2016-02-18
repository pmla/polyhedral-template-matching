#time g++ -std=c++11 -g -O1 -Wall -Wextra main.cpp index_PTM.cpp convex_hull.cpp graph_data.cpp graph_tools.cpp qcprot.c nauty25r9/nauty.c nauty25r9/nautil.c nauty25r9/naugraph.c nauty25r9/schreier.c nauty25r9/naurng.c  -lCGAL -lgmp -frounding-math -Inauty25r9 -Wno-missing-field-initializers -Wno-sign-compare -Wno-write-strings -Wno-unused-parameter

CXX = g++-4.9

CPP_SRC_FILES = ptmmodule.cpp canonical.cpp graph_data.cpp convex_hull.cpp index_PTM.cpp\
	alloy_types.cpp qcprot.cpp deformation_gradient.cpp normalize_vertices.cpp quat.cpp\
	svdpolar/polar_decomposition.cpp\
	#selftest.cpp

HEADER_FILES = $(wildcard *.h) $(wildcard *.hpp)

CPP_OBJECT_FILES = $(CPP_SRC_FILES:%.cpp=%.o)

PYTHONMODULE = ptmmodule.so

PYTHON = python

PYTHONVERSION = $(shell $(PYTHON) -c 'import sys; print "{0[0]}.{0[1]}".format(sys.version_info)')
PYTHONPREFIX = $(shell $(PYTHON) -c 'import sys; print sys.prefix')
PYTHONEXECPREFIX = $(shell $(PYTHON) -c 'import sys; print sys.exec_prefix')
NUMPY_INCLUDE := $(shell $(PYTHON) -c 'import numpy; print numpy.get_include()')/numpy

PYTHONINCLDIR = $(PYTHONPREFIX)/include/python$(PYTHONVERSION)
PYTHONLIBDIR = $(PYTHONEXECPREFIX)/lib/python$(PYTHONVERSION)/config
PYTHONLIB = python$(PYTHONVERSION)

CXXFLAGS = -std=c++11 -fPIC -g -O2 -Wall -Wextra -Wno-missing-field-initializers -Wno-sign-compare -Wno-write-strings -Wno-unused-parameter 

ifeq ($(shell uname),Darwin)
MAKESHARED = -bundle -undefined dynamic_lookup
else
MAKESHARED = -shared
endif

# Rule for linking module
$(PYTHONMODULE): $(CPP_OBJECT_FILES)
	$(CXX) $(MAKESHARED) -fPIC -g -O2 -o $(PYTHONMODULE) $(CPP_OBJECT_FILES) -z defs -L$(PYTHONLIBDIR) -l$(PYTHONLIB)

# Lazy again: all object files depend on all header files
%.o: $(HEADER_FILES)

# Rule for compiling CPP source
%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) -o $@ -I$(PYTHONINCLDIR) -I$(NUMPY_INCLUDE) $<
