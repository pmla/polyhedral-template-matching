CC = gcc
CPP = g++

C_SRC_FILES = canonical.c graph_data.c convex_hull_incremental.c \
	index_ptm.c alloy_types.c deformation_gradient.c \
	normalize_vertices.c reference_templates.c \
	fundamental_mappings.c \
	polar_decomposition.c \
	qcprot/qcprot.c qcprot/quat.c #unittest.c

CPP_SRC_FILES = neighbour_ordering.cpp voronoi/cell.cpp

C_SRC_MODULE_FILE = ptmmodule.c 

HEADER_FILES = alloy_types.h canonical.h convex_hull_incremental.h \
	deformation_gradient.h graph_data.h index_ptm.h \
	normalize_vertices.h reference_templates.h \
	neighbour_ordering.h polar_decomposition.h \
	qcprot/qcprot.h qcprot/quat.h

OBJDIR = .

LIBRARY = libptm.a

C_OBJECT_FILES = $(C_SRC_FILES:%.c=$(OBJDIR)/%.o) 
CPP_OBJECT_FILES = $(CPP_SRC_FILES:%.cpp=$(OBJDIR)/%.o) 
C_OBJECT_MODULE_FILE = $(C_SRC_MODULE_FILE:%.c=$(OBJDIR)/%.o) 

PYTHONMODULE = ptmmodule.so

PYTHON = python

PYTHONVERSION = $(shell $(PYTHON) -c 'import sys; print "{0[0]}.{0[1]}".format(sys.version_info)')
PYTHONPREFIX = $(shell $(PYTHON) -c 'import sys; print sys.prefix')
PYTHONEXECPREFIX = $(shell $(PYTHON) -c 'import sys; print sys.exec_prefix')
NUMPY_INCLUDE := $(shell $(PYTHON) -c 'import numpy; print numpy.get_include()')/numpy

PYTHONINCLDIR = $(PYTHONPREFIX)/include/python$(PYTHONVERSION)
PYTHONLIBDIR = $(PYTHONEXECPREFIX)/lib/python$(PYTHONVERSION)/config
PYTHONLIB = python$(PYTHONVERSION)

CFLAGS = -std=c99 -fPIC -g -O3 -Wall -Wextra -I$(PYTHONINCLDIR) -I$(NUMPY_INCLUDE)
CPPFLAGS = -fPIC -g -O3 -Wall -Wextra -I$(PYTHONINCLDIR) -I$(NUMPY_INCLUDE)

ifeq ($(shell uname),Darwin)
MAKESHARED = -bundle -undefined dynamic_lookup
else
MAKESHARED = -shared
endif

all: $(OBJDIR) $(OBJDIR)/$(PYTHONMODULE)

lib: $(OBJDIR) $(OBJDIR)/$(LIBRARY)

# Rule for linking module
$(OBJDIR)/$(PYTHONMODULE): $(C_OBJECT_MODULE_FILE) $(CPP_OBJECT_MODULE_FILE) $(OBJDIR)/$(LIBRARY)
	$(CPP) $(MAKESHARED) -fPIC -g -O2 -o $@ $^ -L$(PYTHONLIBDIR) -l$(PYTHONLIB) -lm

$(OBJDIR)/$(LIBRARY): $(C_OBJECT_FILES) $(CPP_OBJECT_FILES)
	rm -f $@
	ar -r -s $@ $^

$(OBJDIR):
	mkdir $(OBJDIR)

# Rule for compiling C source
$(OBJDIR)/%.o: %.cpp  $(HEADER_FILES)
	$(CPP) -c $(CPPFLAGS) $(INCLUDES) -o $@ -I$(PYTHONINCLDIR) -I$(NUMPY_INCLUDE) $<

$(OBJDIR)/%.o: %.c  $(HEADER_FILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) -o $@ -I$(PYTHONINCLDIR) -I$(NUMPY_INCLUDE) $<

clean:
	rm -f $(OBJDIR)/*.o $(OBJDIR)/ptmmodule.so

cleanall: clean
	rm -rf build

