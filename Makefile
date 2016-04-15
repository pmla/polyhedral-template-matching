CC = gcc

C_SRC_FILES = canonical.c graph_data.c convex_hull_incremental.c \
	index_PTM.c alloy_types.c qcprot.c deformation_gradient.c \
	normalize_vertices.c quat.c

C_SRC_MODULE_FILE = ptmmodule.c 

C_SRC_SVDPOLAR_FILES = polar_decomposition.c

HEADER_FILES = alloy_types.h canonical.h convex_hull_incremental.h \
	deformation_gradient.h graph_data.h index_PTM.h \
	normalize_vertices.h qcprot.h quat.h reference_templates.h \
	svdpolar/polar_decomposition.h \
	svdpolar/Singular_Value_Decomposition_Givens_QR_Factorization_Kernel.h \
	svdpolar/Singular_Value_Decomposition_Jacobi_Conjugation_Kernel.h \
	svdpolar/Singular_Value_Decomposition_Kernel_Declarations.h \
	svdpolar/Singular_Value_Decomposition_Main_Kernel_Body.h

OBJDIR = .

LIBRARY = libptm.a

C_OBJECT_FILES = $(C_SRC_FILES:%.c=$(OBJDIR)/%.o) 
C_OBJECT_SVDPOLAR_FILES = $(C_SRC_SVDPOLAR_FILES:%.c=$(OBJDIR)/%.o) 
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

CFLAGS = -std=c99 -fPIC -g -O3 -Wall -Wextra

ifeq ($(shell uname),Darwin)
MAKESHARED = -bundle -undefined dynamic_lookup
else
MAKESHARED = -shared
endif

all: $(OBJDIR) $(OBJDIR)/$(PYTHONMODULE)

lib: $(OBJDIR) $(OBJDIR)/$(LIBRARY)

# Rule for linking module
$(OBJDIR)/$(PYTHONMODULE): $(C_OBJECT_MODULE_FILE) $(OBJDIR)/$(LIBRARY)
	$(CC) $(MAKESHARED) -fPIC -g -O2 -o $@ $^ -L$(PYTHONLIBDIR) -l$(PYTHONLIB) -lm

$(OBJDIR)/$(LIBRARY): $(C_OBJECT_FILES) $(C_OBJECT_SVDPOLAR_FILES)
	rm -f $@
	ar -r -s $@ $^

$(OBJDIR):
	mkdir $(OBJDIR)

# Rule for compiling C source
$(OBJDIR)/%.o: svdpolar/%.c  $(HEADER_FILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) -o $@ -I$(PYTHONINCLDIR) -I$(NUMPY_INCLUDE) $<

$(OBJDIR)/%.o: %.c  $(HEADER_FILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) -o $@ -I$(PYTHONINCLDIR) -I$(NUMPY_INCLUDE) $<

clean:
	rm -f $(OBJDIR)/*.o $(OBJDIR)/ptmmodule.so

cleanall: clean
	rm -rf build
