# Generic GNUMakefile

# Just a snippet to stop executing under other make(1) commands
# that won't understand these lines
ifneq (,)
This makefile requires GNU Make.
endif

PROGRAM = benchmark
CPP_FILES = main.cpp graph_data.cpp convex_hull_incremental.cpp \
	index_ptm.cpp alloy_types.cpp deformation_gradient.cpp \
	normalize_vertices.cpp \
	structure_matcher.cpp \
	canonical_coloured.cpp \
	graph_tools.cpp \
	quat.cpp \
	polar.cpp \
	unittest.cpp\
	initialize_data.cpp \
	neighbour_ordering.cpp voronoi/cell.cpp

CPPOBJS := $(patsubst %.cpp, %.o, $(CPP_FILES))
LDFLAGS =
LDLIBS = -lm #-fno-omit-frame-pointer -fsanitize=address

CPP = g++

HEADER_FILES = alloy_types.hpp convex_hull_incremental.hpp \
	deformation_gradient.hpp graph_data.hpp index_ptm.h \
	normalize_vertices.hpp \
	structure_matcher.hpp \
	canonical_coloured.hpp \
	graph_tools.hpp \
	fundamental_mappings.hpp \
	quat.hpp \
	polar.hpp \
	initialize_data.hpp \
	neighbour_ordering.hpp \
	voronoi/cell.hpp

OBJDIR = .

CPP_OBJECT_FILES = $(CPP_SRC_FILES:%.cpp=$(OBJDIR)/%.o) 
C_OBJECT_MODULE_FILE = $(C_SRC_MODULE_FILE:%.c=$(OBJDIR)/%.o) 

#CPPFLAGS = -g -O3 -std=c++11 -Wall -Wextra -Wvla -pedantic #-fno-omit-frame-pointer -fsanitize=address
CPPFLAGS = -g -O3 -Wall -Wextra -Wvla -pedantic #-fno-omit-frame-pointer -fsanitize=address


all: $(PROGRAM)

$(PROGRAM): $(CPPOBJS)
	$(CPP) -c $(CPPFLAGS) $(CPPOBJS)
	$(CPP) $(CPPOBJS) -o $(PROGRAM) $(LDLIBS) $(LDFLAGS)

