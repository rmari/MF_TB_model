#Makefile TB, Romain Mari

CXX = g++
# Compilation Flags
CXXWarnings = -Wall
CXXFLAGS_common = $(CXXWarnings) $(CXXFLAGS_EXTRA) -std=c++11
CXXFLAGS_optimized = -O3 $(CXXFLAGS_common)
CXXFLAGS_debug = -g $(CXXFLAGS_common)

CXXFLAGS_debug += -pg
ifeq ($(CXX),$(filter $(CXX),g++ clang++ gcc))
	CXXFLAGS_debug += -Wextra -Wshadow
endif

INCLUDE_PATH =

OUT_CODE=TBM
SRC = dynamics.cpp main.cpp
LIB = boxing.h configuration.h MersenneTwister.h
OBJ=$(SRC:.cpp=.o)
DATE=$(shell date +%D | sed "s./..g")

all:	CXXFLAGS= $(CXXFLAGS_optimized)
all:	$(OUT_CODE)

debug:	CXXFLAGS= $(CXXFLAGS_debug)
debug:	$(OUT_CODE)

$(OUT_CODE): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ) -lstdc++ -lm

$(OBJ): $(SRC)
	$(CXX) $(CXXFLAGS) $(INCLUDE_PATH) -c $(@:.o=.cpp)

clean:
	rm -f  $(OBJ)

install:
	cp -p $(OUT_CODE) $(install_dir)
