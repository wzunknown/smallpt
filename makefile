ifeq ($(OS),Windows_NT) 
	detected_OS := Windows
else
	detected_OS := $(shell sh -c 'uname 2>/dev/null || echo Unknown')
endif

$(info Checking OS: $(detected_OS))

ifeq ($(detected_OS),Darwin)
	LLVM_DIR=/opt/homebrew/opt/llvm
	CXX=$(LLVM_DIR)/bin/clang++
	CXXFLAGS=-I$(LLVM_DIR)/include -O2 -fopenmp
	LDFLAGS=-L$(LLVM_DIR)/lib
else ifeq ($(detected_OS),Linux)
	CXX=g++
	CXXFLAGS=-O2 -fopenmp
	# LDFLAGS=-fopenmp
endif

EXE=mysmallpt

.PHONY: all
all: $(EXE)

.PHONY: png
png: $(EXE) image.ppm
	convert image.ppm image.png

mysmallpt: main.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

.PHONY: clean
clean:
	-rm -rf *.o mysmallpt *.ppm *.png *.jpg