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
	LDFLAGS=-L$(LLVM_DIR)/lib -lyaml-cpp
else ifeq ($(detected_OS),Linux)
	CXX=g++
	CXXFLAGS=-O2 -fopenmp
	LDFLAGS=-lyaml-cpp
endif

CPPFLAGS=-Iinclude -std=c++17

SRC_DIR=src
OBJ_DIR=obj


SRC=$(wildcard $(SRC_DIR)/*.cpp)
OBJ=$(SRC:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

EXE=mysmallpt

.PHONY: all
all: $(EXE)

.PHONY: png
png: $(EXE) image.ppm
	convert image.ppm image.png

mysmallpt: mysmallpt.o $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

demo: demo.o $(OBJ_DIR)/vector.o $(OBJ_DIR)/sphere.o $(OBJ_DIR)/color.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@ $(LDFLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@ $(LDFLAGS)

$(OBJ_DIR):
	mkdir -p $@

.PHONY: clean
clean:
	-rm -rf *.o mysmallpt *.ppm *.png *.jpg $(OBJ)