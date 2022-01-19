LLVM_DIR=/opt/homebrew/opt/llvm
CXX=$(LLVM_DIR)/bin/clang++
CXXFLAGS=-I$(LLVM_DIR)/include -O2 -fopenmp      # Use this for gcc >= 4.2
LDFLAGS=-L$(LLVM_DIR)/lib

EXE=mysmallpt

all: ${EXE}

mysmallpt: main.cpp
	${CXX} $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# %.o: %.cpp
# 	${CXX} -c $^

.PHONY: clean
clean:
	-rm -rf *.o mysmallpt