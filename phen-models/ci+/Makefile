# flags for speed. Note, that NDEBUG define is VERY important.
OPTFLAGS=-O3 -DNDEBUG 
# activate this for debug
# OPTFLAGS=-O0 -g
LIBS=-lboost_program_options
CPPFLAGS=-I../../numeric_bindings/ -I../../boost-ublasx

all: cip.x
	true
test: cip.x
	./cip.x -O test-cip.txt --input-params 0.1 1.0 1.0 0.1 0.0001 0.0001 0.1 1.0 0.01 0.1 2 50 2000
cip.x: cipModel.o
	g++ -o cip.x cipModel.o -lstdc++ -lm $(OPTFLAGS) -llapack $(LIBS)
cipModel.o: cip-model.cpp solver.hpp solver_impl.hpp
	g++ -std=c++11 -o cipModel.o -c cip-model.cpp $(OPTFLAGS)  $(CPPFLAGS)
clean:
	rm -f cip.x cipModel.o test-cip.txt
