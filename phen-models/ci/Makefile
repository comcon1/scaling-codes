# flags for speed. Note, that NDEBUG define is VERY important.
OPTFLAGS=-O3 -DNDEBUG 
# activate this for debug
# OPTFLAGS=-O0 -g
LIBS=-lboost_program_options
CPPFLAGS=-I../../numeric_bindings/ -I../../boost-ublasx

all: ci.x
	true
test: ci.x
	./ci.x -O test-ci.txt --input-params 0.1 0.1 0.00001 0.001 0.00001 0.1 0.0001 0.00001 2 50 2000
ci.x: ciModel.o
	g++ -o ci.x ciModel.o -lstdc++ -lm $(OPTFLAGS) -llapack $(LIBS)
ciModel.o: ci-model.cpp solver.hpp solver_impl.hpp
	g++ -std=c++11 -o ciModel.o -c ci-model.cpp $(OPTFLAGS)  $(CPPFLAGS)
clean:
	rm -f ci.x ciModel.o test-ci.txt
