# flags for speed. Note, that NDEBUG define is VERY important.
OPTFLAGS=-O3 -DNDEBUG 
# activate this for debug
# OPTFLAGS=-O0 -g
LIBS=-lboost_program_options
CPPFLAGS=-I../../numeric_bindings/ -I../../boost-ublasx

all: er.x
	true
test: er.x
	./er.x -O test-er.txt --input-params 0.1 0.1 0.00001 0.001 0.00001 0.1 0.0001 0.00001 2 50 2000
er.x: erModel.o
	g++ -o er.x erModel.o -lstdc++ -lm $(OPTFLAGS) -llapack $(LIBS)
erModel.o: er-model.cpp solver.hpp solver_impl.hpp
	g++ -std=c++11 -o erModel.o -c er-model.cpp $(OPTFLAGS)  $(CPPFLAGS)
clean:
	rm -f er.x erModel.o test-er.txt
