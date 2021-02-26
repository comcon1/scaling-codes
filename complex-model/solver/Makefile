# flags for speed. Note, that NDEBUG define is VERY important.
OPTFLAGS=-O3 -DNDEBUG 
# activate this for debug
# OPTFLAGS=-O0 -g
LIBS=-lboost_program_options
CPPFLAGS=-I../../numeric_bindings/ -I../../boost-ublasx

all: tw.x
	true
test: tw.x
	./tw.x -O test-tw.txt --input-params 1.0 2000.0 0.001 7.3e-5 4.3 10.0 0.2 5.6 0.0005 3.9 3.0 1e-4 0.0003 350.0 2e-05 6e-05 1.4 --change-d 10.4 1.5
tw.x: twModel.o
	g++ -o tw.x twModel.o -lstdc++ -lm $(OPTFLAGS) -llapack $(LIBS)
twModel.o: tw-model.cpp solver.hpp solver_impl.hpp
	g++ -std=c++11 -o twModel.o -c tw-model.cpp $(OPTFLAGS)  $(CPPFLAGS)
clean:
	rm -f tw.x twModel.o test-tw.txt
