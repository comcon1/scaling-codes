NUM_CORES=4

all: 
	@echo 'This makefile runs scanning procedure.'
	@echo 'First, compile the code in the folder "../solver"'
	@echo 'Next use "make scan" for start.'

scan: logoptim.log
	true

logoptim.log: _current_/1 ../solver/tw.x
	L1M_NTHREADS=$(NUM_CORES) ./optimizer.py

_current_/1:
	mkdir _current_
	touch _current_/1

clean:
	rm -fr logoptim.log _current_ __pycache__
