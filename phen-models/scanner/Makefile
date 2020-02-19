NCORES=4

all: 
	@echo Use one of following targets: erscan ciscan ci+scan

erscan: paramset-er.txt ../er/er.x
	./launcher_1m.py ../er/er.x paramset-er.txt $(NCORES)

paramset-er.txt: mkparamset-er.py
	./mkparamset-er.py

ciscan: paramset-ci.txt ../ci/ci.x
	./launcher_1m.py ../ci/ci.x paramset-ci.txt $(NCORES)

paramset-ci.txt: mkparamset-ci.py
	./mkparamset-ci.py

ci+scan: paramset-cip.txt ../ci+/cip.x
	./launcher_1m.py ../ci+/cip.x paramset-cip.txt $(NCORES)

paramset-cip.txt: mkparamset-cip.py
	./mkparamset-cip.py

clean:
	rm paramset*.txt
	rm -fr Sim_paramset*
