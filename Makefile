CC=g++
LINKS=-lnlopt -lm -lpthread
CFLAGS=-DANMS -g
UTILS=utils/units.o utils/coords.o utils/opt.o utils/vector_ops.o
THERMO=Thermodynamics/yx.o Thermodynamics/wilson.o Thermodynamics/antoine.o
DISTL=Distillation/mc_cabe_thiele.o
OUTPUT_FMTS=*/*.csv *.csv */*.dat *.dat
EXE=tests/opt_test tests/yx_test tests/mt_test tests/opt_test

opt.o: 
	$(CC) $(CFLAGS) -c utils/opt.cpp -o utils/opt.o 

opt: tests/opt_test.o $(UTILS)
	$(CC) -o tests/opt_test tests/opt_test.o $(UTILS) $(LINKS)
	tests/opt_test

yx: tests/yx_test.o $(THERMO) $(UTILS)
	$(CC) -o tests/yx_test tests/yx_test.o $(THERMO) $(UTILS) $(LINKS)
	tests/yx_test

mt: tests/mt_test.o $(DISTL) $(THERMO) $(UTILS)
	$(CC) -o tests/mt_test tests/mt_test.o $(DISTL) $(THERMO) $(UTILS) $(LINKS)
	tests/mt_test

clean:
	rm -f $(DISTL) $(THERMO) $(UTILS) $(OUTPUT_FMTS) $(EXE) tests/*.o