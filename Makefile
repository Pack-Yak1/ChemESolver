CC=g++
LINKS=-lm -lpthread
CFLAGS=-DANMS -g
UTILS=utils/units.o utils/coords.o utils/opt.o utils/vector_ops.o
THERMO=Thermodynamics/yx.o Thermodynamics/wilson.o Thermodynamics/antoine.o
DISTL=Distillation/mc_cabe_thiele.o
OUTPUT_FMTS=*/*.csv *.csv */*.dat *.dat
EXE=tests/opt_test tests/yx_test tests/mt_test tests/opt_test

utils/opt.o: utils/opt.cpp
	$(CC) $(CFLAGS) -c $^ -o utils/opt.o 

opt: $(UTILS) tests/opt_test.o 
	$(CC) -o tests/opt_test tests/opt_test.o $(UTILS) $(LINKS)
	tests/opt_test

yx: $(UTILS) $(THERMO) tests/yx_test.o 
	$(CC) -o tests/yx_test tests/yx_test.o $(THERMO) $(UTILS) $(LINKS)
	tests/yx_test

mt: $(UTILS) $(THERMO) $(DISTL) tests/mt_test.o
	$(CC) -o tests/mt_test tests/mt_test.o $(DISTL) $(THERMO) $(UTILS) $(LINKS)
	tests/mt_test

clean:
	rm -f $(DISTL) $(THERMO) $(UTILS) $(OUTPUT_FMTS) $(EXE) tests/*.o