CC=g++
CFLAGS=-lnlopt -lm -D DEBUG 
UTILS=utils/units.o utils/coords.o utils/opt.o utils/vector_ops.o
THERMO=Thermodynamics/yx.o Thermodynamics/wilson.o Thermodynamics/antoine.o
DISTL=Distillation/mc_cabe_thiele.o
OUTPUT_FMTS=*/*.csv *.csv */*.dat *.dat
EXE=tests/opt_test Thermodynamics/yx tests/mt_test tests/opt_test


opt: tests/opt_test.o $(UTILS)
	$(CC) -o tests/opt_test tests/opt_test.o $(UTILS) $(CFLAGS)
	tests/opt_test

yx: $(THERMO) $(UTILS)
	$(CC) -o Thermodynamics/yx $(THERMO) $(UTILS) $(CFLAGS)
	Thermodynamics/yx

mt: tests/mt_test.o $(DISTL) $(THERMO) $(UTILS)
	$(CC) -o tests/mt_test tests/mt_test.o $(DISTL) $(THERMO) $(UTILS) $(CFLAGS)
	tests/mt_test

clean:
	rm -f $(DISTL) $(THERMO) $(UTILS) $(OUTPUT_FMTS) $(EXE) tests/*.o