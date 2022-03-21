CC=g++
CFLAGS=-lnlopt -lm -D DEBUG
# DEPS=utils/units.cpp

# %.o: %.c $(DEPS)
# 	$(CC) -c -o $@ $< $(CFLAGS)

yx: Thermodynamics/yx.o Thermodynamics/wilson.o Thermodynamics/antoine.o utils/units.o utils/coords.o
	$(CC) -o Thermodynamics/yx Thermodynamics/yx.o Thermodynamics/wilson.o Thermodynamics/antoine.o utils/units.o utils/coords.o $(CFLAGS)
	Thermodynamics/yx

mt: Distillation/mc_cabe_thiele.o Thermodynamics/yx.o Thermodynamics/wilson.o Thermodynamics/antoine.o utils/units.o utils/coords.o
	$(CC) -o Distillation/mc_cabe_thiele Distillation/mc_cabe_thiele.o Thermodynamics/yx.o Thermodynamics/wilson.o Thermodynamics/antoine.o utils/units.o utils/coords.o $(CFLAGS)
	Distillation/mc_cabe_thiele

opt: tests/opt_test.o utils/opt.o utils/vector_ops.o Distillation/mc_cabe_thiele.o Thermodynamics/yx.o Thermodynamics/wilson.o Thermodynamics/antoine.o utils/units.o utils/coords.o
	$(CC) -o tests/opt_test tests/opt_test.o utils/opt.o utils/vector_ops.o Distillation/mc_cabe_thiele.o Thermodynamics/yx.o Thermodynamics/wilson.o Thermodynamics/antoine.o utils/units.o utils/coords.o $(CFLAGS)
	tests/opt_test

clean:
	rm -f Distillation/mc_cabe_thiele.o Thermodynamics/yx.o Thermodynamics/wilson.o Thermodynamics/antoine.o utils/units.o utils/point.o utils/opt.o utils/coords.o utils/vector_ops.o */*.csv *.csv */*.dat *.dat