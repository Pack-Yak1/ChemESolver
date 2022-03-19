CC=g++
CFLAGS=-lnlopt -lm
# DEPS=utils/units.cpp

# %.o: %.c $(DEPS)
# 	$(CC) -c -o $@ $< $(CFLAGS)

yx: Thermodynamics/yx.o Thermodynamics/wilson.o Thermodynamics/antoine.o utils/units.o utils/coords.o
	$(CC) -o Thermodynamics/yx Thermodynamics/yx.o Thermodynamics/wilson.o Thermodynamics/antoine.o utils/units.o utils/coords.o $(CFLAGS)
	Thermodynamics/yx

mt: Distillation/mc_cabe_thiele.o Thermodynamics/yx.o Thermodynamics/wilson.o Thermodynamics/antoine.o utils/units.o utils/coords.o
	$(CC) -o Distillation/mc_cabe_thiele Distillation/mc_cabe_thiele.o Thermodynamics/yx.o Thermodynamics/wilson.o Thermodynamics/antoine.o utils/units.o utils/coords.o $(CFLAGS)
	Distillation/mc_cabe_thiele

opt: utils/opt.o Distillation/mc_cabe_thiele.o Thermodynamics/yx.o Thermodynamics/wilson.o Thermodynamics/antoine.o utils/units.o utils/coords.o
	$(CC) -o utils/opt utils/opt.o Distillation/mc_cabe_thiele.o Thermodynamics/yx.o Thermodynamics/wilson.o Thermodynamics/antoine.o utils/units.o utils/coords.o $(CFLAGS)
	utils/opt

clean:
	rm -f Distillation/mc_cabe_thiele.o Thermodynamics/yx.o Thermodynamics/wilson.o Thermodynamics/antoine.o utils/units.o