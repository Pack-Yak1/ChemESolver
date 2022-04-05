CC=g++
CFLAGS=-pthread #-Wall -Werror -Wextra -g# -DDEBUG
OPTFLAGS=-DANMS
LIBFLAGS=-fPIC

INCLUDE=src/include/
LIB_HEADERS=$(wildcard $(INCLUDE)*.hpp)

LIB_SRC=src/lib/
LIB_BUILD=build/lib/

APP_SRC=src/apps/
APP_BUILD=build/apps/
APP_EXE=bin/

TEST_SRC=src/tests/
TEST_BUILD=build/tests/
TEST_EXE=test_bin/

BINARIES=bin/
LIB_NAME=ChemESolver

# Library modules
UTILS=units coords vector_ops opt
THERMO=yx wilson antoine
DISTL=mc_cabe_thiele
LIB_OBJS=$(addprefix $(LIB_BUILD), $(addsuffix .o, $(UTILS) $(THERMO) $(DISTL)))

# Apps and tests
TESTS=opt_test yx_test mt_test
TEST_OBJS=$(patsubst %, $(TEST_BUILD)%.o, $(TESTS))
TEST_TARGETS=$(patsubst $(TEST_SRC)%.cpp, $(TEST_EXE)%, $(wildcard $(TEST_SRC)*.cpp))
APP_OBJS=$(patsubst $(APP_SRC)%.cpp, $(APP_BUILD)%.o, $(wildcard $(APP_SRC)*.cpp))
APP_TARGETS=$(patsubst $(APP_SRC)%.cpp, $(APP_EXE)%, $(wildcard $(APP_SRC)*.cpp))

# Library static and shared archives
LIBDEPS=$(wildcard build/lib/*.o)
STATICDIR=$(BINARIES)static/
SHAREDDIR=$(BINARIES)shared/
STATICLIB=$(STATICDIR)lib$(LIB_NAME).a
SHAREDLIB=$(SHAREDDIR)lib$(LIB_NAME).so

OUTPUT_FMTS=csv dat txt
OUTPUTS=$(foreach fmt, $(OUTPUT_FMTS), *.$(fmt) */*.$(fmt))

.PHONY: utils.o thermo.o mt.o lib.o clean clean2 tests all all2

# Build apps and libraries for users
all: lib.o apps $(STATICLIB) $(SHAREDLIB)

# Build everything including tests
all2: all tests

# Run all tests in background and log output. Sequential because we test speed
run: tests
	$(foreach test, $(TEST_TARGETS), $(test) > $(test).txt &&)

# Clean object files, executables and libraries
clean:
	rm -f $(LIB_BUILD)*.o $(APP_BUILD)*.o $(TEST_BUILD)*.o $(STATICDIR)*.a $(SHAREDDIR)*.so
	$(patsubst %, find % -type f  ! -name "*.*"  -delete;, $(APP_EXE) $(TEST_EXE))

# Clean output files as well
clean2: clean
	rm -f $(OUTPUTS)

# Rules for building library

$(LIB_BUILD)%.o: $(LIB_SRC)%.cpp $(LIB_HEADERS)
	$(CC) $(CFLAGS) $(OPTFLAGS) -c $< -o $@ -I$(INCLUDE) $(LIBFLAGS)

utils.o: $(addsuffix .o, $(addprefix $(LIB_BUILD), $(UTILS)))

thermo.o: $(addsuffix .o, $(addprefix $(LIB_BUILD), $(THERMO)))

mt.o: $(addsuffix .o, $(addprefix $(LIB_BUILD), $(DISTL)))

lib.o: utils.o thermo.o mt.o

# Rules for building apps and tests

.SECONDARY: $(TEST_OBJS) $(APP_OBJS)

$(APP_BUILD)%.o: $(APP_SRC)%.cpp
	$(CC) $(CFLAGS) -c $^ -o $@ -I$(INCLUDE)

$(APP_EXE)%: $(APP_BUILD)%.o $(LIB_OBJS)
	$(CC) $(CFLAGS) -o $@ $^

$(TEST_BUILD)%.o: $(TEST_SRC)%.cpp
	$(CC) $(CFLAGS) -c $^ -o $@ -I$(INCLUDE)

$(TEST_EXE)%: $(TEST_BUILD)%.o $(LIB_OBJS)
	$(CC) $(CFLAGS) -o $@ $^

tests: $(TEST_TARGETS)

apps: $(APP_TARGETS)

# Rules for building libraries

$(STATICDIR): 
	mkdir $(STATICDIR)

$(STATICLIB): $(STATICDIR)
	ar rcs $@ $(LIBDEPS)

$(SHAREDDIR):
	mkdir $(SHAREDDIR)

$(SHAREDLIB): $(SHAREDDIR)
	$(CC) -shared $(LIBDEPS) -o $@ $(CFLAGS)

# Convenience utilities for development
opt: $(TEST_EXE)opt_test
	test_bin/opt_test

yx: $(TEST_EXE)yx_test
	test_bin/yx_test

mt: $(TEST_EXE)mt_test
	test_bin/mt_test