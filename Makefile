CC=g++
CFLAGS=-pthread# -Wall -Werror -Wextra
OPTFLAGS=-DANMS
LIBFLAGS=-fPIC

INCLUDE=src/include/

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

.PHONY: utils.o thermo.o mt.o clean clean2 lib.o tests all all2 staticlib

all: apps $(STATICLIB) $(SHAREDLIB)

all2: all tests

# Make object files for library

$(LIB_BUILD)%.o: $(LIB_SRC)%.cpp
	$(CC) $(CFLAGS) $(OPTFLAGS) -c $^ -o $@ -I$(INCLUDE) $(LIBFLAGS)

utils.o: $(addsuffix .o, $(addprefix $(LIB_BUILD), $(UTILS)))

thermo.o: $(addsuffix .o, $(addprefix $(LIB_BUILD), $(THERMO)))

mt.o: $(addsuffix .o, $(addprefix $(LIB_BUILD), $(DISTL)))

lib.o: utils.o thermo.o mt.o

# Make app and test object files and executable apps

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

run: tests
	$(foreach test, $(TEST_TARGETS), $(test) > $(test).txt &)

apps: $(APP_TARGETS)

# Clean build and executable files

clean:
	rm -f $(LIB_BUILD)*.o $(APP_BUILD)*.o $(TEST_BUILD)*.o $(STATICDIR)*.a $(SHAREDDIR)*.so
	$(patsubst %, find % -type f  ! -name "*.*"  -delete;, $(APP_EXE) $(TEST_EXE))

# Clean output files as well

clean2: clean
	rm -f $(OUTPUTS)

$(STATICDIR): 
	mkdir $(STATICDIR)

$(STATICLIB): $(STATICDIR) 
	ar rcs $@ $(LIBDEPS)

$(SHAREDDIR):
	mkdir $(SHAREDDIR)

$(SHAREDLIB): $(SHAREDDIR)
	$(CC) -shared $(LIBDEPS) -o $@ $(CFLAGS)

foo:
	g++ -Lbin/shared mt_test.cpp -o test -lChemESolver -I$(INCLUDE)