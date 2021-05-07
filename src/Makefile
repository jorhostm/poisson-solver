IDIR =../include
LDIR =../lib
CC=gcc
CFLAGS= -O3 -I$(IDIR)
FORMAT=

BDIR=../build
TDIR=$(BDIR)/test
SDIR=$(BDIR)/solutions
ODIR=obj

EXAMPLE_LIBS= -lm -lpthread
DEMO_LIBS= -lm

ifeq ($(OS),Windows_NT)
    CFLAGS += -D WIN32
	DEMO_LIBS += -lfreeglut -lglu32 -lopengl32
    FORMAT = .exe
    ifeq ($(PROCESSOR_ARCHITEW6432),AMD64)
        CFLAGS += -D AMD64
    else
        ifeq ($(PROCESSOR_ARCHITECTURE),AMD64)
            CFLAGS += -D AMD64
        endif
        ifeq ($(PROCESSOR_ARCHITECTURE),x86)
            CFLAGS += -D IA32
        endif
    endif
else
    FORMAT = .out
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Linux)
        CFLAGS += -D LINUX
		DEMO_LIBS += -lglut -lGLU -lGL
    endif
    ifeq ($(UNAME_S),Darwin)
        CFLAGS += -D OSX
		DEMO_LIBS += -framework GLUT -framework OpenGL -framework Cocoa
    endif
    UNAME_P := $(shell uname -p)
    ifeq ($(UNAME_P),x86_64)
        CFLAGS += -D AMD64
    endif
    ifneq ($(filter %86,$(UNAME_P)),)
        CFLAGS += -D IA32
    endif
    ifneq ($(filter arm%,$(UNAME_P)),)
        CFLAGS += -D ARM
    endif
endif


_DEPS = poissonsolver.h example.h solver_utils.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_EXAMPLE_OBJ = example.o poissonsolver.o solver_utils.o
EXAMPLE_OBJ = $(patsubst %,$(ODIR)/%,$(_EXAMPLE_OBJ))

_DEMO_OBJ = demo.o
DEMO_OBJ = $(patsubst %,$(ODIR)/%,$(_DEMO_OBJ))

_EXAMPLE = example$(FORMAT)
EXAMPLE = $(patsubst %,$(BDIR)/%,$(_EXAMPLE))


_DEMO = demo$(FORMAT)
DEMO = $(patsubst %,$(BDIR)/%,$(_DEMO))

_TEST = test$(FORMAT)
TEST = $(patsubst %,$(TDIR)/%,$(_TEST))


all: example demo test

example: $(EXAMPLE) $(SDIR)
demo: $(DEMO) $(SDIR)
test: $(TEST)

$(ODIR):
	mkdir -p $(ODIR)

$(SDIR): 
	mkdir -p $(BDIR)/solutions/gnuplot/images
	mkdir -p $(BDIR)/solutions/opengl

$(TDIR): 
	mkdir -p $(TDIR)

$(ODIR)/%.o: %.c $(DEPS) $(ODIR)
	$(CC) -c -o $@ $< $(CFLAGS)

$(EXAMPLE): $(EXAMPLE_OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(EXAMPLE_LIBS)

$(DEMO): $(DEMO_OBJ)
	$(CC)  -o $@ $^ $(CFLAGS) $(DEMO_LIBS)

$(TEST): test.c $(ODIR)/poissonsolver.o $(ODIR)/solver_utils.o $(TDIR)
	$(CC) -o $@  test.c $(ODIR)/poissonsolver.o $(ODIR)/solver_utils.o $(CFLAGS) $(EXAMPLE_LIBS)


clean:
	rm  -rf $(ODIR)
	rm -rf $(BDIR)/test
	rm -rf $(BDIR)/solutions
	rm -rf $(BDIR)/example*
	rm -rf $(BDIR)/test*
	rm -rf $(BDIR)/demo*

.PHONY: clean example demo test