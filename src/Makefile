all: example demo

IDIR =../include
LDIR =../lib
CC=gcc
CFLAGS= -O3 -I$(IDIR)

BDIR=../build
ODIR=obj

EXAMPLE_LIBS= -lm -lpthread
DEMO_LIBS= -lglut -lGL

_DEPS = poissonsolver.h example.h solver_utils.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_EXAMPLE_OBJ = example.o poissonsolver.o solver_utils.o
EXAMPLE_OBJ = $(patsubst %,$(ODIR)/%,$(_EXAMPLE_OBJ))

_DEMO_OBJ = demo.o
DEMO_OBJ = $(patsubst %,$(ODIR)/%,$(_DEMO_OBJ))

_EXAMPLE = example
EXAMPLE = $(patsubst %,$(BDIR)/%,$(_EXAMPLE))

_DEMO = demo
DEMO = $(patsubst %,$(BDIR)/%,$(_DEMO))

$(ODIR):
	mkdir -p $(ODIR)

$(ODIR)/%.o: %.c $(DEPS) $(ODIR)
	$(CC) -c -o $@ $< $(CFLAGS)

$(EXAMPLE): $(EXAMPLE_OBJ)
	mkdir -p $(BDIR)/solutions/gnuplot/images
	mkdir -p $(BDIR)/solutions/opengl
	$(CC) -o $@ $^ $(CFLAGS) $(EXAMPLE_LIBS)

$(DEMO): $(DEMO_OBJ)
	$(CC)  -o $@ $^ $(CFLAGS) $(DEMO_LIBS)

example: $(EXAMPLE)
demo: $(DEMO)



clean:
	rm  -rf $(ODIR)
	
	rm -rf $(BDIR)/solutions
	rm -rf $(BDIR)/example*
	rm -rf $(BDIR)/demo*

.PHONY: clean demo example