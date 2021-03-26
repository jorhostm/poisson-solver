
IDIR =../include
LDIR =../lib/freeglut/x64
CC=gcc
CFLAGS= -O3 -I$(IDIR)

ODIR=obj
BDIR=../build

MAIN_LIBS= -lpthread
DEMO_LIBS= -lfreeglut -lopengl32

_DEPS = poissonsolver.h functions.h solver_utils.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))


_MAIN_OBJ = main.o poissonsolver.o solver_utils.o
MAIN_OBJ = $(patsubst %,$(ODIR)/%,$(_MAIN_OBJ))

_DEMO_OBJ = demo.o
DEMO_OBJ = $(patsubst %,$(ODIR)/%,$(_DEMO_OBJ))

$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(BDIR)/main: $(MAIN_OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(MAIN_LIBS)

demo: $(DEMO_OBJ)
	$(CC)  -o $(BDIR)/$@ $^ $(CFLAGS) $(DEMO_LIBS)