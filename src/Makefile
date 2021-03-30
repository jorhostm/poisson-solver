

IDIR =../include
LDIR =../lib
CC=gcc
CFLAGS= -O3 -I$(IDIR)

BDIR=../build
ODIR=obj

MAIN_LIBS= -lm -lpthread
DEMO_LIBS= -lfreeglut -lopengl32

_DEPS = poissonsolver.h functions.h solver_utils.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_MAIN_OBJ = main.o poissonsolver.o solver_utils.o
MAIN_OBJ = $(patsubst %,$(ODIR)/%,$(_MAIN_OBJ))

_DEMO_OBJ = demo.o
DEMO_OBJ = $(patsubst %,$(ODIR)/%,$(_DEMO_OBJ))

_MAIN = main
MAIN = $(patsubst %,$(BDIR)/%,$(_MAIN))

_DEMO = demo
DEMO = $(patsubst %,$(BDIR)/%,$(_DEMO))


$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(MAIN): $(MAIN_OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(MAIN_LIBS)

$(DEMO): $(DEMO_OBJ)
	$(CC)  -o $@ $^ $(CFLAGS) $(DEMO_LIBS)

main: $(MAIN)
demo: $(DEMO)

clean:
	rm  -f $(MAIN_OBJ) $(DEMO_OBJ) $(MAIN) $(DEMO)
	
	rm -rf $(BDIR)/solutions

.PHONY: clean main demo