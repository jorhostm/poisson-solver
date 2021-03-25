
IDIR =../include
LDIR =../lib/freeglut/x64
CC=gcc
CFLAGS_MAIN= -O3 -I$(IDIR)
CFLAGS_DEMO= -O3 -I$(IDIR) -L$(LDIR) -Wl,--subsystem,windows

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
	$(CC) -c -o $@ $< $(CFLAGS_MAIN)


$(BDIR)/main: $(MAIN_OBJ)
	$(CC) -o $@ $^ $(CFLAGS_MAIN) $(MAIN_LIBS)

$(BDIR)/demo: $(DEMO_OBJ)
	gcc -o $(BDIR)/demo.exe $(DEMO_OBJ) -L"..\lib\freeglut\x64" -lfreeglut -lopengl32 -Wl,--subsystem,windows