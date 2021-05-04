
IDIR =../include
LDIR =../lib
CC=gcc
CFLAGS= -O3 -I$(IDIR)

BDIR=../build
ODIR=obj

EXAMPLE_LIBS= -lm -lpthread
DEMO_LIBS=

ifeq ($(OS),Windows_NT)
    CCFLAGS += -D WIN32
	DEMO_LIBS += -lglut32 -lglu32 -lopengl32
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
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Linux)
        CFLAGS += -D LINUX
		DEMO_LIBS += -lfreeglut -lGLU -lGL
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