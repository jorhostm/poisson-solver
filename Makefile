CC=gcc
CFLAGS=-I.
DEPS = poissonsolver.h functions.h

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

main: main.o poissonsolver.o 
	$(CC) -o main main.o poissonsolver.o -lpthread

demo: demo.c
	$(CC) -o demo demo.c -lfreeglut -lopengl32