
SRC = src
INCL = include
VPATH = src

all : test.o params.o
	gcc params.o test.o -o test -lm

test.o : test.c
	gcc -c $(SRC)/test.c -I$(INCL)

params.o : params.c
	gcc -c $(SRC)/params.c -I$(INCL)

clean : 
	rm -f *.o
