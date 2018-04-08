CC=mpicc
CFLAGS=-Wall -O3

all:monprogramme

monprogramme: main.c
	$(CC) $(CFLAGS) $< bmp_reader.c -o $@ $(LIB)

afficher:
	display copie.bmp

clean:
	rm -f *.o monprogramme *.so copie.bmp
