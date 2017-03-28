default:main

CC = gcc
CFLAGS = -O3 -std=c99
CLIBS = -lm 

main:main.o jacobi.o
	$(CC) $(CLIBS) main.o jacobi.o -o main

main.o:main.c
	$(CC) $(CFLAGS) -c main.c -o main.o

jacobi.o:jacobi.c
	$(CC) $(CFLAGS) -c jacobi.c -o jacobi.o


clean:
	rm -rf ./*.o main

run:main
	./main
