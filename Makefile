main.exe: main.o
	g++ main.o -o main.out

main.o: main.c lud_omp.h
	g++ -c main.c -o main.o