CC = g++
CFLAGS = -Wall -O3 --std=c++11

es04_1.exe : MolDyn_NVE.o
	$(CC) MolDyn_NVE.o -o es04_1.exe
es04_1.o : MolDyn_NVE.cpp MolDyn_NVE.h
	$(CC) -c MolDyn_NVE.cpp -o MolDyn_NVE.o $(CFLAGS)
clean :
	rm *.o *.exe
