CC = g++
CFLAGS = -Wall -O3 --std=c++11

all : es01_1.exe es01_2.exe es01_3.exe
es01_1.exe : es01_1.o random.o
	$(CC) random.o es01_1.o -o es01_1.exe
es01_1.o : es01_1.cpp
	$(CC) -c es01_1.cpp -o es01_1.o $(CFLAGS)
es01_2.exe : es01_2.o random.o
	$(CC) random.o es01_2.o -o es01_2.exe
es01_2.o : es01_2.cpp
	$(CC) -c es01_2.cpp -o es01_2.o $(CFLAGS)
es01_3.exe : es01_3.o random.o
	$(CC) random.o es01_3.o -o es01_3.exe
es01_3.o : es01_3.cpp
	$(CC) -c es01_3.cpp -o es01_3.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o *.exe *.out
