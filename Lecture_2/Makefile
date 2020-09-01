CC = g++
CFLAGS = -Wall -O3 --std=c++11

all : es02_1.exe es02_2.exe
es02_1.exe : es02_1.o random.o
	$(CC) random.o es02_1.o -o es02_1.exe
es02_1.o : es02_1.cpp
	$(CC) -c es02_1.cpp -o es02_1.o $(CFLAGS)
es02_2.exe : es02_2.o random.o
	$(CC) random.o es02_2.o -o es02_2.exe
es02_2.o : es02_2.cpp
	$(CC) -c es02_2.cpp -o es02_2.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o *.exe *.out
