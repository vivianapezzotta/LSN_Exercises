CC = g++
CFLAGS = -Wall -O3 --std=c++11

es05_1.exe : es05_1.o random.o hydrogen.o
	$(CC) random.o hydrogen.o es05_1.o -o es05_1.exe
es05_1.o : es05_1.cpp
	$(CC) -c es05_1.cpp -o es05_1.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
hydrogen.o : hydrogen.cpp hydrogen.h
	$(CC) -c hydrogen.cpp -o hydrogen.o $(CFLAGS)
clean :
	rm *.o *.exe *.out
