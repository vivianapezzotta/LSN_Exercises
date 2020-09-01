CC = g++
CFLAGS = -Wall -O3 --std=c++11

es10_1.exe : es10_1.o random.o class.o
	$(CC) random.o class.o es10_1.o -o es10_1.exe
es10_1.o : es10_1.cpp
	$(CC) -c es10_1.cpp -o es10_1.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
class.o : class.cpp class.h
	$(CC) -c class.cpp -o class.o $(CFLAGS)
clean :
	rm *.o *.exe *.dat
