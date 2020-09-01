CC = g++
CFLAGS = -Wall -O3 --std=c++11

es09.exe : es09.o random.o class.o
	$(CC) random.o class.o es09.o -o es09.exe
es09.o : es09.cpp
	$(CC) -c es09.cpp -o es09.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
class.o : class.cpp class.h
	$(CC) -c class.cpp -o class.o $(CFLAGS)
clean :
	rm *.o *.exe *.dat
