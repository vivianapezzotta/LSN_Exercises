CC = g++
CFLAGS = -Wall -O3 --std=c++11


es03_1.exe : es03_1.o random.o options.o
	$(CC) random.o options.o es03_1.o -o es03_1.exe
es03._1.o : es03_1.cpp
	$(CC) -c es03_1.cpp -o es03_1.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
options.o : options.cpp options.h 
	$(CC) -c options.cpp -o options.o $(CFLAGS)
clean :
	rm *.o *.exe *.out
