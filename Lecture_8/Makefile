CC = g++
CFLAGS = -Wall -O3 --std=c++11

es08.exe : es08.o random.o VMC_QMparticle.o
	$(CC) random.o VMC_QMparticle.o es08.o -o es08.exe
es08.o : es08.cpp
	$(CC) -c es08.cpp -o es08.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
VMC_QMparticle.o : VMC_QMparticle.cpp VMC_QMparticle.h
	$(CC) -c VMC_QMparticle.cpp -o VMC_QMparticle.o $(CFLAGS)
clean :
	rm *.o *.exe *.out
