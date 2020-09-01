#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include "random.h"
#include "options.h"

using namespace std;



int main (){

	int M=1E5;
	int N=100;
	int L=M/N;

	European eur(M, N, 100, 100, 1, 0.1, 0.25);
	
	eur.Call(1);									//Direct Call
	ofstream out;
	out.open("DirectCall.out");							

	for(int i=0; i<N; i++)								//print results
		out << i*L << " " << eur.avg[i] << " " << eur.error[i] << endl;
	out.close();


	
	eur.Put(1);									//Direct Put
	out.open("DirectPut.out");

	for(int i=0; i<N; i++)								//print results
		out << i*L << " " << eur.avg[i] << " " << eur.error[i] << endl;
	out.close();


	eur.Call(100);									//Discrete Call
	out.open("DiscreteCall.out");

	for(int i=0; i<N; i++)								//print results
		out << i*L << " " << eur.avg[i] << " " << eur.error[i] << endl;
	out.close();


	eur.Put(100);									//Discrete Put
	out.open("DiscretePut.out");

	for(int i=0; i<N; i++)								//print results
		out << i*L << " " << eur.avg[i] << " " << eur.error[i] << endl;
	out.close();


return 0;

}

