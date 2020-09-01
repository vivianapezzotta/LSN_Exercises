#include <iostream>
#include <fstream>
#include <cmath>
#include "random.h"

using namespace std;

int main (){

	Random rnd;

	ofstream out;
	 
	int M=1E4;

	//Standard dice
	out.open ("data_unif.out");
	for(int j=0; j<M; j++)				//N=1
		out << rnd.Rannyu(1,6) << endl;
	for(int j=0; j<M; j++){				//N=2
		double sum = 0;
		for(int i=0; i<2; i++)
			sum += rnd.Rannyu(1,6);
		out << sum/2 << endl;
	}
	for(int j=0; j<M; j++){				//N=10
		double sum = 0;
		for(int i=0; i<10; i++)
			sum += rnd.Rannyu(1,6);
		out << sum/10 << endl;
	}
	for(int j=0; j<M; j++){				//N=100
		double sum = 0;
		for(int i=0; i<100; i++)
			sum += rnd.Rannyu(1,6);
		out << sum/100 << endl;
	}
	out.close();

	//Exponential dice (lambda=1)
	out.open ("data_exp.out");
	int lambda=1;
	for(int j=0; j<M; j++)				//N=1
		out << rnd.Exp(lambda) << endl;
	for(int j=0; j<M; j++){				//N=2
		double sum = 0;
		for(int i=0; i<2; i++)
			sum += rnd.Exp(lambda);
		out << sum/2 << endl;
	}
	for(int j=0; j<M; j++){				//N=10
		double sum = 0;
		for(int i=0; i<10; i++)
			sum += rnd.Exp(lambda);
		out << sum/10 << endl;
	}
	for(int j=0; j<M; j++){				//N=100
		double sum = 0;
		for(int i=0; i<100; i++)
			sum += rnd.Exp(lambda);
		out << sum/100 << endl;
	}
	out.close();

	//Lorentzian dice (mu=0, gamma=1)
	out.open ("data_lorentz.out");
	int mu=0;
	int gamma=1;
	for(int j=0; j<M; j++)				//N=1
		out << rnd.Lorentzian(gamma,mu) << endl;
	for(int j=0; j<M; j++){				//N=2
		double sum = 0;
		for(int i=0; i<2; i++)
			sum += rnd.Lorentzian(gamma,mu);
		out << sum/2 << endl;
	}
	for(int j=0; j<M; j++){				//N=10
		double sum = 0;
		for(int i=0; i<10; i++)
			sum += rnd.Lorentzian(gamma,mu);
		out << sum/10 << endl;
	}
	for(int j=0; j<M; j++){				//N=100
		double sum = 0;
		for(int i=0; i<100; i++)
			sum += rnd.Lorentzian(gamma,mu);
		out << sum/100 << endl;
	}

	out.close();

}
