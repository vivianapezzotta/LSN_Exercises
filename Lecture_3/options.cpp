#include <iostream>
#include <cmath>
#include "random.h"
#include "options.h"

using namespace std;

European :: European (int n_simul, int n_blocks, double asset, double strike, double deliv_time, double interest, double vol){

	M = n_simul;
	N = n_blocks;
	S = asset;
	K = strike;
	T = deliv_time;
	r = interest;
	sigma = vol;
}

European :: ~European(){}


void European :: Call(int step){

	avg.clear();
	error.clear();

	int L=M/N;
	
	vector<double> ave(N,0), av2(N,0), sum_prog(N,0), su2_prog(N,0), err_prog(N,0);
	Random rnd;
	
	
	for (int i=0; i<N; i++){											//cycle on N blocks
		double sum = 0;
		for ( int j=0; j<L; j++){										//cycle on L evaluations inside every block
			double s=S;
			for(int t=0; t<(int)step; t++){									//cycle on the number of steps
					
				s *= exp( (r-0.5*sigma*sigma)*((double)T/step) + sigma*rnd.Gauss(0,1)*sqrt((double)T/step) );		//recorsivity for s
			}	
			sum += exp(-r*T)*max(0., s-K);									//considering both the interest (exp(-r*T)) and the profit (max(0, s-K))
		}
		ave[i] = sum/L;
		av2[i] = (ave[i])*(ave[i]);
	}

	for (int i=0; i<N; i++){
		for ( int j=0; j<i+1; j++){
			sum_prog[i] += ave[j];
			su2_prog[i] += av2[j];
		}
		sum_prog[i] /= (i+1);											//cumulative average
		su2_prog[i] /= (i+1);											//cumulative square average
		err_prog[i] = err(sum_prog, su2_prog, i);								//statistical uncertainty
	}

	for (int i=0; i<N; i++){
		avg.push_back(sum_prog[i]);
		error.push_back(err_prog[i]);
	}
}



void European :: Put(int step){

	avg.clear();
	error.clear();

	int L=M/N;
	
	vector<double> ave(N,0), av2(N,0), sum_prog(N,0), su2_prog(N,0), err_prog(N,0);
	Random rnd;
	
	
	for (int i=0; i<N; i++){											//cycle on N blocks
		double sum = 0;
		for ( int j=0; j<L; j++){										//cycle on L evaluations inside every block
			double s=S;
			for(int t=0; t<(int)step; t++){									//cycle on the number of steps
							
				s *= exp( (r-0.5*sigma*sigma)*((double)T/step) + sigma*rnd.Gauss(0,1)*sqrt((double)T/step) );		//recorsivity for s
			}	
			sum += exp(-r*T)*max(0., K-s);									//considering both the interest (exp(-r*T)) and the profit (max(0, s-K))
		}
		ave[i] = sum/L;
		av2[i] = (ave[i])*(ave[i]);
	}

	for (int i=0; i<N; i++){
		for ( int j=0; j<i+1; j++){
			sum_prog[i] += ave[j];
			su2_prog[i] += av2[j];
		}
		sum_prog[i] /= (i+1);											//cumulative average
		su2_prog[i] /= (i+1);											//cumulative square average
		err_prog[i] = err(sum_prog, su2_prog, i);								//statistical uncertainty
	}

	for (int i=0; i<N; i++){
		avg.push_back(sum_prog[i]);
		error.push_back(err_prog[i]);
	}
}


double European :: err (vector<double> AV, vector<double> AV2, int n){
	
	if (n==0)
		return 0;
	else
		return sqrt((AV2[n]-AV[n]*AV[n])/n);
}
