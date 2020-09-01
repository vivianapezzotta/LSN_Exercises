#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include "random.h"

using namespace std;

double err (vector<double>, vector<double>, int);
 
int main (){

	Random rnd;

	int M=1E4;				//n of throws
	int N=100;				//n of blocks
	int L=M/N;
	vector<double> ave(N,0), av2(N,0), sum_prog(N,0), su2_prog(N,0), err_prog(N,0);
	vector<int> x;	
	vector<double> r;			//vector of random numbers - uniform distribution
	
	for(int j=0; j<M; j++)
		r.push_back(rnd.Rannyu());
	for(int i=0; i<N; i++)
		x.push_back(i);


	//1 uniform distribution

	for (int i=0; i<N; i++){				//filling ave and av2
		double sum = 0;
		for ( int j=0; j<L; j++){
			int k = j+i*L;
			sum += M_PI/2*cos(M_PI/2*r[k]);
		}
		ave[i] = sum/L;					//average of the block i
		av2[i] = (ave[i])*(ave[i]);			//square average of the block i
	}

	for (int i=0; i<N; i++){				//filling sum_prog, su2_prog and err_prog
		for ( int j=0; j<i+1; j++){
			sum_prog[i] += ave[j];
			su2_prog[i] += av2[j];
		}
		sum_prog[i] /= (i+1);				//cumulative average
		su2_prog[i] /= (i+1);				//cumulative square average
		err_prog[i] = err(sum_prog, su2_prog, i);	//statistical uncertainty

	}

	ofstream out;						//output data
	out.open ("data_unif.out"); 
	for( int i=0; i<N; i++)
		out << x[i]*L << " " << sum_prog[i] << " " << err_prog[i] << endl;
	out.close();

	//2 importance sampling

	for(int j=0; j<M; j++)
		r[j]=1+sqrt(1-r[j]);				//from the inversion of the cumulative function
	fill(ave.begin(), ave.end(),0);
	fill(av2.begin(), av2.end(),0);
	fill(sum_prog.begin(), sum_prog.end(),0);
	fill(su2_prog.begin(), su2_prog.end(),0);
	fill(err_prog.begin(), err_prog.end(),0);

	for (int i=0; i<N; i++){
		double sum = 0;
		for ( int j=0; j<L; j++){
			int k = j+i*L;
			sum += M_PI/2.*cos(M_PI/2.*r[k])/(2*(1-r[k]));	//evaluation
		}
		ave[i] = sum/L;
		av2[i] = (ave[i])*(ave[i]);
	}

	for (int i=0; i<N; i++){
		for ( int j=0; j<i+1; j++){
			sum_prog[i] += ave[j];
			su2_prog[i] += av2[j];
		}
		sum_prog[i] /= (i+1);				//cumulative average
		su2_prog[i] /= (i+1);				//cumulative square average
		err_prog[i] = err(sum_prog, su2_prog, i);	//statistical uncertainty
	}

	out.open ("data_samp.out"); 
	for( int i=0; i<N; i++)
		out << x[i]*L << " " << sum_prog[i] << " " << err_prog[i] << endl;
	out.close();

	return 0;
}


double err ( vector<double> AV, vector<double> AV2, int n){
	if (n==0)
		return 0;
	else
		return sqrt((AV2[n]-AV[n]*AV[n])/n);
}
