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

	int M=1E4;						//n of throws
	int N=100;						//n of blocks
	int L=M/N;
	vector<double> ave(N,0), av2(N,0), sum_prog(N,0), su2_prog(N,0), err_prog(N,0);
	vector<int> x;	
	vector<double> r;					//vector of random numbers - uniform distribution 
	
	for(int j=0; j<M; j++)
		r.push_back(rnd.Rannyu());
	for(int i=0; i<N; i++)
		x.push_back(i);


	//1 Integral mean value

	for (int i=0; i<N; i++){				//filling ave and av2
		double sum = 0;
		for ( int j=0; j<L; j++){
			int k = j+i*L;
			sum += r[k];
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
	out.open ("data_mean.out"); 
	for( int i=0; i<N; i++)
		out << (x[i]+1)*L << " " << sum_prog[i] << " " << err_prog[i] << endl;
	out.close();


	//2 Variance

	fill(ave.begin(), ave.end(),0);
	fill(av2.begin(), av2.end(),0);
	fill(sum_prog.begin(), sum_prog.end(),0);
	fill(su2_prog.begin(), su2_prog.end(),0);
	fill(err_prog.begin(), err_prog.end(),0);

	for (int i=0; i<N; i++){				//filling ave and av2
		double sum = 0;
		for (int j=0; j<L; j++){
			int k = j+i*L;
			sum += (r[k]-0.5)*(r[k]-0.5);
		}
		ave[i] = sum/L;					//average of the block i
		av2[i] = (ave[i])*(ave[i]);			//square average of the block i
	}

	for (int i=0; i<N; i++){				//filling sum_prog, su2_prog and err_prog
		for (int j=0; j<i+1; j++){
			sum_prog[i] += ave[j];
			su2_prog[i] += av2[j];
		}
		sum_prog[i] /= (i+1);				//cumulative average
		su2_prog[i] /= (i+1);				//cumulative square average
		err_prog[i] = err(sum_prog, su2_prog, i);	//statistical uncertainty
	}

	out.open ("data_var.out");				//output data
	for(int i=0; i<N; i++)
		out << (x[i]+1)*L << " " << sum_prog[i] << " " << err_prog[i] << endl;  
	out.close();


	//3 Pearson's test
	
	int throws = 1E4;	
	int sub = 100;						//number of sub-intervals
	vector<double> ni(sub,0);				//observed in every sub-interval
	vector<double> chij(100,0);				//evaluations of chi

	for(int i=0; i<100; i++){				//cycle on the single evaluation of chi
		r.clear();
		for(int j=0; j<throws; j++)
			r.push_back(rnd.Rannyu());
		sort(r.begin(),r.end());
		fill(ni.begin(), ni.end(),0.);
		int k=0;
		for(int j=0; j<sub; j++)			
			while(r[k]<((j+1.)/sub) and k<throws){	//count the random numbers in every sub-interval
				ni[j]++;
				k++;
			}
		
		for(int l=0; l<sub;l++)				//evaluate chi
			chij[i]+=pow((ni[l])-(throws/sub), 2)/(throws/sub);
	}	

	out.open ("data_chi2.out");				//output data
	for(int i=0; i<100; i++)
		out << i+1 << " " << chij[i] << endl;  
	out.close();

	return 0;
}


double err ( vector<double> AV, vector<double> AV2, int n){
	if (n==0)
		return 0;
	else
		return sqrt((AV2[n]-AV[n]*AV[n])/n);
}
