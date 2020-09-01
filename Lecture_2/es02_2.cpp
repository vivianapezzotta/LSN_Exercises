#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include "random.h"

using namespace std;

double err (double AV, double AV2, int n);

int main (){

	Random rnd;

	int M=1E4;								//n of throws
	int N=100;								//n of blocks
	int L=M/N;
	int steps=100;								//n steps
	

	double dist[steps+1], sum[steps+1];
	double ave[N][steps+1], av2[N][steps+1];				// first index: BLOCK, second index: STEP
	double sum_prog[steps+1], su2_prog[steps+1], err_prog[steps+1];

	//1 Discrete

	for (int k=0; k<N; k++){						//cycle on N blocks
		for (int i=0; i<steps+1; i++){
			sum[i] = 0;
		}
		for ( int j=0; j<L; j++){					//cycle on L simulations inside every block
			int x=0;						//positioning the walker in the origin
			int y=0;
			int z=0;
			
			dist[0]=0;
		
			for(int i=1; i<steps+1; i++){				//moving the walker randomly
				int ran = rnd.Rannyu(1,7);
				if(ran==1)
					x++;
				else if(ran==2)
					x--;
				else if(ran==3)
					y++;
				else if(ran==4)
					y--;
				else if(ran==5)
					z++;
				else
					z--;
				dist[i] = (x*x+y*y+z*z);			//evaluate square distance
			}
			for (int i=0; i<steps+1; i++){
				sum[i] += dist[i];				//sum over the L simulations in every block for every step i
			}
		}
		for( int i=0; i<steps+1; i++){					//evaluate the average square distance after i steps in the k block
			ave[k][i] = sum[i]/L;
			av2[k][i] = (ave[k][i])*(ave[k][i]);
		}
	}

	for ( int j=0; j<N; j++){						//sum on the N blocks for every step i
		for( int i=0; i<steps+1; i++){
			sum_prog[i] += ave[j][i];				//contains square distance for every step (counting all the blocks)
			su2_prog[i] += av2[j][i];
		}
	}
	for( int i=0; i<steps+1; i++){						//for every step i:
		sum_prog[i] /= N;						//cumulative average
		su2_prog[i] /= N;						//cumulative square average
		err_prog[i] = err(sum_prog[i], su2_prog[i], N);			//statistical uncertainty on square distance
	}
	
	ofstream out;								//output data
	out.open ("data_discreteRW.out"); 
	for( int i=0; i<steps+1; i++)						//writing distance and uncertainty derived from propagation of errors
		out << i << " " << sqrt(sum_prog[i]) << " " << 0.5/sqrt(sum_prog[i])*err_prog[i] <<  endl;
	out.close();

	//2 Continuum

	for( int i=0; i<steps+1; i++){						//erase prevoius data and set to 0 the vectors
		for( int k=0; k<N; k++){
			ave[k][i]=0.;
			av2[k][i]=0.;
		}
		sum_prog[i]=0.;
		su2_prog[i]=0.;
		err_prog[i]=0.;
	}

	for (int k=0; k<N; k++){						//cycle on N blocks
		for (int i=0; i<steps+1; i++){
			sum[i] = 0;
		}
		for ( int j=0; j<L; j++){					//cycle on L simulations inside every block
			double x=0.;						//positioning the walker in the origin
			double y=0.;
			double z=0.;
			
			dist[0]=0;

			for(int i=1; i<steps+1; i++){				//moving the walker randomly
				double x_rnd=0;
				double y_rnd=0;
				double z_rnd=0;

				double norm=0;
				do{						//finding a point (randomly extracted in the cube [-1,1]^3
					x_rnd=rnd.Rannyu(-1,1);			//whose distance from the origin is <1
					y_rnd=rnd.Rannyu(-1,1);
					z_rnd=rnd.Rannyu(-1,1);
					norm=sqrt(x_rnd*x_rnd+y_rnd*y_rnd+z_rnd*z_rnd);
				}while(norm>1);

										
				x_rnd=x_rnd/norm;				//normalizing so that the point belongs to the unitary sphere
				y_rnd=y_rnd/norm;
				z_rnd=z_rnd/norm;
				x += x_rnd;					//summing on the steps 
				y += y_rnd;
				z += z_rnd;
				
				dist[i]= (x*x+y*y+z*z);				//evaluate square distance
			}
			for (int i=0; i<steps+1; i++){
				sum[i] += dist[i];				//sum over the L simulations in every block for every step i
			}
		}
		for( int i=0; i<steps+1; i++){					//evaluate the average square distance after i steps in the k block
			ave[k][i] = sum[i]/L;
			av2[k][i] = (ave[k][i])*(ave[k][i]);
		}
	}

	for ( int j=0; j<N; j++){						//sum on the N blocks for every step i
		for( int i=0; i<steps+1; i++){
			sum_prog[i] += ave[j][i];				//contains square distance for every step (counting all the blocks)
			su2_prog[i] += av2[j][i];
		}
	}
	for( int i=0; i<steps+1; i++){
		sum_prog[i] /= N;						//cumulative average
		su2_prog[i] /= N;						//cumulative square average
		err_prog[i] = err(sum_prog[i], su2_prog[i], N);			//statistical uncertainty on square distance
	}
	
						
	out.open ("data_continuumRW.out"); 					//output data
	
	for( int i=0; i<steps+1; i++){						//writing distance and uncertainty derived from propagation of errors
		
		out << i << " " << sqrt(sum_prog[i]) << " " << 0.5/sqrt(sum_prog[i])*err_prog[i] <<  endl;
	}
	out.close();

	return 0;
}

double err (double AV, double AV2, int n){
	if (n==0)
		return 0;
	else
		return sqrt((AV2-AV*AV)/n);
}
