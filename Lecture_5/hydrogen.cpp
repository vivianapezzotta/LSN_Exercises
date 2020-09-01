#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "hydrogen.h"
#include "random.h"


Hydrogen :: Hydrogen(int n_sim, int n_blocks, double x0, double y0, double z0){
	m_M = n_sim;
	m_N = n_blocks;
	m_L = m_M/m_N;
	
	for (int i=0; i<m_M; i++){
		x.push_back(x0);				//sets initial positions
		y.push_back(y0);	
		z.push_back(z0);
	}
	
	for(int i=0; i<m_N; i++){
		ave.push_back(0);
		av2.push_back(0);
		sum_prog.push_back(0);
		su2_prog.push_back(0);
		err_prog.push_back(0);
	}
}

Hydrogen :: ~Hydrogen(){}

void Hydrogen :: Equilibrate( int steps, int state, char distrib ){		//equilibration: 'steps' simulations to let the system equilibrate
										//used also to check on acceptance rate
	m_acc = 0;
	double xnew = 0;
	double ynew = 0;
	double znew = 0;
	Random rnd;
	
	for(int i=0; i<steps; i++){
		
		if( distrib == 'u' ){						//uniform distribution
			
			double xdelta = rnd.Rannyu(-m_delta, m_delta);
			double ydelta = rnd.Rannyu(-m_delta, m_delta);
			double zdelta = rnd.Rannyu(-m_delta, m_delta);
			xnew = x[i] + xdelta;
			ynew = y[i] + ydelta;
			znew = z[i] + zdelta;
		}
		if( distrib == 'g' ){						//gaussian distribution
			
			xnew = rnd.Gauss(x[i], m_delta);
			ynew = rnd.Gauss(y[i], m_delta);
			znew = rnd.Gauss(z[i], m_delta);
		}
		
		double AccRate=0;
		if( state == 100 )
			AccRate = min( 1., Prob100(xnew, ynew, znew)/Prob100(x[i], y[i], z[i]) );	//acceptance rate for Psi{1,0,0}
		if( state == 210 )
			AccRate = min( 1., Prob210(xnew, ynew, znew)/Prob210(x[i], y[i], z[i]) );	//acceptance rate for Psi{2,1,0}
		double r = rnd.Rannyu();
		if ( r <= AccRate ){						//accepted
			m_acc++;						
			x[i+1] = xnew;
			y[i+1] = ynew;
			z[i+1] = znew;
		}
		else{								//rejected
				x[i+1] = x[i];
				y[i+1] = y[i];
				z[i+1] = z[i];
			}
	}
	
	setX( x[steps-1] );		//updates starting position with the result of equilibration
	setY( y[steps-1] );		
	setZ( z[steps-1] );
}

void Hydrogen :: Metropolis( int state, char distrib ){
	fill(ave.begin(), ave.end(),0);
	fill(av2.begin(), av2.end(),0);
	fill(sum_prog.begin(), sum_prog.end(),0);
	fill(su2_prog.begin(), su2_prog.end(),0);
	fill(err_prog.begin(), err_prog.end(),0);
	Random rnd;
	double xnew = 0;
	double ynew = 0;
	double znew = 0;
	
	for (int i=0; i<m_N; i++){
		double sum = 0;
		for ( int j=0; j<m_L; j++){
			int k = i*m_L+j;
			if( distrib == 'u' ){					//uniform distribution
				
				double xdelta = rnd.Rannyu(-m_delta, m_delta);	//uniform sampling
				double ydelta = rnd.Rannyu(-m_delta, m_delta);
				double zdelta = rnd.Rannyu(-m_delta, m_delta);
				xnew = x[k] + xdelta; 				//evaluate new positions
				ynew = y[k] + ydelta; 	
				znew = z[k] + zdelta;
			}
			if( distrib == 'g' ){					//gaussian distribution
				
				xnew = rnd.Gauss(x[k], m_delta);		//gaussian sampling: evaluate new positions
				ynew = rnd.Gauss(y[k], m_delta);
				znew = rnd.Gauss(z[k], m_delta);
			}
 			double AccRate=0;
			if( state == 100 )
 				AccRate = min( 1., Prob100(xnew, ynew, znew)/Prob100(x[k], y[k], z[k]) );	//acceptance rate for Psi{1,0,0}
			if( state == 210 )
 				AccRate = min( 1., Prob210(xnew, ynew, znew)/Prob210(x[k], y[k], z[k]) );	//acceptance rate for Psi{2,1,0}
			double r = rnd.Rannyu();
			if ( r <= AccRate ){									//accepted
				x[k+1] = xnew;
				y[k+1] = ynew;
				z[k+1] = znew;
			}
			else{											//rejected
				x[k+1] = x[k];
				y[k+1] = y[k];
				z[k+1] = z[k];
			}
			sum += sqrt( x[k]*x[k]+y[k]*y[k]+z[k]*z[k]);
		}
		ave[i] = sum/m_L;
		av2[i] = (ave[i])*(ave[i]);
		//cout << ave[i] << "  " << av2[i] << endl;
	}

	for (int i=0; i<m_N; i++){
		for ( int j=0; j<i+1; j++){
			sum_prog[i] += ave[j];
			su2_prog[i] += av2[j];
		}
		sum_prog[i] /= (i+1);				//cumulative average
		su2_prog[i] /= (i+1);				//cumulative square average
		err_prog[i] = err(sum_prog, su2_prog, i);	//statistical uncertainty
		//cout << sum_prog[i] << "  " << 	su2_prog[i] << "  " << err_prog[i] << endl;
	}
}

void Hydrogen :: writeInstantValues( string namefile, int n_points ){		//used for equilibration and autocorrelation
	ofstream out;
	out.open ( namefile ); 
	for( int i=1; i<n_points+1; i++)
		out << i << " " << sqrt(x[i]*x[i]+y[i]*y[i]+z[i]*z[i]) << endl;
	out.close();
}

void Hydrogen :: writeInstantPoints( string namefile, int n_points ){		//used for equilibration
	ofstream out;
	out.open ( namefile );
	for( int i=1; i<n_points+1; i++)
		out << x[i] << " " << y[i] << " " << z[i] << endl;
	out.close();
}
void Hydrogen :: writeProgrBlocks( string namefile ){				//writes cumulative averages (function of n_blocks)
	ofstream out;
	out.open ( namefile ); 
	for( int i=0; i<m_N; i++)
		out << (i+1)*m_L << " " << sum_prog[i] << " " << err_prog[i] << endl;
	out.close();
}

void Hydrogen :: writePoints( string namefile, int step ){
	ofstream out;
	out.open ( namefile );
	for( int i=0; i<m_M; i++){
		if( i%step == 0 ){						//writes a point every "step" points
			out << x[i] << " " << y[i] << " " << z[i] << endl;
		}
	}
	out.close();
}


double Hydrogen :: Prob100 ( double x, double y, double z){			//probability distribution for Psi{1,0,0}
	return exp(-2*sqrt(x*x+y*y+z*z)) / M_PI;
}

double Hydrogen :: Prob210 ( double x, double y, double z){			//probability distribution for Psi{2,1,0}
	return z*z*exp(-sqrt(x*x+y*y+z*z)) / (32*M_PI);
}

double Hydrogen :: err(vector<double> AV, vector<double> AV2, int n){		//error
	if (n==0)
		return 0;
	else
		return sqrt((AV2[n]-AV[n]*AV[n])/n);
}
