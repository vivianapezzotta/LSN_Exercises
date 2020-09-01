#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include "VMC_QMparticle.h"
#include "random.h"

VMC_QMparticle :: VMC_QMparticle(int n_sim, int n_blocks, double x0){								//in: n_sim,  n_blocks, starting position x0
	m_M = n_sim;
	m_N = n_blocks;
	m_L = m_M/m_N;
	x.push_back(x0);
	for(int i=0; i< m_N; i++){
		ave.push_back(0);
		av2.push_back(0);
		sum_prog.push_back(0);
		su2_prog.push_back(0);
		err_prog.push_back(0);
	}
}	
	
VMC_QMparticle :: ~VMC_QMparticle(){};

void VMC_QMparticle :: Optimization( double startMu, double endMu, double startSigma, double endSigma, double var){		//finding optimized parameters mu and sigma for Ground State
	Random rnd;
	double bestMu, bestSigma;
        double bestE=2000;
	for(double mu=startMu; mu<=endMu; mu+=var){
       		for(double sigma=startSigma; sigma<=endSigma; sigma+=var){
       			m_acc = 0;
			double xnew = 0.;
			double x = 0.;
			double sum = 0.;
                	for(int i=0; i<m_M; i++){										//new position
                  		double xdelta = rnd.Rannyu(-m_delta, m_delta);
				xnew = x + xdelta;
                    		double A=min(1., pow(Psi(xnew,mu,sigma),2.)/pow(Psi(x,mu,sigma),2.));
                         	double r = rnd.Rannyu();
                   		if( r <= A ){
					m_acc++;
                        		x = xnew;
                   		}
                    		sum += ( -d2Psi(x,mu,sigma)*0.5+Potential(x)*Psi(x,mu,sigma))/(Psi(x,mu,sigma) );		//computing <H>
                	}
			sum /= m_M;
			cout << mu << "		" << sigma << "		" << sum << "		" << m_acc/(double)m_M << "		" << endl;
                	if(sum < bestE){
                  		bestE = sum;
                    		bestMu = mu; 
				bestSigma = sigma;
                	}
           	}
        }
	cout << endl << "Optimized parameters " << endl << "Mu: "  << bestMu << endl << "Sigma: " << bestSigma << endl <<  "Energy: " << bestE << endl;
} 	
 
void VMC_QMparticle :: Metropolis( double mu, double sigma){
	fill(ave.begin(), ave.end(),0);
	fill(av2.begin(), av2.end(),0);
	fill(sum_prog.begin(), sum_prog.end(),0);
	fill(su2_prog.begin(), su2_prog.end(),0);
	fill(err_prog.begin(), err_prog.end(),0);	
	Random rnd;
	m_acc = 0;
	double xnew = 0.;
	for(int i=0; i<m_N; i++){												//cycle on the blocks
		double sum=0;		
		for(int j=0; j<m_L; j++){											//in every block
			int k = i*m_L+j;
			double xdelta = rnd.Rannyu(-m_delta, m_delta);
			xnew = x[k] + xdelta;
			double A = min( 1., pow(Psi(xnew,mu,sigma),2.)/pow(Psi(x[k],mu,sigma),2.) );				//acceptance rate
			double r = rnd.Rannyu();
                   	if( r <= A ){
				m_acc++;
				x.push_back(xnew);
			}

			else{
				x.push_back(x[k]);
			}
			sum += ( -d2Psi(x[k+1],mu,sigma)*0.5+Potential(x[k+1])*Psi(x[k+1],mu,sigma))/(Psi(x[k+1],mu,sigma) );	//computing <H>
		}
		ave[i] = sum/m_L;			
		av2[i] = (ave[i])*(ave[i]);
	}

	for(int i=0; i<m_N; i++){
		for( int j=0; j<i+1; j++){
			sum_prog[i] += ave[j];
			su2_prog[i] += av2[j];
		}

		sum_prog[i] /= (i+1);            										//Cumulative average				
    		su2_prog[i] /= (i+1);												//Cumulative square average		
		err_prog[i] = err(sum_prog, su2_prog, i);									//statistical uncertainty			
	}
}

void VMC_QMparticle :: writeProgrBlocks( string namefile ){									//writes cumulative averages (function of n_blocks)
	ofstream out;
	out.open ( namefile ); 
	for( int i=0; i<m_N; i++)
		out << (i+1)*m_L << " " << sum_prog[i] << " " << err_prog[i] << endl;
	out.close();
}

void VMC_QMparticle :: writeInstantPoints( string namefile, int n_points ){							//writes the position of the point
	ofstream out;
	out.open ( namefile );
	for( int i=1; i<=n_points+1; i++)
		out << x[i] << endl;
	out.close();
}

double VMC_QMparticle :: Psi(double x, double mu, double sigma){								//groud state wave function
	return exp(-pow(x-mu,2.)/(2.*sigma*sigma)) + exp(-pow(x+mu,2.)/(2.*sigma*sigma));
}

double VMC_QMparticle :: d2Psi(double x, double mu, double sigma){								//compunting second derivative
	double e_neg = exp(-pow(x-mu,2.)/(2.*sigma*sigma))*(pow((x-mu)/(sigma*sigma),2.)-(1./(sigma*sigma)));
	double e_pos = exp(-pow(x+mu,2.)/(2.*sigma*sigma))*(pow((x+mu)/(sigma*sigma),2.)-(1./(sigma*sigma)));
	
	return e_neg + e_pos;
}

double VMC_QMparticle :: Potential(double x){											//potential V considered

	return pow(x,4)-2.5*pow(x,2);

}

double VMC_QMparticle :: err(vector<double> AV, vector<double> AV2, int n){							//statistic error
	if (n==0)
		return 0;
	else
		return sqrt((AV2[n]-AV[n]*AV[n])/n);
}
