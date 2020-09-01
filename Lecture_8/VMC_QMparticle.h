#include <vector>
#include <string>

using namespace std;

#ifndef __VMC_QMparticle__
#define __VMC_QMparticle__

class VMC_QMparticle {

private:
	int m_M;						//n of simulations
	int m_N;						//n of blocks
	int m_L;						//simulations per block
	double m_delta;
	vector<double> x;					//positions
	vector<double> ave, av2, sum_prog, su2_prog, err_prog;
	int m_acc;						//accumulator for acceptance rate
	
protected:


public:

	// constructors
  	VMC_QMparticle(int, int, double);			//in: n_sim,  n_blocks, starting position x0

	// destructor
        ~VMC_QMparticle();

	// methods				
 	void setDelta( double d ) {m_delta = d ;};
	void setX( double X ) {x[0] = X ;};			//sets coordinates
	int GetAcc() { return m_acc; }
	double GetH() { return sum_prog[m_N-1]; };
	void Optimization( double, double, double, double, double );
	void Metropolis( double, double );
	double Psi(double, double, double);
	double d2Psi(double, double, double);
	double Potential(double);
	double err( vector<double>, vector<double>, int);
	void writeProgrBlocks( string );			//writes cumulative averages (function of n_blocks)
	void writeInstantPoints( string, int );			//used for equilibration
};	

#endif // __VMC_QMparticle__

