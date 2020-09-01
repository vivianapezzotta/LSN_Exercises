#include <vector>
#include <string>

using namespace std;

#ifndef __Hydrogen__
#define __Hydrogen__

class Hydrogen {

private:
	int m_M;							//n of simulations
	int m_N;							//n of blocks
	int m_L;							//simulations per block
	double m_delta;
	int m_acc;							//accumulator for acceptance rate
	vector<double> x, y, z;						//positions
	vector<double> ave, av2, sum_prog, su2_prog, err_prog;
	

protected:

	
public:
	// constructor
	Hydrogen(int, int, double, double, double);			//in: n_sim,  n_blocks, starting positions (x0, y0, z0)
	
	// destructor
	~Hydrogen();
	
	// methods
	
	void setX( double X ) {x[0] = X ;};				//set coordinates
	void setY( double Y ) {y[0] = Y ;};
	void setZ( double Z ) {z[0] = Z ;};

	void setDelta( double d ) {m_delta = d ;};
	int getAcc() {return m_acc;};
	void Equilibrate(int, int, char);
	void Metropolis(int, char);
	double Prob100 ( double, double, double);			//probability distribution for Psi{1,0,0}
	double Prob210 ( double, double, double);			//probability distribution for Psi{2,1,0}
	
	void writeInstantValues( string, int );				//used for equilibration
	void writeInstantPoints( string, int );				//used for equilibration

	void writeProgrBlocks( string );				//writes cumulative averages (function of n_blocks)
	void writePoints( string, int );				//writes a point every "step" points
	
	double err( vector<double>, vector<double>, int);
};

#endif // __Hydrogen__
