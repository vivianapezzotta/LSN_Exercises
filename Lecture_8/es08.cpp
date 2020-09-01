#include <iostream>
#include "random.h"
#include "VMC_QMparticle.h"

using namespace std;

int main (){

	int M=1E6;
	int N=100;
   	VMC_QMparticle Part(M, N, 0.5);				//M, N, x0
	
	Part.setDelta( 2.65 );					//delta chosen so that acceptance rate~50%
	
	//Part.Optimization(0.78, 0.82, 0.57, 0.63, 0.004);	//startMu, endMu, startSigma, endSigma, var
	
	Part.Metropolis( 0.808, 0.618 );			//with optimized parameters
	Part.writeProgrBlocks( "Energy.dat" );			//writes energies
	Part.writeInstantPoints( "Points.dat", M );		//writes positions

	return 0;

}
