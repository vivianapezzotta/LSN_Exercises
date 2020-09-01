#include <iostream>
#include "random.h"
#include "hydrogen.h"

using namespace std;
 
int main (){
	
	Hydrogen H(1E6, 100, 0., 0., 0.);		//M, N, x0, y0, z0
	
	int trials = 1E3;				//trials necessary to equilibrate the system
							//and find deltas so that acceptance rate ~ 0.5

//Psi {1,0,0} and uniform distribution

	H.setDelta( 1.27 );								
	H.Equilibrate( trials, 100, 'u');
	//cout << "Acceptance rate: " << H.getAcc() / (double)trials << endl << endl;
	//H.writeInstantValues( "EQvalues100Unif.out", trials );					//used for equilibration
	//H.writeInstantValues( "AutocorrValues100Unif.out", 1E5 );
	//H.writeInstantPoints( "EQpoints100Unif.out", trials );					//used for equilibration
	H.Metropolis( 100, 'u' );
	H.writeProgrBlocks( "mean100Unif.out" );
	H.writePoints( "points100Unif.out", 100 );
	

//Psi {1,0,0} and gaussian distribution

	H.setDelta( 0.76 );								
	H.Equilibrate( trials, 100, 'g');
	//cout << "Acceptance rate: " << H.getAcc() / (double)trials << endl << endl;
	//H.writeInstantValues( "EQvalues100Gauss.out", trials );
	//H.writeInstantPoints( "EQpoints100Gauss.out", trials );
	H.Metropolis( 100, 'g' );
	H.writeProgrBlocks( "mean100Gauss.out" );
	H.writePoints( "points100Gauss.out", 100 );

//Psi {2,1,0} and uniform distribution

	H.setX(0.);
	H.setY(0.);
	H.setZ(2.);
	H.setDelta( 2.95 );								
	H.Equilibrate( trials, 210, 'u');
	//cout << "Acceptance rate: " << H.getAcc() / (double)trials << endl << endl;
	//H.writeInstantValues( "EQvalues210Unif.out", trials );
	//H.writeInstantPoints( "EQpoints210Unif.out", trials );
	H.Metropolis( 210, 'u' );
	H.writeProgrBlocks( "mean210Unif.out" );
	H.writePoints( "points210Unif.out", 100 );
	
//Psi {2,1,0} and gaussian distribution

	H.setDelta( 1.89 );								
	H.Equilibrate( trials, 210, 'g');
	//cout << "Acceptance rate: " << H.getAcc() / (double)trials << endl << endl;
	//H.writeInstantValues( "EQvalues210Gauss.out", trials );
	//H.writeInstantPoints( "EQpoints210Gauss.out", trials );
	H.Metropolis( 210, 'g' );
	H.writeProgrBlocks( "mean210Gauss.out" );
	H.writePoints( "points210Gauss.out", 100 );
	
	
	return 0;
}
