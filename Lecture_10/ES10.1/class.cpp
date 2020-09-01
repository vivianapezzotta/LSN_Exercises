#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>

#include "class.h"
#include "random.h"

using namespace std;

											//INDIVIDUAL
								//**********************************************************//


Individual :: Individual( vector<double> x_cities, vector<double> y_cities, vector<int> path ){
	
	m_Ncities = x_cities.size();
	m_citiesX = x_cities;
	m_citiesY = y_cities;
	m_path = path;
	m_L1 = L1();
	
}

Individual :: ~Individual(){

}

double Individual :: L1(){							//computing the lenght of the path

	double L1 = 0.;
	double dx = 0.;
	double dy = 0.;
	for( int i=0; i<m_Ncities-1; i++ ){
		dx = m_citiesX[m_path[i]] - m_citiesX[m_path[i+1]];
		dy = m_citiesY[m_path[i]] - m_citiesY[m_path[i+1]];
		L1 += sqrt( dx*dx + dy*dy );
	}
	dx = m_citiesX[m_path[m_Ncities-1]] - m_citiesX[m_path[0]];		//adding the distance from the last city to the first one (not computed in the cycle)
	dy = m_citiesY[m_path[m_Ncities-1]] - m_citiesY[m_path[0]];	
	L1 += sqrt( dx*dx + dy*dy );

	return L1;

}
		


void Individual :: check(){							//checks if every individual fulfils the bonds:
										// - first city never changes
	if (m_path[0]!=0){							// - every city is visited only once
		cerr << "La prima città è stata cambiata!" << endl;
		exit(2);
	}

	for( int i=0; i<m_Ncities; i++){
		for( int j=i+1; j<m_Ncities; j++){
			if( m_path[i] == m_path[j] ){
				cerr << "Una città viene visitata più volte!" << endl;
				exit(1);
			}
		}
	}

}


											//GENERATION
						//*****************************************************************************************//


Generation :: Generation( vector<double> x_cities, vector<double> y_cities, double initial_beta ){

	m_Ncities = x_cities.size();
	m_citiesX = x_cities;
	m_citiesY = y_cities;
	Random rnd;
	m_rnd = rnd;
	m_beta = initial_beta;

	
		vector<int> rnd_path = random_path();					//Generating a random path to create a new individual
		Individual new_individual( m_citiesX, m_citiesY, rnd_path );		//Creating a new individual
		new_individual.check();							//Checking that the individual respects the bonds							
		m_individual = new_individual;
	

}

Generation :: ~Generation(){

}

vector<int> Generation :: random_path(){

	vector<int> path;
	for( int i=0; i<m_Ncities; i++){
		path.push_back(i);							//initial settlement for the path: the first city is 0
	}
	for( int i=0; i<1E3; i++){
		path = pair_permutation( path );					//shuffling the path (not changing the first city)
	}

	return path;

}


Individual Generation :: mutate( Individual original ){

	vector<int> mutated_path;							//all the mutations have p(m)=0.07
	mutated_path = original.getPath();
	if( m_rnd.Rannyu() < 0.07 ){
		mutated_path = pair_permutation( mutated_path );
	}
	if( m_rnd.Rannyu() < 0.07 ){
		mutated_path = shift( mutated_path );
	}
	if( m_rnd.Rannyu() < 0.07 ){
		mutated_path = block_permutation( mutated_path );
	}
	if( m_rnd.Rannyu() < 0.07 ){
		mutated_path = inversion( mutated_path );
	}

	Individual mutated( m_citiesX, m_citiesY, mutated_path);
	mutated.check();								//checking after the mutations that the individual respects the bounds

	return mutated;

}

vector<int> Generation :: pair_permutation( vector<int> path ){

	int pos1 = m_rnd.Rannyu( 1, m_Ncities-1 );					//the two positions are selected so that that the first city isn't moved
	int pos2 = m_rnd.Rannyu( pos1+1, m_Ncities );					//and that they are different

	double appo = path[pos1];							//swap
	path[pos1] = path[pos2];
	path[pos2] = appo;

	return path;

}

vector<int> Generation :: shift( vector<int> path ){

	int pos1 = m_rnd.Rannyu( 1, m_Ncities-1 );					//first city does not change and moving the last one is nonsense
	int m = m_rnd.Rannyu( 1, m_Ncities-pos1);					//how many cities are shifted together
	int n = m_rnd.Rannyu( 1, m_Ncities-pos1-m+1);					//how many steps is the shift long

	vector<int> appo;
	for( int i=pos1; i<pos1+m; i++){
		appo.push_back( path[i] );						//saving the positions of the m cities
	}
	for( int i=0; i<n; i++){
		path[pos1+i] = path[pos1+i+m];						//moving n steps of the path to the left
	}
	for( int i=0; i<m; i++){
		path[pos1+n+i] = appo[i];						//adding back the m cities on the right
	}		

	return path;

}

vector<int> Generation :: block_permutation( vector<int> path ){

	int max=0;
	if( m_Ncities%2 == 0 ){								//even	
		max = (m_Ncities-2)/2;							//-1 because the first city can't be moved, then -1 again to have it even and avoid superpositions of the blocks
	}
	else{										//odd
		max = (m_Ncities-1)/2;							//-1 because the first city can't be moved
	}
	int m = m_rnd.Rannyu( 1, max+1 );						//lenght of the blocks to be moved: m <= max to avoid superpositions of the blocks
	int pos1 = m_rnd.Rannyu( 1, m_Ncities-(2*m)+1 );				//position of the first city of the first block to be moved
	int pos2 = m_rnd.Rannyu( pos1+m, m_Ncities-m+1 );				//position of the first city of the second block to be moved
	for( int i=0; i<m; i++){
		double appo = path[pos1+i];
		path[pos1+i] = path[pos2+i];
		path[pos2+i] = appo;
	}

	return path;

}

vector<int> Generation :: inversion( vector<int> path ){

	int pos1 = m_rnd.Rannyu( 1, m_Ncities-1 );					//pos1: first city of the block to invert; can't be the first one, nonsense to be the last one either
	int pos2 = m_rnd.Rannyu( pos1+1, m_Ncities );					//pos2: last city of the block to invert
	int diff = pos2-pos1;
	int i=0;
	do{
		double appo = path[pos1+i];
		path[pos1+i] = path[pos2-i];
		path[pos2-i] = appo;
		i++;
	}
	while( i<diff-1);

	return path;

}

void Generation :: new_generation(){

	Individual mutated = mutate ( m_individual );					//generating a mutating individual
	double oldBoltzmann = Boltzmann( m_individual.getL1() );			//evaluating Boltzmann for the original individual
	double newBoltzmann = Boltzmann( mutated.getL1() );				//evaluating Boltzmann for the mutated individual
	double p = newBoltzmann/oldBoltzmann;
	double r = m_rnd.Rannyu();
	if( r < p ){									//n.b.: always true if L(new)<L(old)
		m_individual = mutated;							//accepting the change with probability p
	}

}

void Generation :: setBeta( double beta ){						//setting a new value for Beta
	
	m_beta = beta;

}

double Generation :: Boltzmann( double L1 ){						//computing Boltzmann weight

	return exp( -m_beta * L1 );

}

void Generation :: stats( string namefile ){

	ofstream fitnessL1;
	fitnessL1.open( namefile, ios::app );
	fitnessL1 << m_individual.getL1() << endl;					//printing the lenght of the best path
	fitnessL1.close();

}

void Generation :: print_best_path( string namefile ){

	ofstream bestPath;								//printing coordinates of the best path
	bestPath.open( namefile );
	vector<int> path = m_individual.getPath();
	for( int i=0; i<m_Ncities; i++){
		bestPath << m_citiesX[path[i]] << "   " << m_citiesY[path[i]] << endl;
	}
	bestPath << m_citiesX[path[0]] << "   " << m_citiesY[path[0]] << endl; 		//adding again the first city to return to hometown at the end
	bestPath.close();

}
