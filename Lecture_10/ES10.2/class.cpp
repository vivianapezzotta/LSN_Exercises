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



Generation :: Generation( int n_individuals, vector<double> x_cities, vector<double> y_cities, Random rnd){

	m_Ncities = x_cities.size();
	m_Nindividuals = n_individuals;
	m_citiesX = x_cities;
	m_citiesY = y_cities;
	m_rnd = rnd;

	for( int i=0; i<m_Nindividuals; i++){
		vector<int> rnd_path = random_path();					//Generating random paths to create the new individuals
		Individual new_individual( m_citiesX, m_citiesY, rnd_path );		//Creating new individuals to create the first generation
		new_individual.check();							//Checking that the individuals respect the bonds
		m_individuals.push_back( new_individual );
	}


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


void Generation :: Sort(){

	vector<double> L1;
	for( int i=0; i<m_Nindividuals; i++){
		L1.push_back( m_individuals[i].getL1() );
	}
	sort( L1.begin(), L1.end() );							//sorting L1s
	vector<Individual> appo;
	for( int i=0; i<m_Nindividuals; i++){
		bool same=false;							//the bool "same" is useful in case two or more paths share the same lenght
		for( int j=0; j<m_Nindividuals; j++){					
			if( L1[i] == m_individuals[j].getL1() and same == false ){	//using the ordered values of L1 to sort the paths (so the individuals)
				appo.push_back( m_individuals[j] );
				same = true;
			}
		}
	}

	m_individuals = appo;								//sorted individuals

}

int Generation :: selectParent(){

	return m_Nindividuals*pow( m_rnd.Rannyu(), 3 );					//using p=3

}

vector<Individual> Generation :: crossover( int pos1, int pos2 ){

	vector<int> parent1 = m_individuals[pos1].getPath();
	vector<int> parent2 = m_individuals[pos2].getPath();
	vector<int> child1;
	vector<int> child2;

	int stop = m_rnd.Rannyu(1, m_Ncities);
	for( int i=0; i<stop; i++){
		child1.push_back( parent1[i] );						//children have the same path as the parents
		child2.push_back( parent2[i] );						//until the stop position
	}	
	for( int i=0; i<m_Ncities; i++){
		for( int j=stop; j<m_Ncities; j++){
			if( parent2[i] == parent1[j] ){					//after the stop, the rest of the cities are ordered
				child1.push_back( parent2[i] );				//following the order in which they are found in the other parent
			}
			if( parent1[i] == parent2[j] ){
				child2.push_back( parent1[i] );
			}
		}
	}
	Individual c1( m_citiesX, m_citiesY, child1 );					//generating the two children and checking if they respect the bounds
	c1.check();
	Individual c2( m_citiesX, m_citiesY, child2 );
	c2.check();
	vector<Individual> appo;
	appo.push_back( c1 );
	appo.push_back( c2 );

	return appo;

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

	Sort();										//individuals are sorted so that the best ones are more likely to be chosen as parents
	vector<Individual> new_gen;
	int pos1 = 0;
	int pos2 = 0;
	for( int i=0; i<m_Nindividuals; i+=2){
		pos1 = selectParent();							//choosing the parents
		pos2 = selectParent();
		while( pos2 == pos1 ){
			pos2 = selectParent();
		}
		if( m_rnd.Rannyu() < 0.7 ){						//crossover with p(c)=0.7
			vector<Individual> children;
			children = crossover( pos1, pos2 );
			vector<Individual> mutated;
			mutated.push_back( mutate( children[0] ) );			//mutating the children
			mutated.push_back( mutate( children[1] ) );
			new_gen.push_back( mutated[0] );
			new_gen.push_back( mutated[1] );
		}
		else{
			vector<Individual> mutated;
			mutated.push_back( mutate( m_individuals[pos1] ) );		//mutating the parents	
			mutated.push_back( mutate( m_individuals[pos2] ) );
			new_gen.push_back( mutated[0] );
			new_gen.push_back( mutated[1] );
		}
	}		
	
	m_individuals = new_gen;							//new generation
	
}


void Generation :: stats( string namefile ){

	ofstream fitnessL1;
	fitnessL1.open( namefile, ios::app );
	Sort();
	double sum=0.;
	int k=0;
	for( int i=0; i<m_Nindividuals/2; i++){						//path lenght averaged on the better half of the individuals
		sum += m_individuals[i].getL1();
		k++;
	}
	sum = sum/k;
	fitnessL1 << m_individuals[0].getL1() << "   " << sum << endl;			//printing the lenght of the best path and the averaged one
	fitnessL1.close();

}

void Generation :: print_best_path( string namefile ){

	ofstream bestPath;								//printing coordinates of the best path
	bestPath.open( namefile );
	Sort();
	vector<int> path = m_individuals[0].getPath();
	for( int i=0; i<m_Ncities; i++){
		bestPath << m_citiesX[path[i]] << "   " << m_citiesY[path[i]] << endl;
	}
	bestPath << m_citiesX[path[0]] << "   " << m_citiesY[path[0]] << endl; 		//adding again the first city to return to hometown at the end
	bestPath.close();

}

vector<int> Generation :: get_path(int pos){
	Sort();
	vector<int> path=m_individuals[pos].getPath();
	return path;
}

void Generation :: set_path(int pos, vector<int> path){
	Sort();
	Individual new_path(m_citiesX, m_citiesY, path);
	m_individuals[pos]=new_path;
	Sort();
}
