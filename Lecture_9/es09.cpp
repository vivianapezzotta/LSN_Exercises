#include "random.h"
#include "class.h"
#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

int main(){

										
	int n_cities = 32;
	int n_individuals = 530;
	int n_generations = 1600;
	Random rnd;
	vector<double> x_cities, y_cities;

	cout << "The traveling salesman problem using Genetic Algorithm" <<endl <<endl;
	cout << "Number of cities: " << n_cities <<endl;
	cout << "Number of generations: " << n_generations <<endl;
	cout << "Number of individuals in each generation: " << n_individuals << endl <<endl;

//*****************************************************************************CIRCLE**************************************************************

	//Positioning cities on a circle with r=1
	cout << "Extracting cities on a circle with r=1..." << endl;
	for(int i=0; i<n_cities; i++){
		double theta = rnd.Rannyu(0, 2*M_PI);
		x_cities.push_back( cos(theta) );
		y_cities.push_back( sin(theta) );
	}
	cout << "Creating Generation 1..." <<endl;
	Generation circle( n_individuals, x_cities, y_cities );			//GENERATION 1: CREATED
	cout << "Generation 1 completed. " <<endl<<endl;

	//Evolution
	cout << "Evolution process in progress..." <<endl;
	for( int i=1; i<=n_generations; i++){
		circle.stats( "CIRCLE_fitness.dat" );
		circle.new_generation();					//EVOLUTION PROCESS
		if(i%100==0) cout<<"Generation "<<i<<endl;
	}

	cout << "Evolution completed. Printing results..." <<endl<<endl;	//PRINT ON FILE

	circle.print_best_path( "CIRCLE_bestpath.dat" );

//******************************************************************************SQUARE************************************************************

	//Positioning cities in a square of side 2
	cout << "Extracting cities in a square of side 2..." << endl;
	for(int i=0; i<n_cities; i++){
		x_cities[i] = rnd.Rannyu(-1, 1);
		y_cities[i] = rnd.Rannyu(-1, 1);
	}
	cout << "Creating Generation 1..." <<endl;
	Generation square( n_individuals, x_cities, y_cities );			//GENERATION 1: CREATED
	cout << "Generation 1 completed. " <<endl<<endl;

	//Evolution
	cout << "Evolution process in progress..." <<endl;
	for( int i=1; i<=n_generations; i++){
		square.stats( "SQUARE_fitness.dat" );
		square.new_generation();					//EVOLUTION PROCESS
		if(i%100==0) cout<<"Generation "<<i<<endl;
	}
	
	cout << "Evolution completed. Printing results..." <<endl<<endl;	//PRINT ON FILE
	
	square.print_best_path( "SQUARE_bestpath.dat" );


	return 0;
}
