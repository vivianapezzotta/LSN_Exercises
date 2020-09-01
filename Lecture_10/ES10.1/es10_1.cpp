#include "random.h"
#include "class.h"

#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

int main(){

	Random rnd;
	int n_cities = 32;
	int n_generations = 1E6;							//number of iterations
	double beta0 = 0.004;								//starting beta -> T0=1/beta0=250
	double scalingT_parameter = 1.015;						//scaling of T
	int cooling_step = 1E3;								//once every 1000 iterations
	vector<double> x_cities, y_cities;

	cout << "The traveling salesman problem using Simulated Annealing" <<endl <<endl;
	cout << "Number of cities to visit: " << n_cities <<endl;
	cout << "Number of generations: " << n_generations <<endl;
	cout <<  endl <<endl;


//*****************************************************************************CIRCLE**************************************************************

	//Positioning cities on a circle with r=1
	cout << "Extracting cities on a circle with r=1..." << endl;
	for(int i=0; i<n_cities; i++){
		double theta = rnd.Rannyu(0, 2*M_PI);
		x_cities.push_back( cos(theta) );
		y_cities.push_back( sin(theta) );
	}
	cout <<"Creating Generation 1..." <<endl;
	double beta = beta0;
	Generation circle( x_cities, y_cities, beta );					//GENERATION 1: CREATED
	cout << "Generation 1 completed. " <<endl<<endl;

	//Evolution
	cout << "Evolution process in progress..." <<endl;
	for( int i=1; i<=n_generations; i++){
		circle.stats( "CIRCLE_fitness.dat" );
		circle.new_generation();						//EVOLUTION PROCESS
		if(i%1000==0) cout<<"Generation "<<i<<endl;
		if(i%cooling_step==0){
			beta *= scalingT_parameter;
			cout <<"Temperature: " << 1./beta << endl;			//1E6/1E3=1000 changes of beta -> 1000 different temperatures
			circle.setBeta( beta );
		}
	}

	cout << "Evolution completed. Printing results..." <<endl<<endl;		//PRINT ON FILE

	circle.print_best_path( "CIRCLE_bestpath.dat" );

//******************************************************************************SQUARE************************************************************

	//Positioning cities in a square of side 2
	cout << "Extracting cities in a square of side 2..." << endl;	
	for(int i=0; i<n_cities; i++){
		x_cities[i] = rnd.Rannyu(-1, 1);
		y_cities[i] = rnd.Rannyu(-1, 1);
	}
	cout <<"Creating Generation 1..." <<endl;
	beta = beta0;
	Generation square( x_cities, y_cities, beta );					//GENERATION 1: CREATED
	cout << "Generation 1 completed. " <<endl<<endl;

	//Evolution
	cout << "Evolution process in progress..." <<endl;
	for( int i=1; i<=n_generations; i++){
		square.stats( "SQUARE_fitness.dat" );
		square.new_generation();						//EVOLUTION PROCESS
		if(i%1000==0) cout<<"Generation "<<i<<endl;
		if(i%cooling_step==0){
			beta *= scalingT_parameter;
			cout <<"Temperature: " << 1./beta << endl;			//1E6/1E3=1000 changes of beta -> 1000 different temperatures
			square.setBeta( beta );
		}
	}
	
	cout << "Evolution completed. Printing results..." <<endl<<endl;		//PRINT ON FILE

	square.print_best_path( "SQUARE_bestpath.dat" );


	return 0;
}
