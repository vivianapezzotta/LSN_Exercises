#include "random.h"
#include "class.h"
#include "mpi.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

int main(int argc, char* argv[]){

	MPI_Init(&argc,&argv);
	int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	Random rnd;
	int seed[4];
	int p[size][2];
	ifstream Primes("Primes");
	if (Primes.is_open()){
		for(int i=0;i<size;i++){
			Primes >> p[i][0] >> p[i][1];
		}
	} else cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes.close();

	ifstream input("seed.in");
	string property;
	if (input.is_open()){
		while ( !input.eof() ){
			input >> property;
			if( property == "RANDOMSEED" ){
				input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
				rnd.SetRandom(seed,p[rank][0],p[rank][1]);
			}
		}
		input.close();
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;

	int n_cities = 32;
	int n_individuals = 530;
	int n_generations = 1600;
	int migration_rate = 10;			//1600/10=160 times they communicate
	vector<double> x_cities, y_cities;
	double* x_c = new double[n_cities];
	double* y_c = new double[n_cities];

	if(rank==0){
		cout << "The traveling salesman problem using Genetic Algorithm in Parallel Programming" <<endl <<endl;
		cout << "Number of cities: " << n_cities <<endl;
		cout << "Number of generations: " << n_generations <<endl;
		cout << "Number of individuals in each generation: " << n_individuals << endl <<endl;
	}

	//Positioning cities on a circle with r=1
	if(rank==0) cout << "Extracting cities on a circle with r=1..." << endl;
	for(int i=0; i<n_cities; i++){
		double theta = rnd.Rannyu(0, 2*M_PI);
		x_c[i] = cos(theta);
		y_c[i] = sin(theta);
	}
	MPI_Bcast( x_c,n_cities,MPI_DOUBLE,0, MPI_COMM_WORLD);
	MPI_Bcast( y_c,n_cities,MPI_DOUBLE,0, MPI_COMM_WORLD);
	if(rank==0) cout << "Creating Generation 1..." <<endl;
	for(int i=0; i<n_cities; i++){
		x_cities.push_back( x_c[i] );
		y_cities.push_back( y_c[i] );
	}
	for(int i=0; i<n_cities; i++){					//generating the cities inside the square
		x_c[i] = rnd.Rannyu(-1, 1);				//the same generated in the other exercises
		y_c[i] = rnd.Rannyu(-1, 1);				
	}
	Generation circle( n_individuals, x_cities, y_cities, rnd );	//GENERATION 1: CREATED
	
	//Evolution
	if(rank==0) cout << "Evolution process in progress..." <<endl;
	for( int i=1; i<=n_generations; i++){
		circle.stats( to_string(rank)+"_CIRCLE_fitness.dat" );
		circle.new_generation();				//EVOLUTION PROCESS
		if(i!=0 and i%migration_rate==0){			//communication
			MPI_Status stat1, stat2, stat3, stat4;
			MPI_Request req1, req2;
			int root=0;
			int itag1=1, itag2=2, itag3=3, itag4=4;
			int target[3];
			if(rank==0){
				target[0]=(int)rnd.Rannyu(1,4);
				do{
					target[1]=(int)rnd.Rannyu(1,4);
				}while(target[1]==0 or target[1]==target[0]);
				do{
					target[2]=(int)rnd.Rannyu(1,4);
				}while(target[2]==0 or target[2]==target[0] or target[2]==target[1]);
			}
			MPI_Bcast(target, 3, MPI_INT, root, MPI_COMM_WORLD);
			
			vector<int> path_exc=circle.get_path(0);
			vector<int> new_path(n_cities);
			int* path1=new int[n_cities];
			int* path2=new int[n_cities];
			for(int c=0;c<n_cities;c++){
				path1[c]=path_exc[c];
			}
			MPI_Barrier(MPI_COMM_WORLD);
			
			if(rank==0){
				MPI_Isend(&path1[0], n_cities, MPI_INTEGER, target[0], itag1, MPI_COMM_WORLD, &req1);
				MPI_Recv(&path2[0], n_cities, MPI_INTEGER, target[0], itag2, MPI_COMM_WORLD, &stat2);
			}
			else if(rank==target[0]){
				MPI_Send(&path1[0], n_cities, MPI_INTEGER, 0, itag2, MPI_COMM_WORLD);
				MPI_Recv(&path2[0], n_cities, MPI_INTEGER, 0, itag1, MPI_COMM_WORLD, &stat1);
			}
			else if(rank==target[1]){
				MPI_Isend(&path1[0], n_cities, MPI_INTEGER, target[2], itag3, MPI_COMM_WORLD, &req2);
				MPI_Recv(&path2[0], n_cities, MPI_INTEGER, target[2], itag4, MPI_COMM_WORLD, &stat4);
			}
			else if(rank==target[2]){
				MPI_Send(&path1[0], n_cities, MPI_INTEGER, target[1], itag4, MPI_COMM_WORLD);
				MPI_Recv(&path2[0], n_cities, MPI_INTEGER, target[1], itag3, MPI_COMM_WORLD, &stat3);
			}
			for(int c=0;c<n_cities;c++){
				new_path[c]=path2[c];
			}
			
			circle.set_path(0, new_path);
			MPI_Barrier(MPI_COMM_WORLD);
			if(rank==0) cout<<" generation "<<i<<" migration completed"<<endl;
		}
	}

	if(rank==0) cout << "Evolution completed. Printing results..." <<endl<<endl;	//PRINT ON FILE

	circle.print_best_path( to_string(rank)+"_CIRCLE_bestpath.dat" );

	//******************************************************************************SQUARE************************************************************

	//Positioning cities in a square of side 2
	if(rank==0) cout << "Extracting cities in a square of side 2..." << endl;

	MPI_Bcast( x_c,n_cities,MPI_DOUBLE,0, MPI_COMM_WORLD);
	MPI_Bcast( y_c,n_cities,MPI_DOUBLE,0, MPI_COMM_WORLD);
	if(rank==0) cout << "Creating Generation 1..." <<endl;
	for(int i=0; i<n_cities; i++){
		x_cities[i] = x_c[i];
		y_cities[i] = y_c[i];
	}
	Generation square( n_individuals, x_cities, y_cities, rnd );		//GENERATION 1: CREATED

	//Evolution
	if(rank==0) cout << "Evolution process in progress..." <<endl;
	for( int i=1; i<=n_generations; i++){
		square.stats( to_string(rank)+"_SQUARE_fitness.dat" );
		square.new_generation();					//EVOLUTION PROCESS
		if(i!=0 and i%migration_rate==0){				//communication
			MPI_Status stat1, stat2, stat3, stat4;
			MPI_Request req1, req2;
			int root=0;
			int itag1=1, itag2=2, itag3=3, itag4=4;
			int target[3];
			if(rank==0){
				target[0]=(int)rnd.Rannyu(1,4);
				do{
					target[1]=(int)rnd.Rannyu(1,4);
				}while(target[1]==0 or target[1]==target[0]);
				do{
					target[2]=(int)rnd.Rannyu(1,4);
				}while(target[2]==0 or target[2]==target[0] or target[2]==target[1]);
			}
			MPI_Bcast(target, 3, MPI_INT, root, MPI_COMM_WORLD);
			
			vector<int> path_exc=square.get_path(0);
			vector<int> new_path(n_cities);
			int* path1=new int[n_cities];
			int* path2=new int[n_cities];
			for(int c=0;c<n_cities;c++){
				path1[c]=path_exc[c];
			}
			MPI_Barrier(MPI_COMM_WORLD);
			
			if(rank==0){
				MPI_Isend(&path1[0], n_cities, MPI_INTEGER, target[0], itag1, MPI_COMM_WORLD, &req1);
				MPI_Recv(&path2[0], n_cities, MPI_INTEGER, target[0], itag2, MPI_COMM_WORLD, &stat2);
			}
			else if(rank==target[0]){
				MPI_Send(&path1[0], n_cities, MPI_INTEGER, 0, itag2, MPI_COMM_WORLD);
				MPI_Recv(&path2[0], n_cities, MPI_INTEGER, 0, itag1, MPI_COMM_WORLD, &stat1);
			}
			else if(rank==target[1]){
				MPI_Isend(&path1[0], n_cities, MPI_INTEGER, target[2], itag3, MPI_COMM_WORLD, &req2);
				MPI_Recv(&path2[0], n_cities, MPI_INTEGER, target[2], itag4, MPI_COMM_WORLD, &stat4);
			}
			else if(rank==target[2]){
				MPI_Send(&path1[0], n_cities, MPI_INTEGER, target[1], itag4, MPI_COMM_WORLD);
				MPI_Recv(&path2[0], n_cities, MPI_INTEGER, target[1], itag3, MPI_COMM_WORLD, &stat3);
			}
			for(int c=0;c<n_cities;c++){
				new_path[c]=path2[c];
			}
			
			square.set_path(0, new_path);
			MPI_Barrier(MPI_COMM_WORLD);
			if(rank==0) cout<<" generation "<<i<<" migration completed"<<endl;
		}
	}
	
	if(rank==0) cout << "Evolution completed. Printing results..." <<endl<<endl;	//PRINT ON FILE

	square.print_best_path( to_string(rank)+"_SQUARE_bestpath.dat" );


	MPI_Finalize();
	return 0;
}
