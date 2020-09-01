#include <vector>
#include <string>

#include "random.h"

using namespace std;

#ifndef __class_h__
#define __class_h__

//*********************************INDIVIDUAL

class Individual{

private:

	int m_Ncities;			
	vector<double> m_citiesX;	
	vector<double> m_citiesY;	
	vector<int> m_path;
	double m_L1;

public:

	//constructors
	Individual(){};
	Individual( vector<double>, vector<double>, vector<int> );
	//destructor
	~Individual();
	
	//methods
	double L1();
	void check();
	double getL1() {return m_L1;};
	vector<int> getPath() {return m_path;};

};

//*******************************GENERATION

class Generation{
	
private:

	int m_Ncities;
	vector<double> m_citiesX;
	vector<double> m_citiesY;
	Individual m_individual;
	double m_beta;
	Random m_rnd;

public:

	//constructor
	Generation( vector<double>, vector<double>, double );
	//destructor
	~Generation();

	//methods
	vector<int> random_path();
	void new_generation();
	Individual mutate ( Individual );
	vector<int> pair_permutation( vector<int> );
	vector<int> shift( vector<int> );
	vector<int> block_permutation( vector<int> );
	vector<int> inversion( vector<int> );
	void stats( string );
	void print_best_path( string );
	void setBeta( double );
	double Boltzmann( double );


};



#endif //__class_h__
