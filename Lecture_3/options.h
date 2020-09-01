#include <vector>

using namespace std;

#ifndef __Options__
#define __Options__


class European {

private:
  
	int M;			//n simulations
	int N;			//n blocks
	double S;		//asset price
	double K;		//strike price
	double T;		//delivery time
	double r;		//risk-free interest rate
	double sigma;		//volatility

protected:

public:

	vector<double> avg;
	vector<double> error;

  // constructors
	European(int, int, double, double, double, double, double);

  // destructor
	~European();

  // methods
	
	void Call(int step);
	void Put(int step);
	double err (vector<double> AV, vector<double> AV2, int n);
	
};

#endif //__Options__
