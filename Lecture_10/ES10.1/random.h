#ifndef __Random__
#define __Random__

class Random {

	private:
  		int m1,m2,m3,m4,l1,l2,l3,l4,n1,n2,n3,n4;

	protected:

	public:
  		Random();	// constructor
        	~Random();     	// destructor
		
		void Init();				//methods
 		void SetRandom(int * , int, int);	
  		void SaveSeed();
 		double Rannyu(void);
  		double Rannyu(double min, double max);
		double Exp(double lambda);
		double Lorentz(double gamma, double mu);
  		double Gauss(double mean, double sigma);
		double Sampling();
};	

#endif // __Random__

