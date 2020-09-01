#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <iomanip>
#include "MolDyn_NVE.h"

using namespace std;

int main(){ 
  Input();             						//Inizialization
  int nconf = 1;
  for(int iblk=1; iblk<=nblk; ++iblk){
	Reset(iblk); 
        for(int istep=1; istep <= nstep; ++istep){		
		Move();					//Move particles with Verlet algorithm
		Measure();				//Properties measurement and block averages
		if(istep%10 == 0){
							
			//ConfXYZ(nconf);			//Write actual configuration in XYZ format 	//Commented to avoid "filesystem full"! 
        		nconf += 1;
		}
	}
	Averages(iblk);
  }
  ConfFinal();     						//Write final configuration to restart
  return 0;
}


void Input(void){ 						//Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;
  double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open("input.dat"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> nblk;
  ReadInput >> old;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps per block = " << nstep << endl;
  cout << "Number of blocks = " << nblk << endl << endl;

  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables

//Measurement of g(r)
  igofr = 4;
  nbins = 100;
  n_props = n_props + nbins;
  bin_size = (box/2.0)/(double)nbins;

//Read initial configuration
  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();

//Prepare initial velocities

  if(old==false){
	cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
 	double sumv[3] = {0.0, 0.0, 0.0};
 	for (int i=0; i<npart; ++i){
    	  vx[i] = rand()/double(RAND_MAX) - 0.5;
   	  vy[i] = rand()/double(RAND_MAX) - 0.5;
    	  vz[i] = rand()/double(RAND_MAX) - 0.5;

    	  sumv[0] += vx[i];
    	  sumv[1] += vy[i];
    	  sumv[2] += vz[i];
  	}
  	for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
 	double sumv2 = 0.0, fs;
  	for (int i=0; i<npart; ++i){
    	  vx[i] = vx[i] - sumv[0];
    	  vy[i] = vy[i] - sumv[1];
    	  vz[i] = vz[i] - sumv[2];

   	  sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
 	}
  	sumv2 /= (double)npart;

 	fs = sqrt(3 * temp / sumv2);   //fs = velocity scale factor 
  	for (int i=0; i<npart; ++i){
   	  vx[i] *= fs;
   	  vy[i] *= fs;
   	  vz[i] *= fs;

   	  xold[i] = Pbc(x[i] - vx[i] * delta);
   	  yold[i] = Pbc(y[i] - vy[i] * delta);
    	  zold[i] = Pbc(z[i] - vz[i] * delta);
	}

  } else if(old==true){
	cout << "Read old configuration from file old.0 " << endl << endl;
	ReadConf.open("old.0");					//Read old configuration
	for (int i=0; i<npart; ++i){
   	  ReadConf >> xold[i] >> yold[i] >> zold[i];
	  xold[i] = xold[i] * box;
    	  yold[i] = yold[i] * box;
    	  zold[i] = zold[i] * box;
	}
	ReadConf.close();
	Move();
	double t=0.;
	for (int i=0; i<npart; ++i) 				
	  t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
	stima_temp = (2./3.) * t/(double)npart;

	for(int i=0; i<npart; ++i){
	  vx[i] *= sqrt(temp/stima_temp);			//Velocity scaling factor
	  vy[i] *= sqrt(temp/stima_temp);
	  vz[i] *= sqrt(temp/stima_temp);
	}

	for(int i=0; i<npart; ++i){
	  xold[i] = Pbc(x[i] - vx[i] * delta);
	  yold[i] = Pbc(y[i] - vy[i] * delta);
	  zold[i] = Pbc(z[i] - vz[i] * delta);
	}	
  }
  return;
}

void Move(void){ 			//Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ 		//Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ 		//Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ 	//Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );    // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8));	 // -Grad_ip V(r)  
      }    
    }
  }
  
  return f;
}

void Measure(){ 			//Properties measurement
  int bin;
  double v, t, vij;
  double dx, dy, dz, dr;
 /* ofstream Epot, Ekin, Etot, Temp;

  Epot.open("output_epot.dat",ios::app);
  Ekin.open("output_ekin.dat",ios::app);
  Temp.open("output_temp.dat",ios::app);
  Etot.open("output_etot.dat",ios::app);*/

  v = 0.0; //reset observables
  t = 0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( xold[i] - xold[j] );
     dy = Pbc( yold[i] - yold[j] );
     dz = Pbc( zold[i] - zold[j] );

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

//Update of the histogram of g(r)
     bin = igofr + (int)(dr/bin_size);
     if(bin < igofr + nbins)
        	blk_av[bin] += 2.0;

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
//Potential energy
       v += vij;
     }
    }          
  }

 /* //Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]); 				
   
  stima_pot = v/(double)npart; 				//Potential energy per particle
  stima_kin = t/(double)npart; 				//Kinetic energy per particle
  stima_temp = (2.0 / 3.0) * t/(double)npart; 		//Temperature
  stima_etot = (t+v)/(double)npart; 				//Total energy per particle

  blk_av[iv] += stima_pot;
  blk_av[ik] += stima_kin;
  blk_av[ie] += stima_etot;
  blk_av[it] += stima_temp;*/

  blk_norm = blk_norm + 1.0;
/*
  Epot << stima_pot  << endl;
  Ekin << stima_kin  << endl;
  Temp << stima_temp << endl;
  Etot << stima_etot << endl;

  Epot.close();
  Ekin.close();
  Temp.close();
  Etot.close();*/

  return;
}


void Reset(int iblk){ 		//Reset block averages
   if(iblk == 1){
       for(int i=0; i<n_props; ++i){
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }
   for(int i=0; i<n_props; ++i)
     	blk_av[i] = 0;
   blk_norm = 0;
}


void ConfFinal(void){ 		//Write final configuration
  ofstream WriteConf, WriteOld;

  cout << "Print final configuration time t to file config.final and time t-dt to file old.final" << endl << endl;
  WriteConf.open("config.final");
  WriteOld.open("old.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
    WriteOld << xold[i]/box << "   " << yold[i]/box << "   " << zold[i]/box << endl;
  }
  WriteConf.close();
  WriteOld.close();
  return;
}


void Averages(int iblk){
  double r, gdir;
  cout << "Block number " << iblk << endl;

//  ofstream AveEpot, AveEkin, AveTemp, AveEtot;
  ofstream Gave;
/*
  glob_av[iv] += blk_av[iv]/blk_norm;
  glob_av[ik] += blk_av[ik]/blk_norm;
  glob_av[ie] += blk_av[ie]/blk_norm;
  glob_av[it] += blk_av[it]/blk_norm;

  glob_av2[iv] += (blk_av[iv]/blk_norm)*(blk_av[iv]/blk_norm);
  glob_av2[ik] += (blk_av[ik]/blk_norm)*(blk_av[ik]/blk_norm);
  glob_av2[ie] += (blk_av[ie]/blk_norm)*(blk_av[ie]/blk_norm);
  glob_av2[it] += (blk_av[it]/blk_norm)*(blk_av[it]/blk_norm);

  AveEpot.open("ave_epot.out",ios::app);
  AveEkin.open("ave_ekin.out",ios::app);
  AveEtot.open("ave_etot.out",ios::app);
  AveTemp.open("ave_temp.out",ios::app);*/
  Gave.open("output_g.out",ios::app);
/*
  AveEpot << iblk << "   " << setprecision(10) << glob_av[iv]/double(iblk) << "   " << Error(glob_av[iv], glob_av2[iv], iblk) << endl;
  AveEkin << iblk << "   " << setprecision(10) << glob_av[ik]/double(iblk) << "   " << Error(glob_av[ik], glob_av2[ik], iblk) << endl;
  AveEtot << iblk << "   " << setprecision(10) << glob_av[ie]/double(iblk) << "   " << Error(glob_av[ie], glob_av2[ie], iblk) << endl;
  AveTemp << iblk << "   " << setprecision(10) << glob_av[it]/double(iblk) << "   " << Error(glob_av[it], glob_av2[it], iblk) << endl;
*/
  for(int i=igofr; i<nbins+igofr; ++i){	
	r = (i-igofr) * bin_size;
     	gdir = (blk_av[i]/blk_norm)/(rho*(double)npart*(4.*M_PI/3.)*(pow(r+bin_size,3)-pow(r,3)));

        glob_av[i] += gdir;
        glob_av2[i] += pow(gdir,2);

	if(iblk == nblk)
        	Gave << iblk << " " << r << " " << glob_av[i]/(double)iblk << " " << setprecision(10) << Error(glob_av[i],glob_av2[i],iblk) << endl;
  }
  cout << "----------------------------" << endl << endl;
 /* 
  AveEpot.close();
  AveEkin.close();
  AveTemp.close();
  AveEtot.close();*/
  Gave.close();
}      
 

void ConfXYZ(int nconf){ 	//Write configuration in .xyz format

  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();

}

double Pbc(double r){  		//Algorithm for periodic boundary conditions with side l=box
    return r - box * rint(r/box);
}

double Error(double sum, double sum2, int iblk){
    if( iblk == 1 ) 
	return 0.0;
    else 
	return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

