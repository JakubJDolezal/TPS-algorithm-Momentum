#include <iostream>
# include <cmath>
# include <cstdlib>
# include <ctime>
# include <iomanip>
#include "Underdamp.hpp"
# include <mpi.h>
#include "TPS.hpp"
#include "SimulatedSystem.hpp"
#include <random>
#include <memory>
#include <vector>
#include <math.h>
#include <fstream>
#include <chrono>
#include <iterator>
#include <algorithm>
#include <stdio.h>
#include <string.h>
using namespace std;

void sampling(int id) {
  
  std::random_device rd;
  int blab=rd();
  std::mt19937 gen(blab);
  double size=7.0711;
  int repeats=200000;
  int M=200;
  int noOfPart=24;
  int timeint=63*5;
  double sizePart=1.;
  double s=-1.10-id*0.05;
  Underdamp hop(size,  gen);
  TPS sampler(M,timeint,gen);
  hop.set(noOfPart,s);
 // sampler.load_stuff(s,hop);
  sampler.create_stuff(hop);
  sampler.sampling_with_ft(repeats,s,hop,0);

}

//Simple Mpi implementation for multiple processes
int main ( int argc, char *argv[] )

{
  int ierr;
  int id;
  int p;

ierr = MPI_Init ( &argc, &argv );

//
//  Get the number of processes.
//
  ierr = MPI_Comm_size ( MPI_COMM_WORLD, &p );
//
//  Determine this processes's rank.
//
  ierr = MPI_Comm_rank ( MPI_COMM_WORLD, &id );



    for (int i=0;i<p;i++)
    {
      if (id == i)
      {
	sampling ( id ); // here run code
      }
    }
//
//  Terminate MPI.
//
  MPI_Finalize ( );
//
//  Terminate.
}

  return 0;
}
