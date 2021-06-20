
#ifndef HOPPER_H
#define HOPPER_H
#include <vector>
#include <array>
#include <random>
#include <chrono> 
#include <complex>


class Underdamp {
protected:
    int noOfPart,chosenColumn, timeInts, ratio ; 
    double dthalf, halfTimeoverMass, precalcSqrt,sumRanIncrem,mass,dt,boxSize, boxSizeS, beta,gamm, alpha, alphasquared,bond, sizeParticle, maxjump, sizex, sizey, posMoves, currentActivity, currentActivity2, div, forceStrength, gammOverMass, BetaovertwoGamm, importedArrayRatio,massdt ; // s is the biasing factor, entropy stores the relative entropy cost of the current run, sizeParticle 1 by default,  size is size of system, over are defined in the usual underdamped langevin way

    std::vector<double> posX, posY, momX, momY; //the positions of the particles as doubles and the momenta
    std::vector<double> virial,dens ; //the sizes of holes
    std::vector<std::vector<double>> relPosX;
    std::vector<std::vector<double>> relPosY; 
    std::vector<std::vector<double>> forceVector; //
    std::vector<double> PredictedForceVector; // our biasing force
    std::vector<double> PredictedForceVectorDer; //
    std::vector<std::vector<double>> distances;
    std::vector<double> bodyForce;
    std::vector<double> stress;  
    std::vector<double> kinEnLoc;
    std::vector<double> potEnLoc;  
    std::vector<std::vector<int>> particlePositionInBoxes;
    std::vector<std::vector<std::vector<int>> > Boxes;
    std::mt19937 & gen;
    std::normal_distribution<double> distribution;
    double runNumber;
    std::vector<double>  forceOnWall;
    int noOfBoxesX, noOfBoxesY, noOfSBoxesX, noOfSBoxesY;
    double Entropy;
 
public:
 Underdamp(double size1, std::mt19937 & gen1)
 :sizeParticle(pow(2.0,1./6.)),
  beta(1.),
  bond(4.),
  sizex(size1+2.0*(sizeParticle-beta/bond)),
  sizey(size1),
  gen(gen1)
  {};

 ~Underdamp()
 {};
 void set( int, double);

 //getters
 

 int get_counter();
 int get_noOfPart();
 double get_entropy();
 std::vector<double>  get_posX();
 std::vector<double>  get_posY();
 std::vector<double>  get_momX();
 std::vector<double>  get_momY();
 std::vector<double>  get_normStress();
 std::vector<double>  get_normkinEnLoc();
 std::vector<double>  get_normPotEn();
 std::vector<double>  get_normDens();
 std::vector<double>  get_normBodyForce();
 std::vector<double> get_force(int , int);
 std::vector<double> closestPoint(std::vector<double> , std::vector<double> ,std::vector<double> );
 bool pointInReg(std::vector<double> , std::vector<double> ,std::vector<double> );

 //setters
 
 void set_posX(std::vector<double>);
 void set_posY(std::vector<double>);
 void set_momX(std::vector<double>);
 void set_momY(std::vector<double>);
 void set_stress();
 
 //update functions
 double repeat_update(int); //update the system *int times
 void update(double &); //normal update that happens like in unbiased system that updates the lth particle
 std::vector<double> bubble_sort(std::vector<double>); //bubble sort stolen from the web, thx web

 

  
 //utility
  bool random_bool();
  int random_int(int);
  int random_int_no_zero(int);
  int random_int_non_neg(int);
  int random_int_ran(int, int);
  double random_real(double);
  double random_real_range(double low, double high);
  double random_real_neg(double);
  
  //calcs
  void Irving_Kirkwood_calc();
  double ovFracLineReg(std::vector<double> , std::vector<double> ,std::vector<double> ,std::vector<double> );
  
  //
  void distancesquared (double&, int, int);
  double potential (int, int);
  void update_force (int, int);
  std::vector<double> get_directional_force(int);
  void update_distances();
  void update_directional_force();
  double mid_point;
  void update_relative_positions(int i, int j);
  void update_state(double &);
  void set_relative_positions();
  void update_boxes();
  void calc_kinetic_Energy(double&);
  void activity(double& );
  void create_stuff();
};


//updates the state of the system after a move
 inline void Underdamp::update_state(double& act)
{
    for (int part1=0;part1<noOfPart;part1++)
    {
      forceVector[part1][0]=0;
      forceVector[part1][1]=0;
      for (int j=-1;j<2;j++)
      {
	for (int k=-1;k<2;k++)
	{
	  if(particlePositionInBoxes[part1][0]+k>=0 and particlePositionInBoxes[part1][0]+k<noOfBoxesX)
	  {
	    if(Boxes[particlePositionInBoxes[part1][0]+k][(particlePositionInBoxes[part1][1]+j+noOfBoxesY)%noOfBoxesY].empty()==false)
	    {
	      for(int l=0;l<Boxes[particlePositionInBoxes[part1][0]+k][(particlePositionInBoxes[part1][1]+j+noOfBoxesY)%noOfBoxesY].size();l++)
	      {
		int i=Boxes[particlePositionInBoxes[part1][0]+k][(particlePositionInBoxes[part1][1]+j+noOfBoxesY)%noOfBoxesY][l];
		if(part1!=i)
		{
		  update_relative_positions(part1,i);
		  distancesquared(distances[part1][i],part1,i);
		  update_force(part1,i);
		  if (distances[i][part1]<sizeParticle*sizeParticle)
		  {
		    if(distances[i][part1]<sizeParticle*sizeParticle/4.)
		    {
		      act+=sizeParticle/2./63./10.;
		    }
		    else
		    {
		      act+=(sizeParticle-sqrt(distances[i][part1]))/63./10.;
		    }
		  }
		}
	      }
	    }
	  }
	}
      }  
    
      //wall force
  //commented out cause channel
      if (posX[part1]>sizex-sizeParticle*1.5)
      {
	forceVector[part1][0]-=((6/pow(sizex-posX[part1],6.0)-12/pow(sizex-posX[part1],12.0))/(sizex-posX[part1])*bond);
	forceOnWall[1]-=((6/pow(sizex-posX[part1],6.0)-12/pow(sizex-posX[part1],12.0))/(sizex-posX[part1])*bond);
      }
      if (posX[part1]<sizeParticle*1.5)
      {
	forceVector[part1][0]+=((6/pow(posX[part1],6.0)-12/pow(posX[part1],12.0))/posX[part1]*bond);
	forceOnWall[0]+=((6/pow(posX[part1],6.0)-12/pow(posX[part1],12.0))/posX[part1]*bond);
      }
  //     if (posY[part1]>size-sizeParticle)
  //     {
  //       forceVector[part1][1]+=(6/pow(posY[part1]-size.,6.0)-12/pow(posY[part1]-size.,12.0))/(size-posY[part1])*4.*bond;;
  //     }
  //     if (posY[part1]<sizeParticle)
  //     {
  //       forceVector[part1][1]-=(6/pow(posY[part1],6.0)-12/pow(posY[part1],12.0))/posY[part1]*4.*bond;
  //     }
    }
}
 
//updates the relative positions
  inline void Underdamp::update_relative_positions(int i, int j)
 {
//   commented out because channel
//   if(posX[j]-posX[i]>size/2)
//   {
//     relPosX[i][j]=posX[j]-posX[i]-size;
//   }
//   else
//   {
//     if(posX[j]-posX[i]<-size/2)
//     {
//       relPosX[i][j]=posX[j]-posX[i]+size;
//     }
//     else
//     {
  relPosX[i][j]=posX[j]-posX[i];
//     }
//   }
  if(posY[j]-posY[i]>sizey/2.)
  {
    relPosY[i][j]=posY[j]-posY[i]-sizey;
  }
  else
  {
    if(posY[j]-posY[i]<-sizey/2.)
    {
      relPosY[i][j]=posY[j]-posY[i]+sizey;
    }
    else
    {
      relPosY[i][j]=posY[j]-posY[i];
    }
  }
 }
 
 
#endif // HOPPER_H
