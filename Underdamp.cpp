#include <cstdlib>
#include "Underdamp.hpp"
#include <iostream>
#include <random>
#include <vector>
#include <math.h>
#include <fstream>
#include <iterator>
#include <unistd.h>
#include <complex>
#include <unistd.h>
#include <string>
#include <stdlib.h> 
#include <functional> 
#include <algorithm>
using namespace std;


const double pi = 3.1415926535897;

double range_round(double thing,double range){
  double thing2=thing;
  if (thing>range/2)
  {
    thing2=thing-range;
  }
  else
  {
    if (thing<-range/2)
    {
      thing2=thing+range;
    }
  }
  return thing2;
}

double range_round_tight(double thing,double range){
  double thing2=thing;
  if (thing>range)
  {
    thing2=thing-range;
  }
  else
  {
    if (thing<0)
    {
      thing2=thing+range;
    }
  }
  return thing2;
}

std::vector<double>  Underdamp::get_posX()
{
  return posX;
}

std::vector<double>  Underdamp::get_posY()
{
  return posY;
}

std::vector<double>  Underdamp::get_momX()
{
  return momX;
}

std::vector<double>  Underdamp::get_momY()
{
  return momY;
}



void Underdamp::set_posX(std::vector<double> setter )
{
  posX=setter;
}

void Underdamp::set_posY(std::vector<double> setter )
{
  posY=setter;
}

void Underdamp::set_momX(std::vector<double> setter )
{
  momX=setter;
}

void Underdamp::set_momY(std::vector<double> setter )
{
  momY=setter;
}

double Underdamp::get_entropy()
{
  return Entropy;
}

int Underdamp::get_noOfPart()
{
  return noOfPart;
}
  
void Underdamp::set(int partNo, double s) 
//this essentially initialises everything, it
/*
partNo-number of particles we want
s- the s-value of if we want to add a guiding force
*/
{
    std::normal_distribution<double> distribution1(0.,1.0);
    double check=2.0*(sizeParticle-beta/bond);
    sizex=sizey+check;
    distribution=distribution1;
    sizeParticle= pow(2.,1./6.);
    noOfPart=partNo;
    Entropy=0;
    dt=1./500.;
    posX.reserve(noOfPart);
    posY.reserve(noOfPart);
    momX.reserve(noOfPart);
    momY.reserve(noOfPart);
    ratio=10;
    boxSize=1.25;

    boxSizeS=boxSize/(double)ratio;
    noOfBoxesX=ceil(sizex/boxSize);
    noOfBoxesY=ceil(sizey/boxSize);
    noOfSBoxesX=noOfBoxesX*ratio;
    noOfSBoxesY=noOfBoxesY*ratio;

    std::vector<std::vector<double> > stuff(noOfPart, (std::vector<double>(2,0)));
    std::vector<std::vector<int> > stuffint(noOfPart, (std::vector<int>(3,0)));
    std::vector<std::vector<std::vector<int>> > stuffint3(noOfBoxesX, (std::vector<std::vector<int>>(noOfBoxesY,std::vector<int>(0,0))));  
    std::vector<std::vector<double> > stuff2(noOfPart, (std::vector<double>(noOfPart,0)));
    distances=stuff2;
    forceVector=stuff;
    relPosX=stuff2;
    relPosY=stuff2;
    
    // we clear old data
    posX.clear();
    posY.clear();
    momX.clear();
    momY.clear();
    //next we distribute the partiles between the walls
    double k=1;
    double l=sizeParticle-beta/bond;
    particlePositionInBoxes=stuffint;
    for(int i=0;i<noOfBoxesX;i++)
    {
      for(int j=0;j<noOfBoxesY;j++)
      {
	stuffint3[i][j].reserve(noOfPart);
      }
    }
    Boxes=stuffint3;
    std::uniform_real_distribution<double> dis(-1, 1);
    for(int i=0;i<noOfPart;i++)
    {
      if(i%(int)(ceil(sqrt((double)noOfPart)))==0 and i!=0)
      {
        l=l+(sizey)/floor(sqrt((double)noOfPart));
        k=1;
      }
      posY.push_back(k);
      posX.push_back(l);
      particlePositionInBoxes[i][0]=(int)(floor(posX[i]/boxSize));
      particlePositionInBoxes[i][1]=(int)(floor(posY[i]/boxSize));
      Boxes[(int)(floor(posX[i]/boxSize))][(int)(floor(posY[i]/boxSize))].push_back(i);
      particlePositionInBoxes[i][2]=Boxes[(int)(floor(posX[i]/boxSize))][(int)(floor(posY[i]/boxSize))].size()-1;
      momX.push_back(dis(gen));
      momY.push_back(dis(gen));
      k=k+((sizey-2.)/ceil(sqrt((double)noOfPart)));
    }
    // Here we set the values of the system such as mass and temperature
    mass=1.;
    bond=4.;
    beta=1.;
    gamm=10.;
    alpha=exp(-gamm/mass*dt);
    alphasquared=exp(-2*gamm/mass*dt);
    precalcSqrt=sqrt((1-alphasquared)/beta);
    dthalf=dt/2.;
    halfTimeoverMass=dthalf/mass;
    gammOverMass=gamm/mass; 
    massdt=mass*dt;
    BetaovertwoGamm=beta/gamm/2;
    //Here we load the force for the s-value we previously obtained, if s==0 we just simulate the system without a guiding force
    if(s!=0)
    {
	    string sStr=to_string(s).substr(0,5);
	    string txt=".txt";
	    {
		 string local="Force/Force";
		 local.append(sStr);
		 local.append(txt);
		 std::ifstream MAct(local);
		 std::vector<double>randomVector {std::istream_iterator<double>(MAct), {}};
		 PredictedForceVector=randomVector;
	   	 std::transform(PredictedForceVector.cbegin(),PredictedForceVector.cend(),PredictedForceVector.begin(),std::negate<double>());
	    }
    }
    else
    {
            std::vector<double>randomVector (std::vector<double>(noOfBoxesX,0));
            PredictedForceVector=randomVector;
    }
    importedArrayRatio=PredictedForceVector.size()/noOfBoxesX/boxSize;
 }
 
//Function for calculating the kinetic energy 
 void Underdamp::calc_kinetic_Energy(double & kin)
 {
   for(int i=0;i<noOfPart;i++)
   {
     kin+=momX[i]*momX[i]*mass/2;
     kin+=momY[i]*momY[i]*mass/2;
   }
 }
 
 //Function for calculating the clustering/activity
 void Underdamp::activity( double & act)
 {
  for (int part1=0;part1<noOfPart;part1++)
  {
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
		if (distances[i][part1]<sizeParticle*sizeParticle)
		{
		  if(distances[i][part1]<sizeParticle*sizeParticle/4)
		  {
		    act+=sizeParticle/2.;
		  }
		  else
		  {
		    act+=sizeParticle-sqrt(distances[i][part1]);
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
 }

   //We calculate the distance between two particles
 void Underdamp::distancesquared(double& dist, int part1, int part2)
{
    double disX=posX[part1]-posX[part2];
    double disY=posY[part1]-posY[part2];
//     Version with no walls
//     double dist=sqrt(pow((range_round(disX,size)),2.)+pow(range_round(disY,size),2.));
//  Version with a channel
    dist =pow(disX,2.)+pow(range_round(disY,sizey),2.);
// Version with walls
//     double dist=sqrt(pow(disX),2.)+pow(disY,2.));

 }
  
   //We update boxes to know where particles are and which other particles they can interact with
void Underdamp::update_boxes()
{
   for (int i=0; i<noOfPart;i++)
   {
     if((int)(floor(posX[i]/boxSize))!=particlePositionInBoxes[i][0] or (int)(floor(posY[i]/boxSize))!=particlePositionInBoxes[i][1])
     {
	Boxes[particlePositionInBoxes[i][0]][particlePositionInBoxes[i][1]][particlePositionInBoxes[i][2]]=Boxes[particlePositionInBoxes[i][0]][particlePositionInBoxes[i][1]][Boxes[particlePositionInBoxes[i][0]][particlePositionInBoxes[i][1]].size()-1];
	particlePositionInBoxes[Boxes[particlePositionInBoxes[i][0]][particlePositionInBoxes[i][1]][Boxes[particlePositionInBoxes[i][0]][particlePositionInBoxes[i][1]].size()-1]][2]=particlePositionInBoxes[i][2];
	Boxes[particlePositionInBoxes[i][0]][particlePositionInBoxes[i][1]].pop_back();
	particlePositionInBoxes[i][0]=(int)(floor(posX[i]/boxSize));
	particlePositionInBoxes[i][1]=(int)(floor(posY[i]/boxSize));
	Boxes[(int)(floor(posX[i]/boxSize))][(int)(floor(posY[i]/boxSize))].push_back(i);
	particlePositionInBoxes[i][2]=Boxes[(int)(floor(posX[i]/boxSize))][(int)(floor(posY[i]/boxSize))].size()-1;
     }
   }
}

 //updates force with WCA forces
 void Underdamp::update_force(int part1, int part2)
 {
    if (distances[part1][part2]<sizeParticle and distances[part1][part2]!=0)
    {
      forceVector[part1][0]-=(6/pow(distances[part1][part2],3.0)-12/pow(distances[part1][part2],6.0))/distances[part1][part2]*bond*relPosX[part1][part2];
      forceVector[part1][1]-=(6/pow(distances[part1][part2],3.0)-12/pow(distances[part1][part2],6.0))/distances[part1][part2]*bond*relPosY[part1][part2];
    }
 }
 



//Here we run the system for "repeats" number of times
 double Underdamp::repeat_update(int repeats)
{
  Entropy=0;
  std::vector<double> init (Boxes.size(),0);
  virial=init;
  std::vector<double> nullVector(2,0);
  forceOnWall=nullVector;
  update_boxes();
//   update_distances();
//   update_directional_force();
  double localAct=0;
  update_state(localAct);
  localAct=0;
  for (int j=0;j<repeats;j++)
  {
    update(localAct);
//     activity(localAct);
//     calc_kinetic_Energy(kinEn);
  }
//   kinEn/=(double)(repeats*noOfPart);
  return localAct;
}

//Here we update the system in a baoab way
void Underdamp::update(double & act)
{ 
 //update velocity 1/2
 //update positiion 1/2
  std::vector<double> oldmomX=momX;
  std::vector<double> oldposX=posX;
  std::vector<std::vector<double>> oldFV=forceVector;
  for (int i=0;i<noOfPart;i++)
  { 
    momX[i]-=(forceVector[i][0]+PredictedForceVector[int(posX[i]*importedArrayRatio)])*dthalf;
    momY[i]-=forceVector[i][1]*dthalf;
  }
  for (int i=0;i<noOfPart;i++)
  { 
    posX[i]+=momX[i]*halfTimeoverMass;
    posY[i]=range_round_tight(posY[i]+momY[i]*halfTimeoverMass,sizey);
  }

  
 //update velocity randomness
  for (int i=0;i<noOfPart;i++)
  {
    momX[i]=alpha*momX[i]+precalcSqrt*distribution(gen);
    momY[i]=alpha*momY[i]+precalcSqrt*distribution(gen);
    posX[i]+=momX[i]*halfTimeoverMass;
    posY[i]=range_round_tight(posY[i]+momY[i]*halfTimeoverMass, sizey);
  }
  //now update all the things that have changed
  update_boxes();
  update_state(act);
  
 //update velocity with final update and deal with KL-divergence for the guiding force
  for (int i=0;i<noOfPart;i++)
  { 
    momX[i]-=(forceVector[i][0]+PredictedForceVector[int(posX[i]*importedArrayRatio)])*dthalf;
    momY[i]-=forceVector[i][1]*dthalf;
    Entropy+=(oldmomX[i]-momX[i]+((forceVector[i][0]+oldFV[i][0])/2.+(PredictedForceVector[int(posX[i]*importedArrayRatio)]+PredictedForceVector[int(oldposX[i]*importedArrayRatio)])/4.)*dt)*BetaovertwoGamm*(PredictedForceVector[int(posX[i]*importedArrayRatio)]+PredictedForceVector[int(oldposX[i]*importedArrayRatio)])/2.;
  }
}



//Function that returns the force
std::vector<double> Underdamp::get_force(int part1, int part2)
 {
	  std::vector<double> res(std::vector<double>(2,0));
	  if (distances[part1][part2]<sizeParticle and distances[part1][part2]!=0)
	  {
	      res[0]-=(6/pow(distances[part1][part2],3.0)-12/pow(distances[part1][part2],6.0))/distances[part1][part2]*bond*relPosX[part1][part2];
	      res[1]-=(6/pow(distances[part1][part2],3.0)-12/pow(distances[part1][part2],6.0))/distances[part1][part2]*bond*relPosY[part1][part2];
	  }
	  return res;
 }
 //Function that returns the potential
double Underdamp::potential(int part1, int part2)
{
    double V=0;
    if (distances[part1][part2]<sizeParticle and distances[part1][part2]!=0)
    {
    	V=(1/pow(distances[part1][part2],6)-1/pow(distances[part1][part2],3))*bond+bond;
    }
    return V;
}

//Utility functions from here dealing with stress calculations from here

  //Function that sets/resets the stress calculation
void Underdamp::set_stress()
{
	std::vector<double>  stuff (noOfSBoxesX*noOfSBoxesY*2,0);
	stress=stuff;
	kinEnLoc=stuff;
	std::vector<double>  stuff2 (noOfSBoxesX*noOfSBoxesY,0);
	potEnLoc=stuff2;
	dens=stuff2;
	bodyForce=stuff2;
	runNumber=0;
}

//returns the normalised stress
std::vector<double>  Underdamp::get_normStress()
{
	for (int i=0;i<stress.size();i++)
 	{
		stress[i]=stress[i]/runNumber/boxSizeS/boxSizeS;
	}
	return stress;	 
}

//returns the normalised body force
std::vector<double>  Underdamp::get_normBodyForce()
{
	for (int i=0;i<bodyForce.size();i++)
 	{
		bodyForce[i]=bodyForce[i]/runNumber/boxSizeS/boxSizeS;
	}
	return bodyForce;	 
}

//returns the normalised kinetic energy
std::vector<double>  Underdamp::get_normkinEnLoc()
{
	for (int i=0;i<kinEnLoc.size();i++)
 	{
		kinEnLoc[i]=kinEnLoc[i]/runNumber/boxSizeS/boxSizeS;
	}
	return kinEnLoc;	 
}

//returns the normalised potential energy
std::vector<double>  Underdamp::get_normPotEn()
{
	for (int i=0;i<potEnLoc.size();i++)
 	{
		potEnLoc[i]=potEnLoc[i]/runNumber/boxSizeS/boxSizeS;
	}
	return potEnLoc;	 
}

//returns the normalised density
std::vector<double>  Underdamp::get_normDens()
{
	for (int i=0;i<potEnLoc.size();i++)
 	{
		dens[i]=dens[i]/runNumber/boxSizeS/boxSizeS;
	}
	return dens;	 
}

//Calculate the irving Kirkwood stress tensor in every box with PBC
void Underdamp::Irving_Kirkwood_calc()
{
 update_boxes();
 runNumber=runNumber+1;
 std::vector<double> stuff (std::vector<double>(2,0));
 for (int i=0;i<noOfBoxesX;i++)
 {
        for (int j=0;j<noOfBoxesY;j++)
        {
                if(Boxes[i][j].empty()==false)
                {

                    for (int k=0;k<Boxes[i][j].size();k++)
                    {
                        int part1=Boxes[i][j][k];
                        int SmallBoxPosXP1=floor(posX[part1]/boxSizeS);
                        int SmallBoxPosYP1=floor(posY[part1]/boxSizeS);
                        dens[SmallBoxPosXP1*noOfBoxesY*ratio+SmallBoxPosYP1]+=1;
		        if (posX[part1]>sizex-1.5*sizeParticle)
		        {
			  double tempForce1=((-6/pow(sizex-posX[part1],6.0)+12/pow(sizex-posX[part1],12.0))/(sizex-posX[part1])*bond);
			  bodyForce[SmallBoxPosXP1*noOfBoxesY*ratio+SmallBoxPosYP1]+=tempForce1;
		        }
		        if (posX[part1]<1.5*sizeParticle)
		        {
			  double tempForce1=((6/pow(posX[part1],6.0)-12/pow(posX[part1],12.0))/posX[part1]*bond);
			  bodyForce[SmallBoxPosXP1*noOfBoxesY*ratio+SmallBoxPosYP1]+=tempForce1;
		        }	
                        // Get the momentum component
                        stress[SmallBoxPosXP1*noOfBoxesY*ratio*2+SmallBoxPosYP1*2]+=momX[part1]*momX[part1];
                        stress[SmallBoxPosXP1*noOfBoxesY*ratio*2+SmallBoxPosYP1*2+1]+=momY[part1]*momY[part1];
                        kinEnLoc[SmallBoxPosXP1*noOfBoxesY*ratio*2+SmallBoxPosYP1*2]+=momX[part1]*momX[part1]/2.;
                        kinEnLoc[SmallBoxPosXP1*noOfBoxesY*ratio*2+SmallBoxPosYP1*2+1]+=momY[part1]*momY[part1]/2.;
                        std::vector<double> positionPart {posX[part1], posY[part1]};
                        for (int i1=-1;i1<2;i1++)
                        {
                            for (int j1=-1;j1<2;j1++)
                            {
                                if((j1!=0) || (i1!=0))
                                {
                                    if(Boxes[(i+i1+noOfBoxesX)%noOfBoxesX][(j+j1+noOfBoxesY)%noOfBoxesY].empty()==false and (i+i1)<noOfBoxesX and (i+i1)>=0)
                                    {
                                        for (int k1=0;k1<Boxes[i+i1][(j+j1+noOfBoxesY)%noOfBoxesY].size();k1++)
                                        {
                                            positionPart[1]=posY[part1];
                                            int part2=Boxes[i+i1][(j+j1+noOfBoxesY)%noOfBoxesY][k1];
                                            update_relative_positions(part1,part2);
                                            distancesquared(distances[part1][part2],part1,part2);
                                            std::vector<double> positionPart2 {posX[part2], posY[part2]};
                                            std::vector<double> forceTemp=get_force(part1, part2);
                                            if(forceTemp[0]!=0 || forceTemp[1]!=0)
                                            {
                                                double potLoc=potential(part1, part2);
                                                int SmallBoxPosXP2=floor(posX[part2]/boxSizeS);
                                                int SmallBoxPosYP2=floor(posY[part2]/boxSizeS);
                                                potEnLoc[SmallBoxPosXP1*noOfBoxesY*ratio+SmallBoxPosYP1]+=potLoc/2.;
                                                potEnLoc[SmallBoxPosXP2*noOfBoxesY*ratio+SmallBoxPosYP2]+=potLoc/2.;
                                                if((j+j1)>=noOfBoxesY)
                                                {
                                                    positionPart2[1]+=sizey;
                                                }
                                                if((j+j1)<0)
                                                {
                                                    positionPart2[1]-=sizey;
                                                }
                                                for(int i2=0;i2<ratio;i2++)
                                                {
                                                    for(int j2=0;j2<ratio;j2++)
                                                    {
                                                        std::vector<double> UR1 {i*boxSize+(i2+1)*boxSizeS, j*boxSize+(j2+1)*boxSizeS};
                                                        std::vector<double> LL1 {i*boxSize+i2*boxSizeS, j*boxSize+j2*boxSizeS};
                                                        double frac=ovFracLineReg( positionPart, positionPart2, LL1, UR1 );
                                                        int number=((i*ratio)+i2)*noOfBoxesY*ratio*2+(j*ratio+j2)*2;
                                                        stress[number]+=frac*relPosX[part1][part2]*forceTemp[0]/2.;
                                                        stress[number+1]+=frac*relPosY[part1][part2]*forceTemp[1]/2.;
                                                    }
                                                }
                                                if((j+j1)>=noOfBoxesY)
                                                {
                                                    positionPart2[1]-=sizey;
                                                    positionPart[1]-=sizey;
                                                }
                                                if((j+j1)<0)
                                                {
                                                    positionPart2[1]+=sizey;
                                                    positionPart[1]+=sizey;
                                                }
                                                for(int i2=0;i2<ratio;i2++)
                                                {
                                                    for(int j2=0;j2<ratio;j2++)
                                                    {
                                                        std::vector<double> UR2 { (i+i1)*boxSize+(i2+1)*boxSizeS, ((j+j1+noOfBoxesY)%noOfBoxesY)*boxSize+(j2+1)*boxSizeS};
                                                        std::vector<double> LL2 { (i+i1)*boxSize+i2*boxSizeS, ((j+j1+noOfBoxesY)%noOfBoxesY)*boxSize+j2*boxSizeS};
                                                        double frac1=ovFracLineReg( positionPart, positionPart2,  LL2, UR2 );
                                                        int number=((i+i1)*ratio+i2)*noOfBoxesY*ratio*2+(((j+j1+noOfBoxesY)%noOfBoxesY)*ratio+j2)*2;
                                                        stress[number]+=frac1*relPosX[part1][part2]*forceTemp[0]/2.;
                                                        stress[number+1]+=frac1*relPosY[part1][part2]*forceTemp[1]/2.;
                                                    }
                                                }
                                                if((j1!=0) && (i1!=0))
                                                {
                                                //First check i1=0 box if it has any contribution
                                                    std::vector<double> UR4 { (i+1)*boxSize, ((j+j1+noOfBoxesY)%noOfBoxesY+1)*boxSize};
                                                    std::vector<double> LL4 { (i)*boxSize, ((j+j1+noOfBoxesY)%noOfBoxesY)*boxSize};
                                                    double frac2=ovFracLineReg( positionPart, positionPart2, LL4, UR4);
                                                    if(frac2!=0)
                                                    {
                                                        for(int i2=0;i2<ratio;i2++)
                                                        {
                                                            for(int j2=0;j2<ratio;j2++)
                                                            {
                                                                std::vector<double> UR6 { i*boxSize+(i2+1)*boxSizeS, ((j+j1+noOfBoxesY)%noOfBoxesY)*boxSize+(j2+1)*boxSizeS};
                                                                std::vector<double> LL6 { i*boxSize+i2*boxSizeS, ((j+j1+noOfBoxesY)%noOfBoxesY)*boxSize+j2*boxSizeS};
                                                                double frac21=ovFracLineReg( positionPart, positionPart2, LL6, UR6);
                                                                int number=(i*ratio+i2)*noOfBoxesY*ratio*2+((j+j1+noOfBoxesY)%noOfBoxesY*ratio+j2)*2;
                                                                stress[number]+=frac21*relPosX[part1][part2]*forceTemp[0]/2.;
                                                                stress[number+1]+=frac21*relPosY[part1][part2]*forceTemp[1]/2.;
                                                            }	
                                                        }	
                                                    }
                                            //Now check j1=0 box if it has any contribution
                                                    if((j+j1)>=noOfBoxesY)
                                                    {
                                                        positionPart2[1]=positionPart2[1]+sizey;
                                                        positionPart[1]=positionPart[1]+sizey;
                                                    }
                                                    if((j+j1)<0)
                                                    {
                                                        positionPart2[1]=positionPart2[1]-sizey;
                                                        positionPart[1]=positionPart[1]-sizey;
                                                    }
                                                    std::vector<double> UR3 { (i+i1+1)*boxSize, (j+1)*boxSize};
                                                    std::vector<double> LL3 { (i+i1)*boxSize, j*boxSize};
                                                    double frac3=ovFracLineReg( positionPart, positionPart2, LL3, UR3 );
                                                    if(frac3!=0)
                                                    {
                                                        for(int i2=0;i2<ratio;i2++)
                                                        {
                                                            for(int j2=0;j2<ratio;j2++)
                                                            {
                                                                std::vector<double> UR31 { (i+i1)*boxSize+(i2+1)*boxSizeS, j*boxSize+(j2+1)*boxSizeS};
                                                                std::vector<double> LL31 { (i+i1)*boxSize+i2*boxSizeS, j*boxSize+j2*boxSizeS};
                                                                double frac31=ovFracLineReg( positionPart, positionPart2,LL31, UR31);
                                                                int number=((i+i1)*ratio+i2)*noOfBoxesY*ratio*2+(j*ratio+j2)*2;
                                                                stress[number]+=frac31*relPosX[part1][part2]*forceTemp[0]/2.;
                                                                stress[number+1]+=frac31*relPosY[part1][part2]*forceTemp[1]/2.;
                                                            }	
                                                        }	
                                                    }
                                                }
                                            }
                                        }
                                            
                                    }
                                }
                                else
                                {
                                    if(Boxes[i][j].size()>k)
                                    {
                                        positionPart[1]=posY[part1];
                                        for (int k1=k+1;k1<Boxes[i][j].size();k1++)
                                        {
                                            
                                            int part2=Boxes[i][j][k1];
                                            update_relative_positions(part1,part2);
                                            distancesquared(distances[part1][part2],part1,part2);
                                            std::vector<double> forceTemp =get_force(part1, part2);
                                            if(forceTemp[0]!=0 || forceTemp[1]!=0)
                                            {
                                                int SmallBoxPosXP2=floor(posX[part2]/boxSizeS);
                                                int SmallBoxPosYP2=floor(posY[part2]/boxSizeS);
                                                double potLoc=potential(part1, part2);
                                                potEnLoc[SmallBoxPosXP1*noOfBoxesY*ratio+SmallBoxPosYP1]+=potLoc/2;
                                                potEnLoc[SmallBoxPosXP2*noOfBoxesY*ratio+SmallBoxPosYP2]+=potLoc/2;
                                                std::vector<double> positionPart2{posX[part2], posY[part2]};
                                                for(int i2=0;i2<ratio;i2++)
                                                {
                                                    for(int j2=0;j2<ratio;j2++)
                                                    {
                                                        std::vector<double> UR4 { (i*ratio+i2+1)*boxSizeS, (j*ratio+j2+1)*boxSizeS};
                                                        std::vector<double> LL4 { (i*ratio+i2)*boxSizeS, (j*ratio+j2)*boxSizeS};
                                                        double fracSB=ovFracLineReg( positionPart, positionPart2, LL4, UR4 );
                                                        if (fracSB!=0)
                                                        {
                                                            int number=(i*ratio+i2)*noOfBoxesY*ratio*2+(j*ratio+j2)*2;
                                                            stress[number]+=forceTemp[0]*relPosX[part1][part2]*fracSB;
                                                            stress[number+1]+=forceTemp[1]*relPosY[part1][part2]*fracSB;	
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }		
                            }
                        }	
                    }
            }
        }
 }
}

//Finds out if point is in square  
    bool Underdamp::pointInReg( std::vector<double> lineA, std::vector<double> regLL, std::vector<double>  regUR ) {
      bool res = 0;
      if ( lineA[0] >= regLL[0] && lineA[1] >= regLL[1]
           && lineA[0] <= regUR[0] && lineA[1] <= regUR[1] )
	{
        res = 1;
	}
      return res;
    }

//For cube border points either lower left, lower right, upper right, upper left we find the closest point to point 
    std::vector<double> Underdamp::closestPoint( std::vector<double> point, std::vector<double> regLL, std::vector<double> regUR ) 
    {
        std::vector<double> res (std::vector<double>(2,0));
	if((point[0]-regUR[0])*(point[0]-regUR[0])<=(point[0]-regLL[0])*(point[0]-regLL[0]))
	{
		if((point[1]-regUR[1])*(point[1]-regUR[1])<=(point[1]-regLL[1])*(point[1]-regLL[1]))
		{
			res[0]=regUR[0];
			res[1]=regUR[1];
		}
		else
		{
			res[0]=regUR[0];
			res[1]=regLL[1];
		}
	}
	else
	{
		if((point[1]-regUR[1])*(point[1]-regUR[1])<=(point[1]-regLL[1])*(point[1]-regLL[1]))
		{
			res[0]=regLL[0];
			res[1]=regUR[1];
		}
		else
		{
			res[0]=regLL[0];
			res[1]=regLL[1];
		}
	}

      return res;
    }
//finds fraction of line between lineA and lineB that is in square defined by regLL and regUR
    double Underdamp::ovFracLineReg( std::vector<double> lineA, std::vector<double> lineB, std::vector<double> regLL, std::vector<double> regUR ) {    
      char diag = 0;
    
      double res = 0.0;
      int    dim = 2;
      
      int insideA = (int)pointInReg( lineA, regLL, regUR);
      int insideB = (int)pointInReg( lineB, regLL, regUR);;
    

      if ( insideA && insideB ) {  // both inside is the easy case
        res = 1.0;
      
      } else  if ( (! insideA) &&  (! insideB) ) {  // both outside

        int nCut = 0;
        int doCut[dim][2];
        double cutVals[dim][2];

        // d for is direction (ie axis) w is for upper/lower edge
        // if d==0 then we are considering a horizontal edge
        for (int d=0;d<dim;d++) {
          for (int w=0;w<2;w++) {
            // assume first that the line does not cut this edge
            doCut[d][w]=0;
            
            // coordinate of edge
            double regVal = ( w==0 ? regLL[d] : regUR[d] );
            double lambda = ( regVal - lineA[d] )/( lineB[d] - lineA[d] );
            // if the points are opp sides of edge
            if ( lambda > 0.0 && lambda < 1.0 ) {
              double cutVal = lineA[1-d] + lambda * ( lineB[1-d] - lineA[1-d] );

              if ( cutVal > regLL[1-d] && cutVal < regUR[1-d] ) {
                doCut  [d][w] = 1;
                cutVals[d][w] = cutVal;
                nCut++;
              }
            }
          }
        }
        
        char resDone = 0;
        
        if ( nCut == 0 ) { res = 0.0; resDone = 1; }
        
        // do we cut two opposite sides?
        if ( ! resDone ) for (int d=0;d<dim;d++) if ( doCut[d][0] && doCut[d][1] ) {
          res = ( regUR[d] - regLL[d] ) / fabs( lineA[d] - lineB[d] );
          resDone = 1;
        }
        
        if ( ! resDone ) {
          // at this point we must be cutting two adjacent sides
          // we always choose the horizontal edges to compute the fraction
          int d=1;
          double cutVal;  // this will be an x coord
          for (int w=0;w<2;w++) if ( doCut[d][w] ) cutVal = cutVals[d][w];

          double dx;
          // now check which vertical wall we cut
          if ( doCut[1-d][0] )  // left
            dx = cutVal - regLL[1-d];
          else                  // right
            dx = regUR[1-d] - cutVal;
           
          
          // this is an independent check that we do based on the vertical distances
          // it is copy-pasted code
          // if the method is correct then we should get the same answer
          if ( 0 ) {
            
            int d=0;
            double cutVal;  // this will be an x coord
            for (int w=0;w<2;w++) if ( doCut[d][w] ) cutVal = cutVals[d][w];

            double dx;
            // now check which vertical wall we cut
            if ( doCut[1-d][0] )  // left
              dx = cutVal - regLL[1-d];
            else                  // right
              dx = regUR[1-d] - cutVal;

            double check = dx/fabs( lineA[1-d] - lineB[1-d] );


          }
        }
        
      
      } else {  // one inside and one outside
        std::vector<double> inPoint  = ( insideA ? lineA : lineB );
        std::vector<double> outPoint = ( insideA ? lineB : lineA );
        
        double vecIO[dim];   // vector from inner point to outer point
        double dist2 = 0.0;
        for (int d=0;d<dim;d++) {
          vecIO[d] = outPoint[d] - inPoint[d];
        }
        //cout << "#vecIO "; Geom::printPoint( vecIO ); cout << endl;
        
        // distances from inner point to relevant edges
        // (ie the edge that may be between the two particles
        //  if the edge is not between then the compare number below is >1 and so its irrelevant)
        double vecEdge[dim];
        for (int d=0;d<dim;d++) {
          if ( vecIO[d] > 0 ) // outer is to right (or above) of inner
            vecEdge[d] = regUR[d] - inPoint[d];
          else
            vecEdge[d] = inPoint[d] - regLL[d];
        }
        //cout << "#edge "; Geom::printPoint( vecEdge ); cout << endl;
        
        double compare[2];
        for (int d=0;d<dim;d++) compare[d] = fabs( vecEdge[d] / vecIO[d] );
        
        // smaller of the two
        res = ( compare[1] < compare[0] ? compare[1] : compare[0]);
      }

      return res;
    }
