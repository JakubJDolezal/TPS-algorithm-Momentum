#include "TPS.hpp"
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
using namespace std;


double TPS::sampling_shifting_higha_force(std::ofstream& pathEntropy,std::ofstream& pathEntChange,std::ofstream& pathActchange, int iter, int cut,  Underdamp& MeasSys) // the main method for sampling the distribution described in the paper with a force, it does so with shifting, erases M-a things and then shifts the distro and generates more
{
  bool direction=boolgen(gen);
  int a=intgen(gen)+M-16;
  double booya=0;
  std::vector<std::vector<double> > stuff2(M, (std::vector<double>(MeasSys.get_noOfPart(),0)));
  std::vector<double> stuff(M,0);
  std::vector<std::vector<double>> newmomXStorage=stuff2;
  std::vector<std::vector<double>> newposYStorage=stuff2;
  std::vector<std::vector<double>> newposXStorage=stuff2;
  std::vector<std::vector<double>> newmomYStorage=stuff2;
  std::vector<double> newA=stuff;
  std::vector<double> newE=stuff;
  newA.reserve(M);
  newE.reserve(M);
  double oldEntropy=0;
  double restOfEntropy=0;
  double newEntropy=0;
  double oldact=0; 	//this is the old activity
  double restOfAct=0;  
  double newact=0; 	//this is the new activity
  if(direction==true) //forward shift where we delete a part in front and add a part in the back
  {
    for (int i=0;i<M-a;i++)
    {
      oldEntropy=oldEntropy+entropyStorage[i];
      oldact=oldact+actBetweenJumps[i];
    }
    int k=0;
    for (int i=M-a;i<M;i++)
    {
      restOfAct=restOfAct+actBetweenJumps[i];
      restOfEntropy=restOfEntropy+entropyStorage[i];
      newposXStorage[k]=posXStorage[i];
      newposYStorage[k]=posYStorage[i];
      newmomXStorage[k]=momXStorage[i];
      newmomYStorage[k]=momYStorage[i];
      newA[k]=actBetweenJumps[i];
      newE[k]=entropyStorage[i];
      k++;
    }
    MeasSys.set_posX(posXStorage[posXStorage.size()-1]);
    MeasSys.set_posY(posYStorage[posYStorage.size()-1]);
    MeasSys.set_momX(momXStorage[momXStorage.size()-1]);
    MeasSys.set_momY(momYStorage[momYStorage.size()-1]);
    for (int i=0;i<M-a;i++)
    {
      double local=MeasSys.repeat_update(timeInts);
      newposXStorage[k]=MeasSys.get_posX();
      newposYStorage[k]=MeasSys.get_posY();
      newmomXStorage[k]=MeasSys.get_momX();
      newmomYStorage[k]=MeasSys.get_momY();
      newA[k]=local;
      newE[k]=MeasSys.get_entropy();
      newact=newact+local;
      newEntropy=newEntropy+MeasSys.get_entropy();
      k++;
    }
    if (force==0)
    {
	newEntropy=0;
	oldEntropy=0;
    }   
    if(exp(newEntropy-oldEntropy+(oldact-newact)*s)>realgen(gen)) // here it sees if it will accept it, if it does the new activity is added to the average and the generated candidate of slices replaces the old slices
    {
      posXStorage=newposXStorage;
      posYStorage=newposYStorage;
      momXStorage=newmomXStorage;
      momYStorage=newmomYStorage;
      booya=newact+restOfAct;
      actBetweenJumps=newA;
      entropyStorage=newE;
    }
    else
    {
      booya=oldact+restOfAct;
    }
  }
  
  else //backward shift where we delete a part in the back and add a part in the front
  {
    MeasSys.set_posX(posXStorage[0]);
    MeasSys.set_posY(posYStorage[0]);
    MeasSys.set_momX(momXStorage[0]);
    MeasSys.set_momY(momYStorage[0]);
    for (int i=1;i<M-a+1;i++) //we then fill it in with the new parts as we generate them 
    {
      double local=MeasSys.repeat_update(timeInts);
      newposXStorage[M-a-i]=MeasSys.get_posX();
      newposYStorage[M-a-i]=MeasSys.get_posY();
      newmomXStorage[M-a-i]=MeasSys.get_momX();
      newmomYStorage[M-a-i]=MeasSys.get_momY();
      newA[M-a-i]=local;
      newE[M-a-i]=MeasSys.get_entropy();
      newact=newact+local;
      newEntropy=newEntropy+MeasSys.get_entropy();
    }
    for (int i=a;i<M;i++)
    {
    oldEntropy=oldEntropy+entropyStorage[i];
    oldact=oldact+actBetweenJumps[i];
    }
    int k=M-a;
    for (int i=0;i<a;i++) //we fill the candidate in from the end with the unchanged parts
    {
      restOfEntropy=restOfEntropy+entropyStorage[i];
      restOfAct=restOfAct+actBetweenJumps[i];
      newposXStorage[k]=posXStorage[i];
      newposYStorage[k]=posYStorage[i];
      newmomXStorage[k]=momXStorage[i];
      newmomYStorage[k]=momYStorage[i];
      newA[k]=actBetweenJumps[i];
      newE[k]=entropyStorage[i];
      k++;
    }
    if (force==0)
    {
	newEntropy=0;
	oldEntropy=0;
    }
    if(exp(newEntropy-oldEntropy+(oldact-newact)*s)>realgen(gen)) // here it sees if it will accept it, if it does the new activity is added to the average and the generated candidate of slices replaces the old slices
    {
      posXStorage=newposXStorage;
      posYStorage=newposYStorage;
      momXStorage=newmomXStorage;
      momYStorage=newmomYStorage;
      booya=newact+restOfAct;
      actBetweenJumps=newA;
      entropyStorage=newE;
    }
    else
    {
      booya=oldact+restOfAct;
    }
  }
  return booya;
}

double TPS::sampling_shifting_higha(std::ofstream& pathActchange, int iter, int cut,  Underdamp& MeasSys) // the main method for sampling the distribution described in the paper, it does so with shifting, erases M-a things and then shifts the distro and generates more
{
  bool direction=boolgen(gen);
  int a=intgen(gen)+M-16;
  double booya=0;
  std::vector<std::vector<double> > stuff2(M, (std::vector<double>(MeasSys.get_noOfPart(),0)));
  std::vector<double> stuff(M,0);
  std::vector<std::vector<double>> newmomXStorage=stuff2;
  std::vector<std::vector<double>> newposYStorage=stuff2;
  std::vector<std::vector<double>> newposXStorage=stuff2;
  std::vector<std::vector<double>> newmomYStorage=stuff2;
  std::vector<double> newA=stuff;
  std::vector<double> newE=stuff;
  std::vector<std::vector<double> > stuff4(M, (std::vector<double>(2,0)));
  std::vector<std::vector<double>> newFOW=stuff4;
  newA.reserve(M);
  double oldEntropy=0;
  double restOfEntropy=0;
  double newEntropy=0;
  double oldact=0; 	//this is the old activity
  double restOfAct=0;  
  double newact=0; 	//this is the new activity
  if(direction==true) //forward shift where we delete a part in front and add a part in the back
  {
    for (int i=0;i<M-a;i++)
    {
      oldEntropy=oldEntropy+entropyStorage[i];
      oldact=oldact+actBetweenJumps[i];
    }
    int k=0;
    for (int i=M-a;i<M;i++)
    {
      restOfAct=restOfAct+actBetweenJumps[i];
      restOfEntropy=restOfEntropy+entropyStorage[i];
      newposXStorage[k]=posXStorage[i];
      newposYStorage[k]=posYStorage[i];
      newmomXStorage[k]=momXStorage[i];
      newmomYStorage[k]=momYStorage[i];
      newA[k]=actBetweenJumps[i];
      newFOW[k][0]=MeasSys.forceOnWall[0];
      newFOW[k][1]=MeasSys.forceOnWall[1];
      newE[k]=entropyStorage[i];
      k++;
    }
    MeasSys.set_posX(posXStorage[posXStorage.size()-1]);
    MeasSys.set_posY(posYStorage[posYStorage.size()-1]);
    MeasSys.set_momX(momXStorage[momXStorage.size()-1]);
    MeasSys.set_momY(momYStorage[momYStorage.size()-1]);
    for (int i=0;i<M-a;i++)
    {
      double local=MeasSys.repeat_update(timeInts);
      newposXStorage[k]=MeasSys.get_posX();
      newposYStorage[k]=MeasSys.get_posY();
      newmomXStorage[k]=MeasSys.get_momX();
      newmomYStorage[k]=MeasSys.get_momY();
      newFOW[k]=forceWallStorage[k];
      newA[k]=local;
      newE[k]=MeasSys.Entropy;
      newact=newact+local;
      newEntropy=newEntropy+MeasSys.Entropy;
      k++;
    }  
    if(exp(newEntropy-oldEntropy+(oldact-newact)*s)>realgen(gen)) // here it sees if it will accept it, if it does the new activity is added to the average and the generated candidate of slices replaces the old slices
    {
      posXStorage=newposXStorage;
      posYStorage=newposYStorage;
      momXStorage=newmomXStorage;
      momYStorage=newmomYStorage;
      booya=newact+restOfAct;
      actBetweenJumps=newA;
      entropyStorage=newE;
      forceWallStorage=newFOW;
    }
    else
    {
      booya=oldact+restOfAct;
    }
  }
  
  else //backward shift where we delete a part in the back and add a part in the front
  {
    MeasSys.set_posX(posXStorage[0]);
    MeasSys.set_posY(posYStorage[0]);
    MeasSys.set_momX(momXStorage[0]);
    MeasSys.set_momY(momYStorage[0]);
    for (int i=1;i<M-a+1;i++) //we then fill it in with the new parts as we generate them 
    {
      double local=MeasSys.repeat_update(timeInts);
      newposXStorage[M-a-i]=MeasSys.get_posX();
      newposYStorage[M-a-i]=MeasSys.get_posY();
      newmomXStorage[M-a-i]=MeasSys.get_momX();
      newmomYStorage[M-a-i]=MeasSys.get_momY();
      newA[M-a-i]=local;
      newE[M-a-i]=MeasSys.Entropy;
      newact=newact+local;
      newEntropy=newEntropy+MeasSys.Entropy;
      newFOW[M-a-i][0]=MeasSys.forceOnWall[0];
      newFOW[M-a-i][1]=MeasSys.forceOnWall[1];
    }
    for (int i=a;i<M;i++)
    {
    oldEntropy=oldEntropy+entropyStorage[i];
    oldact=oldact+actBetweenJumps[i];
    }
    int k=M-a;
    for (int i=0;i<a;i++) //we fill the candidate in from the end with the unchanged parts
    {
      restOfEntropy=restOfEntropy+entropyStorage[i];
      restOfAct=restOfAct+actBetweenJumps[i];
      newposXStorage[k]=posXStorage[i];
      newposYStorage[k]=posYStorage[i];
      newmomXStorage[k]=momXStorage[i];
      newmomYStorage[k]=momYStorage[i];
      newA[k]=actBetweenJumps[i];
      newE[k]=entropyStorage[i];
      newFOW[k]=forceWallStorage[k];
      k++;
    }
    if(exp(newEntropy-oldEntropy+(oldact-newact)*s)>realgen(gen)) // here it sees if it will accept activity is added to the average and the generated candidate of slices replaces the old slices
    {
      posXStorage=newposXStorage;
      posYStorage=newposYStorage;
      momXStorage=newmomXStorage;
      momYStorage=newmomYStorage;
      booya=newact+restOfAct;
      actBetweenJumps=newA;
      entropyStorage=newE;
      forceWallStorage=newFOW;
    }
    else
    {
      booya=oldact+restOfAct;
    }
  }
  return booya;
}

//TPS sampling as popularised by Chandler et al.
void TPS::sampling(int iterations,double svalue,  Underdamp& MeasSys, int id)
{
  s=svalue;
  int cut=4;
  std::vector<double> meanAct(2,0);
  std::vector<double> meanForceOnWall(2,0);
  iterations=iterations/cut;
  string sStr=to_string(s).substr(0,5);
  string txt=to_string(id).substr(0,2)+".txt";
  string local= "Path";
  local.append(sStr);
  local.append(txt);
  std::ofstream path1(local);
  string local4="PropActChange";
  local4.append(sStr);
  local4.append(txt);
  std::ofstream propActChange(local4);  
  string local8= "Kave"; 
  local8.append(sStr);
  local8.append(txt);
  std::ofstream pathKave(local8);
  string local9= "Kend"; 
  local9.append(sStr);
  local9.append(txt);
  std::ofstream pathKend(local9);
  string local8a= "Aave"; 
  local8a.append(sStr);
  local8a.append(txt);
  std::ofstream pathAave(local8a);
  string local9a= "Aend"; 
  local9a.append(sStr);
  local9a.append(txt);
  std::ofstream pathAend(local9a);
  MeasSys.set_stress();
  for (int i=0;i<iterations*cut;i++)
  {
//     double point=sampling_shifting_higha(pathEntropy,pathEntChange,propActChange, i, cut, MeasSys);
    double point=sampling_shifting_higha(propActChange, i, cut, MeasSys);

    if (i%cut==0)
    {
      MeasSys.set_posX(posXStorage[M/2]);
      MeasSys.set_posY(posYStorage[M/2]);
      MeasSys.set_momX(momXStorage[M/2]);
      MeasSys.set_momY(momYStorage[M/2]);
      MeasSys.Irving_Kirkwood_calc();
      meanAct[0]+=actBetweenJumps[M/2];
      meanAct[1]+=actBetweenJumps[M/2]*actBetweenJumps[M/2];
      path1 << point << '\n';
      pathKave << actBetweenJumps[M/2]<< '\n';
      pathKend << actBetweenJumps[0]<< '\n';
      meanForceOnWall[0]=forceWallStorage[M/2][0];
      meanForceOnWall[1]=forceWallStorage[M/2][1];
      for (int k=0;k<momXStorage[0].size();k++)
      {
	storedvelocities.push_back(momXStorage[M/2][k]);
      }
       pathAend << entropyStorage[0]<< '\n';
       pathAave << entropyStorage[M/2]<< '\n';
    }
  }
  pathKave.close();
  pathKend.close();
  pathAave.close();
  pathAend.close();
  path1.close();
  propActChange.close();
  std::vector<double> meanStressVec=MeasSys.get_normStress();
  for(int i=0;i<meanAct.size();i++)
  {
  meanAct[i]/=(double)(iterations);
  }
  meanForceOnWall[0]/=(double)(iterations);
  meanForceOnWall[1]/=(double)(iterations); 
  {
    string local5= "MeanAct";
    local5.append(sStr);
    local5.append(txt);
    std::ofstream meanActOut(local5);
    for (std::vector<double>::const_iterator j = meanAct.begin(); j != meanAct.end(); ++j)
    {
      meanActOut << *j << '\n';
    }
    meanActOut.close();
  }
 {
    string localStress= "MeanStress";
    localStress.append(sStr);
    localStress.append(txt);
    std::ofstream meanStress(localStress);
    for (std::vector<double>::const_iterator j = meanStressVec.begin(); j != meanStressVec.end(); ++j)
    {
      meanStress << *j << '\n';
    }
    meanStress.close();
  }
  std::vector<double> meanPotEn=MeasSys.get_normPotEn();
 {
    string localStress= "MeanPotEn";
    localStress.append(sStr);
    localStress.append(txt);
    std::ofstream meanStress(localStress);
    for (std::vector<double>::const_iterator j = meanPotEn.begin(); j != meanPotEn.end(); ++j)
    {
      meanStress << *j << '\n';
    }
    meanStress.close();
  }
  std::vector<double> meanKinEn=MeasSys.get_normkinEnLoc();
 {
    string localStress= "MeanKinEn";
    localStress.append(sStr);
    localStress.append(txt);
    std::ofstream meanStress(localStress);
    for (std::vector<double>::const_iterator j = meanKinEn.begin(); j != meanKinEn.end(); ++j)
    {
      meanStress << *j << '\n';
    }
    meanStress.close();
  }
  std::vector<double> meanBodyForce=MeasSys.get_normBodyForce();
 {
    string localStress= "MeanBodyForce";
    localStress.append(sStr);
    localStress.append(txt);
    std::ofstream meanStress(localStress);
    for (std::vector<double>::const_iterator j = meanBodyForce.begin(); j != meanBodyForce.end(); ++j)
    {
      meanStress << *j << '\n';
    }
    meanStress.close();
  }
  {
    string local6= "MeanForceOnWall";
    local6.append(sStr);
    local6.append(txt);
    std::ofstream meanVirialOut(local6);
    for (std::vector<double>::const_iterator j = meanForceOnWall.begin(); j != meanForceOnWall.end(); ++j)
    {
      meanVirialOut << *j << '\n';
    }
    meanVirialOut.close();
  }
  std::vector<double> meanDensVec=MeasSys.get_normDens();
  {
    string local6= "MeanDens";
    local6.append(sStr);
    local6.append(txt);
    std::ofstream meanDensityOut(local6);
    for (std::vector<double>::const_iterator j = meanDensVec.begin(); j != meanDensVec.end(); ++j)
    {
      meanDensityOut << *j << '\n';
    }
    meanDensityOut.close();
  }

  std::vector<double> mean(M,0);
  std::vector<double> Correl(M,0);
  local="load/ZAct";
  local.append(sStr);
  local.append(txt);
  string localEnt="load/ZEnt";
  localEnt.append(sStr);
  localEnt.append(txt);
  std::ofstream MAct(local);
  std::ofstream MEnt(localEnt);
  string local3x;
  string local3y;
  string local3ym;
  string local3xm;
  
  for (int i=0;i<M;i++)
  {
    MAct << actBetweenJumps[i] << '\n';
    MEnt << entropyStorage[i] << '\n';
    local3x="load/ZPosX";
    local3x.append(to_string(i));
    local3x.append("s");
    local3x.append(sStr);
    local3x.append(txt);
    local3y="load/ZPosY";
    local3y.append(to_string(i));
    local3y.append("s");
    local3y.append(sStr);
    local3y.append(txt);
    local3ym="load/ZMomY";
    local3ym.append(to_string(i));
    local3ym.append("s");
    local3ym.append(sStr);
    local3ym.append(txt);
    local3xm="load/ZMomX";
    local3xm.append(to_string(i));
    local3xm.append("s");
    local3xm.append(sStr);
    local3xm.append(txt);
    std::ofstream xState(local3x);
    std::ofstream yState(local3y);
    std::ofstream xmState(local3xm);
    std::ofstream ymState(local3ym);
    for (int j=0;j<MeasSys.get_noOfPart();j++)
    {
      xState << posXStorage[i][j] << '\n';
      yState << posYStorage[i][j] << '\n';
      ymState << momYStorage[i][j] << '\n';
      xmState << momXStorage[i][j] << '\n';
    }
    xState.close();
    yState.close();
    ymState.close();
    xmState.close();
  }
  MAct.close();
  MEnt.close();
}

//To load  the previously saved TPS vectors, load must exist
void TPS::load_stuff(double s, Underdamp& MeasSys)
{
    noBoxes=MeasSys.noOfBoxesX;
    string sStr=to_string(s).substr(0,5);
    string txt="0.txt";
    string local="load/ZAct";
    local.append(sStr);
    local.append(txt);
    std::ifstream MAct(local);
    std::vector<double>randomVector {std::istream_iterator<double>(MAct), {}};
    actBetweenJumps=randomVector;
    string localEnt="load/ZEnt";
    localEnt.append(sStr);
    localEnt.append(txt);
    std::ifstream MEnt(localEnt);
    std::vector<double>randomVectorEnt {std::istream_iterator<double>(MEnt), {}};
    entropyStorage=randomVectorEnt;
    MAct.close();
    MEnt.close();
    posXStorage.reserve(M);
    posYStorage.reserve(M);
    momXStorage.reserve(M);
    momYStorage.reserve(M);
    string local3x;
    string local3y;
    string local3ym;
    string local3xm;
    for (int i=0;i<M;i++)
    {      
      local3x="load/ZPosX";
      local3x.append(to_string(i));
      local3x.append("s");
      local3x.append(sStr);
      local3x.append(txt);
      local3y="load/ZPosY";
      local3y.append(to_string(i));
      local3y.append("s");
      local3y.append(sStr);
      local3y.append(txt);
      local3ym="load/ZMomY";
      local3ym.append(to_string(i));
      local3ym.append("s");
      local3ym.append(sStr);
      local3ym.append(txt);
      local3xm="load/ZMomX";
      local3xm.append(to_string(i));
      local3xm.append("s");
      local3xm.append(sStr);
      local3xm.append(txt);
      std::ifstream xState(local3x);
      std::ifstream yState(local3y);
      std::ifstream xmState(local3xm);
      std::ifstream ymState(local3ym);
      std::vector<double>randomVector1 {std::istream_iterator<double>(xState), {}};
      std::vector<double>randomVector1m {std::istream_iterator<double>(xmState), {}};
      std::vector<double>randomVector2 {std::istream_iterator<double>(yState), {}};
      std::vector<double>randomVector2m {std::istream_iterator<double>(ymState), {}};
      posXStorage.push_back(randomVector1);
      momXStorage.push_back(randomVector1m);
      posYStorage.push_back(randomVector2);
      momYStorage.push_back(randomVector2m);
      xState.close();
      yState.close();
      ymState.close();
      xmState.close();
    }
   // std::vector<std::vector<double> > stuff3(M, (std::vector<double>(noBoxes,0)));
   // densityStorage=stuff3;
    //virialStorage=stuff3;
    std::vector<std::vector<double> > stuff4(M, (std::vector<double>(2,0)));
    forceWallStorage=stuff4;
    MeasSys.repeat_update(timeInts);
}


//To create the TPS vectors
void TPS::create_stuff(Underdamp& MeasSys)
{
  noBoxes=MeasSys.noOfBoxesX;
  std::vector<std::vector<double> > stuff3(M, (std::vector<double>(noBoxes,0)));
  densityStorage=stuff3;
  virialStorage=stuff3;
  for (int i=0;i<M;i++)
  {
    actBetweenJumps.push_back(MeasSys.repeat_update(timeInts));
    posXStorage.push_back(MeasSys.get_posX());
    posYStorage.push_back(MeasSys.get_posY());
    momXStorage.push_back(MeasSys.get_momX());
    momYStorage.push_back(MeasSys.get_momY());
    entropyStorage.push_back(MeasSys.Entropy);
    std::vector<std::vector<double> > stuff4(M, (std::vector<double>(2,0)));
    forceWallStorage=stuff4;
  }
}

//If we want to sample the momentum distribution
void TPS::distrbutionGeneration()
{
  std::vector<double> velocities(std::vector<double>(500,0));
  for (int i=0;i<velocities.size();i++)
  {
    velocities[(int)round(storedvelocities[i]/100)]=velocities[(int)round(storedvelocities[i]/100)]+1;
  }
  string sStr=to_string(s).substr(0,4);
  string txt=".txt";
  string local3x="MomentaDistributed";
  local3x.append("s");
  local3x.append(sStr);
  local3x.append(txt);
  std::ofstream velocitiesOut(local3x);
  for (int j=0;j<velocities.size();j++)
  {
    velocitiesOut << velocities[j] << '\n';
  }
}
