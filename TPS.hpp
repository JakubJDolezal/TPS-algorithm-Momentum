#include <vector>
#include <array>
#include <random>
#include <chrono> 
#include <complex>
#include "Underdamp.hpp"

class TPS {
private:
    bool random_bool();
    int random_int(int);
    int random_int_no_zero(int);
    int random_int_non_neg(int);
    int random_int_ran(int, int);
    double random_real(double);
    double random_real_range(double low, double high);
    double random_real_neg(double);
    
protected:
    int timeInts,M;
    int noBoxes;
    double s;
    std::vector<std::vector<double>> posXStorage;
    std::vector<std::vector<double>> posYStorage;
    std::vector<std::vector<double>> momXStorage;
    std::vector<std::vector<double>> momYStorage;
    std::vector<double> storedvelocities;
    std::vector<std::vector<double>> virialStorage;
    std::vector<std::vector<double>> densityStorage;
    std::vector<std::vector<double>> forceWallStorage;
    std::vector<double> acceptance; // this holds the fluctuations of the acceptance rate
    std::vector<double> entropicPath; //holds the fluctuations of the the entropy
    std::vector<double> actBetweenJumps; //the activity of all the segments
    std::vector<double> entropyStorage; //keeps the entropic costs of all the invidual segmensts
    std::vector<double> proposedAct; //all the proposed activities
    std::vector<double> proposedEntropy; //all the proposed entropies
    unsigned int origBlab; //the seed
    double force; // the biasing force
    std::mt19937 & gen;
    std::uniform_real_distribution<double> realgen;
    std::uniform_int_distribution<int>boolgen;
    std::uniform_int_distribution<int> intgen;
 
public:
 TPS(int M1, int timeInterval, std::mt19937 & gen1):
 M(M1),
 timeInts(timeInterval),
 gen(gen1)
 {
   std::uniform_int_distribution<int>boolgen1(0,1);
   boolgen=boolgen1;
   std::uniform_real_distribution<double> realgen1(0, 1);
   realgen=realgen1;
   std::uniform_int_distribution<int> intgen1(0, 10);
   intgen=intgen1;
};
 ~TPS()
 {};
 
 //sampling
  double sampling_shifting_higha_force(std::ofstream& pathEntropy, std::ofstream& pathEntChange, std::ofstream& pathActchange, int iter, int cut, Underdamp&); //performs sampling with deleting front or back M-a sections and then shifting the entire thing forward or backward respectively and regenerating it, a=M-11, assumes no force is present
  double sampling_shifting_higha(std::ofstream&, int, int, Underdamp&); //performs sampling with deleting front or back M-a sections and then shifting the entire thing forward or backward respectively and regenerating it, a=M-11
  void sampling_with_ft(int,double, Underdamp&, int);
  
 //creation/loading
 void create_stuff(Underdamp&);
 void load_stuff(double,Underdamp&);
 void distrbutionGeneration();


  };
  
