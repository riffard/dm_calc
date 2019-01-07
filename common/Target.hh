#ifndef Target_hh
#define Target_hh 1

#include<iostream>

#include "Nucleus.hh"

using namespace std;

class Target{

public:
  Target(string targetName);
  
  ~Target();


  void AddNucleus(double fraction, Nucleus* nucleus);
  void CheckNormalization();

  size_t size(){return fNuclei.size();}
  
  Nucleus* GetComponant(int i, double &fraction){
    fraction = fFraction[i];
    return fNuclei[i];
  }
  
private:

  string fTargetName;
  
  vector<Nucleus*> fNuclei;
  vector<double> fFraction;

  

};

#endif
