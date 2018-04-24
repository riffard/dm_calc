#ifndef NeutrinoCrossSection_hh
#define NeutrinoCrossSection_hh 1

#include <iostream>

#include "Nucleus.hh"


using namespace std;


class NeutrinoCrossSection{

public:

  NeutrinoCrossSection(){
    fCrossSectionName = "";
  }
  virtual ~NeutrinoCrossSection(){}

  
  virtual double GetCrossSection(double E_recoil, double E_nrutrino, Nucleus* nucleus){ return 0;};
  
protected:

  string fCrossSectionName;

  
};


#endif
