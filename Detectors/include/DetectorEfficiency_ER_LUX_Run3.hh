#ifndef DetectorEfficiency_ER_LUX_Run3_hh
#define DetectorEfficiency_ER_LUX_Run3_hh 1

#include "DetectorEfficiency.hh"

#include "TMath.h"
using namespace std;


class DetectorEfficiency_ER_LUX_Run3: public DetectorEfficiency{

public:
  DetectorEfficiency_ER_LUX_Run3(){

    fName = "ER_LUX_Run3";
    
    fMu = 1.24;
    fSigma = 0.43;
  }
  
  virtual ~DetectorEfficiency_ER_LUX_Run3(){}

  virtual double GetEfficiency(double E){
    return 0.5 + 0.5*TMath::Erf( (E - fMu) / sqrt(2) / fSigma);    
  }


private:
  virtual void PrepareGraph(){

    fGrEfficiency = new TGraph();
    for(double E = 0; E < 100; E +=0.01) fGrEfficiency->SetPoint( fGrEfficiency->GetN(), E, GetEfficiency(E));
   
  }


private:

  double fMu, fSigma;;
  

};

#endif
