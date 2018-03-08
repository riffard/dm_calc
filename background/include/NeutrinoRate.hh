#ifndef neutrino_rate_hh
#define neutrino_rate_hh 1


#include <iostream>
#include <vector>

#include "TH1D.h"
#include "THStack.h"
#include "TLegend.h"

#include "Target.hh"
#include "NeutrinoFlux.hh"
#include "NeutrinoCrossSection.hh"
#include "DetectorEfficiency.hh"

using namespace std;



class NeutrinoRate{

public:
  
  NeutrinoRate(Target* target, NeutrinoFlux* neutrino_flux, string required_fluxes, NeutrinoCrossSection* cross_section, DetectorEfficiency* efficiency);
  
  ~NeutrinoRate();

  void GetRate(TH1D* hE, THStack* hs = NULL, TLegend* leg=NULL);
  
private:


  Target* fTarget;
  NeutrinoFlux* fNeutrinoFluxDB;
  NeutrinoCrossSection* fCrossSection;
  DetectorEfficiency* fEfficiency;
  
  string fRequiredFluxes_str;
  vector<string> fRequiredFluxes;
  
};


#endif
