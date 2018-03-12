#ifndef NeutrinoFlux_hh
#define NeutrinoFlux_hh 1

#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "TGraph.h"

#include "Math/Interpolator.h"

using namespace std;

//-------------------------------------------------------
//-------------------------------------------------------
class Flux{

public:

  //-------------------------------------------------------
  //-------------------------------------------------------
  Flux(double scaling){
    fScaling = scaling;
    fInterp = NULL;

    fLine_energy = -1;
    fLine_Amplitude = -1;
    
    mode = "continuous";   
  }

  //-------------------------------------------------------
  //-------------------------------------------------------
  ~Flux(){  if(fInterp) delete fInterp;   }

  //-------------------------------------------------------
  //-------------------------------------------------------
  void Complete(){

    if(fEnergy.size()==1){
      fLine_energy = fEnergy[0];
      fLine_Amplitude = fAmplitude[0];
      gr->SetPoint(gr->GetN(), fLine_energy, 1e-20);
      mode = "line";
      return;
    }
    fInterp = new ROOT::Math::Interpolator(fEnergy, fAmplitude, ROOT::Math::Interpolation::kLINEAR);
  }

  //-------------------------------------------------------
  //-------------------------------------------------------
  void Push_back(double E, double A){
    fEnergy.push_back(E);
    fAmplitude.push_back( A * fScaling);

    gr->SetPoint(gr->GetN(), E, A * fScaling);
  }

  double GetTotalFlux(){return fScaling;} // neutrinos/cm^2/s
  
  //-------------------------------------------------------
  //-------------------------------------------------------
  string GetMode(){return mode;}
  
  //-------------------------------------------------------
  //-------------------------------------------------------
  void GetLine(double &E, double &A){
    E = fLine_energy;
    A = fLine_Amplitude;
  }
    
  //-------------------------------------------------------
  //-------------------------------------------------------
  double Interp(double E){
    if(!fInterp) return -1;
    if(E < fEnergy[0] || E > fEnergy.back()) return 0;
    return fInterp->Eval(E);
  }

  //-------------------------------------------------------
  //-------------------------------------------------------
  TGraph* GetGraph(){ return gr; }


  vector<double>* GetEnergy(){return &fEnergy;};
  vector<double>* GetFlux(){return &fAmplitude;};

  
private:
    
  vector<double> fEnergy;    // in MeV
  vector<double> fAmplitude; // in MeV^-1

  double fLine_energy;
  double fLine_Amplitude;

  double fScaling;           // in neutrinos/cm^2/s
  
  
  ROOT::Math::Interpolator* fInterp;

  TGraph* gr = new TGraph();

  string mode;
  
  
};

//-------------------------------------------------------
//-------------------------------------------------------


//-------------------------------------------------------
//-------------------------------------------------------
class NeutrinoFlux{

public:

  
  //-------------------------------------------------------
  //-------------------------------------------------------
  NeutrinoFlux();
  
  //-------------------------------------------------------
  //-------------------------------------------------------
  ~NeutrinoFlux();
  
  void PrintAvailableData();
  void DrawFluxCollection();

  vector<string> GetFluxCollectionName(){return fFluxFiles;}

  Flux* GetFlux(string key){return fFlux_collection[key];}
  
private:

  string data_base_path;

  vector<string> fFluxFiles;
  vector<double> fScale_constant;

  map<string, Flux*> fFlux_collection;

  
};
//-------------------------------------------------------
//-------------------------------------------------------

#endif
