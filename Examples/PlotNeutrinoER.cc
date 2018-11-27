#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>

#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1D.h"
#include "THStack.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TMinuit.h"
#include "TColor.h"
#include "TLine.h"
#include "TLatex.h"
#include "TSystem.h"

using namespace std;

#include "Target.hh"

#include "HaloModel.hh"

#include "NeutrinoRate.hh"
#include "NeutrinoCrossSection_coherent_NR.hh"
#include "NeutrinoCrossSection_electroweak_ER.hh"
#include "NeutrinoFlux.hh"


int main(int argc, char** argv){

  TApplication* theApp = new TApplication("app", 0, 0);
  
  Target* fTarget;
  NeutrinoFlux* neutrino_fluxes;
  NeutrinoCrossSection* cross_section;
  NeutrinoRate* neutrino_rate;
  
  
  FormFactorDataBase* FormEvenNucleus = new FormFactorDataBase("Helm","Helm");
  FormFactorDataBase* Form129Xe = new FormFactorDataBase("Helm","129Xe_p");
  FormFactorDataBase* Form131Xe = new FormFactorDataBase("Helm","131Xe_p");

  double fraction, Mass, A, Z, Sp, Sn,J;

  fTarget = new Target("Xenon");
  fTarget->AddNucleus(fraction = 1e-3,   new Nucleus("Xe124", Mass = 123.9058958,  A = 124, Z = 54, J = 0,     Sp = 0,      Sn =  0, FormEvenNucleus));
  fTarget->AddNucleus(fraction = 9e-4,   new Nucleus("Xe126", Mass = 125.9042689,  A = 126, Z = 54, J = 0,     Sp = 0,      Sn =  0, FormEvenNucleus));
  fTarget->AddNucleus(fraction = 0.0191, new Nucleus("Xe128", Mass = 127.9035304,  A = 128, Z = 54, J = 0,     Sp = 0,      Sn =  0, FormEvenNucleus));

  fTarget->AddNucleus(fraction = 0.262,  new Nucleus("Xe129", Mass = 128.9047795,  A = 129, Z = 54, J = 0.5,   Sp = 0.010,  Sn =  0.329, Form129Xe));
  fTarget->AddNucleus(fraction = 0.041,  new Nucleus("Xe130", Mass = 129.9035079,  A = 130, Z = 54, J = 0,     Sp = 0,      Sn =  0, FormEvenNucleus));

  fTarget->AddNucleus(fraction = 0.218,  new Nucleus("Xe131", Mass = 130.9050819,  A = 77, Z = 54, J = 3./2., Sp = -0.009, Sn =  -0.272,Form131Xe));

  fTarget->AddNucleus(fraction = 0.269,  new Nucleus("Xe132", Mass = 131.9041545,  A = 132, Z = 54, J = 0,     Sp = 0,      Sn =  0, FormEvenNucleus));
  fTarget->AddNucleus(fraction = 0.104,  new Nucleus("Xe134", Mass = 133.9053945,  A = 134, Z = 54, J = 0,     Sp = 0,      Sn =  0, FormEvenNucleus));
  fTarget->AddNucleus(fraction = 0.089,  new Nucleus("Xe136", Mass = 135.9072195,  A = 136, Z = 54, J = 0,     Sp = 0,      Sn =  0, FormEvenNucleus));

  fTarget->CheckNormalization();

  
  //----------------------------------------------------
  // Neutrino Claculator parameters
  //----------------------------------------------------
  // Neutrino fluxes
  neutrino_fluxes = new NeutrinoFlux();
  // Neutrino cross section
  cross_section = new NeutrinoCrossSection_electroweak_ER;
  //cross_section = new NeutrinoCrossSection_coherent_NR;
  // Neutrino rate calculator
  neutrino_rate = new NeutrinoRate(fTarget, neutrino_fluxes, "All", cross_section);
  
  TH1D* hE = new TH1D("hE_input","Energy input spectrum",1000, 0, 100);
  
  neutrino_rate->GetRate(hE);
  

  TCanvas*c = new TCanvas;
  hE->Draw();

  c->SetLogx();
  c->SetLogy();
  c->Update();
  
  cout<<"Done"<<endl;
  theApp->Run();
  
  return 0; 
}
