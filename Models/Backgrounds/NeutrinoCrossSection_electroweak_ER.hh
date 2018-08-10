#ifndef NeutrinoCrossSection_electroweak_ER_hh
#define NeutrinoCrossSection_electroweak_ER_hh 1

#include "NeutrinoCrossSection.hh"

#include "TMath.h"

class NeutrinoCrossSection_electroweak_ER: public NeutrinoCrossSection{


public:

  NeutrinoCrossSection_electroweak_ER(){

    fCrossSectionName = "electroweak_ER";
    
  }
  
  virtual ~NeutrinoCrossSection_electroweak_ER(){}
  
  
  virtual double GetCrossSection(double E_recoil, double E_neutrino, Nucleus* nucleus ){

    double A_nuc = nucleus->GetA();
    double Z_nuc = nucleus->GetZ();
    double m_nuc = nucleus->GetMGeV()*1e6; // conversion in keV

    // Define constants
    double Gf = 1.17e-17;      // keV^-2 (From PDG)
    double sin2thetaW = 0.231; // From PDG
    double m_e = 511.0;

    // Calculate cross-section for mu and tau neutrinos
    double gv = 2 * sin2thetaW - 0.5;
    double ga = 0.5;
    double dxsecdEr_mutau = (Gf*Gf * m_e) / (2 * TMath::Pi()) * 
        (pow(gv+ga,2) + pow(gv-ga,2) * pow(1 - E_recoil/E_neutrino,2) + (ga*ga - gv*gv) * m_e * E_recoil / pow(E_neutrino,2));

    // Calculate cross-section for electron neutrinos
    gv = 2 * sin2thetaW + 0.5;
    ga = 1.5;
    double dxsecdEr_electron = (Gf*Gf * m_e) / (2 * TMath::Pi()) *
        (pow(gv+ga,2) + pow(gv-ga,2) * pow(1 - E_recoil/E_neutrino,2) + (ga*ga - gv*gv) * m_e * E_recoil / pow(E_neutrino,2));

    // Calculate full cross-section, weighted for survival probability
    double dxsecdEr = 0.55 * dxsecdEr_electron + 0.45 * dxsecdEr_mutau;
    dxsecdEr *= Z_nuc; // Normalize to number of electrons in atom
    dxsecdEr /= 3.88e-16; // cm^2/keV

    if (dxsecdEr < 0) return 0;
    else return dxsecdEr;
    
    return 0;}
  
  
};


#endif
