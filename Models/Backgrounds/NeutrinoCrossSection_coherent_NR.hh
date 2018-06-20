#ifndef NeutrinoCrossSection_coherent_NR_hh
#define NeutrinoCrossSection_coherent_NR_hh 1

#include "NeutrinoCrossSection.hh"

#include "TMath.h"

class NeutrinoCrossSection_coherent_NR: public NeutrinoCrossSection{


public:

  NeutrinoCrossSection_coherent_NR(){

    fCrossSectionName = "coherent_NR";
    
  }
  
  virtual ~NeutrinoCrossSection_coherent_NR(){}
  
  
  virtual double GetCrossSection(double E_recoil, double E_neutrino, Nucleus* nucleus ){

    double A_nuc = nucleus->GetA();
    double Z_nuc = nucleus->GetZ();
    double m_nuc = nucleus->GetMGeV()*1e6; // conversion in keV

    double Gf = 1.17e-17;      // keV^-2 (From PDG)
    double sin2thetaW = 0.231; // From PDG
    double rn0 = 5.79e-6;      // 1/keV = 1.14 fm
    double s = 4.57e-6;        // 1/keV = 0.9 fm
    
    double Qw = A_nuc - Z_nuc - (1 - 4 * sin2thetaW)*Z_nuc;

    double q = sqrt(2 * m_nuc * E_recoil);
    double rn = rn0 * pow(A_nuc, 1./3.);
    double bessel = sin(q*rn) / pow(q*rn,2) - cos(q*rn) / (q*rn);
    double F = 3 * bessel * exp(-pow(q*s,2) / 2) / (q*rn);

    double dxsecdEr = Gf*Gf/(4*TMath::Pi()) * Qw * Qw * m_nuc * (1 - m_nuc*E_recoil/(2*E_neutrino*E_neutrino)) * F * F;
    dxsecdEr = dxsecdEr * 3.88e-16; // cm^2/keV

    if (dxsecdEr < 0) return 0;
    else return dxsecdEr;
    
    

    return 0;}
  
  
};


#endif
