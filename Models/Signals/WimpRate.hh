#ifndef WimpRate_hh
#define WimpRate_hh 1

#include <iostream>
#include <string>

#include "TH1D.h"

#include "Target.hh"
#include "HaloModel.hh"
#include "VDetectorEfficiency.hh"

using namespace std;


class WimpRate{

public:

  WimpRate(string CalcName, Target* target,  HaloModel* halo, VDetectorEfficiency* efficiency = NULL);
  ~WimpRate();
  

  double Rate_SI(double Er, double mChi, double sigma0Si);
  void GetRate(TH1D* h, double mChi, double sigma0Si);
  
private:

  static double Rate_SI(double E_r, double m_wimp,
			double m_nuc, double A_nuc, double Z_nuc,
			double xsec_wn);
 
  
private:


  
  string fCalcName;
  Target* fTarget;
  HaloModel* fHalo;
  VDetectorEfficiency* fEfficiency;  
};



#endif

