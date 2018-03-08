#ifndef wimp_calculator_hh
#define wimp_calculator_hh 1

#include <iostream>
#include <string>

#include "Target.hh"
#include "HaloModel.hh"

using namespace std;


class wimp_calculator{

public:

  wimp_calculator(string CalcName, Target* target,  HaloModel* halo);
  ~wimp_calculator();
  

  double Rate_SI(double Er, double mChi, double sigma0Si);
  
private:

  static double Rate_SI(double E_r, double m_wimp,
			double m_nuc, double A_nuc, double Z_nuc,
			double xsec_wn);
 
  
private:


  
  string fCalcName;
  Target* fTarget;
  HaloModel* fHalo;
  
};



#endif

