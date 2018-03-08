#ifndef HaloModel_h
#define HaloModel_h 1

#include <iostream>


#include "TMath.h"


using namespace std;


class HaloModel{
  
 public:
  HaloModel(double v_earth, double sigma_v, double v_esc, string model = "SMH");
  
  ~HaloModel();

  
  double getHaloIntegrationValue(double v_min);
  
  Double_t getHaloIntegrationValue_TF1(Double_t* x, Double_t* par){ return getHaloIntegrationValue(x[0]);};

  double get_vEarth(){return fv_earth;};

  static double getHaloIntegrationValue_SMH(double v_min, double v_earth, double sigma_v, double v_esc);
  
  static double getHaloIntegrationValue_SMHcutoff(double v_min, double v_earth, double sigma_v, double v_esc);

  
 private:


  double fv_earth, fsigma_v, fv_esc;

  double (*pointeurSurModelHalo)(double,double,double,double);
  
};


#endif
