#ifndef FormFactorDataBase_h
#define FormFactorDataBase_h 1

#include <iostream>
#include <math.h>
#include "TF1.h"

using namespace std;

class FormFactorDataBase{

 public:

  FormFactorDataBase(string FormSIcode, string FormSDcode);
  ~FormFactorDataBase();
  
  double GetFormSD(double q, double A){
    if(pointeurSurFormFactorSD==NULL) return 0;
    else return (*pointeurSurFormFactorSD)(q,A, this);
  };

 
  double GetFormSI(double q, double A){
    if(pointeurSurFormFactorSI==NULL) return 0;
    else return (*pointeurSurFormFactorSI)(q,A, this);
  };


  double  FormFactorSI_TF1(double *x, double* par){

    double mTarget_GeV = par[0];
    double A = par[1];
    
    double q = 7.1668e-3 * sqrt( mTarget_GeV * x[0] );
    
    if(pointeurSurFormFactorSI==NULL) return 0;
    else return (*pointeurSurFormFactorSI)(q,A, this);
  }

   double FormFactorSD_TF1(double *x, double* par){

    double mTarget_GeV = par[0];
    double A = par[1];
    
    double q = 7.1668e-3 * sqrt( mTarget_GeV * x[0] );
    
    if(pointeurSurFormFactorSD==NULL) return 0;
    else return (*pointeurSurFormFactorSD)(q,A, this);
  }
  
  void Set_apan_fpfn(double inapan, double infpfn){ 
    apan = inapan;
    fpfn = infpfn;
    apan_fpfn_defined = true;
  };
  
  bool get_apan_fpfn_defined(){return apan_fpfn_defined;};
  double get_apan(){return apan;};
  double get_fpfn(){return fpfn;};
  

 void DrawFormFactor_SI(){};

 void DrawFormFactor_SD();
 
 private:
  
string FormSIcode,FormSDcode;

double (*pointeurSurFormFactorSI)(double ,double,FormFactorDataBase*); 
double (*pointeurSurFormFactorSD)(double ,double,FormFactorDataBase*); 

  double apan, fpfn;
  bool apan_fpfn_defined;
  
 private:


////////////////////////////////////////////////////////
// F a c t e u r s   d e   f o r m e  S D  1 2 9 X e
////////////////////////////////////////////////////////
static double Form129XeSD_p(double q, double A, FormFactorDataBase* forn);
static double Form129XeSD_n(double q, double A, FormFactorDataBase* form);
static double Form129XeSD(double q, double A,FormFactorDataBase* forn);

////////////////////////////////////////////////////////
// F a c t e u r s   d e   f o r m e  S D  1 2 9 X e
////////////////////////////////////////////////////////
static double Form131XeSD_p(double q, double A, FormFactorDataBase* forn);
static double Form131XeSD_n(double q, double A, FormFactorDataBase* form);
static double Form131XeSD(double q, double A,FormFactorDataBase* forn);

////////////////////////////////////////////////////////
// F a c t e u r s   d e   f o r m e  S I  G e n e r i c
////////////////////////////////////////////////////////

static double Helm(double q, double A,FormFactorDataBase* forn);
static double SolidSphere(double q, double A,FormFactorDataBase* forn);
 static double ThinShell(double q, double A,FormFactorDataBase* forn);
 
////////////////////////////////////////////////////////
// F a c t e u r s   d e   f o r m e  S D  G e n e r i c
////////////////////////////////////////////////////////

 
static double LewinSmithSD(double q, double A,FormFactorDataBase* forn);

};


#endif
