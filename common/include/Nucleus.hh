#ifndef Nucleus_hh
#define Nucleus_hh 1

#include <iostream>
#include <string>

#include "FormFactorDataBase.hh"

using namespace std;

class Nucleus{

public:

  Nucleus(string Name,
	  double Muma,
	  double Amol,
	  double Z,
	  double J,
	  double Sp,
	  double Sn,
	  FormFactorDataBase* FormFactor){

    fName = Name;
    fZ = Z;
    fA = Amol;
    fMuma = Muma;
    fMGeV =  (fMuma * 931.494061/1000.) - Z * 511e-6; 
    fSp = Sp;
    fSn = Sn;
    fJ = J;
    fFormFactor = FormFactor;
  }
  

  ~Nucleus(){};

  double GetZ(){return fZ;}
  double GetA(){return fA;}
  double GetMuma(){return fMuma;}
  double GetMGeV(){return fMGeV;}
  double GetSp(){return fSp;}
  double GetSn(){return fSn;}
  double GetJ(){return fJ;}
  string GetName(){return fName;}

  FormFactorDataBase* GetFormFactor(){return fFormFactor;}

private:

  string fName;
  double fZ, fA, fMuma, fMGeV, fSp, fSn, fJ;

  FormFactorDataBase* fFormFactor;


};



#endif
