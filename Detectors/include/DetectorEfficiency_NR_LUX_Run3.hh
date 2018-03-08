#ifndef DetectorEfficiency_NR_LUX_Run3_hh
#define DetectorEfficiency_NR_LUX_Run3_hh 1

#include "DetectorEfficiency.hh"

#include <vector>

#include "TMath.h"
#include "Math/Interpolator.h"

using namespace std;


class DetectorEfficiency_NR_LUX_Run3: public DetectorEfficiency{

public:
  DetectorEfficiency_NR_LUX_Run3(){

    fName = "NR_LUX_Run3";

    
    fDataFileName = (string)getenv("TOOLS_PATH") + "/data_base/DetectorEfficiency/LUX_NR_Run3.data";


    ifstream file(fDataFileName.c_str());

    if(!file.is_open()){
      cout<<"<Error::DetectorEfficiency_NR_LUX_Run3> Error while opening the file: "<<fDataFileName<<endl;
      exit(1);
    }

    string line;
    getline(file, line);

    double maximum = -1;
    double Energy, Eff;
    while(file >> Energy >> Eff){
      fEnergy.push_back(Energy);
      fEff.push_back(Eff);
      maximum = max(maximum, Eff);
    }

    for(size_t i = 0; i< fEff.size(); ++i) fEff[i] /= maximum;
    
    fInterp = new ROOT::Math::Interpolator(fEnergy, fEff);
    
    
  }
  
  virtual ~DetectorEfficiency_NR_LUX_Run3(){}

  virtual double GetEfficiency(double E){
    if(E < fEnergy[0] || E > fEnergy.back()) return 0;
    return fInterp->Eval(E);
  }


private:
  virtual void PrepareGraph(){

    fGrEfficiency = new TGraph();
    for(size_t i =0; i<fEnergy.size(); ++i) fGrEfficiency->SetPoint( i, fEnergy[i], fEff[i]);
   
  }


private:

  string fDataFileName;
  
  vector<double> fEnergy, fEff;

  ROOT::Math::Interpolator* fInterp;

};

#endif
