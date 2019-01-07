#include "Target.hh"



Target::Target(string targetName){
  fTargetName = targetName;

  fNuclei.clear();
  fFraction.clear();
}


Target::~Target(){

  fFraction.clear();
  for(size_t i =0; i< fNuclei.size(); ++i) delete fNuclei[i];
  fNuclei.clear();

}



void Target::AddNucleus(double fraction, Nucleus* nucleus){

  fFraction.push_back(fraction);
  fNuclei.push_back(nucleus);
  
}


void Target::CheckNormalization(){

  double sum =0;
  for(size_t i=0; i<fFraction.size(); i++)sum += fFraction[i];

  if(sum != 1){
    cout<<"<Warning::Target> Issue in the target "<< fTargetName  <<" abundance:"<<endl;
    cout<<"<Warning::Target> sum != 1 (sum = "<< sum<<")"<<endl;
    cout<<"<Warning::Target> Re-normalization !"<<endl;

    for(size_t i=0; i<fFraction.size(); i++)   fFraction[i] /= sum;

  }
  
}
