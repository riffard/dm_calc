#include "FormFactorDataBase.hh"

#include "TMath.h"

FormFactorDataBase::FormFactorDataBase(string inFormSIcode, string inFormSDcode): FormSIcode(inFormSIcode), FormSDcode(inFormSDcode){

  pointeurSurFormFactorSD = NULL;
  pointeurSurFormFactorSI = NULL;
  apan_fpfn_defined = false;

  cout<<"<FORM FACTOR DATA BASE> SD form factor reference: "<<inFormSDcode<<endl;


  ////////////////////////////////////////////////////////
  // F a c t e u r s   d e   f o r m e  S I
  ////////////////////////////////////////////////////////

  if(inFormSDcode.find("129Xe") != string::npos){

    if(inFormSDcode.find("_p") != string::npos){

      cout<<"<FORM FACTOR DATA BASE> Proton diffusion model form factor choosen for "<<inFormSDcode<<endl;

      pointeurSurFormFactorSD = Form129XeSD_p;

    }else if(inFormSDcode.find("_n") != string::npos){

      cout<<"<FORM FACTOR DATA BASE> Neutron diffusion model form factor choosen for "<<inFormSDcode<<endl;

      pointeurSurFormFactorSD = Form129XeSD_n;

    }else {
      cout<<"<FORM FACTOR DATA BASE> General model form factor choosen for "<<inFormSDcode <<endl;

      pointeurSurFormFactorSD = Form129XeSD;
    }


  }else if(inFormSDcode.find("131Xe") != string::npos){
    if(inFormSDcode.find("_p") != string::npos){

      cout<<"<FORM FACTOR DATA BASE> Proton diffusion model form factor choosen for "<<inFormSDcode<<endl;

      pointeurSurFormFactorSD = Form131XeSD_p;

    }else if(inFormSDcode.find("_n") != string::npos){

      cout<<"<FORM FACTOR DATA BASE> Neutron diffusion model form factor choosen for "<<inFormSDcode<<endl;

      pointeurSurFormFactorSD = Form131XeSD_n;

    }else {
      cout<<"<FORM FACTOR DATA BASE> General model form factor choosen for "<<inFormSDcode <<endl;

      pointeurSurFormFactorSD = Form131XeSD;
    }

  } else if(inFormSDcode.find("Helm") != string::npos){
    pointeurSurFormFactorSD = Helm;

  } else if(inFormSDcode.find("SolidSphere") != string::npos){
    pointeurSurFormFactorSD = SolidSphere;

  } else if(inFormSDcode.find("LewinSmith") != string::npos){
    pointeurSurFormFactorSD = LewinSmithSD;

  }else if(inFormSDcode.find("No") != string::npos){
    cout<<"<FORM FACTOR DATA BASE><WARNING> No SD form factor for the current nucleus"<<endl;
  }else{
    cout<<"<FORM FACTOR DATA BASE><WARNING> No SD form factor found with the ref "<<inFormSDcode<<" !"<<endl;
  }




  ////////////////////////////////////////////////////////
  // F a c t e u r s   d e   f o r m e  S I
  ////////////////////////////////////////////////////////

  cout<<"<FORM FACTOR DATA BASE> SI form factor reference: "<<inFormSIcode<<endl;

  if(inFormSIcode.find("Helm") != string::npos){
    pointeurSurFormFactorSI = Helm;

  }else if(inFormSIcode.find("SolidSphere") != string::npos){
    pointeurSurFormFactorSI = SolidSphere;

  }else if(inFormSIcode.find("ThinShell") != string::npos){
    pointeurSurFormFactorSI = ThinShell;

  }else{

    cout<<"<FORM FACTOR DATA BASE><WARNING> No SI form factor found with the ref "<<inFormSIcode<<" !"<<endl;
  }



}

FormFactorDataBase::~FormFactorDataBase(){

}




////////////////////////////////////////////////
// F a c t e u r s   d e   f o r m e  1 2 9 X e
////////////////////////////////////////////////
double FormFactorDataBase::Form129XeSD_p(double q, double A, FormFactorDataBase* form){

  double b = 2.2853;

  double u = pow(q*b,2)/2.;

  double S00 = exp(-u)*(0.054731*pow(u,0)-0.146897*pow(u,1) + 0.182479*pow(u,2)-0.128112*pow(u,3) + 0.0539978*pow(u,4)-0.0133335*pow(u,5) + 0.00190579*pow(u,6)-0.000148373*pow(u,7) + 5.11732e-06*pow(u,8)-2.06597e-08*pow(u,9));

  //(1b)
  //double S01 = exp(-u)*(-0.102732*pow(u,0) + 0.297105*pow(u,1)-0.387513*pow(u,2) + 0.281816*pow(u,3)-0.122388*pow(u,4) + 0.0317668*pow(u,5)-0.00492337*pow(u,6) + 0.000439836*pow(u,7)-2.02852e-05*pow(u,8) + 3.46755e-07*pow(u,9));

  //(1b+2b)
  double S01 = exp(-u)*(-0.0796645*pow(u,0) + 0.231997*pow(u,1)-0.304198*pow(u,2) + 0.222024*pow(u,3)-0.096693*pow(u,4) + 0.0251835*pow(u,5)-0.00392356*pow(u,6) + 0.000353343*pow(u,7)-1.65058e-05*pow(u,8) + 2.88576e-07*pow(u,9));

  //(1b)
  //double S11 = exp(-u)*(0.048192*pow(u,0)-0.148361*pow(u,1) + 0.202347*pow(u,2)-0.151853*pow(u,3) + 0.0674284*pow(u,4)-0.0179342*pow(u,5) + 0.00286368*pow(u,6)-0.000265795*pow(u,7) + 1.29656e-05*pow(u,8)-2.47418e-07*pow(u,9));

  //(1b+2b)
  double S11 = exp(-u)*(0.02933*pow(u,0)-0.0905396*pow(u,1) + 0.122783*pow(u,2)-0.0912046*pow(u,3) + 0.0401076*pow(u,4)-0.010598*pow(u,5) + 0.00168737*pow(u,6)-0.000156768*pow(u,7) + 7.69202e-06*pow(u,8)-1.48874e-07*pow(u,9));

  double f = S00 + S01 + S11;
  //(1b)
  f /= 0.054731 - 0.102732 + 0.048192;
  //(1b+2b)
  //f /= 0.054731 - 0.0796645 + 0.02933;

  return sqrt(f);
}
double FormFactorDataBase::Form129XeSD_n(double q, double A, FormFactorDataBase* form){

  double b = 2.2853;

  double u = pow(q*b,2)/2.;

  double S00 = exp(-u)*(0.054731*pow(u,0)-0.146897*pow(u,1) + 0.182479*pow(u,2)-0.128112*pow(u,3) + 0.0539978*pow(u,4)-0.0133335*pow(u,5) + 0.00190579*pow(u,6)-0.000148373*pow(u,7) + 5.11732e-06*pow(u,8)-2.06597e-08*pow(u,9));

  //(1b)
  //double S01 = exp(-u)*(-0.102732*pow(u,0) + 0.297105*pow(u,1)-0.387513*pow(u,2) + 0.281816*pow(u,3)-0.122388*pow(u,4) + 0.0317668*pow(u,5)-0.00492337*pow(u,6) + 0.000439836*pow(u,7)-2.02852e-05*pow(u,8) + 3.46755e-07*pow(u,9));

  //(1b+2b)
  double S01 = exp(-u)*(-0.0796645*pow(u,0) + 0.231997*pow(u,1)-0.304198*pow(u,2) + 0.222024*pow(u,3)-0.096693*pow(u,4) + 0.0251835*pow(u,5)-0.00392356*pow(u,6) + 0.000353343*pow(u,7)-1.65058e-05*pow(u,8) + 2.88576e-07*pow(u,9));

  //(1b)
  //double S11 = exp(-u)*(0.048192*pow(u,0)-0.148361*pow(u,1) + 0.202347*pow(u,2)-0.151853*pow(u,3) + 0.0674284*pow(u,4)-0.0179342*pow(u,5) + 0.00286368*pow(u,6)-0.000265795*pow(u,7) + 1.29656e-05*pow(u,8)-2.47418e-07*pow(u,9));

  //(1b+2b)
  double S11 = exp(-u)*(0.02933*pow(u,0)-0.0905396*pow(u,1) + 0.122783*pow(u,2)-0.0912046*pow(u,3) + 0.0401076*pow(u,4)-0.010598*pow(u,5) + 0.00168737*pow(u,6)-0.000156768*pow(u,7) + 7.69202e-06*pow(u,8)-1.48874e-07*pow(u,9));



  double f = S00 - S01 + S11;
  //(1b)
  f /= 0.054731 + 0.102732 + 0.048192;
  //(1b+2b)
  //f /= 0.054731+0.0796645+0.02933;

  return sqrt(f);
}
double FormFactorDataBase::Form129XeSD(double q, double A, FormFactorDataBase* form){

  bool apan_fpfn_defined = form->get_apan_fpfn_defined();


  if(apan_fpfn_defined == false ){
    cout<<"<WARNING> ap/an and fp/fn data not defined. Please use: Set_apan_fpfn(double apan, double fpfn)"<<endl;
    return 0;
  }
  double apan = form->get_apan();

  double b = 2.2853;

  double u = pow(q*b,2)/2.;

  double S00 = exp(-u)*(0.054731*pow(u,0)-0.146897*pow(u,1) + 0.182479*pow(u,2)-0.128112*pow(u,3) + 0.0539978*pow(u,4)-0.0133335*pow(u,5) + 0.00190579*pow(u,6)-0.000148373*pow(u,7) + 5.11732e-06*pow(u,8)-2.06597e-08*pow(u,9));

  //(1b)
  double S01 = exp(-u)*(-0.102732*pow(u,0) + 0.297105*pow(u,1)-0.387513*pow(u,2) + 0.281816*pow(u,3)-0.122388*pow(u,4) + 0.0317668*pow(u,5)-0.00492337*pow(u,6) + 0.000439836*pow(u,7)-2.02852e-05*pow(u,8) + 3.46755e-07*pow(u,9));

  //(1b+2b)
  //double S01 = exp(-u)*(-0.0796645*pow(u,0) + 0.231997*pow(u,1)-0.304198*pow(u,2) + 0.222024*pow(u,3)-0.096693*pow(u,4) + 0.0251835*pow(u,5)-0.00392356*pow(u,6) + 0.000353343*pow(u,7)-1.65058e-05*pow(u,8) + 2.88576e-07*pow(u,9));

  //(1b)
  double S11 = exp(-u)*(0.048192*pow(u,0)-0.148361*pow(u,1) + 0.202347*pow(u,2)-0.151853*pow(u,3) + 0.0674284*pow(u,4)-0.0179342*pow(u,5) + 0.00286368*pow(u,6)-0.000265795*pow(u,7) + 1.29656e-05*pow(u,8)-2.47418e-07*pow(u,9));

  //(1b+2b)
  //double S11 = exp(-u)*(0.02933*pow(u,0)-0.0905396*pow(u,1) + 0.122783*pow(u,2)-0.0912046*pow(u,3) + 0.0401076*pow(u,4)-0.010598*pow(u,5) + 0.00168737*pow(u,6)-0.000156768*pow(u,7) + 7.69202e-06*pow(u,8)-1.48874e-07*pow(u,9));


  double f = pow((apan+1)/apan,2)*S00 + ((apan+1)/apan)*((apan-1)/apan)*S01 + pow((apan-1)/apan,2)*S11;
  //(1b)
  f /= pow((apan+1)/apan,2)*0.054731  -((apan+1)/apan)*((apan-1)/apan)*0.102732 + pow((apan-1)/apan,2)*0.048192;
  //(1b+2b)
  //f /=  pow((apan+1)/apan,2)*0.054731-((apan+1)/apan)*((apan-1)/apan)*0.0796645 + pow((apan-1)/apan,2)*0.02933;

  return f;
}


////////////////////////////////////////////////
// F a c t e u r s   d e   f o r m e  1 3 1 X e
////////////////////////////////////////////////
double FormFactorDataBase::Form131XeSD_p(double q, double A, FormFactorDataBase* form){

  double b = 2.2853;

  double u = pow(q*b,2)/2.;

  double S00 = exp(-u)*(0.0417889*pow(u,0)-0.111171*pow(u,1) + 0.171966*pow(u,2)-0.133219*pow(u,3) + 0.0633805*pow(u,4)-0.0178388*pow(u,5) + 0.00282476*pow(u,6)-0.000231681*pow(u,7) + 7.78223e-06*pow(u,8)-4.49287e-10*pow(u,9));

  //(1b)
  //double S01 = exp(-u)*(-0.0784478*pow(u,0) + 0.230484*pow(u,1)-0.343106*pow(u,2) + 0.263525*pow(u,3)-0.120946*pow(u,4) + 0.0331754*pow(u,5)-0.00528724*pow(u,6) + 0.00046475*pow(u,7)-2.00407e-05*pow(u,8) + 2.90375e-07*pow(u,9));

  //(1b+2b)
  double S01 = exp(-u)*(-0.0608808*pow(u,0) + 0.181473*pow(u,1)-0.272533*pow(u,2) + 0.211776*pow(u,3)-0.0985956*pow(u,4) + 0.027438*pow(u,5)-0.0044424*pow(u,6) + 0.000397619*pow(u,7)-1.74758e-05*pow(u,8) + 2.55979e-07*pow(u,9));

  //(1b)
  //double S11 = exp(-u)*(0.0368132*pow(u,0)-0.118361*pow(u,1) + 0.176773*pow(u,2)-0.137987*pow(u,3) + 0.063821*pow(u,4)-0.0176743*pow(u,5) + 0.00287653*pow(u,6)-0.000263605*pow(u,7) + 1.23239e-05*pow(u,8)-2.17839e-07*pow(u,9));

  //(1b+2b)
  double S11 = exp(-u)*(0.022446*pow(u,0)-0.0733931*pow(u,1) + 0.110509*pow(u,2)-0.0868752*pow(u,3) + 0.0405399*pow(u,4)-0.0113544*pow(u,5) + 0.00187572*pow(u,6)-0.000175285*pow(u,7) + 8.40043e-06*pow(u,8)-1.53632e-07*pow(u,9));

  double f = S00 + S01 + S11;

  //(1b)
  f /= 0.0417889 - 0.0784478 + 0.0368132;
  //(1b+2b)
  //f /= 0.0417889  -0.0608808 + 0.022446;

  return sqrt(f);
}
double FormFactorDataBase::Form131XeSD_n(double q, double A, FormFactorDataBase* form){

  double b = 2.2853;

  double u = pow(q*b,2)/2.;

  double S00 = exp(-u)*(0.0417889*pow(u,0)-0.111171*pow(u,1) + 0.171966*pow(u,2)-0.133219*pow(u,3) + 0.0633805*pow(u,4)-0.0178388*pow(u,5) + 0.00282476*pow(u,6)-0.000231681*pow(u,7) + 7.78223e-06*pow(u,8)-4.49287e-10*pow(u,9));

  //(1b)
  //double S01 = exp(-u)*(-0.0784478*pow(u,0) + 0.230484*pow(u,1)-0.343106*pow(u,2) + 0.263525*pow(u,3)-0.120946*pow(u,4) + 0.0331754*pow(u,5)-0.00528724*pow(u,6) + 0.00046475*pow(u,7)-2.00407e-05*pow(u,8) + 2.90375e-07*pow(u,9));

  //(1b+2b)
  double S01 = exp(-u)*(-0.0608808*pow(u,0) + 0.181473*pow(u,1)-0.272533*pow(u,2) + 0.211776*pow(u,3)-0.0985956*pow(u,4) + 0.027438*pow(u,5)-0.0044424*pow(u,6) + 0.000397619*pow(u,7)-1.74758e-05*pow(u,8) + 2.55979e-07*pow(u,9));

  //(1b)
  //double S11 = exp(-u)*(0.0368132*pow(u,0)-0.118361*pow(u,1) + 0.176773*pow(u,2)-0.137987*pow(u,3) + 0.063821*pow(u,4)-0.0176743*pow(u,5) + 0.00287653*pow(u,6)-0.000263605*pow(u,7) + 1.23239e-05*pow(u,8)-2.17839e-07*pow(u,9));

  //(1b+2b)
  double S11 = exp(-u)*(0.022446*pow(u,0)-0.0733931*pow(u,1) + 0.110509*pow(u,2)-0.0868752*pow(u,3) + 0.0405399*pow(u,4)-0.0113544*pow(u,5) + 0.00187572*pow(u,6)-0.000175285*pow(u,7) + 8.40043e-06*pow(u,8)-1.53632e-07*pow(u,9));


  double f = S00 - S01 + S11;
  //(1b)
  f /= 0.0417889 + 0.0784478 + 0.0368132;
  //(1b+2b)
  //f /= 0.0417889  -0.0608808 + 0.022446;

  return sqrt(f);
}
double FormFactorDataBase::Form131XeSD(double q, double A, FormFactorDataBase* form){

  bool apan_fpfn_defined = form->get_apan_fpfn_defined();

  if(apan_fpfn_defined == false ){
    cout<<"<WARNING> ap/an and fp/fn data not defined. Please use: Set_apan_fpfn(double apan, double fpfn)"<<endl;
    return 0;
  }
  double apan = form->get_apan();

  double b = 2.2853;

  double u = pow(q*b,2)/2.;

  double S00 = exp(-u)*(0.0417889*pow(u,0)-0.111171*pow(u,1) + 0.171966*pow(u,2)-0.133219*pow(u,3) + 0.0633805*pow(u,4)-0.0178388*pow(u,5) + 0.00282476*pow(u,6)-0.000231681*pow(u,7) + 7.78223e-06*pow(u,8)-4.49287e-10*pow(u,9));

  //(1b)
  //double S01 = exp(-u)*(-0.0784478*pow(u,0) + 0.230484*pow(u,1)-0.343106*pow(u,2) + 0.263525*pow(u,3)-0.120946*pow(u,4) + 0.0331754*pow(u,5)-0.00528724*pow(u,6) + 0.00046475*pow(u,7)-2.00407e-05*pow(u,8) + 2.90375e-07*pow(u,9));

  //(1b+2b)
  double S01 = exp(-u)*(-0.0608808*pow(u,0) + 0.181473*pow(u,1)-0.272533*pow(u,2) + 0.211776*pow(u,3)-0.0985956*pow(u,4) + 0.027438*pow(u,5)-0.0044424*pow(u,6) + 0.000397619*pow(u,7)-1.74758e-05*pow(u,8) + 2.55979e-07*pow(u,9));

  //(1b)
  //double S11 = exp(-u)*(0.0368132*pow(u,0)-0.118361*pow(u,1) + 0.176773*pow(u,2)-0.137987*pow(u,3) + 0.063821*pow(u,4)-0.0176743*pow(u,5) + 0.00287653*pow(u,6)-0.000263605*pow(u,7) + 1.23239e-05*pow(u,8)-2.17839e-07*pow(u,9));

  //(1b+2b)
  double S11 = exp(-u)*(0.022446*pow(u,0)-0.0733931*pow(u,1) + 0.110509*pow(u,2)-0.0868752*pow(u,3) + 0.0405399*pow(u,4)-0.0113544*pow(u,5) + 0.00187572*pow(u,6)-0.000175285*pow(u,7) + 8.40043e-06*pow(u,8)-1.53632e-07*pow(u,9));


  double f = pow((apan+1)/apan,2)*S00 + ((apan+1)/apan)*((apan-1)/apan)*S01 + pow((apan-1)/apan,2)*S11;
  //(1b)
  f /= pow((apan+1)/apan,2)*0.0417889 - ((apan+1)/apan)*(apan-1)/apan*0.0784478 + pow((apan-1)/apan,2)*0.0368132;
  //(1b+2b)
  //f /= pow((apan+1)/apan,2)*0.0417889  -((apan+1)/apan)*(apan-1)/apan*0.0608808 + pow((apan-1)/apan,2)*0.022446;


  return sqrt(f);
}



////////////////////////////////////////////////////////
// F a c t e u r s   d e   f o r m e  S I  G e n e r i c
////////////////////////////////////////////////////////

double FormFactorDataBase::Helm(double q, double A,FormFactorDataBase* forn){

  if(q==0)return 1;

  double rn;
  double s = 0.9; // fm

  if(A > 100) {

    double a = 0.52; // fm
    double c = 1.23*pow(A,1./3.) - 0.6; // fm
    double rn2 = pow(c,2) + (7./3.)*pow(TMath::Pi()*a,2) - 5*pow(s,2);
    rn = sqrt(rn2);

  } else {
    rn = 1.25 * pow(A, 1./3.);
  }
  double F = 3;
  F *= exp( - pow( q*s , 2)/2 );
  F *= sin( q * rn ) - q * rn * cos( q * rn );
  F /= pow( q*rn  ,3.);

  return F;

}


double FormFactorDataBase::SolidSphere(double q, double A,FormFactorDataBase* forn){

  if(q==0)return 1;

  double rn;

  double s = 0.9; // fm

  if(A > 100) {

    double a = 0.52; // fm
    double c = 1.23*pow(A,1./3.) - 0.6; // fm

    double rn2 = pow(c,2) + (7./3.)*pow(TMath::Pi()*a,2) - 5*pow(s,2);

    rn = sqrt(rn2);

  } else {
    rn = 1.25 * pow(A, 1./3.);

  }

  double F = 3;

  F *= sin( q * rn ) - q * rn * cos( q * rn );

  F /= pow( q*rn  ,3.);

  return F;


}

double FormFactorDataBase::ThinShell(double q, double A,FormFactorDataBase* forn){

  if(q==0) return 1;

  double rn;

  double s = 0.9; // fm

  if(A > 100) {

    double a = 0.52; // fm
    double c = 1.23*pow(A,1./3.) - 0.6; // fm

    double rn2 = pow(c,2) + (7./3.)*pow(TMath::Pi()*a,2) - 5*pow(s,2);

    rn = sqrt(rn2);

  } else {
    rn = 1.25 * pow(A, 1./3.);

  }

  double F = sin( q * rn ) / ( q * rn );

  return F;


}


double FormFactorDataBase::LewinSmithSD(double q, double A,FormFactorDataBase* forn){

  if(q==0)return 1;

  double rn = 1.0 * pow(A,1./3.);

  double F = sqrt(0.047);

  if( q*rn < 2.55 || q*rn > 4.5)F =  sin( q*rn ) / (q*rn) ;


  return F;
}




///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
