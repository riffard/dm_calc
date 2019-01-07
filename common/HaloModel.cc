#include "HaloModel.hh"





HaloModel::HaloModel(double v_earth, double sigma_v, double v_esc, string model):
  fv_earth(v_earth), fsigma_v(sigma_v),  fv_esc(v_esc)
{

  pointeurSurModelHalo = NULL;

  cout<<"<HaloModelIntegration::Info> Halo model: "<< model<<endl; 
  
  if(model == "SMH"){

    
    pointeurSurModelHalo = getHaloIntegrationValue_SMH;
    
  }else if(model == "SMHcutoff"){
    
    pointeurSurModelHalo = getHaloIntegrationValue_SMHcutoff;
    
  }else{

    cout<<"<HaloModelIntegration::Warning> The chosen halo model: "<< model<< " does not exist !"<<endl; 

    
    pointeurSurModelHalo = getHaloIntegrationValue_SMH;
    
  }
  
  
  
  
}


HaloModel::~HaloModel(){
  
  
}




double HaloModel::getHaloIntegrationValue(double v_min){


  if(pointeurSurModelHalo == NULL){
    cout<<"<HaloModelIntegration::Error> No halo model choosen!"<<endl;
    exit(1);
  }

  return (*pointeurSurModelHalo)(v_min,fv_earth, fsigma_v, fv_esc);
  
}


double HaloModel::getHaloIntegrationValue_SMH(double v_min, double v_earth, double sigma_v, double v_esc){


  //if(v_min > fv_esc) return 0;

  double v_0 = sigma_v*sqrt(2./3.);
  
  double erfmax = erf( (v_min + v_earth) /  v_0  );

  double erfmin = erf( (v_min - v_earth) /  v_0  );
 
    
  double val = (erfmax - erfmin) / (4. * v_earth);
      
  val -= exp( - pow( v_esc / v_0 ,2) ) / ( sqrt( TMath::Pi() ) * v_0 );

  double Norme =  erf(  v_esc  / v_0  );

  Norme -= (2 * v_esc * exp( - pow( v_esc / v_0 ,2)  )  ) / ( sqrt( 2*TMath::Pi() )*sigma_v ) ;
  
  
  val /= Norme;

  val *= 2;

  if(val < 0) return 0;
  
  return val;
}




double HaloModel::getHaloIntegrationValue_SMHcutoff(double v_min, double v_earth, double sigma_v, double v_esc){

  
  double v0bar = sqrt(2./3.)*sigma_v;
  
  double x = v_min / v0bar;

  double y = v_earth / v0bar;

  double z = v_esc / v0bar;

  //double vobar_y = v0bar;
  double vobar_y = v0bar*y;
  
  double N_esc = erf(z) - 2 * z * exp( - pow( z ,2) )/pow(TMath::Pi() , 0.5 );

  if( z < y && x > abs( y - z ) ){
    
    double val = 1/(vobar_y);

    return val;
    
  }else if( z >  y && x < abs(y - z)){

    double val = 1 /( 2 * N_esc * vobar_y ); 

    val *= ( erf( x + y ) - erf( x - y ) - ( 4 / sqrt( TMath::Pi() ) )*y*exp( - pow( z , 2) )  );
    
    return val;
    
  }else if( abs(y - z) < x && x < (y + z)){

    double val = 1 /( 2 * N_esc * vobar_y ); 

    val *= ( erf( z ) - erf( x - y ) - ( 2 / sqrt( TMath::Pi() ) ) * ( y + z - x ) * exp( - pow( z , 2) )  );
   
    return val;
    
  }else if( (y + z) < x){

    return 0;

  }else{

    return -10;
  }
  
  
}
