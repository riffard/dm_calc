#include "WimpRate.hh"





WimpRate::WimpRate(string CalcName, Target* target, HaloModel* halo, DetectorEfficiency* efficiency){

  fCalcName = CalcName;
  fTarget = target;
  fHalo = halo;
  fEfficiency = efficiency;
  cout<<"<Info::WimpRate> Creation of a new WimpRate: "<< CalcName <<endl;

}



WimpRate::~WimpRate(){



  


}


void WimpRate::GetRate(TH1D*h, double mChi, double sigma0Si){

  h->Reset();
  for(int i=1; i<= h->GetNbinsX(); ++i){

    double Er = h->GetXaxis()->GetBinCenter(i);
    double rate = Rate_SI(Er, mChi, sigma0Si);
    h->SetBinContent(i, rate);
  }
  
}


double WimpRate::Rate_SI(double Er, double mChi, double sigma0Si){


  double rate = 0;

  for(size_t i =0 ; i<fTarget->size(); ++i){
  
    double fraction;
    Nucleus* nucleus = fTarget->GetComponant(i,fraction);
    
    double A = nucleus->GetA();
    double Z = nucleus->GetZ();
    double m = nucleus->GetMGeV()*1e6;
    double mChi_eV = mChi*1e6;
    
    rate += fraction*WimpRate::Rate_SI(Er, mChi_eV,
					      m, A, Z,
					      sigma0Si);
    
    
    
  }
  
  if(fEfficiency) rate *= fEfficiency->GetEfficiency(Er);
  

  return rate;
}




double WimpRate::Rate_SI(double E_r, double m_wimp,
				double m_nuc, double A_nuc, double Z_nuc,
				double xsec_wn){
  
  
  
  // Define constants  
  double rho0 = 2.29e-18; // keV^4 = 0.3 GeV/cm^3
  double v0 = 7.33e-4; // 1 = 220 km/s
  double vlab = 7.73e-4; // 1 = 232 km/s
  double vesc = 1.81e-3; // 1 = 544 km/s
  double rn0 = 5.79e-6; // 1/keV = 1.14 fm
  double s = 4.57e-6; // 1/keV = 0.9 fm
  double m_p = 938e3; // keV = 938 MeV
  
  // Calculate nuclear form factor
  double q = sqrt(2 * m_nuc * E_r);
  double rn = rn0 * pow(A_nuc,1./3.);
  double bessel = sin(q*rn) / pow(q*rn,2) - cos(q*rn) / (q*rn);
  double F = 3 * bessel * exp(-pow(q*s,2) / 2) / (q*rn);

  // Calculate phase space factor
  double Nesc = erf(vesc / v0) - (2/sqrt(TMath::Pi())) * (vesc/v0) * exp(-pow(vesc,2) / pow(v0,2));
  double r = 4 * m_wimp * m_nuc / pow(m_wimp + m_nuc,2);
  double E0 = m_wimp * pow(v0, 2) / 2;
  double vmin = v0 * sqrt(E_r / (E0*r));
  double term = (sqrt(TMath::Pi())/4) * (v0/vlab) * (erf((vmin + vlab)/v0) - erf((vmin-vlab)/v0)) - exp(-pow(vesc,2) / pow(v0,2));
  double phase = (1/Nesc) * (2/sqrt(TMath::Pi())) * (1/v0) * term;

  // Calculate nucleus-WIMP xsec
  double m_r_nuc = m_wimp * m_nuc / (m_wimp + m_nuc);
  double m_r_p = m_wimp * m_p / (m_wimp + m_p);
  double  xsec = pow(A_nuc,2) * pow(m_r_nuc / m_r_p,2) * xsec_wn;


  // Calculate differential rate
  double dRdEr = (rho0 * xsec / (2 * m_wimp * pow(m_r_nuc,2))) * pow(F,2) * phase;
  dRdEr = dRdEr * 6.94e76; // Evts/tonne/yr/keV

  if (dRdEr < 0 )dRdEr = 0;
  
  return dRdEr;

}



