#include "NeutrinoRate.hh"

#include <boost/algorithm/string.hpp>

NeutrinoRate::NeutrinoRate(Target* target, NeutrinoFlux* neutrino_flux, string required_fluxes, NeutrinoCrossSection* cross_section){
  
  fTarget = target;
  fNeutrinoFluxDB = neutrino_flux;
  fRequiredFluxes_str = required_fluxes;

  if(required_fluxes.find("All")!= string::npos || required_fluxes.find("all")!= string::npos){
    fRequiredFluxes = neutrino_flux->GetFluxCollectionName();
  } else{
    boost::split(fRequiredFluxes, required_fluxes, boost::is_any_of(",: ;"));  
  }
  fCrossSection = cross_section;

}



NeutrinoRate::~NeutrinoRate(){
  
  
}






void NeutrinoRate::GetRate(TH1D* hE, THStack* hs, TLegend* leg){

  double m_e = 511;


  for(int iflux = 0; iflux < fRequiredFluxes.size(); iflux++){
    
    Flux* flux = fNeutrinoFluxDB->GetFlux(fRequiredFluxes[iflux]);
    
    TH1D* hE_contribution = (TH1D*)hE->Clone(Form("%s_%s", hE->GetName(), fRequiredFluxes[iflux].c_str()));

    hE_contribution->SetLineColor(10+iflux);
    hE_contribution->SetLineWidth(2);
    
    if(hs)hs->Add(hE_contribution);
    if(leg) leg->AddEntry(hE_contribution, Form("%s",fRequiredFluxes[iflux].c_str()), "l");


    
    for(int iE = 1; iE <= hE_contribution->GetNbinsX(); ++iE){
      
      double E_recoil = hE_contribution->GetXaxis()->GetBinCenter(iE);
      double rate = 0;
      
      double fraction = 0;
      for(int incl = 0; incl <fTarget->size(); ++incl){

	Nucleus* nucleus = fTarget->GetComponant(incl, fraction);

	double m_nuc = nucleus->GetMGeV()*1e6;
        double E_nu_min = 0.0;
        if (fCrossSection->GetCrossSectionName() == "coherent_NR")
	    E_nu_min = sqrt( m_nuc * E_recoil / 2); // keV
        else if (fCrossSection->GetCrossSectionName() == "electroweak_ER")
            E_nu_min = 0.5 * (E_recoil + sqrt(E_recoil * (E_recoil + 2 * m_e)));
        else {
            cerr << "Cross-section type not found" << endl;
            exit(EXIT_FAILURE);
        }
	
	double dRdEr  = 0;
	
	if(flux->GetMode() == "line"){
	//----------------------------------------------------
	// Rate calculation for a line
	//----------------------------------------------------	

	  double E_nu, flux_nu;
	  flux->GetLine(E_nu, flux_nu);
	  double E_nu_keV = E_nu*1e3;
	  
	  if (E_nu_keV < E_nu_min){
	    dRdEr = 0;
	  }else{
	      
	    double integral = flux_nu * fCrossSection->GetCrossSection(E_recoil, E_nu_keV, nucleus);
	    dRdEr = (5.61e35 / m_nuc) * integral * 3.16e7;  
	  }
	  
	}else{
	  //----------------------------------------------------
	  // Rate for a continuous spectrum
	  //----------------------------------------------------

	  vector<double> *energies = flux->GetEnergy();
	  vector<double> *amplitude = flux->GetFlux();

	  double integral = 0;
	  for(size_t i=0; i<amplitude->size() -1 ; ++i){

	    if(energies->at(i)*1e3 < E_nu_min )continue;
	    
	    double E_nu = 0.5 * (energies->at(i) + energies->at(i+1));
	    double E_nu_keV = E_nu * 1000;
	    
	    double flux_nu = 0.5 * (amplitude->at(i) + amplitude->at(i+1));
	    
	    double dE_nu = energies->at(i+1) - energies->at(i);
	    //double dE_nu_keV = dE_nu * 1e3;
	    
	    double dxsecdEr = fCrossSection->GetCrossSection(E_recoil, E_nu_keV, nucleus);
	    integral += (flux_nu * dxsecdEr * dE_nu);

	  }

          // Currently, dRdEr is in evts/keV/atom/sec
	  dRdEr = integral * (5.61e35 / m_nuc); // evts/keV/tonne/sec
	  dRdEr = dRdEr * 3.1558e7; // evts/keV/tonne/yr
	  
	  if (dRdEr < 0) dRdEr = 0;

	}
	
	rate += fraction * dRdEr;
       	
      }

      
      hE_contribution->SetBinContent(iE, rate);


	
    }
  

    hE->Add(hE_contribution);
    
    
  }
  

  
  /*
  for(int iE = 0; iE < hE->GetNbinsX(); ++iE){
    double E = hE->GetXaxis()->GetBinCenter(iE);

    //cout<<E<<endl;

    
    
    }*/
  
  
  
}
