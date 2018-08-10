#include "NeutrinoFlux.hh"
#include "style.hh"


#include<fstream>

#include "TAxis.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TLegend.h"

//-------------------------------------------------------
//-------------------------------------------------------
NeutrinoFlux::NeutrinoFlux(){

  
  load_style();

  data_base_path = "/data_base/neutrino_flux";
  if(getenv("DM_CALC_PATH")) data_base_path = (string)getenv("DM_CALC_PATH") + "/data_base/neutrino_flux";
  
  
  string NeutrinoListFile = data_base_path + "/NeutrinoList.txt";
  string ScaleConstantsFile = data_base_path + "/ScaleConstants.txt";

  string line;
  ifstream NeutrinoList(NeutrinoListFile.c_str());
  if(!NeutrinoList.is_open()){
    cout<<"<Error::NeutrinoFlux> Impossible to open "<< NeutrinoListFile<<endl;
    exit(-1);
  }
  while( NeutrinoList >> line) fFluxFiles.push_back(line);

  double val;
  ifstream ScaleConstants(ScaleConstantsFile.c_str());
  while( ScaleConstants >> val) fScale_constant.push_back(val);


  if(fFluxFiles.size() != fScale_constant.size()){
    cout<<"<Error::NeutrinoFlux> Error on files "<<endl;
    exit(100);
  }
  
  for(size_t i=0; i< fFluxFiles.size(); ++i){

    string filename = data_base_path + "/data/" + fFluxFiles[i] + ".txt";

    ifstream file(filename);

    if(!file.is_open()){ 
      cout<<"<Error::NeutrinoFlux> Error while opening: "<< filename <<endl;
      exit(2);
    }

    Flux* flux = new Flux(fScale_constant[i]);
    
    double x,y;
    while(file >> x >> y) flux->Push_back(x,y);
   
    flux->Complete();
    
    fFlux_collection[fFluxFiles[i]] = flux;
  }
 
 
}

//-------------------------------------------------------
//-------------------------------------------------------



//-------------------------------------------------------
//-------------------------------------------------------
NeutrinoFlux::~NeutrinoFlux(){
  
}
//-------------------------------------------------------
//-------------------------------------------------------

void NeutrinoFlux::PrintAvailableData(){

  cout<<endl;
  cout<<"<Information::NeutrinoFlux> These are all the available neutrino flux:"<<endl;
  
  for(auto it : fFlux_collection){
    cout<<"<Information::NeutrinoFlux>    - "<< it.first<<"                : " << it.second->GetTotalFlux()<< "  neutrinos/cm^2/s" <<flush;
    if( it.second->GetMode() == "line") cout<< " (line)"<<flush;
    cout<<endl;
  }
}

//-------------------------------------------------------
//-------------------------------------------------------
void NeutrinoFlux::DrawFluxCollection(){

  TCanvas* c = new TCanvas;
  c->SetLogx();
  c->SetLogy();
  
  
  TMultiGraph* mg = new TMultiGraph();
  TLegend* leg = new TLegend(0.6,0.6,0.95,0.95);


  int color = 10;

  for(auto it : fFlux_collection){
    it.second->GetGraph()->SetLineColor(color);
    it.second->GetGraph()->SetMarkerColor(color);
    it.second->GetGraph()->SetLineWidth(2);
    
    mg->Add(it.second->GetGraph(), "pl");
    leg->AddEntry(it.second->GetGraph(), it.first.c_str(), "l");
    color ++;
  }

  leg->SetNColumns(6);
  
  mg->Draw("a");
  mg->GetXaxis()->SetLimits(2e-2, 1e3);
  mg->GetYaxis()->SetLimits(1e-7, 1e15);
  mg->GetXaxis()->SetRangeUser(2e-2, 1e3);
  mg->GetYaxis()->SetRangeUser(1e-7, 1e15);

  mg->GetXaxis()->SetTitle("Neutrino energy [MeV]");
  mg->GetYaxis()->SetTitle("Flux [neutrinos/cm^2/s/MeV]");
  
  leg->Draw("same");

  c->Update();
  c->Modified();  
  
}
//-------------------------------------------------------
//-------------------------------------------------------
