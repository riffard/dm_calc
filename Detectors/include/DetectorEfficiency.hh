#ifndef DetectorEfficiency_hh
#define DetectorEfficiency_hh 1


#include <iostream>
#include <fstream>

//#include "style.hh"

#include "TAxis.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "Math/Interpolator.h"

using namespace std;


class DetectorEfficiency{

public:
  DetectorEfficiency(){

    //load_style();
    
    fName = "";
    fGrEfficiency = NULL;
  }
  
  virtual ~DetectorEfficiency(){}

  virtual double GetEfficiency(double E){return -1;}

  void Draw(){

    PrepareGraph();

    if(!fGrEfficiency) return;

    fGrEfficiency->SetTitle(fName.c_str());
    fGrEfficiency->GetXaxis()->SetTitle("Energy [keV]");
    fGrEfficiency->GetYaxis()->SetTitle("Efficiency");
    fGrEfficiency->SetLineColor(10);
    fGrEfficiency->SetLineWidth(3);
    
    TCanvas* c = new TCanvas;
    fGrEfficiency->Draw("apl");
    c->Update();
    c->Modified();
    
    
  }

 
  string GetName(){return fName;}

  
private:
  virtual void PrepareGraph(){  }
  
protected:

  string fName;
  TGraph* fGrEfficiency;

  

};



#endif
