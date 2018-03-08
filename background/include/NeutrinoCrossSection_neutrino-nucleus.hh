#ifndef NeutrinoCrossSection_hh
#define NeutrinoCrossSection_hh 1

class NeutrinoCrossSection{


public:

  NeutrinoCrossSection();
  virtual NeutrinoCrossSection();

  
  virtual double GetCrossSection(double E_recoil, double E_nrutrino){ return 0;};
  

};


#endif
