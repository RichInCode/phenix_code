#ifndef plotPhotons_h
#define plotPhotons_h

#include "TGraphErrors.h"
#include "TCanvas.h"

class analyzePhotons;
class TFile;

class plotPhotons
{

 public:

  plotPhotons(analyzePhotons *photo);
  virtual ~plotPhotons(){}

  void drawFGtoBG();
  void drawFG();
  void drawSub();
  void drawPi0Yield();

  void drawInclYield();
  
  void drawRawInclToPi();

  void drawRgamma();
  void drawRgamma_withSys();
  
  void drawDirectYield();
  void drawDirectYield_withSys();

  void drawCorrection();

  void drawCocktailRatio();

  void drawPPG088();
  void drawPPG162();
  void drawPPfit();

  void SaveRgamma(TFile *file);
  void SaveCorrectionUsed(TFile *file);
  void SaveRgammaNumerator(TFile *file);
  void SaveCocktailRatio(TFile *file);
  void SaveDirectYield(TFile *file);
  void SaveInclYield(TFile *file);
  void SavePi0Yield(TFile *file);
  void SaveAll(TFile *file);
  
 private:

  TGraphErrors *p_pi0yield;
  TGraphErrors *p_inclYield;
  TGraphErrors *p_inclToPiYield;
  TGraphErrors *p_directYield;
  TGraphErrors *p_directYieldSys;
  TGraphErrors *p_rgamma;
  TGraphErrors *p_rgammaSys;
  TGraph *p_correction;
  TGraph *p_cocktailratio;

  TCanvas *c_FGtoBG;
  TCanvas *c_FG;
  TCanvas *c_BG;
  TCanvas *c_sub;
  TCanvas *c_pi0yield;
  TCanvas *c_inclyield;
  TCanvas *c_inclToPiYield;
  TCanvas *c_directyield;
  TCanvas *c_rgamma;
  TCanvas *c_correction;
  TCanvas *c_cocktailratio;

  analyzePhotons *myphotons;

  double x[11];
  double xbin[11];

  int pp;

  ClassDef(plotPhotons, 1)

};


#endif /* plotPhotons_h */
