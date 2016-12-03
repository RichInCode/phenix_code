#ifndef analyzePhotons_h
#define analyzePhotons_h

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TFile.h"

class plotPhotons;

class analyzePhotons
{
  friend class plotPhotons;

 public:
  analyzePhotons(int dopp, int alarmLevel);
  virtual ~analyzePhotons(){}

  //=-=-=-=-=-==-=-==-=-=-===-
  // pi0 analysis functions
  int runPi0Analysis(TH3D *somehistFG, TH3D *somehistBG, int normScheme);
  //  int readInFileAndHistos(char *inname);
  int sliceAndDice(TH3D *somehistFG, TH3D *somehistBG);
  int findNormalization_pol2();      // chosen with normScheme = 1
  int findNormalization_pol0();      // chosen with normScheme = 0
  int findNormalization_pol1();      // chosen with normScheme = 5
  int findNormalization_fit();       // chosen wiht normScheme = 4
  int makeRatioHistos();
  int applyNormalization();
  int subtractBG();
  int calculateYield();
  int calculateYield_gaus();

  //=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-
  // inclusive photon analysis functions
  int runInclAnalysis(TH3D *somehistFG);
  int sliceAndDice(TH3D *simehistFG);
  int calculateInclYield();

  //=-=-=-=-=-=-=-=-=-=-=-=
  // calculate R_gamma
  int loadInCorrection(TGraphErrors *someplot);
  int loadInCorrection(TH1D *hho);
  int loadInCorrection(double a, double b, double c);  // use fermi function fit
  int calculateRgammaNum();
  int calculateRgamma(TH1D *somehist);
  
  //=-=-=-=-=-=-=-=-=-=-=-
  // direct photon yield
  int calculateDirectYield(TH1D *somehist, int numberOfEvents);
  int applyAbsoluteNorm(double binwidth, double ptbin, int nevents);

  void setIntegrationLimits(double low, double high)
  {
    low_intLimit = low;
    high_intLimit = high;
  };

 private:
  TH3D *h_FG_mass_vs_pt_vs_pt[5];
  TH3D *h_BG_mass_vs_pt_vs_pt[5];

  TH1D *massFG12[13];   // #gamma e^{+}e^{-} FG
  TH1D *massBG12[13];   // #gamma e^{+}e^{-} mixed BG
  TH1D *normRatio[13];  // FG/BG
  TH1D *normRatio_part[13];  // FG/BG excluding pi0 region
  TH1D *mass_sub[13];   // FG - nBG
  TH1D *massBG12_norm[13]; // n*BG
  TH1D *norm[13];           // n

  TH1D *hist_mass[13];

  TF1 *BG_function[13];
  TF1 *BG_poly[13];
  TF1 *BG_gaus[13];

  char* infilename;

  // some run flag variables
  int cent;
  int normProcedure;
  
  int pp;
  int verbosity;

  double yield[11];
  double yieldErr[11];
  double yield_gaus[11];
  double yieldErr_gaus[11];

  double yieldIncl[11];
  double yieldIncl_err[11];

  double yieldDirect[11];
  double yieldDirect_err[11];

  double rgamma[11];
  double rgamma_err[11];

  double hadronYield_sim[11];
  double simRatio[11];
  double correction[11];
  double correction_err[11];
  double numerator[11];
  double num_err[11];

  double ptindex[11];

  double low_intLimit, high_intLimit;

  TFile *out;

  ClassDef(analyzePhotons, 1)
};


#endif /* analyzePhotons_h */
