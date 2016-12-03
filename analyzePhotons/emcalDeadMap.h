#ifndef emcalDeadMap_h
#define emcalDeadMap_h

// ROOT includes
#include "TH2F.h"

class emcalDeadMap:deadMap
{
 
 public:
  emcalDeadMap() {}   // default constructor
  //emcalDeadMap(const std::string &name, int alarmLevel);   // constructor
  emcalDeadMap(int alarmLevel);   // constructor
  virtual ~emcalDeadMap(){}     // destructor


  // access functions

  virtual int LoadTrueSectorMap(TH2F *histo, int sect); 
  void InitializeWarnSectorMap();
  void ScanSector(int sect);
  void MakeNewDataMap(int sect);
  void FindExpectedPerSectorHits(); 

  void PlotTrueSectorMap(TCanvas *canvas);
  void PlotCutTrueSectorMap(TCanvas *canvas);
  void PlotWarnSectorMap(TCanvas *canvas);

  void SetThreshold(int thresh) { threshold = thresh; }

  const TH2F *GetWarnHistogram(int sect) { return warnSectorMap[sect]; }

 private:


  // encodes the dead map
  //  - a value of 2 indicates a good tower
  //  - a value of -1 indicates a dead tower
  //  - a value of -2 indicates a dead tower neighbor tower
  //  - a value of -3 indicates a hot tower
  //  - a value of -4 indicates a hot neighbor tower
  //  - a value of -5 indicates a sector edge tower
  TH2F *warnSectorMap[8];

  //std::string filename;
  const char *filename;
  
  int verbosity;

  float threshold;

  double sectorMean[8];
  double sectorRMS[8];

  //TH2F *trueSectorMap[8];
  //TH2F *cutTrueSectorMap[8];

  ClassDef(emcalDeadMap, 1)
};

#endif /*emcalDeadMap_h*/
