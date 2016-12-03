#ifndef dchDeadMap_h
#define dchDeadMap_h

#include <map>

//class deadMap;


//typedef dedgePoint edgePoint;

class dchDeadMap: public deadMap
{
 public:
  dchDeadMap() {};
  virtual ~dchDeadMap() {};

  void InitializeWarnSectorMap();
  void ScanSector(int sect);
  void MakeNewDataMap(int sect);
  void FindExpectedPerSectorHits();  


 private:
  
  TH2D *hist_dcMap;

  std::map<int,edgePoint> edgePointHolder;
  std::vector<edgeLine> edgeLineHolder;



};

#endif /*dchDeadMap_h*/
