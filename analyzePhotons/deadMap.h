#ifndef deadMap_h
#define deadMap_h

#include "TH2F.h"

//class TH2F

class edgePoint
{
 public:
  edgePoint() {}
 edgePoint(int idnum, int znum, int ynum): id(idnum), iz(znum), iy(ynum){}
  virtual ~edgePoint() {}

  //copy constructor
  edgePoint(const edgePoint &point) { edgePointPointer = new edgePoint;  *edgePointPointer = *point.edgePointPointer; }

  const int Get_iz() {return iz;}
  const int Get_iy() {return iy;}
  const int Get_id() {return id;}

  void Set_iz(int znum) {iz = znum;}
  void Set_iy(int ynum) {iy = ynum;}
  void Set_id(int idnum) {id = idnum;}

 private:
  int id;
  int iz;
  int iy;

  edgePoint *edgePointPointer;

};


class edgeLine
{
 public:
  edgeLine() {}
 edgeLine(int idvalue, double slopevalue, double interceptvalue):id(idvalue),slope(slopevalue),intercept(interceptvalue) {}
  virtual ~edgeLine() {}

  //copy constructor
  //edgeLine(const edgeLine &line) { edgeLinePointer = new edgeLine;  *edgeLinePointer = *line.edgeLinePointer; }

  const double Get_slope() {return slope;}
  const double Get_intercept() {return intercept;}

  void Set_slope(int value) {slope = value;}
  void Set_intercept(int value) {intercept = value;}


 private:
  int id;
  double slope;
  double intercept;

  edgeLine *edgeLinePointer;
};

class deadMap
{
 public:
  deadMap() {}
  virtual ~deadMap() {}

  int LoadTrueSectorMap(TH2F *histo, int sect, char *name);
  virtual void ScanSector() {}

 protected:
  //edgePoint *edgep;
  TH2F *trueSectorMap[8];
  TH2F *cutTrueSectorMap[8];

};


#endif /*deadMap_h*/
