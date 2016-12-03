#include "deadMap.h"
#include "dchDeadMap.h"

#include <iostream>
#include <map>

/*
void dchDeadMap::BreakupMapIntoSectors()
{
  float armEdge = 1.5;
  float chargeEdge = 0.;

  int binx = trueSectorMap[0]->GetXaxis()->FindBin(armEdge);
  int biny = trueSectorMap[0]->GetYaxis()->FindBin(chargeEdge);

  // west arm, positive charge sector
  for(int i=1; i<=binx; i++)
    {
      for(int j=1; j<=biny; j++)
	{
	  trueSectorMap
	}
    }
}
*/

void dchDeadMap::ScanSector(int sect)
{
  
}

/*
int dchDeadMap::LoadTrueSectorMap(TH2F *histo, int sect, char *)
{

  char inputhistoname[256];
  char outputhistoname[256];
 
  sprintf(inputhistoname,"hist_modMap_%d",sect);
  sprintf(outputhistoname,"trueSectorMap_%d",sect);
  trueSectorMap[sect] = (TH2F*)histo->Clone(outputhistoname);

  sprintf(outputhistoname,"cutTrueSectorMap_%d",sect);
  cutTrueSectorMap[sect] = (TH2F*)histo->Clone(outputhistoname);
  //  cutTrueSectorMap[sect]->Reset();
 

  return 1;
}
*/
