#include "deadMap.h"



int deadMap::LoadTrueSectorMap(TH2F *histo, int sect, char *name)
{

  char inputhistoname[256];
  char outputhistoname[256];
 
  //sprintf(inputhistoname,"hist_modMap_%d",sect);
  sprintf(inputhistoname,name,sect);
  sprintf(outputhistoname,"trueSectorMap_%d",sect);
  trueSectorMap[sect] = (TH2F*)histo->Clone(outputhistoname);

  sprintf(outputhistoname,"cutTrueSectorMap_%d",sect);
  cutTrueSectorMap[sect] = (TH2F*)histo->Clone(outputhistoname);
  //  cutTrueSectorMap[sect]->Reset();
 

  return 1;
}

//ClassImp(deadMap);
