
//ROOT includes
#include "TCanvas.h"

//C++ stdlib
#include <iostream>
#include <string>
#include <cmath>

#include "deadMap.h"
#include "emcalDeadMap.h"

using namespace std;

ClassImp(emcalDeadMap);

//emcalDeadMap::emcalDeadMap(const string &name, int alarmlevel = 0):filename(name.c_str()),
//						       verbosity(alarmlevel)

emcalDeadMap::emcalDeadMap(int alarmlevel = 0): verbosity(alarmlevel),
						threshold(0.)
{
  
  // LoadTrueSectorMap();
  //InitializeWarnSectorMap();

}


int emcalDeadMap::LoadTrueSectorMap(TH2F *histo, int sect)
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


void emcalDeadMap::InitializeWarnSectorMap()
{

  for(int isect=0; isect<8; isect++)
    {
      char outputhistoname[256];
      sprintf(outputhistoname,"warnSectorMap_%d",isect);
      warnSectorMap[isect] = (TH2F*)trueSectorMap[isect]->Clone(outputhistoname);
      warnSectorMap[isect]->Reset();
  
      int nbinz = warnSectorMap[isect]->GetNbinsX();
      int nbiny = warnSectorMap[isect]->GetNbinsY();
      
      for(int iz=1; iz<=nbinz; iz++)
	{
	  for(int iy=1; iy<=nbiny; iy++)
	    {
	      // assume everything is alive
	      warnSectorMap[isect]->SetBinContent(iz, iy, 2);
	    }
	}
    }
}

void emcalDeadMap::FindExpectedPerSectorHits()
{

  //  cout << "expected hits per sector = ";

  for(int isect=0; isect<8; isect++)
    {
      int nbinz = cutTrueSectorMap[isect]->GetNbinsX();
      int nbiny = cutTrueSectorMap[isect]->GetNbinsY();

      /*
      sectorMean[isect] = cutTrueSectorMap[isect]->Integral()/(nbinz*nbiny);
      //sectorRMS[isect] = sqrt(trueSectorMap[isect]->GetEntries())/(nbinz*nbiny);
      sectorRMS[isect] = sqrt(sectorMean[isect]);
      */

      double bincounter = 0.;


      // calculate mean
      for(int iz=1; iz<=nbinz; iz++)
	{
	  for(int iy=1; iy<=nbiny; iy++)
	    {
	      double value = cutTrueSectorMap[isect]->GetBinContent(iz, iy);
	      //cout << value << endl;
	      
	      if(value > 0)
		{
		  sectorMean[isect] += value;
		  bincounter++;
		}
	    }
	}

      //  cout << "pre-mean = " << sectorMean[isect] << endl;
      //cout << "bincounter = " << bincounter << endl;

      sectorMean[isect] = sectorMean[isect]/bincounter;

      // calculate sigma
      for(int iz=1; iz<=nbinz; iz++)
	{
	  for(int iy=1; iy<=nbiny; iy++)
	    {
	      double value = cutTrueSectorMap[isect]->GetBinContent(iz, iy);
	      if(value > 0)
		{
		  sectorRMS[isect] += pow(value-sectorMean[isect],2);
		}
	    }
	}
      
      sectorRMS[isect] = sqrt(sectorRMS[isect]/(double)bincounter);

      cout << sectorMean[isect] << " +/- " << sectorRMS[isect] << endl;
    }
}

void emcalDeadMap::ScanSector(int sect)
{
  int nbinz = trueSectorMap[sect]->GetNbinsX();
  int nbiny = trueSectorMap[sect]->GetNbinsY();

  for(int iz=1; iz<=nbinz; iz++)
    {
      for(int iy=1; iy<=nbiny; iy++)
	{
	  float value = trueSectorMap[sect]->GetBinContent(iz, iy);
	  float diff = value - sectorMean[sect];
	  
	  // dead
	  if(diff <= -threshold*sectorRMS[sect])
	    {
	      warnSectorMap[sect]->SetBinContent(iz, iy, -1);
	      
	      // dead 3x3 neighbors
	      for(int i=-1; i<=1; i++)
		{
		  for(int j=-1; j<=1; j++)
		    {
		      int iztemp = iz+i;
		      int iytemp = iy+j;

		      if(!(j==0 && i==0))
			warnSectorMap[sect]->SetBinContent(iztemp, iytemp, -2);
		    }
		}
	    }
	  // hot
	  if(diff >= threshold*sectorRMS[sect])
	    {
	      warnSectorMap[sect]->SetBinContent(iz, iy, -3);
	      
	      // hot 3x3 neighbors
	      for(int i=-1; i<=1; i++)
		{
		  for(int j=-1; j<=1; j++)
		    {
		      int iztemp = iz+i;
		      int iytemp = iy+j;

		      if(!(j==0 && i==0))
			warnSectorMap[sect]->SetBinContent(iztemp, iytemp, -4); 
		    }
		}
	    }  // else
	}  // for (iy=0; iy<nbinsy; iy++)
    }  // for (iz=0; iz<nbinsy; iz++)

  //return 1;
}

void emcalDeadMap::MakeNewDataMap(int sect)
{
  int nbinz = cutTrueSectorMap[sect]->GetNbinsX();
  int nbiny = cutTrueSectorMap[sect]->GetNbinsY();

  //cutTrueSectorMap[sect]->Reset();

  for(int iz=1; iz<=nbinz; iz++)
    {
      for(int iy=1; iy<=nbiny; iy++)
	{
	  float warnlevel = warnSectorMap[sect]->GetBinContent(iz, iy);
	  //if( !(warnlevel==-1 || warnlevel==-3))
	  if(warnlevel == 2)
	    {
	      float value = trueSectorMap[sect]->GetBinContent(iz, iy);
	      cutTrueSectorMap[sect]->SetBinContent(iz, iy, value);
	    }
	  else
	    {
	      cutTrueSectorMap[sect]->SetBinContent(iz, iy, 0.);
	    }
	}
    }
}

void emcalDeadMap::PlotTrueSectorMap(TCanvas *canvas)
{

  for(int isect=0; isect<8; isect++)
    {
      canvas->cd(isect+1);
      trueSectorMap[isect]->Draw("colz");
    }

}

void emcalDeadMap::PlotCutTrueSectorMap(TCanvas *canvas)
{

  for(int isect=0; isect<8; isect++)
    {
      canvas->cd(isect+1);
      cutTrueSectorMap[isect]->Draw("colz");
    }

}

void emcalDeadMap::PlotWarnSectorMap(TCanvas *canvas)
{
  for(int isect=0; isect<8; isect++)
    {
      canvas->cd(isect+1);
      warnSectorMap[isect]->SetContour(6);
      warnSectorMap[isect]->Draw("colz");
    }
}
