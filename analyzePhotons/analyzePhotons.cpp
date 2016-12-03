// ROOT includes
#include "Riostream.h"
#include "TVirtualFitter.h"
#include "TF1.h"
#include "TMath.h"

//PHENIX includes

// std includes
#include <iostream>
#include <math.h>

#include "analyzePhotons.h"

ClassImp(analyzePhotons);

analyzePhotons::analyzePhotons(int dopp, int alarmLevel):pp(dopp), verbosity(alarmLevel), low_intLimit(0.1), high_intLimit(0.18)
{
  if(pp)
    {
      ptindex[0] = 0.6;
      ptindex[1] = 1.0;
      ptindex[2] = 1.6;
      ptindex[3] = 3.5;
    }
  else
    {
      ptindex[0] = 0.5;
      ptindex[1] = 0.7;
      ptindex[2] = 0.9;
      ptindex[3] = 1.1;
      ptindex[4] = 1.3;
      ptindex[5] = 1.5;
      ptindex[6] = 1.7;
      ptindex[7] = 1.9;
      ptindex[8] = 2.25;
      ptindex[9] = 3.;
      ptindex[10] = 4.25;
    }

  char name[256];
  
  for(int i=0; i<13; i++)
    {
      sprintf(name, "BG_gaus_%d", i);
      BG_gaus[i] = new TF1(name, "gaus");
      sprintf(name, "BG_poly_%d", i);
      BG_poly[i] = new TF1(name, "pol3");
      sprintf(name, "BG_function_%d", i);
      BG_function[i] = new TF1(name, "gaus(0) + pol3(3)"); 
    }
}

int analyzePhotons::runPi0Analysis(TH3D *somehistFG, TH3D *somehistBG, int normScheme = 1)
{

  normProcedure = normScheme;

  // if( !readInFileAndHistos(infilename) )
  //   {
  //     std::cout << "There is an error in reading file, stop and check" << std::endl;

  //     return 0;
  //   }

  sliceAndDice(somehistFG, somehistBG);
   
  makeRatioHistos();
   
  if(normProcedure == 1)
    findNormalization_pol2();
  else if(normProcedure == 0)
    findNormalization_pol0();
  else if(normProcedure == 5)
    findNormalization_pol1();
  else if(normProcedure == 4)
    findNormalization_fit();
  else
    {
      std::cout << "you picked an unknown normScheme..stop and check" << std::endl;
      return 0;
    }
    

  applyNormalization();

  subtractBG();

  calculateYield();

  calculateYield_gaus();

  //writePi0ResultsToFile();

  return 1;
}

int analyzePhotons::runInclAnalysis(TH3D *someHistFG)
{

  sliceAndDice(someHistFG);

  calculateInclYield();

  return 1;
}

/*
int analyzePhotons::readInFileAndHistos(char* inname)
{
  //TFile *inputfile = new TFile(inname);
  std::cout << inname << std::endl;
  TFile inputfile(inname);
  std::cout << "yep" << std::endl;
  inputfile.ls();

  // if( !inputfile )
  //   return 0;
  
  // FG
  h_FG_mass_vs_pt_vs_pt[0] = (TH3D*)inputfile.Get("ATM_mass_vs_ATM_pT_vs_ph_pT_FG12_0")->Clone("h_FG_mass_vs_pt_vs_pt_0");
  h_FG_mass_vs_pt_vs_pt[1] = (TH3D*)inputfile.Get("ATM_mass_vs_ATM_pT_vs_ph_pT_FG12_1")->Clone("h_FG_mass_vs_pt_vs_pt_1");
  h_FG_mass_vs_pt_vs_pt[2] = (TH3D*)inputfile.Get("ATM_mass_vs_ATM_pT_vs_ph_pT_FG12_2")->Clone("h_FG_mass_vs_pt_vs_pt_2");
  h_FG_mass_vs_pt_vs_pt[3] = (TH3D*)inputfile.Get("ATM_mass_vs_ATM_pT_vs_ph_pT_FG12_3")->Clone("h_FG_mass_vs_pt_vs_pt_3");
  h_FG_mass_vs_pt_vs_pt[4] = (TH3D*)inputfile.Get("ATM_mass_vs_ATM_pt_vs_pt_FG_MB")->Clone("h_FG_mass_vs_pt_vs_pt_4");

  // BG
  h_BG_mass_vs_pt_vs_pt[0] = (TH3D*)inputfile.Get("ATM_mass_vs_ATM_pT_vs_ph_pT_BG12_0")->Clone("h_BG_mass_vs_pt_vs_pt_0");
  h_BG_mass_vs_pt_vs_pt[1] = (TH3D*)inputfile.Get("ATM_mass_vs_ATM_pT_vs_ph_pT_BG12_1")->Clone("h_BG_mass_vs_pt_vs_pt_1");
  h_BG_mass_vs_pt_vs_pt[2] = (TH3D*)inputfile.Get("ATM_mass_vs_ATM_pT_vs_ph_pT_BG12_2")->Clone("h_BG_mass_vs_pt_vs_pt_2");
  h_BG_mass_vs_pt_vs_pt[3] = (TH3D*)inputfile.Get("ATM_mass_vs_ATM_pT_vs_ph_pT_BG12_3")->Clone("h_BG_mass_vs_pt_vs_pt_3");
  h_BG_mass_vs_pt_vs_pt[4] = (TH3D*)inputfile.Get("ATM_mass_vs_ATM_pt_vs_pt_BG_MB")->Clone("h_BG_mass_vs_pt_vs_pt_4");

  return 1;
}
*/
int analyzePhotons::sliceAndDice(TH3D *somehistFG, TH3D *somehistBG)
{
 
  char name_FG[256];
  char name_BG[256];

  double pT[13];
  double pTBin[13];

  int begin = 0;
  int end = 0;

  for(int ipt=0; ipt<13; ipt++)
    {
      sprintf(name_FG, "massFG12_%d", ipt);
      sprintf(name_BG, "massBG12_%d", ipt);
      
      switch (pp)
	{

	case 1:
	  
	  if(ipt < 2)  	  // 200 MeV pT bins  (these are empty histogram bins)
	    {
	      begin = ipt*2+1;
	      end = ipt*2+2;
	      
	      pT[ipt] = ipt*0.2+0.1;
	      pTBin[ipt] = 0.1;
	    }

	  else if(ipt < 4)  // 400 MeV pT bins
	    {
	      begin = end+1;
	      end = begin+3;

	      pT[ipt] = ipt*0.4+0.2;
	      pTBin[ipt] = 0.2;
	    }
	  else if(ipt < 5)  // 800 MeV pT bins
	    {
	      begin = end+1;
	      end = begin+7;
	      
	      pT[ipt] = 1.6;
	      pTBin[ipt] = 0.4;
	    }
	  else if(ipt < 6)  // all the rest to 5GeV
	    {
	      begin = end+1;
	      end = begin+29;

	      pT[ipt] = 3.5;
	      pTBin[ipt] = 1.5;
	    }
	  
	  break;
	  
	default:
	  if(ipt < 10)    // 200 MeV pT bins
	    {
	      
	      begin = ipt*2+1;
	      end = ipt*2+2;
	      
	      pT[ipt] = ipt*0.2+0.1;
	      pTBin[ipt] = 0.1;
	    }
	  
	  else if(ipt<11)              // 500 MeV pT bins
	    {
	      
	      begin = end+1;
	      end = end+5;
	      
	      if(ipt == 10)
		{
		  pT[ipt] = 2.25;
		  pTBin[ipt] = 0.25;
		}
	      else
		{
		  pT[ipt] = pT[ipt-1]+0.5;
		  pTBin[ipt] = 0.25;
		}
	    }
	  
	  else if(ipt < 12)    // 1 GeV pT bin
	    {
	      begin = end+1;
	      end = end+10;
	      
	      pT[ipt] = 3.;
	      pTBin[ipt] = 0.50;
	    }
	  
	  else           // 1.5 GeV bin
	    {
	      
	      begin = end+1;
	      end = end+15;
	      
	      pT[ipt] = 4.25;
	      pTBin[ipt] = 0.75;
	      
	    }
	  break;
	}

      if(verbosity)
	{
	  cout << "pi0 bin start: " << begin << ", bin end: " << end << " for bin " << ipt << endl;
	  cout << "binwidth = " << somehistFG->GetYaxis()->GetBinWidth(1) << endl;;
	}

      if(pp && ipt>5)
	break; 

      char rangename[256];
      sprintf(rangename, "%1.2f < p^{ee}_{T} < %1.2f GeV/c", pT[ipt]-pTBin[ipt], pT[ipt]+pTBin[ipt]);

      somehistFG->GetYaxis()->SetRange(begin, end);
      //      somehistFG->GetZaxis()->SetRange(1,6);
      massFG12[ipt] = (TH1D*)somehistFG->Project3D("x")->Clone(name_FG);
      massFG12[ipt]->SetTitle(rangename);
      //h_FG_mass_vs_pt_vs_pt[centindex]->GetYaxis()->SetRange(begin, end);
      //massFG12[ipt] = (TH1D*)h_FG_mass_vs_pt_vs_pt[centindex]->Project3D("x")->Clone(name_FG);
      massFG12[ipt]->Sumw2();

      somehistBG->GetYaxis()->SetRange(begin, end);
      //somehistBG->GetZaxis()->SetRange(1,6);
      massBG12[ipt] = (TH1D*)somehistBG->Project3D("x")->Clone(name_BG);
      //h_BG_mass_vs_pt_vs_pt[centindex]->GetYaxis()->SetRange(begin, end);
      //massBG12[ipt] = (TH1D*)h_BG_mass_vs_pt_vs_pt[centindex]->Project3D("x")->Clone(name_BG);
      //massBG12[ipt]->Sumw2();


      // make nice bins
      massFG12[ipt]->Rebin(5);
      massBG12[ipt]->Rebin(5);
      /*
      
      massFG12[ipt]->Rebin(4); 
      massBG12[ipt]->Rebin(4); 
      
      if(ipt < 10)
	{
	  massFG12[ipt]->Rebin(2);
	  massBG12[ipt]->Rebin(2);
	}
      else
	{
	  massFG12[ipt]->Rebin(3);
	  massBG12[ipt]->Rebin(3);
	}
      */
    }
  
  return 1;
}


int analyzePhotons::sliceAndDice(TH3D *somehistFG)
{
 
  char name_FG[256];

  double pT[13];
  double pTBin[13];

  int begin = 0;
  int end = 0;

  somehistFG->GetZaxis()->SetRange(1,1);

  for(int ipt=0; ipt<13; ipt++)
    {
      sprintf(name_FG, "hist_mass_%d", ipt);
      

      switch (pp)
	{
	  
	case 1:
	  if(ipt < 2)  // 200 MeV pT bins (empty throw away bins)
	    {
	      begin = ipt*10+1;
	      end = begin+9;
	      
	      pT[ipt] = ipt*0.2+0.1;
	      pTBin[ipt] = 0.1;
	    }
	  else if(ipt < 4)  // 400 MeV pT bins
	    {
	      begin = end+1;
	      end = begin+19;

	      pT[ipt] = ipt*0.4 + 0.2;
	      pTBin[ipt] = 0.2;
	    }
	  else if(ipt < 5)  // 800 MeV pT bins
	    {
	      begin = end + 1;
	      end = begin + 39;
	      
	      pT[ipt] = 1.6;
	      pTBin[ipt] = 0.4;
	    }
	  else if(ipt < 6)
	    {
	      begin = end+1;
	      end = begin+149;
	      
	      pT[ipt] = 3.5;
	      pTBin[ipt] = 1.5;
	  }
	  
	  break;
	  
	default:
	  if(ipt < 10)    // 200 MeV pT bins
	    {
	      
	      begin = ipt*10+1;
	      end = begin+9;
	      
	      pT[ipt] = ipt*0.2+0.1;
	      pTBin[ipt] = 0.1;
	    }
	  
	  else if(ipt<11)              // 500 MeV pT bins
	    {
	      
	      begin = end+1;
	      end = end+25;
	      
	      if(ipt == 10)
		{
		  pT[ipt] = 2.25;
		  pTBin[ipt] = 0.25;
		}
	      else
		{
		  pT[ipt] = pT[ipt-1]+0.5;
		  pTBin[ipt] = 0.25;
		}
	    }
	  
	  else if(ipt < 12)    // 1 GeV pT bin
	    {
	      begin = end+1;
	      end = end+50;
	      
	      pT[ipt] = 3.;
	      pTBin[ipt] = 0.50;
	    }
	  
	  else           // 1.5 GeV bin
	    {
	      
	      begin = end+1;
	      end = end+75;
	      
	      pT[ipt] = 4.25;
	      pTBin[ipt] = 0.75;
	      
	    }
	  break;
	}
  
      if(verbosity)
	{
	  cout << "binwidth = " << somehistFG->GetYaxis()->GetBinWidth(1) << " for number of bins = " << somehistFG->GetNbinsY() << endl;
	  cout << "inclusive bin start: " << begin << ", bin end: " << end << endl;
	}
  
      if(pp && ipt>5)
	break;
      
      somehistFG->GetYaxis()->SetRange(begin, end);
      hist_mass[ipt] = (TH1D*)somehistFG->Project3D("x")->Clone(name_FG);
      
    }
  
  return 1;
}


int analyzePhotons::makeRatioHistos()
{
  
  char name_norm[256];

  double binwidth;
  int exclude1, exclude2;
  
  for(int ipt=0; ipt<13; ipt++)
    {

      if(pp && ipt>5)
	break;

      sprintf(name_norm, "norm_%d", ipt);
      norm[ipt] = (TH1D*)massFG12[ipt]->Clone(name_norm);
      norm[ipt]->Reset();

      sprintf(name_norm, "normRatio_%d", ipt);
      normRatio[ipt] = (TH1D*)massFG12[ipt]->Clone(name_norm);
      normRatio[ipt]->Divide(massBG12[ipt]);

      sprintf(name_norm, "normRatio_part_%d", ipt);
      normRatio_part[ipt] = (TH1D*)massFG12[ipt]->Clone(name_norm);
      normRatio_part[ipt]->Reset();

      sprintf(name_norm, "massBG12_norm_%d", ipt);
      massBG12_norm[ipt] = (TH1D*)massBG12[ipt]->Clone(name_norm);

      binwidth = normRatio[ipt]->GetBinWidth(1);

      // define exclusion region in pi0 peak
      exclude1 = floor( 0.1/binwidth )+1;
      exclude2 = floor( 0.2/binwidth );

      // fill in the histogram with empty bins in the pi0 region
      //  used for the normalization calcualtion
      for(int inbin=0; inbin<exclude1; inbin++)
	{
	  normRatio_part[ipt]->SetBinContent( inbin, normRatio[ipt]->GetBinContent(inbin) );
	  normRatio_part[ipt]->SetBinError( inbin, normRatio[ipt]->GetBinError(inbin) );
	}
      for(int inbin=exclude2+1; inbin<=normRatio_part[ipt]->GetNbinsX(); inbin++)
	{
	  normRatio_part[ipt]->SetBinContent( inbin, normRatio[ipt]->GetBinContent(inbin) );
	  normRatio_part[ipt]->SetBinError( inbin, normRatio[ipt]->GetBinError(inbin) );
	}

      
    }
  
  return 1;
}

int analyzePhotons::findNormalization_pol2()
{

  for(int ipt=2; ipt<13; ipt++)
    {

      if(pp && ipt>5)
	break;

      //TF1 *mypol2 = new TF1("mypol2","[0] + [1]*x + [2]*x*x",0.,0.5);

      //mypol2->SetParameter(2, 0.5);

      // something funny with this bin...
      if(pp && ipt==3)
	normRatio_part[ipt]->SetBinError(4,normRatio_part[ipt]->GetBinError(3));

      normRatio_part[ipt]->Fit("pol2","S","",0.,0.5);

      TVirtualFitter *fitter = TVirtualFitter::GetFitter();

      if( fitter )
	{
	  fitter->GetConfidenceIntervals(norm[ipt], 0.66); // fills norm histograms with a histogram of the fit result with errors
	}
      else
	return 0;
    }


  return 1;
}

int analyzePhotons::findNormalization_pol1()
{

  for(int ipt=2; ipt<13; ipt++)
    {
      if(pp && ipt>5)
	break;
      
      normRatio_part[ipt]->Fit("pol1","S","",0.01,0.5);

      TVirtualFitter *fitter = TVirtualFitter::GetFitter();

      if( fitter )
	{
	  fitter->GetConfidenceIntervals(norm[ipt], 0.66); // fills norm histograms with a histogram of the fit result with errors
	}
      else
	return 0;
    }

  return 1;
}

int analyzePhotons::findNormalization_pol0()
{

  for(int ipt=2; ipt<13; ipt++)
    {
      if(pp && ipt>5)
	break;
      
      normRatio_part[ipt]->Fit("pol0","S","",0.01,0.5);

      TVirtualFitter *fitter = TVirtualFitter::GetFitter();

      if( fitter )
	{
	  fitter->GetConfidenceIntervals(norm[ipt], 0.66); // fills norm histograms with a histogram of the fit result with errors
	}
      else
	return 0;
    }

  return 1;
}

int analyzePhotons::findNormalization_fit()
{
  for(int ipt=0; ipt<13; ipt++)
    {

      if(pp && ipt>5)
	break;

      BG_function[ipt]->SetParameter(0,1e9);
      BG_function[ipt]->SetParameter(1,0.14);
      BG_function[ipt]->SetParameter(2,0.02);

      //      massFG12[ipt]->Rebin(2);

      massFG12[ipt]->Fit(BG_function[ipt],"","");

      BG_gaus[ipt]->SetParameter(0,BG_function[ipt]->GetParameter(0));
      BG_gaus[ipt]->SetParameter(1,BG_function[ipt]->GetParameter(1));
      BG_gaus[ipt]->SetParameter(2,BG_function[ipt]->GetParameter(2));

      cout << "parameter = " << BG_gaus[ipt]->GetParameter(0) << endl;


      //      BG_poly[ipt]->SetParameters(BG_function[ipt]->GetParameter(3), BG_function[ipt]->GetParameter(4), BG_function[ipt]->GetParameter(5), BG_function[ipt]->GetParameter(6));
    }

  return 1;
}

int analyzePhotons::applyNormalization()
{

  for(int ipt=2; ipt<13; ipt++)
    {
      if(pp && ipt>5)
	break;

      massBG12_norm[ipt]->Multiply(norm[ipt]);
    }
  
  return 1;

}

int analyzePhotons::subtractBG()
{

  char hname[256];

  for(int ipt=2; ipt<13; ipt++)
    {

      if(pp && ipt>5)
	break;

      sprintf(hname, "mass_sub_%d", ipt);
      mass_sub[ipt] = (TH1D*)massFG12[ipt]->Clone(hname);
 
      //if(ipt < 5)
      //mass_sub[ipt]->Add(massBG12_norm[ipt], -1);

      // do the subtraction the long way to avoid empty bin issue
      //  creating negative yields
      mass_sub[ipt]->Reset();
      //      mass_sub[ipt]->Sumw2();
      for(int ibin=1; ibin<=mass_sub[ipt]->GetNbinsX(); ibin++)
	{
	  double valueFG = massFG12[ipt]->GetBinContent(ibin);
	  double valueErrFG = massFG12[ipt]->GetBinError(ibin);
	  double valueBG = massBG12_norm[ipt]->GetBinContent(ibin);
	  double valueErrBG = massBG12_norm[ipt]->GetBinError(ibin);
	  double value;
	  double valueErr;
	  if(valueFG > 0)
	  {
	    value = valueFG - valueBG;
	    valueErr = sqrt(valueErrFG*valueErrFG + valueErrBG*valueErrBG);

	    mass_sub[ipt]->SetBinContent(ibin, value);
	    mass_sub[ipt]->SetBinError(ibin, valueErr);
	  }
	}

      /*
      char rangename[256];
      sprintf(rangename, "%f < p^{ee}_{T} < %f", pT[ipt]-pTBin[ipt], pT[ipt]+pTBin[ipt]);
      mass_sub[ipt]->SetTitle(rangename);
      */
    }
  
  return 1;
}

int analyzePhotons::calculateInclYield()
{

  for(int ipt=2; ipt<13; ipt++)
    {

      if(pp && ipt>5)
	break;

      yieldIncl[ipt-2] = hist_mass[ipt]->Integral();
      yieldIncl_err[ipt-2] = sqrt(yieldIncl[ipt-2]);
           

      if(pp)
	{
	  if(ipt < 4)
	    {
	      yieldIncl[ipt-2] = yieldIncl[ipt-2]/0.4;
	      yieldIncl_err[ipt-2] = yieldIncl_err[ipt-2]/0.4;
	    }
	  else if(ipt < 5)
	    {
	      yieldIncl[ipt-2] = yieldIncl[ipt-2]/0.8;
	      yieldIncl_err[ipt-2] = yieldIncl_err[ipt-2]/0.8;
	    }
	  else
	    {
	      yieldIncl[ipt-2] = yieldIncl[ipt-2]/3.;
	      yieldIncl_err[ipt-2] = yieldIncl_err[ipt-2]/3.;
	    }
	}
      else
	{
	  if(ipt < 10)
	    {
	      yieldIncl[ipt-2] = yieldIncl[ipt-2]/0.2;
	      yieldIncl_err[ipt-2] = yieldIncl_err[ipt-2]/0.2;
	    }
	  else if(ipt < 11)
	    {
	      yieldIncl[ipt-2] = yieldIncl[ipt-2]/0.5;
	      yieldIncl_err[ipt-2] = yieldIncl_err[ipt-2]/0.5;
	    }
	  else if(ipt < 12)
	    {
	      yieldIncl[ipt-2] = yieldIncl[ipt-2]/1.;
	      yieldIncl_err[ipt-2] = yieldIncl_err[ipt-2]/1.;
	    }
	  else
	    {
	      yieldIncl[ipt-2] = yieldIncl[ipt-2]/1.5;
	      yieldIncl_err[ipt-2] = yieldIncl_err[ipt-2]/1.5;
	    }
	}
      
    }

  return 1;
}

int analyzePhotons::calculateYield()
{
  int integral_start, integral_end;

  double negativeYield[11];
 
  for(int ipt=2; ipt<13; ipt++)
    {

      if(pp && ipt>5)
	break;

      integral_start = floor( low_intLimit/mass_sub[ipt]->GetBinWidth(1) );
      integral_end = floor( high_intLimit/mass_sub[ipt]->GetBinWidth(1) );  
      
      cout << "normProcedure = " << normProcedure << endl;
      
      if(pp && ipt>5)
	break;

      if(normProcedure != 4)
	{
	  
	  
	  yield[ipt-2] = 0.;
	  yieldErr[ipt-2] = 0.;
	  negativeYield[ipt-2] = 0.;
	  
	  //	  yield[ipt-2] = mass_sub[ipt]->IntegralAndError( integral_start, integral_end, yieldErr[ipt-2] );
	  for(int j=integral_start; j<=integral_end; j++)
	    {
	      if( mass_sub[ipt]->GetBinContent(j) >= 0)
		{
		  yield[ipt-2]+= mass_sub[ipt]->GetBinContent(j);
		  yieldErr[ipt-2]+= mass_sub[ipt]->GetBinError(j)*mass_sub[ipt]->GetBinError(j);
	      
		}
	      else
		negativeYield[ipt-2]+=mass_sub[ipt]->GetBinContent(j);
	    }
	  
	  yieldErr[ipt-2] = sqrt(yieldErr[ipt-2]);
	  
	  cout << "yield = " << yield[ipt-2] << " negative yield = " << negativeYield[ipt-2] << endl;
	  
	}
      else
	{
	  
	  yield[ipt-2] = BG_gaus[ipt]->Integral( low_intLimit, high_intLimit )/massFG12[ipt]->GetXaxis()->GetBinWidth(1);
	  yieldErr[ipt-2] = 0.;
	}
      
      cout << "yield = " << yield[ipt-2] << endl;
    }
  
  return 1;
}

int analyzePhotons::calculateYield_gaus()
{
  TF1 *mygaus = new TF1("mygaus","gaus",0.,1);
  mygaus->SetLineWidth(1);
  mygaus->SetLineColor(kRed);

  for(int ipt=2; ipt<13; ipt++)
    {

      if(pp && ipt>5)
	break;

      mass_sub[ipt]->Fit("mygaus","R","",0.05,0.3);

      yield_gaus[ipt] = mygaus->Integral(0.100, 0.180)/mass_sub[ipt]->GetBinWidth(1);
      yieldErr_gaus[ipt] = mygaus->IntegralError(0.100,0.180)/mass_sub[ipt]->GetBinWidth(1);
    }

  return 1;
}


int analyzePhotons::calculateRgammaNum()
{

  double cov, sigma_r2;
  
  for(int ipt=0; ipt<11; ipt++)
    {

      if(pp && ipt>3)
	break;

      if(yield[ipt]==0) continue;

      numerator[ipt] = yieldIncl[ipt]/yield[ipt];
      
      // for the error propagation
      cov = yieldIncl_err[ipt]*yieldErr[ipt]/numerator[ipt];
      sigma_r2 = pow(yieldIncl[ipt]/yield[ipt],2)*( pow(yieldIncl_err[ipt]/yieldIncl[ipt],2) + pow(yieldErr[ipt]/yield[ipt],2) - 2.*cov/(yieldIncl[ipt]*yield[ipt]) );
      num_err[ipt] = sqrt(sigma_r2);
    }

  return 1;
}

int analyzePhotons::calculateRgamma(TH1D *somehist)
{

  for(int ipt=0; ipt < 11; ipt++)
    { 

      if(pp && ipt>3)
	break;

      simRatio[ipt] = somehist->GetBinContent(ipt+1);

      rgamma[ipt] = correction[ipt]*numerator[ipt]/simRatio[ipt];
      rgamma_err[ipt] = rgamma[ipt]*num_err[ipt]/numerator[ipt];
    }
  

  return 1;
}

int analyzePhotons::calculateDirectYield(TH1D *somehist, int numberOfEvents)
{

  cout << "direct yield points:" << endl;
  
  for(int ipt=0; ipt<11; ipt++)
    {

      if(pp && ipt>3)
	break;

      int desiredBin = floor(ptindex[ipt]/somehist->GetBinWidth(1))+1;
      hadronYield_sim[ipt] = somehist->GetBinContent(desiredBin);
      
      //applyAbsoluteNorm(somehist->GetBinWidth(1), ptindex[ipt], numberOfEvents);

      yieldDirect[ipt] = hadronYield_sim[ipt]*(rgamma[ipt] - 1.);  // factor of 1000 to put in correct cross-section units
      ////yieldDirect_err[ipt] = yieldDirect[ipt]*(rgamma_err[ipt]/rgamma[ipt]);
      yieldDirect_err[ipt] = hadronYield_sim[ipt]*rgamma_err[ipt];

      cout << yieldDirect[ipt] << ", " << hadronYield_sim[ipt] << endl;

      //yieldDirect[ipt] = hadronYield_sim[ipt]*(rgamma[ipt] - 1.)*1000.;  // factor of 1000 to put in correct cross-section units
      //yieldDirect_err[ipt] = hadronYield_sim[ipt]*rgamma_err[ipt]*1000.;
    }

  return 1;
}

int analyzePhotons::applyAbsoluteNorm(double binwidth, double ptbin, int nevents)
{
  for(int ipt=0; ipt<11; ipt++)
    {
      if(pp && ipt>3)
	break;
      hadronYield_sim[ipt] = hadronYield_sim[ipt]/(2.*TMath::Pi()*ptbin*binwidth*nevents);
    }
  
  return 1;
}


int analyzePhotons::loadInCorrection(TGraphErrors *someplot)
{

  TH1D *hho = new TH1D("hho","hho",1000,0.,5);
  
  double *y = someplot->GetY();
  double *y_e = someplot->GetEY();
  
  for(int i=0; i<someplot->GetN(); i++)
    {
      hho->SetBinContent(i+1, y[i]);
      hho->SetBinError(i+1,y_e[i]);
    }

  for(int ipt=0; ipt<11; ipt++)
    {
      if(pp && ipt>3)
	break;

      int desiredBin = floor(ptindex[ipt]/hho->GetBinWidth(1))+1;

      correction[ipt] = hho->GetBinContent(desiredBin);
      correction_err[ipt] = hho->GetBinError(desiredBin);

      cout << "+_+_+_+_+_+_" << endl;
      cout << correction[ipt] << endl;
    }

  //hack for now
  if(pp)
    {
      correction[0] = 0.15;
      correction[1] = 0.18;
      correction[2] = 0.22;
    }

  return 1;
}

int analyzePhotons::loadInCorrection(TH1D *hho)
{

  for(int ipt=0; ipt<11; ipt++)
    {
      if(pp && ipt>3)
	break;

      int desiredBin = floor(ptindex[ipt]/hho->GetBinWidth(1))+1;

      correction[ipt] = hho->GetBinContent(desiredBin);
      correction_err[ipt] = hho->GetBinError(desiredBin);
    }

  return 1;
}

int analyzePhotons::loadInCorrection(double a, double b, double c)
{

  TF1 *corr_from_fit = new TF1("corr_from_fit","[1]/([2] + exp([0]*x))", 0.4, 5.);   // fermi func

  corr_from_fit->SetParameters(a, b, c);
  
  for(int ipt=0; ipt<11; ipt++)
    {

      if(pp && ipt>3)
	break;

      correction[ipt] = corr_from_fit->Eval(ptindex[ipt]);
      correction_err[ipt] = 0.;  // not sure what to do with this here at the moment, just set to zero for now
    }

  return 1;
  
}
