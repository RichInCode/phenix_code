#include "plotPhotons.h"
#include "analyzePhotons.h"

#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TLine.h"
#include "TArrow.h"
#include "TF1.h"

#include <cmath>
#include <vector>
#include <iostream>

ClassImp(plotPhotons);

plotPhotons::plotPhotons(analyzePhotons *photo): c_FGtoBG(0),
						 c_FG(0),
						 c_BG(0),
						 c_sub(0),
						 c_pi0yield(0),
						 c_inclyield(0),
						 c_inclToPiYield(0),
						 c_directyield(0),
						 c_rgamma(0),
						 c_correction(0),
						 c_cocktailratio(0)
{

  myphotons = photo;

  pp = myphotons->pp;

  if(pp)
    {
      x[0] = 0.6;
      x[1] = 1.;
      x[2] = 1.6;
      x[3] = 3.5;
    }
  else
    {
      x[0] = 0.5;
      x[1] = 0.7;
      x[2] = 0.9;
      x[3] = 1.1;
      x[4] = 1.3;
      x[5] = 1.5;
      x[6] = 1.7;
      x[7] = 1.9;
      x[8] = 2.25;
      x[9] = 3.;
      x[10] = 4.25;
    }

  for(int i=0; i<11; i++)
    xbin[i]=0.;

  
}

void plotPhotons::drawFGtoBG()
{ 
  c_FGtoBG = new TCanvas("c_FGtoBG","c_FGtoBG");
  // c_FGtoBG->Divide(3,4);
  c_FGtoBG->Divide(3,3);

  for(int i=1; i<9; i++)
    {
      if(pp && i>4)
	break;

      c_FGtoBG->cd(i+1);
      myphotons->normRatio_part[i+1]->GetXaxis()->SetTitle("e^{+}e^{-}#gamma mass");
      myphotons->normRatio_part[i+1]->GetYaxis()->SetTitle("counts");
      myphotons->normRatio_part[i+1]->Draw();
    }

  c_FGtoBG->SaveAs("FGtoBG.pdf");
}

void plotPhotons::drawFG()
{
  c_FG = new TCanvas("c_FG","c_FG");
  //c_FG->Divide(3,4);
  c_FG->Divide(3,3);

  gStyle->SetOptStat(0);

  for(int i=1; i<9; i++)
    {
      if(pp && i>4)
	break;

      c_FG->cd(i+1);
      myphotons->massFG12[i+1]->GetXaxis()->SetTitle("e^{+}e^{-}#gamma mass");
      myphotons->massFG12[i+1]->GetYaxis()->SetTitle("counts");
      myphotons->massFG12[i+1]->Draw();
      myphotons->massBG12_norm[i+1]->SetLineColor(kRed);
      myphotons->massBG12_norm[i+1]->Draw("same");
    }

  c_FG->SaveAs("FG.pdf");
  
}

void plotPhotons::drawSub()
{
  c_sub = new TCanvas("c_sub","c_sub");
  // c_sub->Divide(3,4);
  c_sub->Divide(3,3);

  for(int i=1; i<9; i++)
    {
      if(pp && i>4)
	break;

      c_sub->cd(i+1);
      myphotons->mass_sub[i+1]->GetXaxis()->SetTitle("e^{+}e^{-}#gamma mass");
      myphotons->mass_sub[i+1]->GetYaxis()->SetTitle("counts");
      myphotons->mass_sub[i+1]->Draw();
    }

  c_sub->SaveAs("FGsub.pdf");

}

void plotPhotons::drawPi0Yield()
{
  c_pi0yield = new TCanvas("c_pi0yield","c_pi0yield");

  // divide out the bin width
  for(int ipt=0; ipt<11;ipt++)
    {
      if(pp && ipt>3)
	break;
      
      if(pp)
	{
	  if(ipt<2)
	    {
	      myphotons->yield[ipt] = myphotons->yield[ipt]/0.4;
	      myphotons->yieldErr[ipt] = myphotons->yieldErr[ipt]/0.4;
	    }
	  else if(ipt<3)
	    {
	      myphotons->yield[ipt] = myphotons->yield[ipt]/0.8;
	      myphotons->yieldErr[ipt] = myphotons->yieldErr[ipt]/0.8;
	    }
	  else
	    {
	      myphotons->yield[ipt] = myphotons->yield[ipt]/3.;
	      myphotons->yieldErr[ipt] = myphotons->yieldErr[ipt]/3.;
	    }
	}
      else
	{
	  if(ipt<8)
	    {
	      myphotons->yield[ipt] = myphotons->yield[ipt]/0.2;
	      myphotons->yieldErr[ipt] = myphotons->yieldErr[ipt]/0.2;
	    }
	  else if(ipt<9)
	    {
	      myphotons->yield[ipt] = myphotons->yield[ipt]/0.5;
	      myphotons->yieldErr[ipt] = myphotons->yieldErr[ipt]/0.5;
	    }
	  else if(ipt<10)
	    {
	      myphotons->yield[ipt] = myphotons->yield[ipt]/1.;
	      myphotons->yieldErr[ipt] = myphotons->yieldErr[ipt]/1.;
	    }
	  else
	    {
	      myphotons->yield[ipt] = myphotons->yield[ipt]/1.5;
	      myphotons->yieldErr[ipt] = myphotons->yieldErr[ipt]/1.5;
	    }
	}

      xbin[ipt]=0.;
    }
  

  p_pi0yield = new TGraphErrors(11, x, myphotons->yield, xbin, myphotons->yieldErr);
  p_pi0yield->SetMarkerStyle(20);
  p_pi0yield->SetMarkerColor(kBlack);
  p_pi0yield->GetYaxis()->SetTitle("dN/dp_{T}");
  p_pi0yield->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  p_pi0yield->SetTitle("Raw Pion Tagged Photon Yield");
  p_pi0yield->SetName("p_pi0yield");

  if(pp)
    {
      for(int i=10; i>3; i--)
	p_pi0yield->RemovePoint(i);
    }

  c_pi0yield->cd()->SetLogy();
  p_pi0yield->GetXaxis()->SetRangeUser(0,5);
  p_pi0yield->Draw("ap");


  c_pi0yield->SaveAs("pi0yield.pdf");
}

void plotPhotons::drawInclYield()
{
  c_inclyield = new TCanvas("c_inclyield","c_inclyield");

  p_inclYield = new TGraphErrors(11, x, myphotons->yieldIncl, xbin, myphotons->yieldIncl_err);
  p_inclYield->SetName("p_inclYield");
  p_inclYield->SetMarkerStyle(20);
  p_inclYield->SetMarkerColor(kBlack);
  p_inclYield->GetYaxis()->SetTitle("dN/dp_{T}");
  p_inclYield->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  p_inclYield->SetTitle("Raw Inclusive Photon Yield");
  
  if(pp)
    {
      for(int i=10; i>3; i--)
	p_inclYield->RemovePoint(i);
    }

  c_inclyield->cd()->SetLogy();
  p_inclYield->GetXaxis()->SetRangeUser(0,5.);
  p_inclYield->SetName("p_inclyield");
  p_inclYield->Draw("ap");

  c_inclyield->SaveAs("inclyield.pdf");

}


void plotPhotons::drawRawInclToPi()
{
  c_inclToPiYield = new TCanvas("c_inclToPiYield","c_inclToPiYield");

  if(pp)
    {
      double points[4], points_err[4];
      double shortx[4], shortx_bin[4];
      
      for(int i=0; i<4; i++)
	{
	  points[i] = myphotons->numerator[i];
	  points_err[i] = myphotons->num_err[i];

	  shortx[i] = x[i];
	  shortx_bin[i] = xbin[i];
	}

      p_inclToPiYield = new TGraphErrors(4, shortx, points, shortx_bin, points_err);

    }
  else
    p_inclToPiYield = new TGraphErrors(11, x, myphotons->numerator, xbin, myphotons->num_err);

  p_inclToPiYield->SetName("p_inclToPiYield");

  p_inclToPiYield->SetMarkerStyle(20);
  p_inclToPiYield->SetMarkerColor(kBlack);
  p_inclToPiYield->GetYaxis()->SetTitle("incl/pi0 yield (raw)");
  p_inclToPiYield->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  p_inclToPiYield->SetTitle("Inclusive to #pi^{0} tagged Ratio");


  /*
  if(pp)
    {
      double points[4], points_err[4];
      double shortx[4], shortx_bin[4];

      
      
      std::cout << "here brah?" << std::endl;
      for(int i=10; i>3; i--)
	p_inclToPiYield->RemovePoint(i);
    }
  */
  

  std::cout << "points: " << p_inclToPiYield->GetN() << std::endl;
  
  for(int i=0; i<11; i++)
    {
      std::cout << i << std::endl;
      std::cout << myphotons->numerator[i] << " +/- " << myphotons->num_err[i] << std::endl;
      std::cout << x[i] << " +/- " << xbin[i] << std::endl;
    }

  c_inclToPiYield->cd();
  p_inclToPiYield->GetXaxis()->SetRangeUser(0,5.);
  p_inclToPiYield->Draw("ap");

  c_inclToPiYield->SaveAs("inclToPiYield.pdf");
}

void plotPhotons::drawDirectYield()
{

  p_directYield = new TGraphErrors(11, x, myphotons->yieldDirect, xbin, myphotons->yieldDirect_err);
  p_directYield->SetName("p_directYield");
  p_directYield->SetMarkerStyle(20);
  p_directYield->SetMarkerColor(kBlack);
  p_directYield->GetYaxis()->SetTitle("1/nevents 1/2#pi p_{T} d^{2}N/dydp_{T}");
  p_directYield->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  p_directYield->SetTitle("Direct Photon Yield");
  p_directYield->GetYaxis()->SetRangeUser(1e-7,10);

  if(pp)
    {
      for(int i=10; i>3; i--)
	p_directYield->RemovePoint(i);
    }
  
  p_directYield->GetXaxis()->SetRangeUser(0.,5.);
  if(!c_directyield)
    {
      c_directyield = new TCanvas("c_directyield","c_directyield");
      c_directyield->cd()->SetLogy();
      p_directYield->Draw("ap");
    }
  else
    {
      c_directyield->cd()->SetLogy();
      p_directYield->Draw("psame");
    }

  c_directyield->SaveAs("directYield.pdf");

  
}

void plotPhotons::drawDirectYield_withSys()
{
  
  /**************
systematic on cocktail: 10%

   **************/

  double sys[11];
  double xsys[11];
  
  double perr[11];

  std::vector<TArrow> varrow;
  std::vector<TLine> vline;
  
  for(int i=0; i<11; i++)
    {
      sys[i] = myphotons->yieldDirect[i]*sqrt( (0.04*0.04 + 0.1*0.1) );
      xsys[i] = 0.1;
      //perr[i] = fabs(myphotons->yieldDirect[i])*1.04*myphotons->rgamma_err[i]/myphotons->rgamma[i];
      perr[i] = myphotons->yieldDirect_err[i];
      //perr[i] = myphotons->hadronYield_sim[i]*1.04*myphotons->rgamma_err[i]/myphotons->rgamma[i];

      std::cout << myphotons->rgamma_err[i]/myphotons->rgamma[i] << std::endl;
      std::cout << myphotons->yieldDirect[i] << " +/- " << perr[i] << std::endl;
      
      if(myphotons->yieldDirect[i] < 0.)
	{
	  //std::cout << "preparing errors:" << std::endl;
	  //std::cout << x[i] << " " << myphotons->yieldDirect[i]+perr[i] << std::endl; 

	  varrow.push_back( TArrow(x[i],
				   1e-5,
				   x[i],
				   myphotons->yieldDirect[i]+perr[i]) );
	  vline.push_back( TLine(x[i]-0.1,
				 myphotons->yieldDirect[i]+perr[i],
				 x[i]+0.1,
				 myphotons->yieldDirect[i]+perr[i]) );
	}
    }

  p_directYield = new TGraphErrors(11, x, myphotons->yieldDirect, xbin, perr);
  
  p_directYieldSys = new TGraphErrors(11, x, myphotons->yieldDirect, xsys, sys);
  p_directYieldSys->SetMarkerStyle(20);
  p_directYieldSys->SetMarkerColor(kBlack);
  p_directYieldSys->GetYaxis()->SetTitle("1/nevents 1/2#pi p_{T} d^{2}N/dydp_{T}");
  p_directYieldSys->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  p_directYieldSys->SetTitle("Direct Photon Yield");
  p_directYieldSys->GetYaxis()->SetRangeUser(1e-7,10);

  c_directyield->cd()->SetLogy();
  p_directYieldSys->GetXaxis()->SetRangeUser(0.,2.);
  p_directYieldSys->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  p_directYieldSys->GetYaxis()->SetTitle("#gamma^{direct}");
  p_directYieldSys->SetTitle("Direct Photon Yield");
  p_directYieldSys->SetFillStyle(1001);
  p_directYieldSys->SetFillColor(14);
  p_directYieldSys->Draw("a2");
  p_directYieldSys->Draw("pX");
  p_directYield->Draw("psame");

  for(unsigned int i=0; i<varrow.size(); i++)
    {
      varrow[i].Draw("same");
     
    }

  for(unsigned int i=0; i<vline.size(); i++)
    {
      vline[i].Draw("same");
     
    }

}

void plotPhotons::drawRgamma()
{
  c_rgamma = new TCanvas("c_rgamma","c_rgamma");

  p_rgamma = new TGraphErrors(11, x, myphotons->rgamma, xbin, myphotons->rgamma_err);
  p_rgamma->SetName("p_rgamma");
  p_rgamma->SetMarkerStyle(20);
  p_rgamma->SetMarkerColor(kBlack);
  p_rgamma->GetYaxis()->SetTitle("dN/dp_{T}");
  p_rgamma->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  p_rgamma->SetTitle("Direct Photon Yield");

  if(pp)
    {
      for(int i=10; i>3; i--)
	p_rgamma->RemovePoint(i);
    }
  
  TH1D *dum = new TH1D("dum","dum",1,0,2);
  dum->GetYaxis()->SetRangeUser(0.8, 2.5);
  dum->GetYaxis()->SetTitle("R_{#gamma} = #gamma^{incl}/#gamma^{hadron}");
  dum->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  dum->SetTitle("");
  dum->Draw();
  //  c_rgamma->cd()->SetLogy();
  p_rgamma->Draw("psame");

  c_rgamma->SaveAs("rgamma.pdf");

  /*
  TLegend *leg1 = new TLegend(0.6,0.6,0.9,0.9);
  leg1->AddEntry(p_rgamma,"external conversions","p");
  leg1->AddEntry(p_ppg086pp,"ppg086","p");
  leg1->SetFillColor(0);
  leg1->Draw("same");
  */
}

void plotPhotons::drawRgamma_withSys()
{

  double sys[11];
  
  double xsys[11];

  double perr[11];
  
  for(int i=0; i<11; i++)
    {
      sys[i] = 0.04*myphotons->rgamma[i];
      xsys[i] = 0.1;
      perr[i] = 1.04*myphotons->rgamma_err[i];
    }
  
  // add in the type A systematic
  p_rgamma = new TGraphErrors(11, x, myphotons->rgamma, xbin, perr);

  p_rgammaSys = new TGraphErrors(11, x, myphotons->rgamma, xsys, sys);
  p_rgammaSys->SetName("p_rgammaSys");
  p_rgammaSys->SetMarkerStyle(20);
  p_rgammaSys->SetMarkerColor(kBlack);
  p_rgammaSys->SetFillStyle(1001);
  p_rgammaSys->SetFillColor(14);
  p_rgammaSys->GetYaxis()->SetTitle("#gamma^{incl}/#gamma^{hadron}");
  p_rgammaSys->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  p_rgammaSys->SetTitle("R_{#gamma}");
  p_rgammaSys->GetYaxis()->SetRangeUser(0.5,2.5);
  p_rgammaSys->GetXaxis()->SetRangeUser(0.,2.);

  c_rgamma->cd();

  p_rgammaSys->Draw("a2");
  p_rgammaSys->Draw("pX");
  p_rgamma->Draw("psame");

  TLine *aline = new TLine(0,1,2,1);
  aline->SetLineStyle(2);
  aline->SetLineColor(kBlack);
  aline->Draw("same");
  
}

void plotPhotons::drawPPG088()
{
  double ppg086_pt[6] = {1.18449, 1.69126, 2.19773, 2.70364, 3.3488, 4.37494};
  double ppg086_ptbin[6] = {0., 0., 0., 0., 0., 0.};
  double ppg086_r[6] = {-0.00107015, 0.0259054, 0.0337655, 0.0389738, 0.0307654, 0.0352367};
  double ppg086_r_err[6] = {0.00392405, 0.00603416, 0.00989513, 0.0158787, 0.0200106, 0.0445122};
  double ppg086_r_sys[6] = {0.0116212, 0.0123477, 0.0130569, 0.0150306, 0.0148773, 0.0214418};

  double ppg086_Rgamma[6], ppg086_Rgamma_err[6], ppg086_Rgamma_sys[6];

  for(int i=0; i<6; i++)
    {
      ppg086_Rgamma[i] = 1./(1. - ppg086_r[i]);
      //      ppg086_Rgamma_err[i] = ppg086_Rgamma[i]*ppg086_r_err[i]/ppg086_r[i];
      //ppg086_Rgamma_sys[i] = ppg086_Rgamma[i]*ppg086_r_sys[i]/ppg086_r[i];
      ppg086_Rgamma_err[i] = (1./(1. - (ppg086_r[i] + ppg086_r_err[i]))) - ppg086_Rgamma[i];
      ppg086_Rgamma_sys[i] = (1./(1. - (ppg086_r[i] + ppg086_r_sys[i]))) - ppg086_Rgamma[i];
    }

  TGraph *p_ppg086pp = new TGraphErrors(6, ppg086_pt, ppg086_Rgamma, ppg086_ptbin, ppg086_Rgamma_err);
  p_ppg086pp->SetMarkerStyle(24);
  p_ppg086pp->SetMarkerColor(kRed);
  p_ppg086pp->SetLineColor(kRed);

  c_rgamma->cd();
  p_ppg086pp->Draw("psame");


  // direct yield
  double ptd[6] = {1.18449, 1.69126, 2.19773, 2.70364, 3.3488, 4.37494};
  double yd[6] = {-0.000231666, 0.000669745, 0.000163268, 4.65505e-5, 7.78764e-6, 1.2448e-6};
  double yst[6] = {0.000848211, 0.000156004, 4.78465e-5, 1.89656e-5, 5.06527e-6, 1.57247e-6};
  //double ysys[6] = {0.00251232, 0.000339391, 6.91015e-5, 1.9658e-5, 3.99712e-6, 7.87164e-7};

  TGraphErrors *p_ppg088pp_d = new TGraphErrors(6, ptd, yd, 0, yst);
  p_ppg088pp_d->SetMarkerStyle(24);
  p_ppg086pp->SetMarkerColor(kRed);
  p_ppg086pp->SetLineColor(kRed);
  
  c_directyield->cd();
  p_ppg088pp_d->Draw("psame");

  
}

void plotPhotons::drawPPG162()
{

  c_directyield = new TCanvas("c_directyield","c_directyield");

  // just draw the 0-20%
  double pt[11] = {0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.25, 3., 4.25};
  double y[11] = {4.162295, 3.863694, 1.531712, 0.500764, 0.248402, 0.105543, 0.071944, 0.028950, 0.007928, 0.001125, 0.000025};
  double y_stat[11] = {4.162295, 0.908128, 0.313389, 0.119182, 0.053970, 0.025248, 0.014094, 0.007115, 0.001929, 0.000274, 0.000033};
  double y_sys[11] = {4.963615, 1.412782, 0.499280, 0.184227, 0.077102, 0.032901, 0.016161, 0.007324, 0.002126, 0.000242, 0.000012};

  double ptbin[11];

  // scale the value down to pp value by ncol;
  for(int i=0; i<11; i++)
    {
      ptbin[i] = 0.1;
      y[i] = y[i]/(770.6/42.);
      y_stat[i] = y_stat[i]/(770.6/42.);
      y_sys[i] = y_sys[i]/(770.6/42.);
    }

  TGraphErrors *p_162 = new TGraphErrors(11, pt, y, 0, y_stat);
  TGraphErrors *p_162_sys = new TGraphErrors(11, pt, y, ptbin, y_sys);

  p_162_sys->SetFillStyle(1001);
  p_162_sys->SetFillColor(11);

  c_directyield->cd();
  p_162_sys->Draw("a2");
  p_162_sys->Draw("px");
  p_162->SetMarkerStyle(20);
  p_162->SetMarkerColor(kBlue);
  p_162->SetLineColor(kBlue);
  p_162->Draw("psame");
  
}

void plotPhotons::drawPPfit()
{
  /*
  TF1 *func = new TF1("func","[0]*pow(1. + x*x/[1],[2])",0,2);
  func->SetParameters(8.3e-3, 2.26, -3.45);
  func->SetLineColor(kRed);
  func->SetLineWidth(2);

  TF1 *func_high = new TF1("func_high","[0]*pow(1. + x*x/[1],[2])",0,2);
  func_high->SetParameters(func->GetParameter(0)+7.5e-3, 2.26+0.78, -3.45+0.08);

  TF1 *func_low = new TF1("func_low","[0]*pow(1. + x*x/[1],[2])",0,2);
  func_low->SetParameters(func->GetParameter(0)-7.5e-3, 2.26-0.78, -3.45-0.08);
  */
  
  TFile *ppin = new TFile("/phenix/u/workarea/rpetti/analysisRun9/macros/pp_data_fit_histogram.root");
  TH1D *h_powerfit111_2 = (TH1D*)ppin->Get("h_powerfit111_2");
  h_powerfit111_2->SetLineColor(kBlue);
  h_powerfit111_2->SetMarkerColor(kBlue);

  c_directyield->cd();
  /*
  func->Draw("same");
  func_high->Draw("same");
  func_low->Draw("same");
  */

  h_powerfit111_2->Draw("same");
}

void plotPhotons::drawCorrection()
{
  c_correction = new TCanvas("c_correction","c_correction");

  p_correction = new TGraph(11, x, myphotons->correction);
  p_correction->SetName("p_correction");
  p_correction->SetMarkerStyle(20);
  p_correction->SetMarkerColor(kBlack);
  p_correction->GetYaxis()->SetTitle("#varepsilon f");
  p_correction->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  p_correction->SetTitle("#pi^{0} tagging efficiency");

  if(pp)
    {
      for(int i=10; i>3; i--)
	p_correction->RemovePoint(i);
    }
  
  p_correction->GetXaxis()->SetRangeUser(0.,2.);
  p_correction->Draw("ap");
}

void plotPhotons::drawCocktailRatio()
{
  c_cocktailratio = new TCanvas("c_cocktailratio","c_cocktailratio");

  p_cocktailratio = new TGraph(11, x, myphotons->simRatio);
  p_cocktailratio->SetName("p_cocktailratio");
  p_cocktailratio->SetMarkerStyle(20);
  p_cocktailratio->SetMarkerColor(kBlack);
  p_cocktailratio->GetXaxis()->SetRangeUser(0,2);
  p_cocktailratio->GetYaxis()->SetTitle("#gamma^{hadron}/#gamma^{#pi^{0}}");
  p_cocktailratio->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  p_cocktailratio->SetTitle("Cocktail Ratio");

  if(pp)
    {
      for(int i=10; i>3; i--)
	p_cocktailratio->RemovePoint(i);
    }
  
  p_cocktailratio->GetXaxis()->SetRangeUser(0.,5.);
  p_cocktailratio->Draw("ap");

  c_cocktailratio->SaveAs("cocktailRatio.pdf");
}

void plotPhotons::SaveRgamma(TFile *file)
{
  file->cd();
  
  p_rgamma->Write();
}

void plotPhotons::SaveCorrectionUsed(TFile *file)
{
  file->cd();

  p_correction->Write();
}

void plotPhotons::SaveRgammaNumerator(TFile *file)
{
  file->cd();

  p_inclToPiYield->Write();
}

void plotPhotons::SaveCocktailRatio(TFile *file)
{
  file->cd();

  p_cocktailratio->Write();
}

void plotPhotons::SaveDirectYield(TFile *file)
{
  file->cd();
  
  p_directYield->Write();
}

void plotPhotons::SaveInclYield(TFile *file)
{
  file->cd();
  
  p_inclYield->Write();
}

void plotPhotons::SavePi0Yield(TFile *file)
{
  file->cd();

  p_pi0yield->Write();
}

void plotPhotons::SaveAll(TFile *file)
{
  SaveRgamma(file);
  SaveCorrectionUsed(file);
  SaveRgammaNumerator(file);
  SaveCocktailRatio(file);
  SaveDirectYield(file);
  SaveInclYield(file);
  SavePi0Yield(file);
}
