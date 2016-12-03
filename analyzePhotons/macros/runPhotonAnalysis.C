void runPhotonAnalysis(char *piinfile = "/phenix/plhf/rpetti/taxi/Run7AuAu200/Run7_Conversions__anatrain_taxi294/sum_total/postCabana_pi0_histos_allHBD_taxi294.root", char *inclinfile = "/phenix/plhf/rpetti/taxi/Run7AuAu200/Run7_Conversions__anatrain_taxi294/sum_total/inclusive_histos_allHBD_taxi294.root", char* corrinfile = "/phenix/u/workarea/rpetti/fromCip/o_gp_sim_pl_0.root", char* simratioinfile = "/phenix/u/workarea/rpetti/fromBen/withParents/cocktail_ratio_0.root", char* hadroninfile = "/phenix/u/workarea/rpetti/fromBen/withParents/cent0.root", int centindex = 0, int normSchemeIndex = 1)
{
  //=-=-=-=-=-=-==-==-
  // load the library
  gSystem->Load("libanalyzePhotons.so");

  //readInFileAndHistos("/direct/phenix+plhf/rpetti/taxi/Run7AuAu200/Run7_Conversions__anatrain_taxi294/sum_total/postCabana_pi0_histos_allHBD_taxi294.root");

  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  // read in the file and grab the histograms
  TFile *inputfile = new TFile(piinfile);

  inputfile->ls();

  TH3D *h_FG_mass_vs_pt_vs_pt[5];
  TH3D *h_BG_mass_vs_pt_vs_pt[5];

  // FG
  h_FG_mass_vs_pt_vs_pt[0] = (TH3D*)inputfile->Get("ATM_mass_vs_ATM_pT_vs_ph_pT_FG12_0")->Clone("h_FG_mass_vs_pt_vs_pt_0");
  h_FG_mass_vs_pt_vs_pt[1] = (TH3D*)inputfile->Get("ATM_mass_vs_ATM_pT_vs_ph_pT_FG12_1")->Clone("h_FG_mass_vs_pt_vs_pt_1");
  h_FG_mass_vs_pt_vs_pt[2] = (TH3D*)inputfile->Get("ATM_mass_vs_ATM_pT_vs_ph_pT_FG12_2")->Clone("h_FG_mass_vs_pt_vs_pt_2");
  h_FG_mass_vs_pt_vs_pt[3] = (TH3D*)inputfile->Get("ATM_mass_vs_ATM_pT_vs_ph_pT_FG12_3")->Clone("h_FG_mass_vs_pt_vs_pt_3");
  h_FG_mass_vs_pt_vs_pt[4] = (TH3D*)inputfile->Get("ATM_mass_vs_ATM_pt_vs_pt_FG_MB")->Clone("h_FG_mass_vs_pt_vs_pt_4");

  // BG
  h_BG_mass_vs_pt_vs_pt[0] = (TH3D*)inputfile->Get("ATM_mass_vs_ATM_pT_vs_ph_pT_BG12_0")->Clone("h_BG_mass_vs_pt_vs_pt_0");
  h_BG_mass_vs_pt_vs_pt[1] = (TH3D*)inputfile->Get("ATM_mass_vs_ATM_pT_vs_ph_pT_BG12_1")->Clone("h_BG_mass_vs_pt_vs_pt_1");
  h_BG_mass_vs_pt_vs_pt[2] = (TH3D*)inputfile->Get("ATM_mass_vs_ATM_pT_vs_ph_pT_BG12_2")->Clone("h_BG_mass_vs_pt_vs_pt_2");
  h_BG_mass_vs_pt_vs_pt[3] = (TH3D*)inputfile->Get("ATM_mass_vs_ATM_pT_vs_ph_pT_BG12_3")->Clone("h_BG_mass_vs_pt_vs_pt_3");
  h_BG_mass_vs_pt_vs_pt[4] = (TH3D*)inputfile->Get("ATM_mass_vs_ATM_pt_vs_pt_BG_MB")->Clone("h_BG_mass_vs_pt_vs_pt_4");


  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  // pass the relevant histogram to the analyzer
  analyzePhotons *photon = new analyzePhotons(0,0);


  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  // run the pi0 extraction analysis
  photon->runPi0Analysis(h_FG_mass_vs_pt_vs_pt[centindex], h_BG_mass_vs_pt_vs_pt[centindex], normSchemeIndex);
  

  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-
  // plot the results of the analysis
  plotPhotons *plotter = new plotPhotons(photon);
 
  plotter->drawFG();
  plotter->drawFGtoBG();
  plotter->drawSub();
  plotter->drawPi0Yield();
  

  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-==-
  // read in the inclusive photon file and grab the histograms
  TFile *inputfile2 = new TFile(inclinfile);

  TH3D *h_convertmass_vs_pt_vs_cut[5];
  
  h_convertmass_vs_pt_vs_cut[0] = (TH3D*)inputfile2->Get("recal_mass_vs_recal_pair_pt_vs_cut_0");
  h_convertmass_vs_pt_vs_cut[1] = (TH3D*)inputfile2->Get("recal_mass_vs_recal_pair_pt_vs_cut_1");
  h_convertmass_vs_pt_vs_cut[2] = (TH3D*)inputfile2->Get("recal_mass_vs_recal_pair_pt_vs_cut_2");
  h_convertmass_vs_pt_vs_cut[3] = (TH3D*)inputfile2->Get("recal_mass_vs_recal_pair_pt_vs_cut_3");
  
  // make the min. bias
  h_convertmass_vs_pt_vs_cut[4] = (TH3D*)h_convertmass_vs_pt_vs_cut[0]->Clone("h_convertmass_vs_pt_vs_cut_4");
  h_convertmass_vs_pt_vs_cut[4]->Add(h_convertmass_vs_pt_vs_cut[1]);
  h_convertmass_vs_pt_vs_cut[4]->Add(h_convertmass_vs_pt_vs_cut[2]);
  h_convertmass_vs_pt_vs_cut[4]->Add(h_convertmass_vs_pt_vs_cut[3]);


  photon->runInclAnalysis(h_convertmass_vs_pt_vs_cut[centindex]);
  plotter->drawInclYield();

  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  // read in hadron sim information
  TFile *insimratio = new TFile(simratioinfile);
  TH1D *h_ratio = insimratio->Get("ratio_rebin");
  TFile *inhadron = new TFile(hadroninfile);
  
  TH1D *h_hadrons = inhadron->Get("inc");
  TH1D *h_pi = inhadron->Get("pi0");
  int npievents = h_pi->GetEntries()/2.;

  //=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-
  // read in the correction and put it in memory
  TFile *incorr = new TFile(corrinfile);
  TGraphErrors *p_corr = incorr->Get("ho");
  photon->loadInCorrection(p_corr);

  //-=-=-=-=-=-=-=-=-==-
  // R_gamma analysis
  photon->calculateRgamma(h_ratio);
  plotter->drawRgamma();


  //=-=-=-=-=-=-=-==-=-=-=-=-=
  // direct photon analysis
  photon->calculateDirectYield(h_hadrons, npievents);
  plotter->drawDirectYield();

}
