using namespace RooFit;
using namespace RooStats ;

static const Int_t NCAT = 4;  // chiara
Int_t MINmass= 130;
Int_t MAXmass= 1000;
std::string filePOSTfix="";
double signalScaler=1.00;

void AddSigData(RooWorkspace*, Float_t);
void AddBkgData(RooWorkspace*, Float_t);

void SigModelResponseFcnFit(RooWorkspace*);
void SigModelFitConvBW(RooWorkspace*, Float_t, Double_t);



TPaveText* get_labelCMS( int legendQuadrant = 0 , std::string year="2012", bool sim=false) {

  if( legendQuadrant!=0 && legendQuadrant!=1 && legendQuadrant!=2 && legendQuadrant!=3 ) {
    std::cout << "WARNING! Legend quadrant '" << legendQuadrant << "' not yet implemented for CMS label. Using 2." << std::endl;
    legendQuadrant = 2;
  }

  float x1, y1, x2, y2;
  if( legendQuadrant==1 ) {
    x1 = 0.63;
    y1 = 0.83;
    x2 = 0.8;
    y2 = 0.87;
  } else if( legendQuadrant==0 ) {
    x1 = 0.175;
    y1 = 0.953;
    x2 = 0.6;
    y2 = 0.975;
  } else if( legendQuadrant==3 ) {
    x1 = 0.25;
    y1 = 0.2;
    x2 = 0.42;
  }
 
  TPaveText* cmslabel = new TPaveText( x1, y1, x2, y2, "brNDC" );
  cmslabel->SetFillColor(kWhite);
  cmslabel->SetTextSize(0.038);
  if( legendQuadrant==0 ) cmslabel->SetTextAlign(11);
  cmslabel->SetTextSize(0.038);
  cmslabel->SetTextFont(42);
  
    std::string leftText;
   
     
    if (sim)  leftText = "CMS Simulation"; 
    else {
     leftText = "CMS Preliminary, 19.5 fb^{-1}";
    }
    cmslabel->AddText(leftText.c_str());
    return cmslabel;

}




TPaveText* get_labelSqrt( int legendQuadrant ) {

  if( legendQuadrant!=0 && legendQuadrant!=1 && legendQuadrant!=2 && legendQuadrant!=3 ) {
    std::cout << "WARNING! Legend quadrant '" << legendQuadrant << "' not yet implemented for Sqrt label. Using 2." << std::endl;
    legendQuadrant = 2;
  }


  float x1, y1, x2, y2;
  if( legendQuadrant==1 ) {
    x1 = 0.63;
    y1 = 0.78;
    x2 = 0.8;
    y2 = 0.82;
  } else if( legendQuadrant==2 ) {
    x1 = 0.25;
    y1 = 0.78;
    x2 = 0.42;
    y2 = 0.82;
  } else if( legendQuadrant==3 ) {
    x1 = 0.25;
    y1 = 0.16;
    x2 = 0.42;
    y2 = 0.2;
  } else if( legendQuadrant==0 ) {
    x1 = 0.65;
    y1 = 0.953;
    x2 = 0.87;
    y2 = 0.975;
  }


  TPaveText* label_sqrt = new TPaveText(x1,y1,x2,y2, "brNDC");
  label_sqrt->SetFillColor(kWhite);
  label_sqrt->SetTextSize(0.038);
  label_sqrt->SetTextFont(42);
  label_sqrt->SetTextAlign(31); 
  label_sqrt->AddText("#sqrt{s} = 8 TeV");
  return label_sqrt;

}



RooArgSet* defineVariables() {

  // define variables of the input ntuple 
  RooRealVar* PhotonsMass  = new RooRealVar("PhotonsMass", "M(gg)",MINmass, MAXmass,"GeV");
  RooRealVar* PhotonsMassTrue  = new RooRealVar("PhotonsMassTrue", "M(gg)",130, 1000,"GeV");
  RooRealVar* PhotonsMass_scaleUP  = new RooRealVar("PhotonsMass_scaleUP", "M(gg)",130, 1000,"GeV");
  RooRealVar* PhotonsMass_scaleDN  = new RooRealVar("PhotonsMass_scaleDN", "M(gg)",130, 1000,"GeV");

  RooRealVar* ph1_eta = new RooRealVar("ph1_eta", "eta(g1)",-10,10,"");
  RooRealVar* ph2_eta = new RooRealVar("ph2_eta", "eta(g2)",-10,10,"");
  RooRealVar* ph1_r9 = new RooRealVar("ph1_r9", "R9(g1)",-10,10,"");
  RooRealVar* ph2_r9 = new RooRealVar("ph2_r9", "R9(g2)",-10,10,"");
  RooRealVar* evweight = new RooRealVar("evweight","weightings",0,1000,"");
  RooRealVar* nvtx = new RooRealVar("nvtx", "nvtx", 0, 50, "");
  RooArgSet* ntplVars = new RooArgSet(*PhotonsMass,*PhotonsMassTrue, *ph1_eta, *ph2_eta, *ph1_r9, *ph2_r9, *evweight, *nvtx);
  
  return ntplVars;
}

RooArgSet* defineVariablesM250() {

  // define variables of the input ntuple //livia
  RooRealVar* PhotonsMass  = new RooRealVar("PhotonsMass", "M(gg)",130, 1000,"GeV");
  RooRealVar* PhotonsMassTrue  = new RooRealVar("PhotonsMassTrue", "M(gg)",130, 1000,"GeV");
  RooRealVar* PhotonsMass_scaleUP  = new RooRealVar("PhotonsMass_scaleUP", "M(gg)",130, 1000,"GeV");
  RooRealVar* PhotonsMass_scaleDN  = new RooRealVar("PhotonsMass_scaleDN", "M(gg)",130, 1000,"GeV");
  RooRealVar* ph1_eta = new RooRealVar("ph1_eta", "eta(g1)",-10,10,"");
  RooRealVar* ph2_eta = new RooRealVar("ph2_eta", "eta(g2)",-10,10,"");
  RooRealVar* ph1_r9 = new RooRealVar("ph1_r9", "R9(g1)",-10,10,"");
  RooRealVar* ph2_r9 = new RooRealVar("ph2_r9", "R9(g2)",-10,10,"");
  RooRealVar* evweight = new RooRealVar("evweight","Reweightings",0,1000,"");
  RooRealVar* nvtx = new RooRealVar("nvtx", "nvtx", 0, 50, "");
  //  RooRealVar* btagCategory = new RooRealVar("btagCategory","event category",0.9,2.1,"") ;
  //  RooRealVar* newweight = new RooRealVar("newweight","Reweightings",0., 1000,"");
  RooArgSet* ntplVars = new RooArgSet(*PhotonsMass,*PhotonsMassTrue,*PhotonsMass_scaleUP, *PhotonsMass_scaleDN, *ph1_eta, *ph2_eta, *ph1_r9, *ph2_r9, *evweight);
  return ntplVars;
}

RooArgSet* defineVariables_newWeight() {

  
  RooRealVar* PhotonsMass  = new RooRealVar("PhotonsMass", "M(gg)",MINmass, MAXmass,"GeV");
  RooRealVar* ph1_eta = new RooRealVar("ph1_eta", "eta(g1)",-10,10,"");
  RooRealVar* ph2_eta = new RooRealVar("ph2_eta", "eta(g2)",-10,10,"");
  RooRealVar* ph1_r9 = new RooRealVar("ph1_r9", "R9(g1)",-10,10,"");
  RooRealVar* ph2_r9 = new RooRealVar("ph2_r9", "R9(g2)",-10,10,"");
  RooRealVar* nvtx = new RooRealVar("nvtx", "nvtx", 0, 50, "");
  RooRealVar* newweight = new RooRealVar("newweight","Reweightings",0., 1000,"");
  RooArgSet* ntplVars = new RooArgSet(*PhotonsMass, *ph1_eta, *ph2_eta, *ph1_r9, *ph2_r9,*newweight,   *nvtx);
  
  return ntplVars;
}



void runfits(const Float_t mass=150, Bool_t dobands = false, Float_t width=0.1) {

  //******************************************************************//
  //  Running mode  corresponds to the following cases
  //         - full run set:
  //         - create signal and background data sets 
  //         - make and fit signal and background  models 
  //         - write signal and background workspaces in root files
  //         - write data card
  //*******************************************************************//

  TString fileBaseName("HighMass-hgg");    
  
  TString card_name("HighMass-hgg_models_Bkg_8TeV_test_SYST.rs");
  HLFactory hlf("HLFactory", card_name, false);
  RooWorkspace* w = hlf.GetWs();
 

  // Luminosity:
  Float_t Lum = 19500.0;  
  RooRealVar lumi("lumi","lumi",Lum);
  w->import(lumi); 
  
  MINmass=130;
  MAXmass=1000;
  
  std::cout<<"  MIN MASS: "<<MINmass<<"   MAX MASS: "<<MAXmass<<std::endl;
  

  Double_t MMIN = MINmass; 
  Double_t MMAX = MAXmass; 
  w->var("PhotonsMass")->setMin(MMIN);
  w->var("PhotonsMass")->setMax(MMAX);


  //  w->Print("v");
  
  cout << endl; cout << "Now AddSigData" << endl;
  // AddSigData(w, mass, 0);   
 
  cout << endl; cout << "Now SigModelFit" << endl;
  
  bool isEscale=false;
  bool isEsmear=true;
  if(isEscale){
  //ENERGY SCALE
  std::cout<<"----------------------------------- STANDARD --------------------------------------------"<<std::endl;
  RooPlot* p;
  p_noscale = SigModelResponseFcnFit(w, mass,0);
  std::cout<<"----------------------------------- SCALE UP --------------------------------------------"<<std::endl;
  RooPlot* p_scaleUP;
  p_scaleUP=SigModelResponseFcnFit(w, mass,1);
  std::cout<<"----------------------------------- SSCALE DN --------------------------------------------"<<std::endl;
  RooPlot* p_scaleDN;
  p_scaleDN=SigModelResponseFcnFit(w, mass,2);
  
  TCanvas* c1 = new TCanvas("c1", "c1", 1);
  c1->cd();
  TPaveText* label_cms = get_labelCMS(0, "2012", true);
  TPaveText* label_sqrt = get_labelSqrt(0);
  p_noscale->Draw();
  p_noscale->GetXaxis()->SetTitle("#frac{m_{#gamma#gamma}}{m_{H}} -1");
  p_scaleUP->Draw("same");
  p_scaleDN->Draw("same");
  label_cms->Draw("same");
  label_sqrt->Draw("same");
  TH1F* h_noscale= new TH1F("h_noscale", "h", 600, -0.5, 0.5);
  h_noscale->SetMarkerColor(kGray+2);
  h_noscale->SetLineColor(kBlack);
  TH1F* h_scaleup= new TH1F("h_scaleUP", "h", 600, -0.5, 0.5);
  h_scaleUP->SetMarkerColor(kBlue-2);
  h_scaleUP->SetLineColor(kBlue);
  TH1F* h_scaledn= new TH1F("h_scaleDN", "h", 600, -0.5, 0.5);
  h_scaleDN->SetMarkerColor(kRed-2);
  h_scaleDN->SetLineColor(kRed);

  TLegend* legmc = new TLegend(0.6, 0.6, 0.85, 0.89, "", "brNDC");
  legmc->SetTextSize(0.0206044);  
  legmc->SetTextFont(42);
  legmc->SetBorderSize(0);
  legmc->SetFillStyle(0);
  legmc->SetHeader("CAT: 2");
  legmc->AddEntry(h_noscale, "Standard", "PL");
  legmc->AddEntry(h_scaleup, "scale UP", "PL");
  legmc->AddEntry(h_scaledn, "scale DOWN", "PL");
  legmc->Draw("same");
  c1->SaveAs("systPlot/SignalShape_scale_cat3.png");
  c1->SaveAs("systPlot/SignalShape_scale_cat3.pdf");
  }

  if(isEsmear){
  //ENERGY SMEAR
  std::cout<<"----------------------------------- STANDARD --------------------------------------------"<<std::endl;
  RooPlot* p_nosmear;
  p_nosmear = SigModelResponseFcnFitRES(w, mass,0);
  std::cout<<"----------------------------------- SMEAR UP --------------------------------------------"<<std::endl;
  RooPlot* p_smearUP;
  p_smearUP=SigModelResponseFcnFitRES(w, mass,1);
  std::cout<<"----------------------------------- SMEAR DN --------------------------------------------"<<std::endl;
  RooPlot* p_smearDN;
  p_smearDN=SigModelResponseFcnFitRES(w, mass,2);
  
  TCanvas* c1 = new TCanvas("c1", "c1", 1);
  c1->cd();
  TPaveText* label_cms = get_labelCMS(0, "2012", true);
  TPaveText* label_sqrt = get_labelSqrt(0);
  p_smearDN->Draw();
  p_smearDN->GetXaxis()->SetTitle("#frac{m_{#gamma#gamma}}{m_{H}} -1");
  p_nosmear->Draw("same");
  p_smearUP->Draw("same");
  label_cms->Draw("same");
  label_sqrt->Draw("same");
  TH1F* h_nosmear= new TH1F("h_nosmear", "h", 600, -0.5, 0.5);
  h_nosmear->SetMarkerColor(kGray+2);
  h_nosmear->SetLineColor(kBlack);
  TH1F* h_smearup= new TH1F("h_smearUP", "h", 600, -0.5, 0.5);
  h_smearUP->SetMarkerColor(kBlue-2);
  h_smearUP->SetLineColor(kBlue);
  TH1F* h_smeardn= new TH1F("h_smearDN", "h", 600, -0.5, 0.5);
  h_smearDN->SetMarkerColor(kRed-2);
  h_smearDN->SetLineColor(kRed);

  TLegend* legmc = new TLegend(0.6, 0.6, 0.85, 0.89, "", "brNDC");
  legmc->SetTextSize(0.0206044);  
  legmc->SetTextFont(42);
  legmc->SetBorderSize(0);
  legmc->SetFillStyle(0);
  legmc->SetHeader("CAT: 3");
  legmc->AddEntry(h_nosmear, "Standard", "PL");
  legmc->AddEntry(h_smearup, "smear UP", "PL");
  legmc->AddEntry(h_smeardn, "smear DOWN", "PL");
  legmc->Draw("same");
  c1->SaveAs("systPlot/SignalShape_smear_cat3.png");
  c1->SaveAs("systPlot/SignalShape_smear_cat3.pdf");
  }


  w->var("PhotonsMass")->setMin(MINmass);
  w->var("PhotonsMass")->setMax(MAXmass);
 
  // SigModelFitConvBW(w, mass, width);      
  
  return;
}

// Signal Data Set
void AddSigData(RooWorkspace* w, Float_t mass, int Sys) {

  Int_t ncat = NCAT;
  TString inDir = "";

  Float_t MASS(mass);

  
  // Variables
  RooArgSet* ntplVars = defineVariables();
  
  int iMass = abs(mass);      
  TFile sigFile1("histograms_CMS-HGG_19032013.root");  
  
  //chain summing up all production modes
  TChain* sigTree1  = new TChain();
  sigTree1->Add(TString::Format("histograms_CMS-HGG_19032013_WithSyst.root/sigTree1", iMass));
  /* sigTree1->Add(TString::Format("histograms_CMS-HGG_19032013.root/ggh_m%d_8TeV", iMass));
  sigTree1->Add(TString::Format("histograms_CMS-HGG_19032013.root/vbf_m%d_8TeV", iMass));
  sigTree1->Add(TString::Format("histograms_CMS-HGG_19032013.root/wzh_m%d_8TeV", iMass));
  sigTree1->Add(TString::Format("histograms_CMS-HGG_19032013.root/tth_m%d_8TeV", iMass));*/

  std::cout<<sigTree1->GetEntries()<<std::endl;
  //sigTree1->SetTitle("sigTree1");
  //sigTree1->SetName("sigTree1");

  Double_t NevTOT;
  Double_t Nev[NCAT];
  // common preselection cut

 std::string suffix;
  if(Sys==0)suffix = "";
  if(Sys==1)suffix = "_scaleUP";
  if(Sys==2)suffix = "_scaleDN";

  TString mainCut= TString::Format(("PhotonsMass"+suffix+">=130 && PhotonsMass"+suffix+"<=1000").c_str(), mass, mass);   // livia
 
  
  
  // Create signal dataset composed with different productions, the weight is already applied in our ntuples
  RooDataSet sigWeighted("sigWeighted","dataset",sigTree1,*ntplVars,mainCut,"evweight");
  cout << endl;
  cout << "sigWeighted" << endl;
  sigWeighted.Print("v");
  cout << "---- nX:  " << sigWeighted.sumEntries() << endl; 
  NevTOT = sigWeighted.sumEntries();
  
  // apply a common preselection cut; split in categories
  cout << endl;
  RooDataSet* signal[NCAT];
  for (int c=0; c<ncat; ++c) {

  
    
    // 1)  prime 4 cat livia
    if (c==0) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var(("PhotonsMass"+suffix).c_str()),mainCut+TString::Format("&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 )"));
    if (c==1) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var(("PhotonsMass"+suffix).c_str()),mainCut+TString::Format("&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9<0.94 || ph2_r9<0.94 ) "));
    if (c==2) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var(("PhotonsMass"+suffix).c_str()),mainCut+TString::Format("&& (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 && ph2_r9>0.94 )"));
    if (c==3) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var(("PhotonsMass"+suffix).c_str()),mainCut+TString::Format("&& (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9<0.94 || ph2_r9<0.94 ) "));
   


    w->import(*signal[c],Rename(TString::Format(("SigWeight"+suffix+"_cat%d").c_str(),c)));
    
    cout << "cat " << c << ", signal[c]: " << endl;
    signal[c]->Print("v");
    cout << "---- for category " << c << ", nX for signal[c]:  " << signal[c]->sumEntries() << endl; 
    cout << endl;

    Nev[c] = signal[c]->sumEntries();
  }

  // Create full weighted signal data set without categorization
  RooDataSet* signalAll = (RooDataSet*) sigWeighted.reduce(*w->var(("PhotonsMass"+suffix).c_str()),mainCut);
  w->import(*signalAll, Rename(("SigWeight"+suffix).c_str()));
  cout << "now signalAll" << endl;
  signalAll->Print("v");
  cout << "---- nX for signalAll:  " << signalAll->sumEntries() << endl; 
  cout << endl;

  std::cout<<" NevTOT: "<<NevTOT<<std::endl;
  Double_t NevTOTEB = Nev[0]+Nev[1];
  Double_t NevTOTEE = Nev[2]+Nev[3];

  for(int i = 0; i<NCAT; i++) std::cout<<" Nev cat"<<i<<" : "<<Nev[i]<<" Eff cat"<<i<<" : "<<Nev[i]/NevTOT<<std::endl;
  
  std::cout<<"*************R9 Migration*****************"<<std::endl;
  std::cout<<" cat0: "<<Nev[0]/NevTOTEB<<std::endl;
  std::cout<<" cat1: "<<Nev[1]/NevTOTEB<<std::endl;
  std::cout<<" cat2: "<<Nev[2]/NevTOTEE<<std::endl;
  std::cout<<" cat3: "<<Nev[3]/NevTOTEE<<std::endl;

}




RooPlot*  SigModelResponseFcnFit(RooWorkspace* w, Float_t mass, int Sys) {

  
  //chain summing up all production modes
  TChain* sigTree1  = new TChain();
  sigTree1->Add("histograms_CMS-HGG_19032013_WithSyst_M250.root/sigTree1");
  /* sigTree1->Add("histograms_CMS-HGG_19032013.root/ggh_m250_8TeV");
  sigTree1->Add("histograms_CMS-HGG_19032013.root/vbf_m250_8TeV");
  sigTree1->Add("histograms_CMS-HGG_19032013.root/wzh_m250_8TeV");
  sigTree1->Add("histograms_CMS-HGG_19032013.root/tth_m250_8TeV");*/
  sigTree1->SetTitle("sigTree1");
  sigTree1->SetName("sigTree1");
  std::string suffix;
  if(Sys==0)suffix = "";
  if(Sys==1)suffix = "_scaleUP";
  if(Sys==2)suffix = "_scaleDN";


  // Variables
  RooArgSet* ntplVars = defineVariablesM250();
  TString mainCut1 = TString::Format(("PhotonsMass"+suffix+" > 130 && PhotonsMass"+suffix+"<1000").c_str());   // livia
  RooDataSet sigWeighted("sigWeighted","dataset",sigTree1,*ntplVars,mainCut1,"evweight");

   
  RooRealVar *mH = new RooRealVar("MH", "MH", MINmass, MAXmass);
  mH->setVal(mass);
  mH->setConstant();
  w->import(*mH);

  RooFormulaVar *massReduced_formula     = new RooFormulaVar("massReduced_formula","","@0/@1 -1",RooArgList(*w->var(("PhotonsMass"+suffix).c_str()),*w->var("PhotonsMassTrue")));
  RooRealVar* massReduced = (RooRealVar*) sigWeighted.addColumn(*massReduced_formula);
  massReduced->SetName("massReduced");
  massReduced->SetTitle("massReduced");
  w->import(*massReduced);  
  massReduced->setRange(-0.5, 0.5);

 // common preselection cut
  TString mainCut = TString::Format("massReduced>-0.5 && massReduced <0.5");   // livia
  
  
  RooDataSet* signal[NCAT];
  RooCBShape* ResponseCBpos[NCAT];
  RooCBShape* ResponseCBneg[NCAT];
  RooGaussian* ResponseGauss[NCAT];
  RooAddPdf* ResponseAddGauss[NCAT];
  RooAddPdf* ResponseAdd[NCAT];
  
  TCanvas* c1 = new TCanvas("c1", "c1", 1);
  c1->cd();
  TPaveText* label_cms = get_labelCMS(0, "2012", true);
  TPaveText* label_sqrt = get_labelSqrt(0);
  
  for(int c = 1; c<NCAT; c++){
    
    TLatex *lat  = new TLatex(0.65,0.9,TString::Format("Cat: %d", c));  
    lat->SetTextSize(0.038);
    lat->SetTextAlign(11);
    lat->SetTextFont(42); 
    lat->SetNDC();

   // 1)  prime 4 cat livia
   if (c==0) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("massReduced"),mainCut+TString::Format("&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 )"));
   if (c==1) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("massReduced"),mainCut+TString::Format("&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9<0.94 || ph2_r9<0.94 ) "));
   if (c==2) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("massReduced"),mainCut+TString::Format("&& (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 && ph2_r9>0.94 )"));
   if (c==3) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("massReduced"),mainCut+TString::Format("&& (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9<0.94 || ph2_r9<0.94 ) "));

   // w->import(*signal[c],Rename(TString::Format("SigWeight_cat%d",c)));
   //add cb neg +pos

   //cb pos                                                                                                                     
   RooFormulaVar CBpos_mean(TString::Format("ReducedMass_CBpos_sig_mean_cat%d",c),"", "@0", *w->var(TString::Format("ReducedMass_sig_mean_cat%d",c)));
   RooFormulaVar CBpos_sigma(TString::Format("ReducedMass_CBpos_sig_sigma_cat%d",c), "", "@0", *w->var(TString::Format("ReducedMass_sig_sigma_cat%d",c)));
   (*w->var(TString::Format("ReducedMass_sig_sigma_cat%d",c))).setConstant();
   RooFormulaVar CBpos_alphaCB(TString::Format("ReducedMass_CBpos_sig_alphaCB_cat%d",c),"", "@0", *w->var( TString::Format("ReducedMass_sig_alphaCBpos_cat%d",c)));
   (*w->var(TString::Format("ReducedMass_sig_alphaCBpos_cat%d",c))).setConstant();
   RooFormulaVar CBpos_n(TString::Format("ReducedMass_CBpos_sig_n_cat%d",c),"", "@0", *w->var( TString::Format("ReducedMass_sig_Npos_cat%d",c)));
   (*w->var(TString::Format("ReducedMass_sig_Npos_cat%d",c))).setConstant();
   //cb neg
   RooFormulaVar CBneg_n(TString::Format("ReducedMass_CBneg_sig_n_cat%d",c),"", "@0", *w->var( TString::Format("ReducedMass_sig_Nneg_cat%d",c)));
   (*w->var(TString::Format("ReducedMass_sig_Nneg_cat%d",c))).setConstant();
   RooFormulaVar CBneg_alphaCB(TString::Format("ReducedMass_CBneg_sig_alphaCB_cat%d",c),"", "@0", *w->var( TString::Format("ReducedMass_sig_alphaCBneg_cat%d",c)));
   (*w->var(TString::Format("ReducedMass_sig_alphaCBneg_cat%d",c))).setConstant();
   
   ResponseCBpos[c] =  new RooCBShape(TString::Format("ResponseCBpos_cat%d",c),TString::Format("ResponseCBpos_cat%d",c) , *massReduced, CBpos_mean, CBpos_sigma,  CBpos_alphaCB, CBpos_n) ;
   

   ResponseCBneg[c] =  new RooCBShape(TString::Format("ResponseCBneg_cat%d",c),TString::Format("ResponseCBneg_cat%d",c) , *massReduced, CBpos_mean, CBpos_sigma,  CBneg_alphaCB, CBneg_n) ;
   

   
   RooFormulaVar CB_frac(TString::Format("ReducedMass_CBpos_sig_frac_cat%d",c),"", "@0", *w->var(TString::Format("ReducedMass_sig_frac_cat%d",c)));
   (*w->var(TString::Format("ReducedMass_sig_frac_cat%d",c))).setConstant();
   w->import(CB_frac);  
   ResponseAdd[c]= new RooAddPdf(TString::Format("ResponseAddPdf_cat%d",c),TString::Format("ResponseAddPdf_cat%d",c) , RooArgList(*ResponseCBneg[c], *ResponseCBpos[c]), CB_frac);
   w->import(*ResponseAdd[c]);
   
  
  
  
   RooFitResult* fitresults = (RooFitResult* ) ResponseAdd[c]->fitTo(*signal[c],SumW2Error(kTRUE), Range(-1, 1), RooFit::Save(kTRUE));
   std::cout<<TString::Format("******************************** Signal Fit results CB+CB  mass %f cat %d***********************************", mass, c)<<std::endl;
    fitresults->Print("V");
   



    std::cout<<" &"<<(*w->var( TString::Format("ReducedMass_sig_Nneg_cat%d",c))).getVal()<<std::endl;
    std::cout<<"               "<<std::endl;
    std::cout<<" &"<<(*w->var( TString::Format("ReducedMass_sig_Npos_cat%d",c))).getVal()<<std::endl;
    std::cout<<"               "<<std::endl;
    std::cout<<" & "<<(*w->var( TString::Format("ReducedMass_sig_alphaCBneg_cat%d",c))).getVal()<<std::endl;
    std::cout<<"               "<<std::endl;
    std::cout<<" &"<<(*w->var( TString::Format("ReducedMass_sig_alphaCBpos_cat%d",c))).getVal()<<std::endl;
    std::cout<<"               "<<std::endl;
    std::cout<<" &"<<(*w->var( TString::Format("ReducedMass_sig_frac_cat%d",c))).getVal()<<std::endl;
    std::cout<<"               "<<std::endl;
    std::cout<<" &"<<(*w->var( TString::Format("ReducedMass_sig_mean_cat%d",c))).getVal()<<std::endl;
    std::cout<<"               "<<std::endl;
    std::cout<<" &"<<(*w->var( TString::Format("ReducedMass_sig_sigma_cat%d",c))).getVal()<<std::endl;
    std::cout<<"               "<<std::endl;

    RooPlot* plotReturn = massReduced->frame(Range(-0.12, 0.12),Title("Mass Reduced"), Bins(60));
    if(Sys==0)  signal[c]->plotOn(plotReturn, MarkerColor(kGray+2), LineColor(kGray+2));
    if(Sys==1)  signal[c]->plotOn(plotReturn, MarkerColor(kBlue-2), LineColor(kBlue-2));
    if(Sys==2)  signal[c]->plotOn(plotReturn, MarkerColor(kRed-2), LineColor(kRed-2));

    if(Sys==0)  ResponseAdd[c]->plotOn(plotReturn, LineColor(kBlack));
    if(Sys==1)  ResponseAdd[c]->plotOn(plotReturn, LineColor(kBlue));
    if(Sys==2)  ResponseAdd[c]->plotOn(plotReturn, LineColor(kRed));


    RooPlot* plotG = massReduced->frame(Range(-0.12, 0.12),Title("Mass Reduced"), Bins(60));
    signal[c]->plotOn(plotG);
   

    ResponseAdd[c]->plotOn(plotG, LineColor(kBlue));
    ResponseAdd[c]->plotOn(plotG,Components(TString::Format("ResponseCBneg_cat%d",c)), LineColor(kGreen), LineStyle(kDashed));
    ResponseAdd[c]->plotOn(plotG,Components(TString::Format("ResponseCBpos_cat%d",c)), LineColor(kRed), LineStyle(kDashed));
  
    plotG->GetYaxis()->SetRangeUser(0.0001,plotG->GetMaximum()*10 );
    plotG->GetXaxis()->SetTitle("#Delta m");
    plotG->GetXaxis()->SetTitleFont(42);
    plotG->GetXaxis()->SetTitleSize(0.05);
    TLegend* legmc = new TLegend(0.6, 0.6, 0.85, 0.89, "", "brNDC");
    legmc->SetTextSize(0.0206044);  
    legmc->SetTextFont(42);
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    
    legmc->AddEntry(plotG->getObject(0),"m_{#gamma#gamma} = 250 GeV","LPE");    
    

    legmc->AddEntry(plotG->getObject(1),"Sum of two CB ","L");
    legmc->AddEntry(plotG->getObject(2),"CB 1","L");   
    legmc->AddEntry(plotG->getObject(3),"CB 2","L");   
    plotG->Draw();
    
    lat->Draw("same");
    label_cms->Draw("same");
    label_sqrt->Draw("same"); 
    legmc->Draw("same");
    
    
    c1->SetLogy();
    
    c1->SaveAs(TString::Format("plots/responseFcnFitCBCB_cat%d_LOG.png",c)); 
    c1->SaveAs(TString::Format("plots/responseFcnFitCBCB_cat%d_LOG.pdf",c)); 
    c1->SaveAs(TString::Format("plots/responseFcnFitCBCB_cat%d_LOG.eps",c)); 
    
    
    plotG->GetYaxis()->SetRangeUser(0.0001,plotG->GetMaximum()*0.12 );
    lat->Draw("same");
    label_cms->Draw("same");
    label_sqrt->Draw("same"); 
    legmc->Draw("same");
    c1->SetLogy(0);  

    c1->SaveAs(TString::Format("plots/responseFcnFitCBCB_cat%d.png",c)); 
    c1->SaveAs(TString::Format("plots/responseFcnFitCBCB_cat%d.pdf",c)); 
    

  }
  

  return plotReturn;

}





RooPlot*  SigModelResponseFcnFitRES(RooWorkspace* w, Float_t mass, int Sys) {

  
  //chain summing up all production modes
  TChain* sigTree1  = new TChain();
  // sigTree1->Add("histograms_CMS-HGG_19032013_WithSyst_M250.root/sigTree1");

  if(Sys==0){
  sigTree1->Add("histograms_CMS-HGG_NOSYST.root/ggh_m250_8TeV");
  sigTree1->Add("histograms_CMS-HGG_NOSYST.root/vbf_m250_8TeV");
  sigTree1->Add("histograms_CMS-HGG_NOSYST.root/wzh_m250_8TeV");
  sigTree1->Add("histograms_CMS-HGG_NOSYST.root/tth_m250_8TeV");
  }else if (Sys==1){
  sigTree1->Add("histograms_CMS-HGG_E_RES_UP.root/ggh_m250_8TeV");
  sigTree1->Add("histograms_CMS-HGG_E_RES_UP.root/vbf_m250_8TeV");
  sigTree1->Add("histograms_CMS-HGG_E_RES_UP.root/wzh_m250_8TeV");
  sigTree1->Add("histograms_CMS-HGG_E_RES_UP.root/tth_m250_8TeV");
  }else if (Sys==2){
  sigTree1->Add("histograms_CMS-HGG_E_RES_DN.root/ggh_m250_8TeV");
  sigTree1->Add("histograms_CMS-HGG_E_RES_DN.root/vbf_m250_8TeV");
  sigTree1->Add("histograms_CMS-HGG_E_RES_DN.root/wzh_m250_8TeV");
  sigTree1->Add("histograms_CMS-HGG_E_RES_DN.root/tth_m250_8TeV");
  }

  sigTree1->SetTitle("sigTree1");
  sigTree1->SetName("sigTree1");
  std::string suffix;
  if(Sys==0)suffix = "";
  if(Sys==1)suffix = "_smearUP";
  if(Sys==2)suffix = "_smearDN";


  // Variables
  RooArgSet* ntplVars = defineVariablesM250();
  TString mainCut1 = TString::Format("PhotonsMass > 130 && PhotonsMass<1000");   // livia
  RooDataSet sigWeighted("sigWeighted","dataset",sigTree1,*ntplVars,mainCut1,"evweight");

   
  RooRealVar *mH = new RooRealVar("MH", "MH", MINmass, MAXmass);
  mH->setVal(mass);
  mH->setConstant();
  w->import(*mH);

  // RooFormulaVar *massReduced_formula     = new RooFormulaVar("massReduced_formula","","@0/@1-1",RooArgList(*w->var("PhotonsMass"),*w->var("PhotonsMassTrue")));
  RooFormulaVar *massReduced_formula     = new RooFormulaVar("massReduced_formula","","@0/@1-1",RooArgList(*w->var("PhotonsMass"),*w->var("PhotonsMassTrue")));
  RooRealVar* massReduced = (RooRealVar*) sigWeighted.addColumn(*massReduced_formula);
  massReduced->SetName("massReduced");
  massReduced->SetTitle("massReduced");
  w->import(*massReduced);  
  massReduced->setRange(-0.5, 0.5);

 // common preselection cut
  TString mainCut = TString::Format("massReduced>-0.5 && massReduced <0.5");   // livia
 
  
  
  RooDataSet* signal[NCAT];
  RooCBShape* ResponseCBpos[NCAT];
  RooCBShape* ResponseCBneg[NCAT];
  RooGaussian* ResponseGauss[NCAT];
  RooAddPdf* ResponseAddGauss[NCAT];
  RooAddPdf* ResponseAdd[NCAT];
  
  TCanvas* c1 = new TCanvas("c1", "c1", 1);
  c1->cd();
  TPaveText* label_cms = get_labelCMS(0, "2012", true);
  TPaveText* label_sqrt = get_labelSqrt(0);
  
  for(int c = 3; c<NCAT; c++){
    
    TLatex *lat  = new TLatex(0.65,0.9,TString::Format("Cat: %d", c));  
    lat->SetTextSize(0.038);
    lat->SetTextAlign(11);
    lat->SetTextFont(42); 
    lat->SetNDC();

   // 1)  prime 4 cat livia
   if (c==0) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("massReduced"),mainCut+TString::Format("&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 )"));
   if (c==1) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("massReduced"),mainCut+TString::Format("&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9<0.94 || ph2_r9<0.94 ) "));
   if (c==2) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("massReduced"),mainCut+TString::Format("&& (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 && ph2_r9>0.94 )"));
   if (c==3) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("massReduced"),mainCut+TString::Format("&& (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9<0.94 || ph2_r9<0.94 ) "));

   // w->import(*signal[c],Rename(TString::Format("SigWeight_cat%d",c)));
   //add cb neg +pos

   //cb pos                                                                                                                     
   RooFormulaVar CBpos_mean(TString::Format("ReducedMass_CBpos_sig_mean_cat%d",c),"", "@0", *w->var(TString::Format("ReducedMass_sig_mean_cat%d",c)));
   (*w->var(TString::Format("ReducedMass_sig_mean_cat%d",c))).setConstant();
   RooFormulaVar CBpos_sigma(TString::Format("ReducedMass_CBpos_sig_sigma_cat%d",c), "", "@0", *w->var(TString::Format("ReducedMass_sig_sigma_cat%d",c)));
   //
   RooFormulaVar CBpos_alphaCB(TString::Format("ReducedMass_CBpos_sig_alphaCB_cat%d",c),"", "@0", *w->var( TString::Format("ReducedMass_sig_alphaCBpos_cat%d",c)));
   // (*w->var(TString::Format("ReducedMass_sig_alphaCBpos_cat%d",c))).setConstant();
   RooFormulaVar CBpos_n(TString::Format("ReducedMass_CBpos_sig_n_cat%d",c),"", "@0", *w->var( TString::Format("ReducedMass_sig_Npos_cat%d",c)));
   (*w->var(TString::Format("ReducedMass_sig_Npos_cat%d",c))).setConstant();
   //cb neg
   RooFormulaVar CBneg_n(TString::Format("ReducedMass_CBneg_sig_n_cat%d",c),"", "@0", *w->var( TString::Format("ReducedMass_sig_Nneg_cat%d",c)));
   (*w->var(TString::Format("ReducedMass_sig_Nneg_cat%d",c))).setConstant();
   RooFormulaVar CBneg_alphaCB(TString::Format("ReducedMass_CBneg_sig_alphaCB_cat%d",c),"", "@0", *w->var( TString::Format("ReducedMass_sig_alphaCBneg_cat%d",c)));
   // (*w->var(TString::Format("ReducedMass_sig_alphaCBneg_cat%d",c))).setConstant();
   
   ResponseCBpos[c] =  new RooCBShape(TString::Format("ResponseCBpos_cat%d",c),TString::Format("ResponseCBpos_cat%d",c) , *massReduced, CBpos_mean, CBpos_sigma,  CBpos_alphaCB, CBpos_n) ;
   

   ResponseCBneg[c] =  new RooCBShape(TString::Format("ResponseCBneg_cat%d",c),TString::Format("ResponseCBneg_cat%d",c) , *massReduced, CBpos_mean, CBpos_sigma,  CBneg_alphaCB, CBneg_n) ;
   

   
   RooFormulaVar CB_frac(TString::Format("ReducedMass_CBpos_sig_frac_cat%d",c),"", "@0", *w->var(TString::Format("ReducedMass_sig_frac_cat%d",c)));
   (*w->var(TString::Format("ReducedMass_sig_frac_cat%d",c))).setConstant();
   w->import(CB_frac);  
   ResponseAdd[c]= new RooAddPdf(TString::Format("ResponseAddPdf_cat%d",c),TString::Format("ResponseAddPdf_cat%d",c) , RooArgList(*ResponseCBneg[c], *ResponseCBpos[c]), CB_frac);
   w->import(*ResponseAdd[c]);
   
  
  
  
   RooFitResult* fitresults = (RooFitResult* ) ResponseAdd[c]->fitTo(*signal[c],SumW2Error(kTRUE), Range(-1, 1), RooFit::Save(kTRUE));
   std::cout<<TString::Format("******************************** Signal Fit results CB+CB  mass %f cat %d***********************************", mass, c)<<std::endl;
    fitresults->Print("V");
   



    std::cout<<" &"<<(*w->var( TString::Format("ReducedMass_sig_Nneg_cat%d",c))).getVal()<<std::endl;
    std::cout<<"               "<<std::endl;
    std::cout<<" &"<<(*w->var( TString::Format("ReducedMass_sig_Npos_cat%d",c))).getVal()<<std::endl;
    std::cout<<"               "<<std::endl;
    std::cout<<" & "<<(*w->var( TString::Format("ReducedMass_sig_alphaCBneg_cat%d",c))).getVal()<<std::endl;
    std::cout<<"               "<<std::endl;
    std::cout<<" &"<<(*w->var( TString::Format("ReducedMass_sig_alphaCBpos_cat%d",c))).getVal()<<std::endl;
    std::cout<<"               "<<std::endl;
    std::cout<<" &"<<(*w->var( TString::Format("ReducedMass_sig_frac_cat%d",c))).getVal()<<std::endl;
    std::cout<<"               "<<std::endl;
    std::cout<<" &"<<(*w->var( TString::Format("ReducedMass_sig_mean_cat%d",c))).getVal()<<std::endl;
    std::cout<<"               "<<std::endl;
    std::cout<<" &"<<(*w->var( TString::Format("ReducedMass_sig_sigma_cat%d",c))).getVal()<<std::endl;
    std::cout<<"               "<<std::endl;

    RooPlot* plotReturn = massReduced->frame(Range(-0.12, 0.12),Title("Mass Reduced"), Bins(60));
    if(Sys==0)  signal[c]->plotOn(plotReturn, MarkerColor(kGray+2), LineColor(kGray+2));
    if(Sys==1)  signal[c]->plotOn(plotReturn, MarkerColor(kBlue-2), LineColor(kBlue-2));
    if(Sys==2)  signal[c]->plotOn(plotReturn, MarkerColor(kRed-2), LineColor(kRed-2));

    if(Sys==0)  ResponseAdd[c]->plotOn(plotReturn, LineColor(kBlack));
    if(Sys==1)  ResponseAdd[c]->plotOn(plotReturn, LineColor(kBlue));
    if(Sys==2)  ResponseAdd[c]->plotOn(plotReturn, LineColor(kRed));


    RooPlot* plotG = massReduced->frame(Range(-0.12, 0.12),Title("Mass Reduced"), Bins(60));
    signal[c]->plotOn(plotG);
   

    ResponseAdd[c]->plotOn(plotG, LineColor(kBlue));
    ResponseAdd[c]->plotOn(plotG,Components(TString::Format("ResponseCBneg_cat%d",c)), LineColor(kGreen), LineStyle(kDashed));
    ResponseAdd[c]->plotOn(plotG,Components(TString::Format("ResponseCBpos_cat%d",c)), LineColor(kRed), LineStyle(kDashed));
  
    plotG->GetYaxis()->SetRangeUser(0.0001,plotG->GetMaximum()*10 );
    plotG->GetXaxis()->SetTitle("#Delta m");
    plotG->GetXaxis()->SetTitleFont(42);
    plotG->GetXaxis()->SetTitleSize(0.05);
    TLegend* legmc = new TLegend(0.6, 0.6, 0.85, 0.89, "", "brNDC");
    legmc->SetTextSize(0.0206044);  
    legmc->SetTextFont(42);
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    
    legmc->AddEntry(plotG->getObject(0),"m_{#gamma#gamma} = 250 GeV","LPE");    
    

    legmc->AddEntry(plotG->getObject(1),"Sum of two CB ","L");
    legmc->AddEntry(plotG->getObject(2),"CB 1","L");   
    legmc->AddEntry(plotG->getObject(3),"CB 2","L");   
    plotG->Draw();
    
    lat->Draw("same");
    label_cms->Draw("same");
    label_sqrt->Draw("same"); 
    legmc->Draw("same");
    
    
    c1->SetLogy();
    
    c1->SaveAs(TString::Format("plots/responseFcnFitCBCB_cat%d_LOG.png",c)); 
    c1->SaveAs(TString::Format("plots/responseFcnFitCBCB_cat%d_LOG.pdf",c)); 
    c1->SaveAs(TString::Format("plots/responseFcnFitCBCB_cat%d_LOG.eps",c)); 
    
    
    plotG->GetYaxis()->SetRangeUser(0.0001,plotG->GetMaximum()*0.12 );
    lat->Draw("same");
    label_cms->Draw("same");
    label_sqrt->Draw("same"); 
    legmc->Draw("same");
    c1->SetLogy(0);  

    c1->SaveAs(TString::Format("plots/responseFcnFitCBCB_cat%d.png",c)); 
    c1->SaveAs(TString::Format("plots/responseFcnFitCBCB_cat%d.pdf",c)); 
    
  
  }
  

  return plotReturn;

}





// Fit signal with model gauss pdfs
void SigModelFitConvBW(RooWorkspace* w, Float_t mass, Double_t width) {

  Int_t ncat = NCAT;
  RooDataSet* sigToFit[NCAT];
  
  
  Float_t MASS(mass);  
  Float_t minMassFit(mass*0.8);
  Float_t maxMassFit(mass*1.2); 
  RooRealVar* PhotonsMass = w->var("PhotonsMass"); 

  std::cout<<"----------------------------------------------------------------------------------------"<<std::endl;
  //  PhotonsMass->setRange("sigrange",minMassFit-20,maxMassFit+20); 
  TPaveText* label_cms = get_labelCMS(0, "2012", true);
  TPaveText* label_sqrt = get_labelSqrt(0);
 
  // Fit Signal 
  for (int c=0; c<ncat; ++c) {
    cout << "---------- Category = " << c << endl;
    //  sigToFit[c] = (RooDataSet*) w->data(TString::Format("SigWeight_cat%d",c));
    //  w->import(*sigToFit[c]);
 
    //cb
    RooFormulaVar CBpos_mean_draw(TString::Format("CBpos_mean_draw_cat%d",c),"","@0+@1",RooArgList(*w->var(TString::Format("ReducedMass_sig_mean_cat%d",c)),*w->var("MH")));
    RooFormulaVar CBpos_mean(TString::Format("CBpos_mean_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_mean_cat%d",c)));
    RooFormulaVar CBpos_sigma(TString::Format("CBpos_sigma_cat%d",c),"","@0*@1",RooArgList(*w->var(TString::Format("ReducedMass_sig_sigma_cat%d",c)),*w->var("MH")) );
    RooFormulaVar CBpos_alphaCB(TString::Format("CBpos_alphaCB_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_alphaCBpos_cat%d",c)) );
    RooFormulaVar CBneg_alphaCB(TString::Format("CBneg_alphaCB_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_alphaCBneg_cat%d",c)) );
    RooFormulaVar CBpos_n(TString::Format("CBpos_n_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_Npos_cat%d",c)) );
    RooFormulaVar CBneg_n(TString::Format("CBneg_n_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_Nneg_cat%d",c)) );
    RooFormulaVar CBpos_frac(TString::Format("CBpos_frac_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_frac_cat%d",c)) );
    
    RooCBShape ResCBpos_draw(TString::Format("ResCBpos_draw_cat%d",c),TString::Format("ResCBpos_draw_cat%d",c) , *PhotonsMass, CBpos_mean_draw, CBpos_sigma, CBpos_alphaCB, CBpos_n) ;
    RooCBShape ResCBneg_draw(TString::Format("ResCBneg_draw_cat%d",c),TString::Format("ResCBneg_draw_cat%d",c) , *PhotonsMass, CBpos_mean_draw, CBpos_sigma, CBneg_alphaCB, CBneg_n) ;
    RooCBShape ResCBpos(TString::Format("ResCBpos_cat%d",c),TString::Format("ResCBpos_cat%d",c) , *PhotonsMass, CBpos_mean, CBpos_sigma,CBpos_alphaCB, CBpos_n) ;
    RooCBShape ResCBneg(TString::Format("ResCBneg_cat%d",c),TString::Format("ResCBneg_cat%d",c) , *PhotonsMass, CBpos_mean, CBpos_sigma,CBneg_alphaCB, CBneg_n) ;
    PhotonsMass->setBins(40000, "cache");  
    //add CB pos + CB neg
    RooAddPdf ResAddPdf(TString::Format("ResAddPdf_cat%d",c),TString::Format("ResAddPdf_cat%d",c) , RooArgList(ResCBneg, ResCBpos), CBpos_frac);
    RooAddPdf ResAddPdf_draw(TString::Format("ResAddPdf_draw_cat%d",c),TString::Format("ResAddPdf_draw_cat%d",c) , RooArgList(ResCBneg_draw, ResCBpos_draw), CBpos_frac);


    
    //BW
    RooFormulaVar meanBW(TString::Format("meanBW_cat%d",c),"","@0",*w->var("MH"));  
    RooRealVar sigmaBW_var(TString::Format("sigmaBW_var_cat%d",c), TString::Format("sigmaBW_var_cat%d",c), width);
    std::cout<<" width:--------> "<<width<<std::endl;
    sigmaBW_var.setConstant();
    w->import(sigmaBW_var);
    
    
    RooFormulaVar* sigmaBW;
    if(width<1)sigmaBW = new RooFormulaVar(TString::Format("sigmaBW_cat%d",c),"","@0",*w->var(TString::Format("sigmaBW_var_cat%d",c))); 
    else if(width==2) sigmaBW = new RooFormulaVar(TString::Format("sigmaBW_cat%d",c),"","@0*0.02",*w->var("MH"));   
    else if(width==5) sigmaBW = new RooFormulaVar(TString::Format("sigmaBW_cat%d",c),"","@0*0.05",*w->var("MH"));   
    else if(width==7) sigmaBW = new RooFormulaVar(TString::Format("sigmaBW_cat%d",c),"","@0*0.07",*w->var("MH"));   
    else if(width==10) sigmaBW = new RooFormulaVar(TString::Format("sigmaBW_cat%d",c),"","@0*0.10",*w->var("MH"));   
    RooBreitWigner SigModelBW(TString::Format("SigModelBW_cat%d",c),TString::Format("SigModelBW_cat%d",c), *PhotonsMass, meanBW, *sigmaBW);

  
    //CONV 
    RooFFTConvPdf  ConvolutedRes_CB(TString::Format("PhotonsMassSig_cat%d",c),TString::Format("PhotonsMassSig_cat%d",c), *PhotonsMass,SigModelBW, ResAddPdf);
    //RooFitResult* fitresults_CB = (RooFitResult* ) ConvolutedRes_CB.fitTo(*sigToFit[c],RooFit::Range("sigrange"), RooFit::Save(kTRUE));
    // std::cout<<TString::Format("******************************** Signal Fit results CB mass %f cat %d***********************************", mass, c)<<std::endl;
    // fitresults_CB->Print("V");
    w->import(ConvolutedRes_CB);
    // std::cout<<".............> "<<c<<std::endl;
    
    RooHistFunc* rooFunc_norm = getRooHistFunc(c,w->var("MH") );
    w->import(*rooFunc_norm);

    std::cout<<"SIG NORM ----->"<<rooFunc_norm->getVal(*w->var("MH"));
    // w->Print("V");

    if(width < 0.&& mass <400 && mass!=350.){ //if i want to plot the fit
      sigToFit[c] = (RooDataSet*) w->data(TString::Format("SigWeight_cat%d",c));

      //  RooFitResult* fitresults_CB = (RooFitResult* ) ConvolutedRes_CB.fitTo(*sigToFit[c], RooFit::Save(kTRUE));
      //  fitresults_CB->Print("V");
      
      RooPlot* plotOnlyResPdf = w->var("PhotonsMass")->frame(Range(minMassFit-20,maxMassFit+20),Bins(100));
      sigToFit[c]->plotOn(plotOnlyResPdf, LineColor(kRed), LineStyle(kDashed));
      double max = plotOnlyResPdf->GetMaximum();
      
      // Plot to verify everything is ok
      RooPlot* plotPhotonsMassAll = w->var("PhotonsMass")->frame(Range(minMassFit-20,maxMassFit+20),Bins(100));
      sigToFit[c]->plotOn(plotPhotonsMassAll);
      SigModelBW.plotOn(plotPhotonsMassAll, LineColor(kGreen), LineStyle(kDashed));
      //  ResAddPdf_draw.plotOn(plotPhotonsMassAll, LineColor(kRed), LineStyle(kDashed));
      ConvolutedRes_CB.plotOn(plotPhotonsMassAll, LineColor(kBlue));

      
      TCanvas* c1 = new TCanvas("c1","PhotonsMass",0,0,800,800);
      c1->cd(1);
      
      plotPhotonsMassAll->Draw();  
      plotPhotonsMassAll->GetYaxis()->SetRangeUser(0.0001, max*1.2);
      TLegend *legmc = new TLegend(0.5491457,0.75,0.801457,0.9340659, TString::Format("Category %d",c), "brNDC");
      legmc->AddEntry(plotPhotonsMassAll->getObject(0),"Simulation","LPE");
      legmc->AddEntry(plotPhotonsMassAll->getObject(1),"BW","L");
      //legmc->AddEntry(plotPhotonsMassAll->getObject(2)," CB + CB ","L");
      legmc->AddEntry(plotPhotonsMassAll->getObject(2),"BW #otimes Resolution","L");
      legmc->SetTextSize(0.0206044);
      legmc->SetTextFont(42);
      legmc->SetBorderSize(0);
      legmc->SetFillStyle(0);
      legmc->Draw();
      
      label_cms->Draw("same");
      label_sqrt->Draw("same");
      

      
      int massI(mass);
      c1->SetLogy();
     
      plotPhotonsMassAll->GetXaxis()->SetTitle("m_{#gamma #gamma}[GeV]");
     
      c1->SetLogy(0);
      c1->SaveAs("plots/prelimSignalCBCconvBW"+TString::Format("_M%d_cat%d.png",massI, c));
      c1->SaveAs("plots/prelimSignalCBCconvBW"+TString::Format("_M%d_cat%d.pdf",massI, c));
      
      c1->SetLogy();
      plotPhotonsMassAll->GetYaxis()->SetRangeUser(0.0001,max*10. );
     
      c1->SaveAs("plots/prelimSignalCBCconvBW"+TString::Format("_M%d_cat%d_LOG.png",massI,c));
      c1->SaveAs("plots/prelimSignalCBCconvBW"+TString::Format("_M%d_cat%d_LOG.pdf",massI,c));
      c1->SaveAs("plots/prelimSignalCBCconvBW"+TString::Format("_M%d_cat%d_LOG.eps",massI,c));


      /*      
      //plot signal model at different widths
      bool plotW = false;
      if(plotW && c==0){
	RooRealVar var_01("var_w01", "var_w01", 0.1);
	var_01.setConstant();	
	RooFormulaVar sigmaBW_01("w01", "w01","@0", var_01);     
	RooBreitWigner SiBW_01("sigBW_01","sigBW_01" , *PhotonsMass, meanBW, sigmaBW_01);
	RooFFTConvPdf  ConvolutedRes_01("conv01", "conv01", *PhotonsMass,SiBW_01, ResAddPdf);

	RooRealVar var_3("var_w3", "var_w3",3);
	var_3.setConstant();	
	RooFormulaVar sigmaBW_3("w3", "w3","@0",  var_3);     
	RooBreitWigner SiBW_3("sigBW_3","sigBW_3" , *PhotonsMass, meanBW, sigmaBW_3);
	RooFFTConvPdf  ConvolutedRes_3("conv3", "conv3", *PhotonsMass,SiBW_3, ResAddPdf);

	RooRealVar var_6("var_w6", "var_w6", 6);
	var_6.setConstant();	
	RooFormulaVar sigmaBW_6("w6", "w6","@0", var_6);     
	RooBreitWigner SiBW_6("sigBW_6","sigBW_6" , *PhotonsMass, meanBW, sigmaBW_6);
	RooFFTConvPdf  ConvolutedRes_6("conv6", "conv6", *PhotonsMass,SiBW_6, ResAddPdf);

	RooRealVar var_10("var_w10", "var_w10", 10);
	var_10.setConstant();	
	RooFormulaVar sigmaBW_10("w10", "w10","@0", var_10);     
	RooBreitWigner SiBW_10("sigBW_10","sigBW_10" , *PhotonsMass, meanBW, sigmaBW_10);
	RooFFTConvPdf  ConvolutedRes_10("conv10", "conv10", *PhotonsMass,SiBW_10, ResAddPdf);

	RooRealVar var_15("var_w15", "var_w15",15);
	var_15.setConstant();	
	RooFormulaVar sigmaBW_15("w15", "w15","@0",  var_15);     
	RooBreitWigner SiBW_15("sigBW_15","sigBW_15" , *PhotonsMass, meanBW, sigmaBW_15);
	RooFFTConvPdf  ConvolutedRes_15("conv15", "conv15", *PhotonsMass,SiBW_15, ResAddPdf);

	RooPlot* plotWidths = w->var("PhotonsMass")->frame(Range(minMassFit-20,maxMassFit+20),Bins(100));
	ConvolutedRes_15.plotOn( plotWidths, LineColor(kAzure+3));
	ConvolutedRes_10.plotOn( plotWidths, LineColor(kAzure+2));
	ConvolutedRes_6.plotOn( plotWidths, LineColor(kAzure+1));
	ConvolutedRes_3.plotOn( plotWidths, LineColor(kViolet+1));
	ConvolutedRes_01.plotOn( plotWidths, LineColor(kViolet-9));
	plotWidths->Draw();

	label_cms->Draw("same");
	label_sqrt->Draw("same");
      
	TLegend* leg = new TLegend(0.598851,0.6044755,0.84253,0.928252,"", "brNDC");
  
	leg->SetBorderSize(0.);
	leg->SetFillColor(kWhite);
	leg->SetTextFont(42);
	plotWidths->GetYaxis()->SetRangeUser(0.001, 1.);
	plotWidths->GetXaxis()->SetTitle("m_{#gamma #gamma}[GeV]");
	plotWidths->GetYaxis()->SetTitle(" ");
	leg->AddEntry(plotWidths->getObject(0), "Width = 15 GeV", "L");
	leg->AddEntry(plotWidths->getObject(1), "Width = 10 GeV", "L");
	leg->AddEntry(plotWidths->getObject(2),"Width = 6 GeV", "L");
	leg->AddEntry(plotWidths->getObject(3),"Width = 3 GeV", "L");
	leg->AddEntry(plotWidths->getObject(4), "Width = 0.1 GeV", "L");
	leg->Draw("same");

	c1->SaveAs("plots/SignalModels_differentWidths.png");
	c1->SaveAs("plots/SignalModels_differentWidths.pdf");
      
      }


      */

    }
    
    // IMPORTANT: fix all pdf parameters to constant
    
    w->defineSet(TString::Format("ConvolutedPdfParam_cat%d",c),RooArgSet( *w->var(TString::Format("ReducedMass_sig_sigma_cat%d",c)), 
									  *w->var(TString::Format("ReducedMass_sig_alphaCBpos_cat%d",c)),
									  *w->var(TString::Format("ReducedMass_sig_alphaCBneg_cat%d",c)),
									  *w->var(TString::Format("ReducedMass_sig_Npos_cat%d",c)),
									  *w->var(TString::Format("ReducedMass_sig_Nneg_cat%d",c)),	   
									  *w->var(TString::Format("ReducedMass_sig_frac_cat%d",c)),  
									  *w->var(TString::Format("ReducedMass_sig_mean_cat%d",c)),
									  
									  *w->var(TString::Format("sigmaBW_var_cat%d",c))));
    
    SetConstantParams(w->set(TString::Format("ConvolutedPdfParam_cat%d",c)));
    
    //w->Print("V");
    
  }

}

