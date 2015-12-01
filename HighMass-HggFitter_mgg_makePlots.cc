
 
 /** \macro H2GGFitter.cc
 *
 * The analysis root trees produced in a simple format 
 *
 *     TFile file(filename,"RECREATE", "X->jj input tree for unbinned maximum-likelihood fit");
 *     TTree* outTree  = new TTree("XTojj","X->jj input tree for unbinned maximum-likelihood fit");
 *     Float_t mass;
 *     Int_t CAT3;
 *     Float_t weight;
 *
 *     outTree->Branch("mass",&mass,"mass/F");
 *     outTree->Branch("weight",&weight,"weight/F");
 *     outTree->Branch("CAT4",&CAT4,"CAT4/I");
 *     {
 *       .............
 *       outTree->Fill();
 *     }
 *
 *     file.Write();
 *     file.Close();
 *     delete outTree;
 *
 * are used as input files. They have to be produced for 
 * data and Monte Carlo signal and background data sets 
 * after all analysis selections to be applied.  *
 */
// Loading:  .L HH2ggbbFitter.cc
// Running:  runfits("hgg120-shapes-combined-Unbinned.root")  
//                
#include "CMS_lumi.C"
using namespace RooFit;
using namespace RooStats ;

static const Int_t NCAT = 4;  // chiara

std::string filePOSTfix="";
double signalScaler=1.00;

void AddSigData(RooWorkspace*, Float_t);
void AddBkgData(RooWorkspace*, Float_t);
void MakeRooKeysPDFMCBkg(RooWorkspace*, Float_t);
void SigModelFitGauss(RooWorkspace*, Float_t);
void SigModelFitCBC(RooWorkspace*, Float_t);
RooFitResult*  BkgModelFitBernstein(RooWorkspace*, Bool_t, Float_t, bool);
RooFitResult*  BkgModelFitExpo(RooWorkspace*, Bool_t, Float_t, bool);
RooFitResult*  BkgModelFitDiJetFunc(RooWorkspace*, Bool_t, Float_t, bool);
RooFitResult*  BkgModelFitExpolFunc(RooWorkspace*, Bool_t, Float_t, bool);
void MakePlots(RooWorkspace*, Float_t, RooFitResult*, bool);
void MakeSigWS(RooWorkspace* w, const char* filename);
void MakeBkgWS(RooWorkspace* w, const char* filename);
void SetConstantParams(const RooArgSet* params);
void MakeParameterTrendvsMassGauss();
void MakeParameterTrendvsMassCBC();
void MakePlotMassDataMC(RooWorkspace*, Float_t);
void MakePlotNVTXDataMC(RooWorkspace*, Float_t);



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
  } else if( legendQuadrant==2 ) {
    x1 =  0.25;
    y1 = 0.83;
    x2 =  0.42;
    y2 = 0.87;
  } else if( legendQuadrant==3 ) {
    x1 = 0.25;
    y1 = 0.2;
    x2 = 0.42;
    y2 = 0.24;
  } else if( legendQuadrant==0 ) {
    x1 = 0.175;
    y1 = 0.953;
    x2 = 0.6;
    y2 = 0.975;
  }

  
  TPaveText* cmslabel = new TPaveText( x1, y1, x2, y2, "brNDC" );
  cmslabel->SetFillColor(kWhite);
  cmslabel->SetTextSize(0.038);
  if( legendQuadrant==0 ) cmslabel->SetTextAlign(11);
  cmslabel->SetTextSize(0.038);
  cmslabel->SetTextFont(42);
 
    std::string leftText;
   
     
    if (sim)  leftText = "CMS Simulation"; //cwr ->remove 2011
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
  label_sqrt->SetTextAlign(31); // align right
  label_sqrt->AddText("#sqrt{s} = 8 TeV");
  return label_sqrt;

}








RooArgSet* defineVariables() {

  // define variables of the input ntuple //livia
  RooRealVar* PhotonsMass  = new RooRealVar("PhotonsMass", "M(gg)",100, 2000,"GeV");
  RooRealVar* PhotonsMassTrue  = new RooRealVar("PhotonsMassTrue", "M(gg)",100, 2000,"GeV");
  RooRealVar* ph1_eta = new RooRealVar("ph1_eta", "eta(g1)",-10,10,"");
  RooRealVar* ph2_eta = new RooRealVar("ph2_eta", "eta(g2)",-10,10,"");
  RooRealVar* ph1_r9 = new RooRealVar("ph1_r9", "R9(g1)",-10,10,"");
  RooRealVar* ph2_r9 = new RooRealVar("ph2_r9", "R9(g2)",-10,10,"");
  RooRealVar* evweight = new RooRealVar("evweight","Reweightings",0,10,"");
  RooRealVar* nvtx = new RooRealVar("nvtx", "nvtx", 0, 50, "");
  //  RooRealVar* btagCategory = new RooRealVar("btagCategory","event category",0.9,2.1,"") ;
  
  RooArgSet* ntplVars = new RooArgSet(*PhotonsMass, *PhotonsMassTrue, *ph1_eta, *ph2_eta, *ph1_r9, *ph2_r9,  *evweight, *nvtx);
  
  return ntplVars;
}

void runfits(const Float_t mass=150, Bool_t dobands = false) {

  //******************************************************************//
  //  Running mode  corresponds to the following cases
  //         - full run set:
  //         - create signal and background data sets 
  //         - make and fit signal and background  models 
  //         - write signal and background workspaces in root files
  //         - write data card
  //*******************************************************************//
  cout << "Now plot Data vs MC bkg" << endl;
  TString fileBaeName(TString::Format("HighMass-hgg.m%.1f", mass));    
  TString fileBkgName(TString::Format("HighMass-hgg.inputbkg_8TeV", mass));
  
  TString card_name("HighMass-hgg_models_Bkg_8TeV_test.rs");
  HLFactory hlf("HLFactory", card_name, false);
  RooWorkspace* w = hlf.GetWs();
 
  Double_t MMIN = 100.; //130 //livia   //100  //180   //270
  Double_t MMAX = 800.; //450      //livia //200 // 300  //2000
  w->var("PhotonsMass")->setMin(MMIN);
  w->var("PhotonsMass")->setMax(MMAX);

 
  cout << "Now plot Data vs MC bkg" << endl;
  
  //MakePlotNVTXDataMC(w, mass);
  // PlotVtxRecoEff();

  // Make plots for data and fit results
  cout << endl; cout << "Preparing final plots" << endl;
  //MakePlots(w, mass, fitresults_bern, blind);
  
  //compare MC private with Moriond13 NON RD
  /* PlotSigShapeMCofficial(w, 0); 
   PlotSigShapeMCofficial(w, 1); 
  PlotSigShapeMCofficial(w, 2); 
  PlotSigShapeMCofficial(w, 3); 
  */
  // MakePlotNVTXDataMC(w, mass); //MAKE FIG 2
  MakePlotMassDataMC(w, mass); //MAKE FIG 2
  //MakeMeanSigShapeTrend(w); //MAKE FIG6
  // MakeFitSIGMASigShapeTrend(); //MAKE FIG8
  
  /*PlotSigShapeCorrected(w, 0); 
  PlotSigShapeCorrected(w, 1); 
  PlotSigShapeCorrected(w, 2); 
  PlotSigShapeCorrected(w, 3); 
  */ 
  
 

  return;
}


void MakeMeanSigShapeTrend(RooWorkspace* w){


  Float_t sig0[7];
  Float_t sig1[7];
  Float_t sig2[7];
  Float_t sig3[7];

  Float_t sig0rms[7];
  Float_t sig1rms[7];
  Float_t sig2rms[7];
  Float_t sig3rms[7];

  Float_t sig0alphaPOS[7];
  Float_t sig1alphaPOS[7];
  Float_t sig2alphaPOS[7];
  Float_t sig3alphaPOS[7];

  Float_t sig0alphaNEG[7];
  Float_t sig1alphaNEG[7];
  Float_t sig2alphaNEG[7];
  Float_t sig3alphaNEG[7];

  Float_t sig0nPOS[7];
  Float_t sig1nPOS[7];
  Float_t sig2nPOS[7];
  Float_t sig3nPOS[7];

  Float_t sig0nNEG[7];
  Float_t sig1nNEG[7];
  Float_t sig2nNEG[7];
  Float_t sig3nNEG[7];

  Float_t sig0frac[7];
  Float_t sig1frac[7];
  Float_t sig2frac[7];
  Float_t sig3frac[7];

  Float_t sig0Err[7];
  Float_t sig1Err[7];
  Float_t sig2Err[7];
  Float_t sig3Err[7];

  Float_t sig0rmsErr[7];
  Float_t sig1rmsErr[7];
  Float_t sig2rmsErr[7];
  Float_t sig3rmsErr[7];

  Float_t sig0alphaPOSErr[7];
  Float_t sig1alphaPOSErr[7];
  Float_t sig2alphaPOSErr[7];
  Float_t sig3alphaPOSErr[7];

  Float_t sig0alphaNEGErr[7];
  Float_t sig1alphaNEGErr[7];
  Float_t sig2alphaNEGErr[7];
  Float_t sig3alphaNEGErr[7];

  Float_t sig0nPOSErr[7];
  Float_t sig1nPOSErr[7];
  Float_t sig2nPOSErr[7];
  Float_t sig3nPOSErr[7];

  Float_t sig0nNEGErr[7];
  Float_t sig1nNEGErr[7];
  Float_t sig2nNEGErr[7];
  Float_t sig3nNEGErr[7];

  Float_t sig0fracErr[7];
  Float_t sig1fracErr[7];
  Float_t sig2fracErr[7];
  Float_t sig3fracErr[7];

  Int_t nmass = 7;
  Float_t mass[7] = {150, 200, 250, 300, 400, 600, 800};
  Float_t massErr[7] = {0, 0, 0, 0, 0, 0, 0};

  PlotSigShape(w, 0,sig0,sig0rms, sig0Err,sig0rmsErr ); 
  PlotSigShape(w, 1,sig1,sig1rms, sig1Err,sig1rmsErr); 
  PlotSigShape(w, 2,sig2,sig2rms, sig2Err,sig2rmsErr); 
  PlotSigShape(w, 3,sig3,sig3rms, sig3Err,sig3rmsErr); 


  //plot mean
  TGraphErrors* g0 = new TGraphErrors(7, mass, sig0, massErr, sig0Err);
  TGraphErrors* g1 = new TGraphErrors(7, mass, sig1, massErr, sig1Err);
  TGraphErrors* g2 = new TGraphErrors(7, mass, sig2, massErr, sig2Err);
  TGraphErrors* g3 = new TGraphErrors(7, mass, sig3, massErr, sig3Err);
 

  g0->SetMarkerColor(kPink-8);
  g1->SetMarkerColor(kOrange+1);
  g2->SetMarkerColor(kSpring-8);
  g3->SetMarkerColor(kBlue-7);

  g0->SetLineColor(kPink-8);
  g1->SetLineColor(kOrange+1);
  g2->SetLineColor(kSpring-8);
  g3->SetLineColor(kBlue-7);


  TCanvas* c1 = new TCanvas("c1", "c1", 1);
  c1->cd();

  TMultiGraph* mg= new TMultiGraph();
  mg->Add(g0);
  mg->Add(g1);
  mg->Add(g2);
  mg->Add(g3);
  mg->Draw("APL");
  mg->GetYaxis()->SetTitle("#mu of #frac{#Delta m}{m}");
  mg->GetXaxis()->SetTitle("m_{X} [GeV]");
  mg->GetYaxis()->SetRangeUser(-0.02, 0.02);
  TLegend* legmc = new TLegend(0.62, 0.6, 0.87, 0.89, "", "brNDC");
  legmc->AddEntry(g0,"Cat 0","LP");
  legmc->AddEntry(g1,"Cat 1","LP");
  legmc->AddEntry(g2,"Cat 2","LP");
  legmc->AddEntry(g3,"Cat 3","LP");
  legmc->SetTextFont(42);
  legmc->SetBorderSize(0);
  legmc->SetFillStyle(0);
  legmc->Draw();


  TPaveText* label_cms = get_labelCMS(0, "2012", true);
  TPaveText* label_sqrt = get_labelSqrt(0);
  label_cms->Draw("same");
  label_sqrt->Draw("same");
 
  c1->SetLogy(0);
  c1->SaveAs("plots/MeanSigShape.png");
  c1->SaveAs("plots/MeanSigShape.pdf");
  c1->SaveAs("~/www/plotsNota/MeanSigShape.png");
  c1->SaveAs("~/www/plotsNota/MeanSigShape.pdf");

  c1->Clear();

  //plot RMS

  Float_t sig0Ratio[7];
  Float_t sig1Ratio[7];
  Float_t sig2Ratio[7];
  Float_t sig3Ratio[7];

  for(int i  = 0; i< 7; i++){

    sig0Ratio[i] = sig0rms[i]/sig0rms[0];
    sig1Ratio[i] = sig1rms[i]/sig1rms[0];
    sig2Ratio[i] = sig2rms[i]/sig2rms[0];
    sig3Ratio[i] = sig3rms[i]/sig3rms[0];

  }

  TGraphErrors* g0rms = new TGraphErrors(7, mass, sig0rms, massErr, sig0rmsErr);
  TGraphErrors* g1rms = new TGraphErrors(7, mass, sig1rms, massErr, sig1rmsErr);
  TGraphErrors* g2rms = new TGraphErrors(7, mass, sig2rms, massErr, sig2rmsErr);
  TGraphErrors* g3rms = new TGraphErrors(7, mass, sig3rms, massErr, sig3rmsErr);
 

  g0rms->SetMarkerColor(kPink-8);
  g1rms->SetMarkerColor(kOrange+1);
  g2rms->SetMarkerColor(kSpring-8);
  g3rms->SetMarkerColor(kBlue-7);

  g0rms->SetLineColor(kPink-8);
  g1rms->SetLineColor(kOrange+1);
  g2rms->SetLineColor(kSpring-8);
  g3rms->SetLineColor(kBlue-7);
 
  TMultiGraph* mgrms= new TMultiGraph();
  mgrms->Add(g0rms);
  mgrms->Add(g1rms);
  mgrms->Add(g2rms);
  mgrms->Add(g3rms);
  mgrms->Draw("APL");
  mgrms->GetYaxis()->SetTitle("#sigma_{CB} of #frac{#Delta m}{m}");
  mgrms->GetXaxis()->SetTitle("m_{X} [GeV]");
  mgrms->GetYaxis()->SetRangeUser(0., 0.05);
  
  legmc->Draw();


  TPaveText* label_cms = get_labelCMS(0, "2012", true);
  TPaveText* label_sqrt = get_labelSqrt(0);
  label_cms->Draw("same");
  label_sqrt->Draw("same");
 
  c1->SetLogy(0);
  c1->SaveAs("plots/RMSSigShape.png");
  c1->SaveAs("plots/RMSSigShape.pdf");
  c1->SaveAs("~/www/plotsNota/RMSSigShape.png");
  c1->SaveAs("~/www/plotsNota/RMSSigShape.pdf");


  /*
  c1->Clear();

  //plot alpha POS
  TGraphErrors* g0alphaPOS = new TGraphErrors(7, mass, sig0alphaPOS, massErr, sig0alphaPOSErr);
  TGraphErrors* g1alphaPOS = new TGraphErrors(7, mass, sig1alphaPOS, massErr, sig0alphaPOSErr);
  TGraphErrors* g2alphaPOS = new TGraphErrors(7, mass, sig2alphaPOS, massErr, sig0alphaPOSErr);
  TGraphErrors* g3alphaPOS = new TGraphErrors(7, mass, sig3alphaPOS, massErr, sig0alphaPOSErr);
 

  g0alphaPOS->SetMarkerColor(kPink-8);
  g1alphaPOS->SetMarkerColor(kOrange+1);
  g2alphaPOS->SetMarkerColor(kSpring-8);
  g3alphaPOS->SetMarkerColor(kBlue-7);

  g0alphaPOS->SetLineColor(kPink-8);
  g1alphaPOS->SetLineColor(kOrange+1);
  g2alphaPOS->SetLineColor(kSpring-8);
  g3alphaPOS->SetLineColor(kBlue-7);
 
  TMultiGraph* mgalphaPOS= new TMultiGraph();
  mgalphaPOS->Add(g0alphaPOS);
  mgalphaPOS->Add(g1alphaPOS);
  mgalphaPOS->Add(g2alphaPOS);
  mgalphaPOS->Add(g3alphaPOS);
  mgalphaPOS->Draw("APL");
  mgalphaPOS->GetYaxis()->SetTitle("#alpha^{POS}_{CB} of #frac{#Delta m}{m}");
  mgalphaPOS->GetXaxis()->SetTitle("m_{X} [GeV]");
  mgalphaPOS->GetYaxis()->SetRangeUser(0., 0.05);
  
  legmc->Draw();


  TPaveText* label_cms = get_labelCMS(0, "2012", true);
  TPaveText* label_sqrt = get_labelSqrt(0);
  label_cms->Draw("same");
  label_sqrt->Draw("same");
 
  c1->SetLogy(0);
  c1->SaveAs("plots/ALPHAPOSSigShape.png");
  c1->SaveAs("plots/ALPHAPOSSigShape.pdf");
  c1->SaveAs("~/www/plotsNota/ALPHAPOSSigShape.png");
  c1->SaveAs("~/www/plotsNota/ALPHAPOSSigShape.pdf");

  c1->Clear();

  //plot alpha NEG

  TGraphErrors* g0alphaNEG = new TGraphErrors(7, mass, sig0alphaNEG, massErr, sig0alphaNEGErr);
  TGraphErrors* g1alphaNEG = new TGraphErrors(7, mass, sig1alphaNEG, massErr, sig0alphaNEGErr);
  TGraphErrors* g2alphaNEG = new TGraphErrors(7, mass, sig2alphaNEG, massErr, sig0alphaNEGErr);
  TGraphErrors* g3alphaNEG = new TGraphErrors(7, mass, sig3alphaNEG, massErr, sig0alphaNEGErr);
 

  g0alphaNEG->SetMarkerColor(kPink-8);
  g1alphaNEG->SetMarkerColor(kOrange+1);
  g2alphaNEG->SetMarkerColor(kSpring-8);
  g3alphaNEG->SetMarkerColor(kBlue-7);

  g0alphaNEG->SetLineColor(kPink-8);
  g1alphaNEG->SetLineColor(kOrange+1);
  g2alphaNEG->SetLineColor(kSpring-8);
  g3alphaNEG->SetLineColor(kBlue-7);
 
  TMultiGraph* mgalphaNEG= new TMultiGraph();
  mgalphaNEG->Add(g0alphaNEG);
  mgalphaNEG->Add(g1alphaNEG);
  mgalphaNEG->Add(g2alphaNEG);
  mgalphaNEG->Add(g3alphaNEG);
  mgalphaNEG->Draw("APL");
  mgalphaNEG->GetYaxis()->SetTitle("#alpha^{NEG}_{CB} of #frac{#Delta m}{m}");
  mgalphaNEG->GetXaxis()->SetTitle("m_{H} [GeV]");
  mgalphaNEG->GetYaxis()->SetRangeUser(0., 0.05);
  
  legmc->Draw();


  TPaveText* label_cms = get_labelCMS(0, "2012", true);
  TPaveText* label_sqrt = get_labelSqrt(0);
  label_cms->Draw("same");
  label_sqrt->Draw("same");
 
  c1->SetLogy(0);
  c1->SaveAs("plots/ALPHANEGSigShape.png");
  c1->SaveAs("plots/ALPHANEGSigShape.pdf");
  c1->SaveAs("~/www/plotsNota/ALPHANEGSigShape.png");
  c1->SaveAs("~/www/plotsNota/ALPHANEGSigShape.pdf");


  c1->Clear();

  //plot n POS

  TGraphErrors* g0nPOS = new TGraphErrors(7, mass, sig0nPOS, massErr, sig0nPOSErr);
  TGraphErrors* g1nPOS = new TGraphErrors(7, mass, sig1nPOS, massErr, sig0nPOSErr);
  TGraphErrors* g2nPOS = new TGraphErrors(7, mass, sig2nPOS, massErr, sig0nPOSErr);
  TGraphErrors* g3nPOS = new TGraphErrors(7, mass, sig3nPOS, massErr, sig0nPOSErr);
 

  g0nPOS->SetMarkerColor(kPink-8);
  g1nPOS->SetMarkerColor(kOrange+1);
  g2nPOS->SetMarkerColor(kSpring-8);
  g3nPOS->SetMarkerColor(kBlue-7);

  g0nPOS->SetLineColor(kPink-8);
  g1nPOS->SetLineColor(kOrange+1);
  g2nPOS->SetLineColor(kSpring-8);
  g3nPOS->SetLineColor(kBlue-7);
 
  TMultiGraph* mgnPOS= new TMultiGraph();
  mgnPOS->Add(g0nPOS);
  mgnPOS->Add(g1nPOS);
  mgnPOS->Add(g2nPOS);
  mgnPOS->Add(g3nPOS);
  mgnPOS->Draw("APL");
  mgnPOS->GetYaxis()->SetTitle("#n^{POS}_{CB} of #frac{#Delta m}{m}");
  mgnPOS->GetXaxis()->SetTitle("m_{H} [GeV]");
  mgnPOS->GetYaxis()->SetRangeUser(0., 0.05);
  
  legmc->Draw();


  TPaveText* label_cms = get_labelCMS(0, "2012", true);
  TPaveText* label_sqrt = get_labelSqrt(0);
  label_cms->Draw("same");
  label_sqrt->Draw("same");
 
  c1->SetLogy(0);
  c1->SaveAs("plots/NPOSSigShape.png");
  c1->SaveAs("plots/NPOSSigShape.pdf");
  c1->SaveAs("~/www/plotsNota/NPOSSigShape.png");
  c1->SaveAs("~/www/plotsNota/NPOSSigShape.pdf");


  c1->Clear();

  //plot n NEG

  TGraphErrors* g0nNEG = new TGraphErrors(7, mass, sig0nNEG, massErr, sig0nNEGErr);
  TGraphErrors* g1nNEG = new TGraphErrors(7, mass, sig1nNEG, massErr, sig0nNEGErr);
  TGraphErrors* g2nNEG = new TGraphErrors(7, mass, sig2nNEG, massErr, sig0nNEGErr);
  TGraphErrors* g3nNEG = new TGraphErrors(7, mass, sig3nNEG, massErr, sig0nNEGErr);
 

  g0nNEG->SetMarkerColor(kPink-8);
  g1nNEG->SetMarkerColor(kOrange+1);
  g2nNEG->SetMarkerColor(kSpring-8);
  g3nNEG->SetMarkerColor(kBlue-7);

  g0nNEG->SetLineColor(kPink-8);
  g1nNEG->SetLineColor(kOrange+1);
  g2nNEG->SetLineColor(kSpring-8);
  g3nNEG->SetLineColor(kBlue-7);
 
  TMultiGraph* mgnNEG= new TMultiGraph();
  mgnNEG->Add(g0nNEG);
  mgnNEG->Add(g1nNEG);
  mgnNEG->Add(g2nNEG);
  mgnNEG->Add(g3nNEG);
  mgnNEG->Draw("APL");
  mgnNEG->GetYaxis()->SetTitle("#n^{NEG}_{CB} of #frac{#Delta m}{m}");
  mgnNEG->GetXaxis()->SetTitle("m_{H} [GeV]");
  mgnNEG->GetYaxis()->SetRangeUser(0., 0.05);
  
  legmc->Draw();


  TPaveText* label_cms = get_labelCMS(0, "2012", true);
  TPaveText* label_sqrt = get_labelSqrt(0);
  label_cms->Draw("same");
  label_sqrt->Draw("same");
 
  c1->SetLogy(0);
  c1->SaveAs("plots/NNEGSigShape.png");
  c1->SaveAs("plots/NNEGSigShape.pdf");
  c1->SaveAs("~/www/plotsNota/NNEGSigShape.png");
  c1->SaveAs("~/www/plotsNota/NNEGSigShape.pdf");


  c1->Clear();

  //plot n NEG

  TGraphErrors* g0frac = new TGraphErrors(7, mass, sig0frac, massErr, sig0fracErr);
  TGraphErrors* g1frac = new TGraphErrors(7, mass, sig1frac, massErr, sig0fracErr);
  TGraphErrors* g2frac = new TGraphErrors(7, mass, sig2frac, massErr, sig0fracErr);
  TGraphErrors* g3frac = new TGraphErrors(7, mass, sig3frac, massErr, sig0fracErr);
 

  g0frac->SetMarkerColor(kPink-8);
  g1frac->SetMarkerColor(kOrange+1);
  g2frac->SetMarkerColor(kSpring-8);
  g3frac->SetMarkerColor(kBlue-7);

  g0frac->SetLineColor(kPink-8);
  g1frac->SetLineColor(kOrange+1);
  g2frac->SetLineColor(kSpring-8);
  g3frac->SetLineColor(kBlue-7);
 
  TMultiGraph* mgfrac= new TMultiGraph();
  mgfrac->Add(g0frac);
  mgfrac->Add(g1frac);
  mgfrac->Add(g2frac);
  mgfrac->Add(g3frac);
  mgfrac->Draw("APL");
  mgfrac->GetYaxis()->SetTitle("frac_{CB} of #frac{#Delta m}{m}");
  mgfrac->GetXaxis()->SetTitle("m_{H} [GeV]");
  mgfrac->GetYaxis()->SetRangeUser(0., 0.05);
  
  legmc->Draw();


  TPaveText* label_cms = get_labelCMS(0, "2012", true);
  TPaveText* label_sqrt = get_labelSqrt(0);
  label_cms->Draw("same");
  label_sqrt->Draw("same");
 
  c1->SetLogy(0);
  c1->SaveAs("plots/fracSigShape.png");
  c1->SaveAs("plots/fracSigShape.pdf");
  c1->SaveAs("~/www/plotsNota/fracSigShape.png");
  c1->SaveAs("~/www/plotsNota/fracSigShape.pdf");


*/

  std::cout<<"----- MEAN -------"<<std::endl;
  for (int i = 0; i< 7; i++)std::cout<< sig0[i]<<"     ";
  std::cout<<"   "<<std::endl;
  for (int i = 0; i< 7; i++)std::cout<< sig1[i]<<"     ";
  std::cout<<"   "<<std::endl;
  for (int i = 0; i< 7; i++)std::cout<< sig2[i]<<"     ";
  std::cout<<"   "<<std::endl;
  for (int i = 0; i< 7; i++)std::cout<< sig3[i]<<"     ";
  std::cout<<"   "<<std::endl;
  std::cout<<"----- SIGMA -------"<<std::endl;
  for (int i = 0; i< 7; i++)std::cout<< sig0rms[i]<<"     ";
  std::cout<<"   "<<std::endl;
  for (int i = 0; i< 7; i++)std::cout<< sig1rms[i]<<"     ";
  std::cout<<"   "<<std::endl;
  for (int i = 0; i< 7; i++)std::cout<< sig2rms[i]<<"     ";
  std::cout<<"   "<<std::endl;
  for (int i = 0; i< 7; i++)std::cout<< sig3rms[i]<<"     ";
  std::cout<<"   "<<std::endl;

  std::cout<<"----- MEAN -------"<<std::endl;
  for (int i = 0; i< 7; i++)std::cout<< sig0Err[i]<<"     ";
  std::cout<<"   "<<std::endl;
  for (int i = 0; i< 7; i++)std::cout<< sig1Err[i]<<"     ";
  std::cout<<"   "<<std::endl;
  for (int i = 0; i< 7; i++)std::cout<< sig2Err[i]<<"     ";
  std::cout<<"   "<<std::endl;
  for (int i = 0; i< 7; i++)std::cout<< sig3Err[i]<<"     ";
  std::cout<<"   "<<std::endl;
  std::cout<<"----- SIGMA -------"<<std::endl;
  for (int i = 0; i< 7; i++)std::cout<< sig0rmsErr[i]<<"     ";
  std::cout<<"   "<<std::endl;
  for (int i = 0; i< 7; i++)std::cout<< sig1rmsErr[i]<<"     ";
  std::cout<<"   "<<std::endl;
  for (int i = 0; i< 7; i++)std::cout<< sig2rmsErr[i]<<"     ";
  std::cout<<"   "<<std::endl;
  for (int i = 0; i< 7; i++)std::cout<< sig3rmsErr[i]<<"     ";
  std::cout<<"   "<<std::endl;

  /*  double corr0[7];
  double corr1[7];
  double corr2[7];
  double corr3[7];

  for (int i = 0; i< 7; i++){

    corr0[i] = sqrt(sig0rms[0]*sig0rms[0] - sig0rms[i]*sig0rms[i]);
    corr1[i] = sqrt(sig1rms[0]*sig1rms[0] - sig1rms[i]*sig1rms[i]);
    corr2[i] = sqrt(sig2rms[0]*sig2rms[0] - sig2rms[i]*sig2rms[i]);
    corr3[i] = sqrt(sig3rms[0]*sig3rms[0] - sig3rms[i]*sig3rms[i]);
  }


  //plot correction
 c1->Clear();

  //plot RMS
  TGraph* g0corr = new TGraph(7, mass,corr0 );
  TGraph* g1corr = new TGraph(7, mass, corr1);
  TGraph* g2corr = new TGraph(7, mass,corr2 );
  TGraph* g3corr = new TGraph(7, mass, corr3);
 

  g0corr->SetMarkerColor(kPink-8);
  g1corr->SetMarkerColor(kOrange+1);
  g2corr->SetMarkerColor(kSpring-8);
  g3corr->SetMarkerColor(kBlue-7);

  g0corr->SetLineColor(kPink-8);
  g1corr->SetLineColor(kOrange+1);
  g2corr->SetLineColor(kSpring-8);
  g3corr->SetLineColor(kBlue-7);
 
  TMultiGraph* mgrms= new TMultiGraph();
  mgcorr->Add(g0rms);
  mgcorr->Add(g1rms);
  mgcorr->Add(g2rms);
  mgcorr->Add(g3rms);
  mgcorr->Draw("APL");
  mgcorr->GetYaxis()->SetTitle("Smearing of #frac{#Delta m}{m}");
  mgcorr->GetXaxis()->SetTitle("m_{H} [GeV]");
  mgcorr->GetYaxis()->SetRangeUser(0., 0.01);
  
  legmc->Draw();


  TPaveText* label_cms = get_labelCMS(0, "2012", true);
  TPaveText* label_sqrt = get_labelSqrt(0);
  label_cms->Draw("same");
  label_sqrt->Draw("same");
 
  c1->SetLogy(0);
  c1->SaveAs("plots/RMSCORRSigShape.png");
  c1->SaveAs("plots/RMSCORRSigShape.pdf");
  c1->SaveAs("~/www/plotsNota/RMSCORRSigShape.png");
  c1->SaveAs("~/www/plotsNota/RMSCORRSigShape.pdf");*/

}



void MakeFitSIGMASigShapeTrend(){


  Float_t sig0rms[7]={0.00970353 ,   0.00963705 ,   0.00911126 ,   0.00880483 ,   0.00875701 ,   0.00842481 ,   0.00827415  };      
  Float_t sig1rms[7]={0.0120229 ,   0.0112211 ,   0.0111983 ,   0.0105888 ,   0.0101585 ,   0.0101493 ,   0.00987028        };
  Float_t sig2rms[7]={0.0185319 ,   0.017388 ,   0.0168306 ,   0.0168686 ,   0.0158064 ,   0.0152674 ,   0.0149237          };
  Float_t sig3rms[7]={0.0202395 ,   0.0188697 ,   0.0182573 ,   0.0178627 ,   0.0170726 ,   0.0165754 ,   0.0165886           };


  Float_t sig0rmsErr[7]={0.000171741 ,   0.000153105 ,   0.00031357 ,   0.000131645 ,   0.000112986 ,   9.29385e-05 ,   8.3895e-05     };    
  Float_t sig1rmsErr[7]={0.000190229 ,   0.000169746 ,   0.000427034 ,   0.000178806 ,   0.000164291 ,   0.000171058 ,   0.000178646   };  
  Float_t sig2rmsErr[7]={0.000444902 ,   0.000364538 ,   0.000766693 ,   0.000335629 ,   0.000277103 ,   0.000255246 ,   0.00026249    };  
  Float_t sig3rmsErr[7]={0.000401814 ,   0.00035561 ,   0.000846053 ,   0.000383018 ,   0.000365251 ,   0.000397828 ,   0.000447633   	};

 



  
  Float_t sig0Ratio[7];
  Float_t sig1Ratio[7];
  Float_t sig2Ratio[7];
  Float_t sig3Ratio[7];

  
  Float_t sig0RatioErr[7];
  Float_t sig1RatioErr[7];
  Float_t sig2RatioErr[7];
  Float_t sig3RatioErr[7];
  
  for(int i  = 0; i< 7; i++){
    
    sig0Ratio[i] = sig0rms[i]/sig0rms[0];
    sig1Ratio[i] = sig1rms[i]/sig1rms[0];
    sig2Ratio[i] = sig2rms[i]/sig2rms[0];
    sig3Ratio[i] = sig3rms[i]/sig3rms[0];

    sig0RatioErr[i] = sqrt(pow((1./sig0rms[0]),2)*sig0rmsErr[i]*sig0rmsErr[i]+pow((sig0rms[i]/sig0rms[0]/sig0rms[0]),2)*sig0rmsErr[0]*sig0rmsErr[0]);
    sig1RatioErr[i] = sqrt(pow((1./sig1rms[0]),2)*sig1rmsErr[i]*sig1rmsErr[i]+pow((sig1rms[i]/sig1rms[0]/sig1rms[0]),2)*sig1rmsErr[0]*sig1rmsErr[0]);
    sig2RatioErr[i] = sqrt(pow((1./sig2rms[0]),2)*sig2rmsErr[i]*sig2rmsErr[i]+pow((sig2rms[i]/sig2rms[0]/sig2rms[0]),2)*sig2rmsErr[0]*sig2rmsErr[0]);
    sig3RatioErr[i] = sqrt(pow((1./sig3rms[0]),2)*sig3rmsErr[i]*sig3rmsErr[i]+pow((sig3rms[i]/sig3rms[0]/sig3rms[0]),2)*sig3rmsErr[0]*sig3rmsErr[0]);
   
  }


  Int_t nmass = 7;
  Float_t mass[7] = {150, 200, 250, 300, 400, 600, 850};
  Float_t massErr[7] = {0, 0, 0, 0, 0, 0, 0};


 
  TCanvas* c1 = new TCanvas("c1", "c1", 1);
  c1->cd();


  TPaveText* label_cms = get_labelCMS(0, "2012", true);
  TPaveText* label_sqrt = get_labelSqrt(0);
  label_cms->Draw("same");
  label_sqrt->Draw("same");
 
 
  //plot RMS
  TGraphErrors* g0rms = new TGraphErrors(7, mass, sig0Ratio,massErr, sig0RatioErr);
  TGraphErrors* g1rms = new TGraphErrors(7, mass, sig1Ratio,massErr, sig1RatioErr);
  TGraphErrors* g2rms = new TGraphErrors(7, mass, sig2Ratio,massErr, sig2RatioErr);
  TGraphErrors* g3rms = new TGraphErrors(7, mass, sig3Ratio,massErr, sig3RatioErr);
 

  g0rms->SetMarkerColor(kPink-8);
  g1rms->SetMarkerColor(kOrange+1);
  g2rms->SetMarkerColor(kSpring-8);
  g3rms->SetMarkerColor(kBlue-7);

  g0rms->SetLineColor(kPink-8);
  g1rms->SetLineColor(kOrange+1);
  g2rms->SetLineColor(kSpring-8);
  g3rms->SetLineColor(kBlue-7);
 

  TLegend* legmc = new TLegend(0.62, 0.6, 0.87, 0.89, "", "brNDC");
 
  legmc->SetTextFont(42);
  legmc->SetBorderSize(0);
  legmc->SetFillStyle(0);

  TF1* f0 = new TF1("f0",  "pol2", 150, 850);
  // f0->SetLineStyle(kDashed);
  f0->SetLineColor(kPink-8);
  f0->SetLineWidth(2);
  f0->SetParLimits(0, 0.008, 0.012);
  f0->SetParLimits(2, 0., 100);
  TF1* f1 = new TF1("f1",  "pol2", 150, 850);
  // f1->SetLineStyle(kDashed);
  f1->SetLineColor(kOrange+1);
  f1->SetLineWidth(2);
  f1->SetParLimits(0, 0.01, 0.014);
  f1->SetParLimits(2, 0., 100);
  TF1* f2 = new TF1("f2",  "pol2", 150, 850);
  // f2->SetLineStyle(kDashed);
  f2->SetLineColor(kSpring-8);
  f2->SetLineWidth(2);
  f2->SetParLimits(0, 0.018, 0.02);
  f2->SetParLimits(2, 0., 100);
  TF1* f3 = new TF1("f3",  "pol2", 150, 850);
  // f3->SetLineStyle(kDashed);
  f3->SetLineColor(kBlue-7);
  f3->SetLineWidth(2);
  f3->SetParLimits(0, 0.018, 0.03);
  f3->SetParLimits(2, 0., 100);


  g0rms->Fit("f0");
  g1rms->Fit("f1");
  g2rms->Fit("f2");
  g3rms->Fit("f3");

  /* TMultiGraph* mgrms= new TMultiGraph();
  mgrms->Add(g0rms);
  mgrms->Add(g1rms);
  mgrms->Add(g2rms);
  mgrms->Add(g3rms);
  mgrms->Draw("APL");
  mgrms->GetYaxis()->SetTitle("#sigma_{CB}/#sigma_{CB}_{@150} of #frac{#Delta m}{m}");
  mgrms->GetXaxis()->SetTitle("m_{X} [GeV]");
  mgrms->GetYaxis()->SetRangeUser(0.6, 1.4);
  
*/

  f0->GetYaxis()->SetTitle("#sigma_{CB}/#sigma_{CB}_{@150} of #frac{#Delta m}{m}");
  f0->GetXaxis()->SetTitle("m_{X} [GeV]");
  f0->GetYaxis()->SetRangeUser(0.6, 1.4);

  f0->Draw("L");
  f1->Draw("Lsame");
  f2->Draw("Lsame");
  f3->Draw("Lsame");
  legmc->AddEntry(f0,"Cat 0","L");
  legmc->AddEntry(f1,"Cat 1","L");
  legmc->AddEntry(f2,"Cat 2","L");
  legmc->AddEntry(f3,"Cat 3","L");

  legmc->Draw();


  TPaveText* label_cms = get_labelCMS(0, "2012", true);
  TPaveText* label_sqrt = get_labelSqrt(0);
  // label_cms->Draw("same");
  //label_sqrt->Draw("same");
  gPad->Update();
  gStyle->SetOptFit(0011);
  /* TPaveStats* sb0=(TPaveStats*)(g0rms->GetListOfFunctions()->FindObject("stats"));
  sb0->SetOptFit(11);
  sb0->SetX1NDC(.6);
  sb0->SetX2NDC(.85);
  sb0->SetY1NDC(.45);
  sb0->SetY2NDC(.56);
  sb0->SetTextColor(kPink-8);
  gPad->Modified();
  TPaveStats* sb1=(TPaveStats*)(g1rms->GetListOfFunctions()->FindObject("stats"));
  sb1->SetOptFit(11);
  sb1->SetX1NDC(.6);
  sb1->SetX2NDC(.85);
  sb1->SetY1NDC(.56);
  sb1->SetY2NDC(.67);
  sb1->SetTextColor(kOrange+1);
  gPad->Modified();
  TPaveStats* sb2=(TPaveStats*)(g2rms->GetListOfFunctions()->FindObject("stats"));
  sb2->SetOptFit(11);
  sb2->SetX1NDC(.6);
  sb2->SetX2NDC(.85);
  sb2->SetY1NDC(.67);
  sb2->SetY2NDC(.78);
  sb2->SetTextColor(kSpring-8);
  gPad->Modified();

  TPaveStats* sb3=(TPaveStats*)(g3rms->GetListOfFunctions()->FindObject("stats"));
  sb3->SetOptFit(11);
  sb3->SetX1NDC(.6);
  sb3->SetX2NDC(.85.);
  sb3->SetY1NDC(.78);
  sb3->SetY2NDC(.89);
  sb3->SetTextColor(kBlue-7);
  gPad->Modified();
  */

  int iPos=11 ;
  CMS_lumi( c1,true,iPos );
  c1->SetLogy(0);
  c1->SaveAs("plots/RMSFITSigShape.png");
  c1->SaveAs("plots/RMSFITSigShape.pdf");
  c1->SaveAs("~/www/plotsNota/RMSFITSigShape.png");
  c1->SaveAs("~/www/plotsNota/RMSFITSigShape.pdf");
  c1->SaveAs("~/www/plotsPAS/RMSFITSigShape.png");
  c1->SaveAs("~/www/plotsPAS/RMSFITSigShape.pdf");
  c1->SaveAs("~/www/plotsPAS/RMSFITSigShape.C");
  c1->SaveAs("~/www/plotsPAS/RMSFITSigShape.root");

 
  /* TFile* f = new TFile("sigShapeCorrections.root", "RECREATE");
  f->cd();
  f0->Write();
  f1->Write();
  f2->Write();
  f3->Write();
  f->Write();
  f->Close();
*/

}




void MakePlots(RooWorkspace* w, Float_t mass, RooFitResult* fitresults, bool blind) {

  Int_t ncat = NCAT;
  
  cout << endl; cout << "Retreive everything:" << endl; 
  w->Print();

  // retrieve data sets from the workspace
  RooDataSet* dataAll   = (RooDataSet*) w->data("Data");
  RooDataSet* signalAll = (RooDataSet*) w->data("SigWeight");
  
  // maximum 9 cat...
  RooDataSet* data[9];  
  RooDataSet* signal[9];
  RooAbsPdf*  PhotonsMassGaussSig[9];
  RooAbsPdf*  PhotonsMassGaussSig_bis[9];
  RooAbsPdf*  PhotonsMassSig[9];
  RooExtendPdf*  PhotonsMassBkg[9];  
  
  for (int c=0; c<ncat; ++c) {
    data[c]   = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    signal[c] = (RooDataSet*) w->data(TString::Format("SigWeight_cat%d",c));
    PhotonsMassGaussSig[c]  = (RooAbsPdf*)  w->pdf(TString::Format("PhotonsMassGaussSig_cat%d",c));
    PhotonsMassGaussSig_bis[c]     = (RooAbsPdf*)  w->pdf(TString::Format("PhotonsMassGaussSig_cat%d_bis",c));
    PhotonsMassSig[c]       = (RooAbsPdf*)  w->pdf("PhotonsMassSig"+TString::Format("_cat%d",c));
    PhotonsMassBkg[c]       = (RooExtendPdf*)  w->pdf(TString::Format("PhotonsMassBkg_cat%d",c));
  }
  
  // retrieve mass observable from the workspace
  RooRealVar* PhotonsMass = w->var("PhotonsMass");  
  PhotonsMass->setUnit("GeV");
  
  // retrieve pdfs after the fits
  RooAbsPdf* PhotonsMassGaussSigAll = w->pdf("PhotonsMassGaussSig");
  RooAbsPdf* PhotonsMassGaussSigAll_bis    = w->pdf("PhotonsMassGaussSig_bis");
  RooAbsPdf* PhotonsMassSigAll      = w->pdf("PhotonsMassSig");
  RooAbsPdf* PhotonsMassBkgAll      = w->pdf("PhotonsMassBkgAll");

  
 Float_t minMassFit, maxMassFit;
 
 if(mass<200){ //m150
    minMassFit = 100;
    maxMassFit = 200;
  } else if(mass>=200 && mass <500){//m200 and m250
    minMassFit = 180;
    maxMassFit = 500;
  }else{//m300 m400
    minMassFit = 280;
    maxMassFit = 500;
  }


  Float_t MASS(mass);
  int iMass = abs(mass);
  
  Int_t nBinsMass(0.2*mass);  // chiara


  /*  // ****************************SIG ALL
  cout << endl; cout << "Progress plotting: signal" << endl;     

  RooPlot* plotPhotonsMassAll = PhotonsMass->frame(Range(minMassFit,maxMassFit),Bins(nBinsMass));
  signalAll->plotOn(plotPhotonsMassAll);
  
  gStyle->SetOptTitle(0);
  PhotonsMassSigAll->plotOn(plotPhotonsMassAll);
  PhotonsMassSigAll->plotOn(plotPhotonsMassAll,Components("PhotonsMassGaussSig"),LineStyle(kDashed),LineColor(kGreen));
  PhotonsMassSigAll->plotOn(plotPhotonsMassAll,Components("PhotonsMassGaussSig_bis"),LineStyle(kDashed),LineColor(kRed));
  //PhotonsMassSigAll->paramOn(plotPhotonsMassAll, ShowConstants(true), Layout(0.65,0.75,0.85), Format("NEU",AutoPrecision(2)));
  //plotPhotonsMassAll->getAttText()->SetTextSize(0.03);
  
  TCanvas* c1 = new TCanvas("c1","PhotonsMass",0,0,500,500);
  c1->cd(1);
  plotPhotonsMassAll->Draw();  
  c1->SaveAs("plots/sigmodel_"+TString::Format("%d.png", iMass));
  c1->SaveAs("plots/sigmodel_"+TString::Format("%d.root", iMass));
*/

  // ****************************SIG CAT
  cout << endl; cout << "Progress plotting: signal per categories" << endl;     
  
  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
    
  RooPlot* plotPhotonsMass[9];
  for (int c=0; c<ncat; ++c) {
    plotPhotonsMass[c] = PhotonsMass->frame(Range(mass-30,mass+30),Bins(nBinsMass));
    signal[c]->plotOn(plotPhotonsMass[c],LineColor(kWhite),MarkerColor(kWhite));    
    PhotonsMassSig[c]  ->plotOn(plotPhotonsMass[c]);
    PhotonsMassSig[c]  ->plotOn(plotPhotonsMass[c],Components("PhotonsMassGaussSig"+TString::Format("_cat%d",c)),LineStyle(kDashed),LineColor(kGreen));
    PhotonsMassSig[c]  ->plotOn(plotPhotonsMass[c],Components("PhotonsMassGaussSig"+TString::Format("_cat%d_bis",c)),LineStyle(kDashed),LineColor(kRed));
    //PhotonsMassSig[c]  ->paramOn(plotPhotonsMass[c], ShowConstants(true), Layout(0.65,0.75,0.85), Format("NEU",AutoPrecision(2)));
    signal[c]  ->plotOn(plotPhotonsMass[c]);

    //  TCanvas* dummy = new TCanvas("dummy", "dummy",0, 0, 400, 400);
    // TH1F *hist = new TH1F("hist", "hist", 400, minMassFit, maxMassFit);
 
    plotPhotonsMass[c]->SetTitle("");      
    plotPhotonsMass[c]->SetMinimum(0.0);
    plotPhotonsMass[c]->SetMaximum(1.40*plotPhotonsMass[c]->GetMaximum());
    plotPhotonsMass[c]->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV)");

    TCanvas* ctmp = new TCanvas("ctmp","PhotonsMass Background Categories",0,0,500,500);
    plotPhotonsMass[c]->Draw();  
    plotPhotonsMass[c]->SetAxisRange(0.000001,plotPhotonsMass[c]->GetMaximum()*1.2,"Y");
    plotPhotonsMass[c]->Draw("SAME");  
    TLegend *legmc = new TLegend(0.62,0.75,0.92,0.9);
    legmc->AddEntry(plotPhotonsMass[c]->getObject(5),"Simulation","LPE");
    legmc->AddEntry(plotPhotonsMass[c]->getObject(1),"Parametric Model","L");
    legmc->AddEntry(plotPhotonsMass[c]->getObject(3),"Gaussian 2 component","L");
    legmc->AddEntry(plotPhotonsMass[c]->getObject(2),"Gaussian 1 component","L");
    
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->Draw();    
    
    TLatex *lat  = new TLatex(mass-20.,0.91*plotPhotonsMass[c]->GetMaximum(),"#scale[1.0]{CMS Preliminary}");
    lat->Draw();
    
    ctmp->SaveAs("plots/sigmodel_"+TString::Format("%d_cat%d.root", iMass, c));
    ctmp->SaveAs("plots/sigmodel_"+TString::Format("%d_cat%d.png", iMass, c));
    ctmp->SetLogy();
    ctmp->SaveAs("plots/sigmodel_"+TString::Format("%d_cat%d_LOG.png", iMass, c));
  }
 


  // ****************************BKG CAT LOG
  cout << endl; cout << "Progress plotting: background" << endl;     
  TCanvas* c4 = new TCanvas("c4","PhotonsMass Background Categories",0,0,400,400);
  RooPlot* plotPhotonsMassBkg[9];
  for (int c=0; c<ncat; ++c) {
    plotPhotonsMassBkg[c] = PhotonsMass->frame(Range(minMassFit, maxMassFit),nBinsMass);

    data[c]->plotOn(plotPhotonsMassBkg[c],RooFit::Invisible());    
    PhotonsMassBkg[c]->plotOn(plotPhotonsMassBkg[c],LineColor(kBlue),Range("fitrange"),NormRange("fitrange")); 
   
     if( blind ) {
      PhotonsMass->setRange("unblind_up",mass+10.,maxMassFit);
      PhotonsMass->setRange("unblind_down",minMassFit,mass-10.);
      data[c]->plotOn(plotPhotonsMassBkg[c], RooFit::CutRange("unblind_up"));    
      data[c]->plotOn(plotPhotonsMassBkg[c], RooFit::CutRange("unblind_down"));    
    } else {
      data[c]->plotOn(plotPhotonsMassBkg[c], Range("fitrange"));    
      } 
      
    plotPhotonsMassBkg[c]->Draw();  
    gPad->SetLogy(1);
    plotPhotonsMassBkg[c]->SetAxisRange(1.5,plotPhotonsMassBkg[c]->GetMaximum()*1.5,"Y");
    c4->SaveAs("plots/backgrounds_log_"+TString::Format("cat%d.png", c));
    c4->SaveAs("plots/backgrounds_log_"+TString::Format("cat%d.root", c));
  }
  
  // ****************************BKG LIN
  TCanvas* c5 = new TCanvas("c5","PhotonsMass Background Categories",0,0,400,400);
  RooPlot* plotPhotonsMassBkg[9];
  for (int c=0; c<ncat; ++c) {
    plotPhotonsMassBkg[c] = PhotonsMass->frame(Range(minMassFit, maxMassFit),nBinsMass);
    data[c]->plotOn(plotPhotonsMassBkg[c],RooFit::Invisible());    
    PhotonsMassBkg[c]->plotOn(plotPhotonsMassBkg[c],LineColor(kBlue),Range("fitrange"),NormRange("fitrange")); 
   
       if( blind ) {
      PhotonsMass->setRange("unblind_up",mass+10.,maxMassFit);
      PhotonsMass->setRange("unblind_down",minMassFit,mass-10.);
      data[c]->plotOn(plotPhotonsMassBkg[c], RooFit::CutRange("unblind_up"));    
      data[c]->plotOn(plotPhotonsMassBkg[c], RooFit::CutRange("unblind_down"));    
    } else {
      data[c]->plotOn(plotPhotonsMassBkg[c], Range("fitrange"));    
      } 
    

       plotPhotonsMassBkg[c]->SetAxisRange(0.001,plotPhotonsMassBkg[c]->GetMaximum()*1.5,"Y");
       plotPhotonsMassBkg[c]->Draw();  
       c5->SaveAs("plots/backgrounds_"+TString::Format("cat%d.png", c));
       c5->SaveAs("plots/backgrounds_"+TString::Format("cat%d.root", c));
  }
}




Double_t effSigma(TH1 *hist) {
  
  TAxis *xaxis = hist->GetXaxis();
  Int_t nb = xaxis->GetNbins();
  if(nb < 10) {
    std::cout << "effsigma: Not a valid histo. nbins = " << nb << std::endl;
    return 0.;
  }
  
  Double_t bwid = xaxis->GetBinWidth(1);
  if(bwid == 0) {
    std::cout << "effsigma: Not a valid histo. bwid = " << bwid << std::endl;
    return 0.;
  }
  Double_t xmax = xaxis->GetXmax();
  Double_t xmin = xaxis->GetXmin();
  Double_t ave = hist->GetMean();
  Double_t rms = hist->GetRMS();

  Double_t total=0.;
  for(Int_t i=0; i<nb+2; i++) {
    total+=hist->GetBinContent(i);
  }
  if(total < 100.) {
    std::cout << "effsigma: Too few entries " << total << std::endl;
    return 0.;
  }
  Int_t ierr=0;
  Int_t ismin=999;

  Double_t rlim=0.683*total;
  Int_t nrms=rms/(bwid);    // Set scan size to +/- rms
  if(nrms > nb/10) nrms=nb/10; // Could be tuned...

  Double_t widmin=9999999.;
  for(Int_t iscan=-nrms;iscan<nrms+1;iscan++) { // Scan window centre
    Int_t ibm=(ave-xmin)/bwid+1+iscan;
    Double_t x=(ibm-0.5)*bwid+xmin;
    Double_t xj=x;
    Double_t xk=x;
    Int_t jbm=ibm;
    Int_t kbm=ibm;
    Double_t bin=hist->GetBinContent(ibm);
    total=bin;
    for(Int_t j=1;j<nb;j++){
      if(jbm < nb) {
        jbm++;
        xj+=bwid;
        bin=hist->GetBinContent(jbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
      if(kbm > 0) {
        kbm--;
        xk-=bwid;
        bin=hist->GetBinContent(kbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
    }
    Double_t dxf=(total-rlim)*bwid/bin;
    Double_t wid=(xj-xk+bwid-dxf)*0.5;
    if(wid < widmin) {
      widmin=wid;
      ismin=iscan;
    }
  }
  if(ismin == nrms || ismin == -nrms) ierr=3;
  if(ierr != 0) std::cout << "effsigma: Error of type " << ierr << std::endl;

  return widmin;
}

Double_t effSigma(TH1 *hist) {

  TAxis *xaxis = hist->GetXaxis();
  Int_t nb = xaxis->GetNbins();
  if(nb < 10) {
    std::cout << "effsigma: Not a valid histo. nbins = " << nb << std::endl;
    return 0.;
  }

  Double_t bwid = xaxis->GetBinWidth(1);
  if(bwid == 0) {
    std::cout << "effsigma: Not a valid histo. bwid = " << bwid << std::endl;
    return 0.;
  }
  Double_t xmax = xaxis->GetXmax();
  Double_t xmin = xaxis->GetXmin();
  Double_t ave = hist->GetMean();
  Double_t rms = hist->GetRMS();

  Double_t total=0.;
  for(Int_t i=0; i<nb+2; i++) {
    total+=hist->GetBinContent(i);
  }
  if(total < 100.) {
    std::cout << "effsigma: Too few entries " << total << std::endl;
    return 0.;
  }
  Int_t ierr=0;
  Int_t ismin=999;

  Double_t rlim=0.683*total;
  Int_t nrms=rms/(bwid);    // Set scan size to +/- rms
  if(nrms > nb/10) nrms=nb/10; // Could be tuned...

  Double_t widmin=9999999.;
  for(Int_t iscan=-nrms;iscan<nrms+1;iscan++) { // Scan window centre
    Int_t ibm=(ave-xmin)/bwid+1+iscan;
    Double_t x=(ibm-0.5)*bwid+xmin;
    Double_t xj=x;
    Double_t xk=x;
    Int_t jbm=ibm;
    Int_t kbm=ibm;
    Double_t bin=hist->GetBinContent(ibm);
    total=bin;
    for(Int_t j=1;j<nb;j++){
      if(jbm < nb) {
        jbm++;
        xj+=bwid;
        bin=hist->GetBinContent(jbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
      if(kbm > 0) {
        kbm--;
        xk-=bwid;
        bin=hist->GetBinContent(kbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
    }
    Double_t dxf=(total-rlim)*bwid/bin;
    Double_t wid=(xj-xk+bwid-dxf)*0.5;
    if(wid < widmin) {
      widmin=wid;
      ismin=iscan;
    }
  }
  if(ismin == nrms || ismin == -nrms) ierr=3;
  if(ierr != 0) std::cout << "effsigma: Error of type " << ierr << std::endl;

  return widmin;
}

void R2JJFitter(double mass, std::string postfix="")
{
    filePOSTfix=postfix;
    if(postfix!="")
    {
      // for optimization studies
      MMIN=1000;
      if(mass==1000)
         signalScaler=0.034246;
      if((mass>1000)&&(mass<500))
         signalScaler=0.02469;
      if(mass==2000)
         signalScaler=2.0;
    };
    runfits(mass, 1);
    if(postfix!="")
    {
      // for optimization studies
      MMIN=1000;
      if(mass==1000)
         signalScaler=0.033500;
      if((mass>1000)&&(mass<2000))
         signalScaler=0.02016;
      if(mass==2000)
         signalScaler=2.22222;
    };
    runfits(mass, 0);
    runfits(mass, 2);
}



void MakeParameterTrendvsMassGauss(){

  TString wsDir = "workspaces/"+filePOSTfix;
  Int_t ncat = NCAT;
  Float_t masses[5];
  masses[0] = 150.0;
  masses[1] = 200.0;
  masses[2] = 250.0;
  masses[3] = 300.0;
  masses[4] = 400.0.;

  Float_t masses_err[5];
  masses_err[0] = 0.;
  masses_err[1] = 0.;
  masses_err[2] = 0.;
  masses_err[3] = 0.;
  masses_err[4] = 0.;


  TFile* file[5];
  for(imass = 0; imass<5;imass++) file[imass]= new TFile(wsDir+TString::Format("HighMass-hgg.m%.1f_8TeV.inputsig.root",masses[imass]));

  RooWorkspace* ws[5];
  for(imass = 0; imass<5;imass++) ws[imass] = (RooWorkspace*) file[imass]->Get("w_all");

 
  Float_t m0[5];
  Float_t m1[5];
  Float_t sigma0[5];
  Float_t sigma1[5];
  Float_t frac[5];
  Float_t sig[5];
  
  Float_t m0_err[5];
  Float_t m1_err[5];
  Float_t sigma0_err[5];
  Float_t sigma1_err[5];
  Float_t frac_err[5];
  Float_t sig_err[5];

  TCanvas* c1 = new TCanvas("c1", "c1", 1);
  c1->cd();
  TPaveText* label_cat[5];
  TPaveText* label_cms = get_labelCMS(0, "2012", true);
  TPaveText* label_sqrt = get_labelSqrt(0);
    

  for (int c=0; c<ncat; ++c) { //per ogni categoria

    for(imass = 0; imass<5;imass++){//guardo tutte le masse
      m0[imass]= ws[imass]->var(TString::Format("PhotonsMass_sig_m0_cat%d",c))->getVal()/masses[imass];
      //  m1[imass]= ws[imass]->var(TString::Format("PhotonsMass_sig_m1_cat%d",c))->getVal()/masses[imass];
      sigma0[imass]= ws[imass]->var(TString::Format("PhotonsMass_sig_sigma0_cat%d",c))->getVal()/masses[imass];
      sigma1[imass]= ws[imass]->var(TString::Format("PhotonsMass_sig_sigma1_cat%d",c))->getVal()/masses[imass];
      frac[imass]= ws[imass]->var(TString::Format("PhotonsMass_sig_frac_cat%d",c))->getVal();
      //sig[imass]= ws[imass]->var(TString::Format("nsig_cat%d",c))->getVal();
      m0_err[imass]= ws[imass]->var(TString::Format("PhotonsMass_sig_m0_cat%d",c))->getError()/masses[imass];
      //  m1_err[imass]= ws[imass]->var(TString::Format("PhotonsMass_sig_m1_cat%d",c))->getError()/masses[imass];
      sigma0_err[imass]= ws[imass]->var(TString::Format("PhotonsMass_sig_sigma0_cat%d",c))->getError()/masses[imass];
      sigma1_err[imass]= ws[imass]->var(TString::Format("PhotonsMass_sig_sigma1_cat%d",c))->getError()/masses[imass];
      frac_err[imass]= ws[imass]->var(TString::Format("PhotonsMass_sig_frac_cat%d",c))->getError();  
      é//sig_err[imass]= ws[imass]->var(TString::Format("nsig_cat%d",c))->getError();   
      std::cout<<masses[imass]<<"    "<<m0[imass]<<"    "<<m0_err[imass]<<std::endl;
  }

 

  
    label_cat[c] = new TPaveText(0.6, 0.8, 0.7, 0.9, "brNDC" );
    label_cat[c]->SetFillColor(kWhite);
    label_cat[c]->SetBorderSize(1.5);
    label_cat[c]->SetTextSize(0.038);
    label_cat[c]->SetTextAlign(11);
    label_cat[c]->SetTextFont(42);
    label_cat[c]->AddText(TString::Format("Cat %d", c));

    TGraphErrors* m0Graph = new TGraphErrors(5,masses, m0, masses_err, m0_err);
    m0Graph->Draw("APE");
    m0Graph->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    m0Graph->GetYaxis()->SetRangeUser(0.8, 1.2);
    m0Graph->GetYaxis()->SetTitle("m0 / m_{#gamma#gamma}");
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    label_cat[c]->Draw("same");
    c1->SaveAs(TString::Format("plots/m0vsmasses_%d.png", c));

    /*   TGraphErrors* m1Graph = new TGraphErrors(5,masses, m1, masses_err, m1_err);
    m1Graph->Draw("APE");
    m1Graph->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    m1Graph->GetYaxis()->SetTitle("m1 [GeV]");
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    label_cat[c]->Draw("same");
    c1->SaveAs(TString::Format("plots/m1vsmasses_%d.png", c));*/

    TGraphErrors* sigma0Graph = new TGraphErrors(5,masses, sigma0, masses_err, sigma0_err);
    sigma0Graph->Draw("APE");
    sigma0Graph->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    sigma0Graph->GetYaxis()->SetTitle("#sigma0 / m_{#gamma#gamma}");
    sigma0Graph->GetYaxis()->SetRangeUser(0., 0.07);
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    label_cat[c]->Draw("same");
    c1->SaveAs(TString::Format("plots/sigma0vsmasses_%d.png", c));

    TGraphErrors* sigma1Graph = new TGraphErrors(5,masses, sigma1, masses_err, sigma1_err);
    sigma1Graph->Draw("APE");
    sigma1Graph->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    sigma1Graph->GetYaxis()->SetTitle("#sigma1 / m_{#gamma#gamma}");
    sigma1Graph->GetYaxis()->SetRangeUser(0., 0.07);
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    label_cat[c]->Draw("same");
    c1->SaveAs(TString::Format("plots/sigma1vsmasses_%d.png", c));


    TGraphErrors* fracGraph = new TGraphErrors(5,masses, frac, masses_err, frac_err);
    fracGraph->Draw("APE");
    fracGraph->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    fracGraph->GetYaxis()->SetTitle("%Gauss0 ");
    fracGraph->GetYaxis()->SetRangeUser(0., 2.);
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    label_cat[c]->Draw("same");
    c1->SaveAs(TString::Format("plots/fracvsmasses_%d.png", c));

    /*   TGraphErrors* sigGraph = new TGraphErrors(5,masses, sig, masses_err, sig_err);
    sigGraph->Draw("APE");
    sigGraph->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    sigGraph->GetYaxis()->SetTitle("Signal yield ");
    sigGraph->GetYaxis()->SetRangeUser(0.001, 100.);
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    label_cat[c]->Draw("same");
    c1->SetLogy();
    c1->SaveAs(TString::Format("plots/sigYieldvsmasses_%d.png", c));*/




}




}






void MakeParameterTrendvsMassCBC(){
  gSystem->SetIncludePath("-I/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.32.00-cms5/include");
  gSystem->Load("libRooFit");
  gROOT->ProcessLine(".x /afs/cern.ch/work/s/soffi/CMSSW_5_2_3-globe/src/HiggsAnalysis/HiggsTo2photons/h2gglobe/ChiaraFitLimits/RooCBCrujffPdf.cxx+");
  TString wsDir = "workspaces/"+filePOSTfix;
  Int_t ncat = NCAT;
  Float_t masses[5];
  masses[0] = 150.0;
  masses[1] = 200.0;
  masses[2] = 250.0;
  masses[3] = 300.0;
  masses[4] = 400.0.;

  Float_t masses_err[5];
  masses_err[0] = 0.;
  masses_err[1] = 0.;
  masses_err[2] = 0.;
  masses_err[3] = 0.;
  masses_err[4] = 0.;


  TFile* file[5];
  for(imass = 0; imass<5;imass++) file[imass]= new TFile(wsDir+TString::Format("HighMass-hgg.m%.1f_8TeV.inputsig.root",masses[imass]));

  RooWorkspace* ws[5];
  for(imass = 0; imass<5;imass++) ws[imass] = (RooWorkspace*) file[imass]->Get("w_all");

 
  Float_t mean[5];
  Float_t sigma[5];
  Float_t alphaC[5];
  Float_t alphaCB[5];
  Float_t n[5];
  Float_t sig[5];

  Float_t mean_err[5];
  Float_t sigma_err[5];
  Float_t alphaC_err[5];
  Float_t alphaCB_err[5];
  Float_t n_err[5];
  Float_t sig_err[5];
  
  

  TCanvas* c1 = new TCanvas("c1", "c1", 1);
  c1->cd();
  TPaveText* label_cat[5];
  TPaveText* label_cms = get_labelCMS(0, "2012", true);
  TPaveText* label_sqrt = get_labelSqrt(0);
    

  for (int c=0; c<ncat; ++c) { //per ogni categoria

    for(imass = 0; imass<5;imass++){//guardo tutte le masse
      mean[imass]= ws[imass]->var(TString::Format("PhotonsMass_sig_mean_cat%d",c))->getVal()/masses[imass];
      sigma[imass]= ws[imass]->var(TString::Format("PhotonsMass_sig_sigma_cat%d",c))->getVal()/masses[imass];
      alphaC[imass]= ws[imass]->var(TString::Format("PhotonsMass_sig_alphaC_cat%d",c))->getVal();
      alphaCB[imass]= ws[imass]->var(TString::Format("PhotonsMass_sig_alphaCB_cat%d",c))->getVal();
      n[imass]= ws[imass]->var(TString::Format("PhotonsMass_sig_n_cat%d",c))->getVal();
      //  sig[imass]= ws[imass]->var(TString::Format("nsigCBC_cat%d",c))->getVal();
      mean_err[imass]= ws[imass]->var(TString::Format("PhotonsMass_sig_mean_cat%d",c))->getError()/masses[imass];
      sigma_err[imass]= ws[imass]->var(TString::Format("PhotonsMass_sig_sigma_cat%d",c))->getError()/masses[imass];
      alphaC_err[imass]= ws[imass]->var(TString::Format("PhotonsMass_sig_alphaC_cat%d",c))->getError();
      alphaCB_err[imass]= ws[imass]->var(TString::Format("PhotonsMass_sig_alphaCB_cat%d",c))->getError();
      n_err[imass]= ws[imass]->var(TString::Format("PhotonsMass_sig_n_cat%d",c))->getError(); 
      //  sig_err[imass]= ws[imass]->var(TString::Format("nsigCBC_cat%d",c))->getError();    
      std::cout<<masses[imass]<<"    "<<mean[imass]<<"    "<<mean_err[imass]<<std::endl;
  }

 

  
    label_cat[c] = new TPaveText(0.6, 0.8, 0.7, 0.9, "brNDC" );
    label_cat[c]->SetFillColor(kWhite);
    label_cat[c]->SetBorderSize(1.5);
    label_cat[c]->SetTextSize(0.038);
    label_cat[c]->SetTextAlign(11);
    label_cat[c]->SetTextFont(42);
    label_cat[c]->AddText(TString::Format("Cat %d", c));

    TGraphErrors* meanGraph = new TGraphErrors(5,masses, mean, masses_err, mean_err);
    meanGraph->Draw("APE");
    meanGraph->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    meanGraph->GetYaxis()->SetRangeUser(0.8, 1.2);
    meanGraph->GetYaxis()->SetTitle("mean / m_{#gamma#gamma}");
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    label_cat[c]->Draw("same");
    c1->SaveAs(TString::Format("plots/meanvsmasses_%d.png", c));

    TGraphErrors* sigmaGraph = new TGraphErrors(5,masses, sigma, masses_err, sigma_err);
    sigmaGraph->Draw("APE");
    sigmaGraph->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    sigmaGraph->GetYaxis()->SetTitle("sigma / m_{#gamma#gamma}");
    sigmaGraph->GetYaxis()->SetRangeUser(0.,0.1);
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    label_cat[c]->Draw("same");
    c1->SaveAs(TString::Format("plots/sigmavsmasses_%d.png", c));

    TGraphErrors* alphaCGraph = new TGraphErrors(5,masses, alphaC, masses_err, alphaC_err);
    alphaCGraph->Draw("APE");
    alphaCGraph->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    alphaCGraph->GetYaxis()->SetTitle("#alpha_{C}");
    alphaCGraph->GetYaxis()->SetRangeUser(0., 0.4);
    //Graph->GetYaxis()->SetRangeUser(0., 12.);
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    label_cat[c]->Draw("same");
    c1->SaveAs(TString::Format("plots/alphaCvsmasses_%d.png", c));

    TGraphErrors* alphaCBGraph = new TGraphErrors(5,masses, alphaCB, masses_err, alphaCB_err);
    alphaCBGraph->Draw("APE");
    alphaCBGraph->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    alphaCBGraph->GetYaxis()->SetTitle("#alpha_{CB}");
    alphaCBGraph->GetYaxis()->SetRangeUser(0., 3.);
    //Graph->GetYaxis()->SetRangeUser(0., 20.);
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    label_cat[c]->Draw("same");
    c1->SaveAs(TString::Format("plots/alphaCBvsmasses_%d.png", c));


    TGraphErrors* nGraph = new TGraphErrors(5,masses, n, masses_err, n_err);
    nGraph->Draw("APE");
    nGraph->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    nGraph->GetYaxis()->SetTitle("n ");
    nGraph->GetYaxis()->SetRangeUser(0., 100.);
    //nGraph->GetYaxis()->SetRangeUser(0., 2.);
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    label_cat[c]->Draw("same");
    c1->SaveAs(TString::Format("plots/nvsmasses_%d.png", c));

    /*
    TGraphErrors* sigGraph = new TGraphErrors(5,masses, sig, masses_err, sig_err);
    sigGraph->Draw("APE");
    sigGraph->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    sigGraph->GetYaxis()->SetTitle("Signal yield ");
    sigGraph->GetYaxis()->SetRangeUser(0.001, 100.);
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    label_cat[c]->Draw("same");
    c1->SetLogy();
    c1->SaveAs(TString::Format("plots/sigYieldvsmasses_%d.png", c));
    
    */


}




}



void PlotSigShape(RooWorkspace* w, Int_t c, Float_t* mean, Float_t* rms,Float_t* meanErr, Float_t* rmsErr) {

  Int_t nmass = 7;
  Int_t masses[7] = {150, 200, 250, 300, 400, 600, 800};

  TString inDir = "";

  // Variables
  RooArgSet* ntplVars = defineVariables();
  
  RooRealVar* PhotonsMass = w->var("PhotonsMass");
  //  TFile sigFile1("histograms_CMS-HGG_11032014_MC.root");   //ggh prod mode tree livia
  TFile sigFile1("histograms_CMS-HGG_08052014_MC.root");   //ggh prod mode tree livia
  //  TFile sigFile2("histograms_CMS-HGG_HighMass_600800.root"); //high mass
  // TFile sigFile2("histograms_CMS-HGG_SM_RD_M150.root"); 
  TCanvas* c1 = new TCanvas("c1", "c1", 1);
  c1->cd();
  RooPlot* plot;


  RooRealVar* rmsvar = new RooRealVar("rmsvar", "rmsvar", 0.2);
  RooRealVar* meanvar = new RooRealVar("meanvar", "meanvar", 0.);
  TChain* sigTree1;
  RooDataSet* signal;
 
  for(int iMass=0;iMass<7;iMass++){ //per ogni massa

    
  //chain summing up all production modes
    sigTree1  = new TChain();
    
 
    sigTree1->Add(TString::Format("histograms_CMS-HGG_08052014_MC.root/ggh_m%d_8TeV", masses[iMass]));
    /*  sigTree1->Add(TString::Format("histograms_CMS-HGG_17042014_MC.root/vbf_m%d_8TeV", masses[iMass]));
    sigTree1->Add(TString::Format("histograms_CMS-HGG_17042014_MC.root/wzh_m%d_8TeV", masses[iMass]));
    sigTree1->Add(TString::Format("histograms_CMS-HGG_17042014_MC.root/tth_m%d_8TeV", masses[iMass]));*/


    sigTree1->SetTitle("sigTree1");
    sigTree1->SetName("sigTree1");
  // common preselection cut
  TString mainCut = TString::Format("PhotonsMass>=(%d*0.7) && PhotonsMass<=(%d*1.3)", masses[iMass], masses[iMass]);   // livia

  RooDataSet sigWeighted("sigWeighted","dataset",sigTree1,*ntplVars,mainCut,"evweight");
  cout << endl;
  cout << "sigWeighted" << endl;
  sigWeighted.Print("v");
  cout << "---- nX:  " << sigWeighted.sumEntries() << endl; 
  
  // apply a common preselection cut; split in categories
  cout << endl;
  
  RooFormulaVar *massReduced_formula     = new RooFormulaVar("massReduced_formula","","@0/@1-1",RooArgList(*w->var("PhotonsMass"),*w->var("PhotonsMassTrue")));

  RooRealVar* massReduced = (RooRealVar*) sigWeighted.addColumn(*massReduced_formula);
  massReduced->SetName("massReduced");
  massReduced->SetTitle("massReduced");
  if(iMass==0) w->import(*massReduced);   

  // 1)  prime 4 cat livia
  if (c==0) signal = (RooDataSet*) sigWeighted.reduce(*w->var("massReduced"),mainCut+TString::Format("&& (abs(ph1_eta)<1.5 && abs(ph2_eta)<1.5) && (ph1_r9>0.94 && ph2_r9>0.94 )"));
  if (c==1) signal = (RooDataSet*) sigWeighted.reduce(*w->var("massReduced"),mainCut+TString::Format("&& (abs(ph1_eta)<1.5 && abs(ph2_eta)<1.5) && (ph1_r9<0.94 || ph2_r9<0.94 ) "));
  if (c==2) signal = (RooDataSet*) sigWeighted.reduce(*w->var("massReduced"),mainCut+TString::Format("&& (abs(ph1_eta)>1.5 || abs(ph2_eta)>1.5) && (ph1_r9>0.94 && ph2_r9>0.94 )"));
  if (c==3) signal = (RooDataSet*) sigWeighted.reduce(*w->var("massReduced"),mainCut+TString::Format("&& (abs(ph1_eta)>1.5 || abs(ph2_eta)>1.5) && (ph1_r9<0.94 || ph2_r9<0.94 ) "));
  
  
  //fit the resp fcn
  RooCBShape* ResponseCBpos;
  RooCBShape* ResponseCBneg;
  RooAddPdf* ResponseAdd;
  //cb pos                                                                                                                     
   RooFormulaVar CBpos_mean(TString::Format("ReducedMass_CBpos_sig_mean_cat%d",c),"", "@0", *w->var(TString::Format("ReducedMass_sig_mean_cat%d",c)));
   RooFormulaVar CBpos_sigma(TString::Format("ReducedMass_CBpos_sig_sigma_cat%d",c), "", "@0", *w->var(TString::Format("ReducedMass_sig_sigma_cat%d",c)));
   RooFormulaVar CBpos_alphaCB(TString::Format("ReducedMass_CBpos_sig_alphaCB_cat%d",c),"", "@0", *w->var( TString::Format("ReducedMass_sig_alphaCBpos_cat%d",c)));
   RooFormulaVar CBpos_n(TString::Format("ReducedMass_CBpos_sig_n_cat%d",c),"", "@0", *w->var( TString::Format("ReducedMass_sig_Npos_cat%d",c)));
   //cb neg
   RooFormulaVar CBneg_n(TString::Format("ReducedMass_CBneg_sig_n_cat%d",c),"", "@0", *w->var( TString::Format("ReducedMass_sig_Nneg_cat%d",c)));
   RooFormulaVar CBneg_alphaCB(TString::Format("ReducedMass_CBneg_sig_alphaCB_cat%d",c),"", "@0", *w->var( TString::Format("ReducedMass_sig_alphaCBneg_cat%d",c)));
   ResponseCBpos =  new RooCBShape(TString::Format("ResponseCBpos_cat%d",c),TString::Format("ResponseCBpos_cat%d",c) , *w->var("massReduced"), CBpos_mean, CBpos_sigma,  CBpos_alphaCB, CBpos_n) ;
  
   ResponseCBneg =  new RooCBShape(TString::Format("ResponseCBneg_cat%d",c),TString::Format("ResponseCBneg_cat%d",c) , *w->var("massReduced"), CBpos_mean, CBpos_sigma,  CBneg_alphaCB, CBneg_n) ;
   
   
   RooFormulaVar CB_frac(TString::Format("ReducedMass_CBpos_sig_frac_cat%d",c),"", "@0", *w->var(TString::Format("ReducedMass_sig_frac_cat%d",c)));
   w->import(CB_frac); 

   ResponseAdd= new RooAddPdf(TString::Format("ResponseAddPdf_cat%d",c),TString::Format("ResponseAddPdf_cat%d",c) , RooArgList(*ResponseCBneg, *ResponseCBpos), CB_frac);
   RooArgSet* params;
   params = ResponseAdd->getParameters(*massReduced);
   params->readFromFile("/afs/cern.ch/work/s/soffi/CMSSW611-Analysis/src/h2gglobe/ChiaraFitLimits/initParSigM150.txt");
   params->writeToStream(std::cout, kFALSE);
  
  RooFitResult* fitresults = (RooFitResult* ) ResponseAdd->fitTo(*signal,SumW2Error(kTRUE),Range(-1, 1),  RooFit::Save(kTRUE));
   //  std::cout<<TString::Format("******************************** Signal Fit results CB+Gauss  mass %f cat %d***********************************", mass, c)<<std::endl;
   std::cout<<TString::Format("******************************** Signal Fit results CB+CB  mass %f cat %d***********************************", iMass, c)<<std::endl;
    fitresults->Print("V");
    params->writeToStream(std::cout, kFALSE);


    mean[iMass] = CBpos_mean->getVal();
    rms[iMass] = CBpos_sigma->getVal();
    meanErr[iMass] = CBpos_mean->getPropagatedError(*fitresults);
    rmsErr[iMass] = CBpos_sigma->getPropagatedError(*fitresults);
    /*  alphaPOS[iMass] = CBpos_alphaCB->getVal();
    alphaPOSErr[iMass] = CBpos_alphaCB->getPropagatedError(*fitresults);
    alphaNEG[iMass] = CBneg_alphaCB->getVal();
    alphaNEGErr[iMass] = CBneg_alphaCB->getPropagatedError(*fitresults);
    nPOS[iMass] = CBpos_n->getVal();
    nPOSErr[iMass] = CBpos_n->getPropagatedError(*fitresults);
    nNEG[iMass] = CBneg_n->getVal();
    nNEGErr[iMass] = CBneg_n->getPropagatedError(*fitresults);
    frac[iMass] = CB_frac->getVal();
    fracErr[iMass] = CB_frac->getPropagatedError(*fitresults);
    */
 
  if(iMass==0) plot = massReduced->frame(Range(-0.12, 0.12), Bins(20), Title("Mass Reduced"));

  if(iMass==0) signal->plotOn(plot, DrawOption("C"), MarkerSize(0), MarkerColor(kWhite), LineWidth(2),LineColor(kBlack), XErrorSize(0), DataError(RooAbsData::None), Rescale(1./signal->sumEntries()));
   if(iMass==2) signal->plotOn(plot, DrawOption("C"), MarkerSize(0), MarkerColor(kWhite), LineWidth(2),LineColor(kPink-8), XErrorSize(0), DataError(RooAbsData::None), Rescale(1./signal->sumEntries()));
   if(iMass==4) signal->plotOn(plot, DrawOption("C"), MarkerSize(0), MarkerColor(kWhite), LineWidth(2),LineColor(kBlue+3), XErrorSize(0), DataError(RooAbsData::None), Rescale(1./signal->sumEntries()));
   if(iMass==5) signal->plotOn(plot, DrawOption("C"), MarkerSize(0), MarkerColor(kWhite), LineWidth(2),LineColor(kOrange-3), XErrorSize(0), DataError(RooAbsData::None), Rescale(1./signal->sumEntries()));
   if(iMass==6) signal->plotOn(plot, DrawOption("C"), MarkerSize(0), MarkerColor(kWhite), LineWidth(2),LineColor(kSpring-8), XErrorSize(0), DataError(RooAbsData::None), Rescale(1./signal->sumEntries()));

  if(iMass==0) plot->Draw();
  else plot->Draw("same");

  /*  meanvar = signal->meanVar(*w->var("massReduced"));
  mean[iMass] = meanvar->getVal();
  meanErr[iMass] = meanvar->getError();
  rmsvar = signal->rmsVar(*w->var("massReduced"));
  rms[iMass] = rmsvar->getVal();
  rmsErr[iMass] = rmsvar->getError();*/

  }
  plot->SetAxisRange(0.0001,0.8,"Y");
  plot->GetXaxis()->SetTitle("#frac{#Delta m}{m} ");
  plot->GetYaxis()->SetTitle("Normalized To Unity");
  plot->GetXaxis()->SetTitleFont(42);
  plot->GetXaxis()->SetTitleSize(0.05);

  TLegend* legmc = new TLegend(0.62, 0.6, 0.87, 0.89, "", "brNDC");
  legmc->AddEntry(plot->getObject(0),"m_{X} = 150 GeV","L");
  legmc->AddEntry(plot->getObject(1),"m_{X} = 250 GeV","L");
  legmc->AddEntry(plot->getObject(2),"m_{X} = 400 GeV","L");
  legmc->AddEntry(plot->getObject(3),"m_{X} = 600 GeV","L");
  legmc->AddEntry(plot->getObject(4),"m_{X} = 800 GeV","L");
  legmc->SetTextSize(0.0206044);  
  legmc->SetTextFont(42);
  legmc->SetBorderSize(0);
  legmc->SetFillStyle(0);
  legmc->Draw(); 

  TLatex *lat  = new TLatex(0.65,0.9,TString::Format("Cat: %d", c));
  
  lat->SetTextSize(0.038);
  lat->SetTextAlign(11);
  lat->SetTextFont(42);
  lat->Draw("same");
  lat->SetNDC();



 

  TPaveText* label_cms = get_labelCMS(0, "2012", true);
  TPaveText* label_sqrt = get_labelSqrt(0);
  label_cms->Draw("same");
  label_sqrt->Draw("same");

 
  c1->SetLogy(0);
  c1->SaveAs(TString::Format("plots/massesShape_cat%d.png",c));
  c1->SaveAs(TString::Format("plots/massesShape_cat%d.pdf",c));
  c1->SaveAs(TString::Format("~/www/plotsNota/massesShape_cat%d.png",c));
  c1->SaveAs(TString::Format("~/www/plotsNota/massesShape_cat%d.pdf",c));
  if(c==0){
    c1->SaveAs(TString::Format("~/www/plotsPAS/massesShape_cat%d.png",c));
    c1->SaveAs(TString::Format("~/www/plotsPAS/massesShape_cat%d.pdf",c));
    c1->SaveAs(TString::Format("~/www/plotsPAS/massesShape_cat%d.C",c));
    c1->SaveAs(TString::Format("~/www/plotsPAS/massesShape_cat%d.root",c));
  }
  plot->SetAxisRange(0.0001,0.8,"Y");
  c1->SetLogy(1);
  c1->SaveAs(TString::Format("plots/massesShape_cat%d_LOG.png",c));
  c1->SaveAs(TString::Format("plots/massesShape_cat%d_LOG.pdf",c));
  if(c==0){
  c1->SaveAs(TString::Format("~/www/plotsPAS/massesShape_cat%d_LOG.png",c));
  c1->SaveAs(TString::Format("~/www/plotsPAS/massesShape_cat%d_LOG.pdf",c));
  c1->SaveAs(TString::Format("~/www/plotsPAS/massesShape_cat%d_LOG.C",c));
  c1->SaveAs(TString::Format("~/www/plotsPAS/massesShape_cat%d_LOG.root",c));
  }
  /*  std::cout<<" MEAN  "<<std::endl;
  for (int i = 0; i< 7; i++)std::cout<< mean[i]<<"     ";
  std::cout<<"   "<<std::endl;
  for (int i = 0; i< 7; i++)std::cout<< meanErr[i]<<"     ";
  std::cout<<" RMS  "<<std::endl;
  for (int i = 0; i< 7; i++)std::cout<< rms[i]<<"     ";
  std::cout<<"   "<<std::endl;
  for (int i = 0; i< 7; i++)std::cout<< rmsErr[i]<<"     ";
  std::cout<<" a POS  "<<std::endl;
  for (int i = 0; i< 7; i++)std::cout<< alphaPOS[i]<<"     ";
  std::cout<<"   "<<std::endl;
  for (int i = 0; i< 7; i++)std::cout<< alphaPOSErr[i]<<"     ";
  std::cout<<"  a NEG "<<std::endl;
  for (int i = 0; i< 7; i++)std::cout<< alphaNEG[i]<<"     ";
  std::cout<<"   "<<std::endl;
  for (int i = 0; i< 7; i++)std::cout<< alphaNEGErr[i]<<"     ";
  std::cout<<"  n POS "<<std::endl;
  for (int i = 0; i< 7; i++)std::cout<< nPOS[i]<<"     ";
  std::cout<<"   "<<std::endl;
  for (int i = 0; i< 7; i++)std::cout<< nPOSErr[i]<<"     ";
  std::cout<<"  n NEG "<<std::endl;
  for (int i = 0; i< 7; i++)std::cout<< nNEG[i]<<"     ";
  std::cout<<"   "<<std::endl;
  for (int i = 0; i< 7; i++)std::cout<< nNEGErr[i]<<"     ";
  std::cout<<"  frac "<<std::endl;
  for (int i = 0; i< 7; i++)std::cout<< frac[i]<<"     ";
  std::cout<<"   "<<std::endl;
  for (int i = 0; i< 7; i++)std::cout<< fracErr[i]<<"     ";
  std::cout<<"   "<<std::endl;
 */
}







void PlotSigShapeCorrected(RooWorkspace* w, Int_t c) {

  Int_t nmass = 7;
  Int_t masses[7] = {150, 200, 250, 300, 400, 600, 800};

  TString inDir = "";

  // Variables
  RooArgSet* ntplVars = defineVariables();
  
  RooRealVar* PhotonsMass = w->var("PhotonsMass");
  TFile sigFile1("histograms_CMS-HGG_08052014_MC.root");   //ggh prod mode tree livia
  
 


  RooRealVar* rmsvar = new RooRealVar("rmsvar", "rmsvar", 0.2);
  RooRealVar* meanvar = new RooRealVar("meanvar", "meanvar", 0.);
  TChain* sigTree1;
  RooDataSet* signal;

  TH1F* h[7];
  RooPlot* plot;
  TFile* f = new TFile("sigShapeCorrections.root", "READ");
  TF1* fcn = (TF1*)f->Get(TString::Format("f%d",c));

  for(int iMass=0;iMass<7;iMass++){ //per ogni massa
    // if(iMass==1 || iMass ==3)continue;
  //chain summing up all production modes
    sigTree1  = new TChain();
    
 
    sigTree1->Add(TString::Format("histograms_CMS-HGG_08052014_MC.root/ggh_m%d_8TeV", masses[iMass]));
    /*   sigTree1->Add(TString::Format("histograms_CMS-HGG_17042014_MC.root/vbf_m%d_8TeV", masses[iMass]));
    sigTree1->Add(TString::Format("histograms_CMS-HGG_17042014_MC.root/wzh_m%d_8TeV", masses[iMass]));
    sigTree1->Add(TString::Format("histograms_CMS-HGG_17042014_MC.root/tth_m%d_8TeV", masses[iMass]));*/
    

    sigTree1->SetTitle("sigTree1");
    sigTree1->SetName("sigTree1");
  // common preselection cut
  TString mainCut = TString::Format("PhotonsMass>=(%d*0.7) && PhotonsMass<=(%d*1.3)", masses[iMass], masses[iMass]);   // livia

  RooDataSet sigWeighted("sigWeighted","dataset",sigTree1,*ntplVars,mainCut,"evweight");
  cout << endl;
  cout << "sigWeighted" << endl;
  sigWeighted.Print("v");
  cout << "---- nX:  " << sigWeighted.sumEntries() << endl; 
  
  // apply a common preselection cut; split in categories
  cout << endl;
  

  Float_t corr = fcn->Eval(masses[iMass]);
  if(iMass==0)corr =1;
  std::cout<<"Mass: "<<masses[iMass]<<" corr: "<<corr<<std::endl;
  RooRealVar rooCorr ("rooCorr", "rooCorr", corr, "");
  rooCorr->setConstant();

  RooFormulaVar *massReduced_formula     = new RooFormulaVar("massReduced_formula","","(@0/@1-1)/@2",RooArgList(*w->var("PhotonsMass"),*w->var("PhotonsMassTrue"), rooCorr));

  RooRealVar* massReduced = (RooRealVar*) sigWeighted.addColumn(*massReduced_formula);
  massReduced->SetName("massReduced");
  massReduced->SetTitle("massReduced");
  if(iMass==0) w->import(*massReduced);   

  // 1)  prime 4 cat livia
  if (c==0) signal = (RooDataSet*) sigWeighted.reduce(*w->var("massReduced"),mainCut+TString::Format("&& (abs(ph1_eta)<1.5 && abs(ph2_eta)<1.5) && (ph1_r9>0.94 && ph2_r9>0.94 )"));
  if (c==1) signal = (RooDataSet*) sigWeighted.reduce(*w->var("massReduced"),mainCut+TString::Format("&& (abs(ph1_eta)<1.5 && abs(ph2_eta)<1.5) && (ph1_r9<0.94 || ph2_r9<0.94 ) "));
  if (c==2) signal = (RooDataSet*) sigWeighted.reduce(*w->var("massReduced"),mainCut+TString::Format("&& (abs(ph1_eta)>1.5 || abs(ph2_eta)>1.5) && (ph1_r9>0.94 && ph2_r9>0.94 )"));
  if (c==3) signal = (RooDataSet*) sigWeighted.reduce(*w->var("massReduced"),mainCut+TString::Format("&& (abs(ph1_eta)>1.5 || abs(ph2_eta)>1.5) && (ph1_r9<0.94 || ph2_r9<0.94 ) "));
  


  h[iMass]=new TH1F(TString::Format("h_M%d_CAT%d", masses[iMass], c), TString::Format("h_M%d_CAT%d", masses[iMass], c), 20, -0.12, 0.12);
  // h[iMass]= (TH1F*)signal->createHistogram("Hist",*w->var("massReduced"),Binning(60));
  h[iMass]->SetName(TString::Format("h_M%d_CAT%d", masses[iMass], c));
  h[iMass]->Sumw2();

  std::cout<<"("+mainCut+"&& (abs(ph1_eta)<1.5 && abs(ph2_eta)<1.5) && (ph1_r9>0.94 && ph2_r9>0.94 ))*evweight"<<std::endl;
  if (c==0)sigTree1->Draw(TString::Format("(PhotonsMass/PhotonsMassTrue-1)/%f>>h_M%d_CAT0",corr, masses[iMass]),"("+mainCut+TString::Format("&& (abs(ph1_eta)<1.5 && abs(ph2_eta)<1.5) && (ph1_r9>0.94 && ph2_r9>0.94 ))*evweight") );
  if (c==1)sigTree1->Draw(TString::Format("(PhotonsMass/PhotonsMassTrue-1)/%f>>h_M%d_CAT1",corr, masses[iMass]),"("+mainCut+TString::Format("&& (abs(ph1_eta)<1.5 && abs(ph2_eta)<1.5) && (ph1_r9<0.94 || ph2_r9<0.94 ) )*evweight") );
  if (c==2)sigTree1->Draw(TString::Format("(PhotonsMass/PhotonsMassTrue-1)/%f>>h_M%d_CAT2",corr, masses[iMass]),"("+mainCut+TString::Format("&& (abs(ph1_eta)>1.5 || abs(ph2_eta)>1.5) && (ph1_r9>0.94 && ph2_r9>0.94 ))*evweight") );
  if (c==3)sigTree1->Draw(TString::Format("(PhotonsMass/PhotonsMassTrue-1)/%f>>h_M%d_CAT3",corr, masses[iMass]),"("+mainCut+TString::Format("&& (abs(ph1_eta)>1.5 || abs(ph2_eta)>1.5) && (ph1_r9<0.94 || ph2_r9<0.94 ))*evweight") );


  std::cout<<signal->sumEntries()<<"   "<<h[iMass]->Integral()<<std::endl;


  if(iMass==0) plot = massReduced->frame(Range(-0.12, 0.12), Bins(20), Title("Mass Reduced"));

  if(iMass==0) signal->plotOn(plot, DrawOption("C"), MarkerSize(0), MarkerColor(kWhite), LineWidth(2),LineColor(kBlack), XErrorSize(0), DataError(RooAbsData::None), Rescale(1./signal->sumEntries()));  
  if(iMass==2) signal->plotOn(plot, DrawOption("C"), MarkerSize(0), MarkerColor(kWhite), LineWidth(2),LineColor(kPink-8), XErrorSize(0), DataError(RooAbsData::None), Rescale(1./signal->sumEntries()));
  if(iMass==4) signal->plotOn(plot, DrawOption("C"), MarkerSize(0), MarkerColor(kWhite), LineWidth(2),LineColor(kBlue+3), XErrorSize(0), DataError(RooAbsData::None), Rescale(1./signal->sumEntries()));
  if(iMass==5) signal->plotOn(plot, DrawOption("C"), MarkerSize(0), MarkerColor(kWhite), LineWidth(2),LineColor(kOrange-3), XErrorSize(0), DataError(RooAbsData::None), Rescale(1./signal->sumEntries()));
  if(iMass==6) signal->plotOn(plot, DrawOption("C"), MarkerSize(0), MarkerColor(kWhite), LineWidth(2),LineColor(kSpring-8), XErrorSize(0), DataError(RooAbsData::None), Rescale(1./signal->sumEntries()));


 

  }

  
  TCanvas* c1 = new TCanvas("c1", "c1", 1);
  c1->cd();
  //-------pad 1-------//
  TPad * pad1 = new TPad("pad1", "pad1",0.01,0.13,0.75,1.);  
  
    
  pad1->SetRightMargin(0.1);
  
  pad1->SetLogy();
  pad1->Draw();
  pad1->cd();
  
  pad1->Clear();
  
  plot->Draw();
  plot->SetAxisRange(0.0001,0.4,"Y");
  // plot->GetXaxis()->SetTitle("#frac{#Delta m}{m} ");
  plot->GetYaxis()->SetTitle("Normalized To Unity");
  plot->GetXaxis()->SetTitleFont(42);
  plot->GetXaxis()->SetTitleSize(0.05);

  h[0]->Scale(1./h[0]->Integral());
  h[2]->Scale(1./h[2]->Integral());
  h[4]->Scale(1./h[4]->Integral());
  h[5]->Scale(1./h[5]->Integral());
  h[6]->Scale(1./h[6]->Integral());

  h[0]->SetLineColor(kBlack);
  h[2]->SetLineColor(kPink-8);
  h[4]->SetLineColor(kBlue+3);
  h[5]->SetLineColor(kOrange-3);
  h[6]->SetLineColor(kSpring-8);
  /*  h[0]->SetFillColor(kBlack);
  h[2]->SetFillColor(kPink-8);
  h[4]->SetFillColor(kBlue+3);
  h[5]->SetFillColor(kOrange-3);
  h[6]->SetFillColor(kSpring-8);
*/

  TLegend* legmc = new TLegend(0.62, 0.6, 0.87, 0.89, "", "brNDC");
  legmc->AddEntry(plot->getObject(0),"m_{#gamma#gamma} = 150 GeV","L");
  legmc->AddEntry(plot->getObject(1),"m_{#gamma#gamma} = 250 GeV","L");
  legmc->AddEntry(plot->getObject(2),"m_{#gamma#gamma} = 400 GeV","L");
  legmc->AddEntry(plot->getObject(3),"m_{#gamma#gamma} = 600 GeV","L");
  legmc->AddEntry(plot->getObject(4),"m_{#gamma#gamma} = 800 GeV","L");
  legmc->SetTextSize(0.0206044);  
  legmc->SetTextFont(42);
  legmc->SetBorderSize(0);
  legmc->SetFillStyle(0);
  legmc->Draw(); 

  TLatex *lat  = new TLatex(0.65,0.9,TString::Format("Cat: %d", c));
  
  lat->SetTextSize(0.038);
  lat->SetTextAlign(11);
  lat->SetTextFont(42);
  lat->Draw("same");
  lat->SetNDC();



  TPaveText* label_cms = get_labelCMS(0, "2012", true);
  TPaveText* label_sqrt = get_labelSqrt(0);
  label_cms->Draw("same");
  label_sqrt->Draw("same");

  c1->cd();
  //-------pad 2------//
    TPad * pad2 = new TPad("pad2", "pad2",0.01,0.001,0.75,0.2);
    pad2->SetGrid();
    
    //pad2->SetTopMargin(0.01);
    pad2->SetBottomMargin(0.4);
    pad2->SetRightMargin(0.1);
    pad2->Draw();
    pad2->cd();

    h[0]->GetXaxis()->SetRangeUser(-0.12, 0.12);
    h[0]->GetXaxis()->SetTitle("#frac{#Delta m}{m}");
    h[0]->GetYaxis()->SetRangeUser(0., 2.);
    for(int i = 0; i<7; i++){
      std::cout<<" mass: "<<masses[i]<<std::endl;
      if(i==1 || i ==3)continue;
      TH1F* hratio = (TH1F*) h[i]->Clone();
      hratio->Divide(h[0]);
      for(int j = 0; j<20; j++){
	double error = sqrt(pow(1./h[0]->GetBinContent(j),2)*(h[i]->GetBinError(j))*(h[i]->GetBinError(j))+pow(h[i]->GetBinContent(j)/h[0]->GetBinContent(j)/h[0]->GetBinContent(j),2)*h[0]->GetBinError(j)*h[0]->GetBinError(j));
	hratio->SetBinError(j, error);
	std::cout<<"Mass: "<<h[i]->GetBinCenter(j)<<" i: "<<hratio->GetBinContent(j)<<" Err: "<<hratio->GetBinError(j)<<std::endl;
      }
      // hratio->SetFillStyle(3004);
      hratio->SetMarkerSize(0);
      if(i==0){	
	hratio->Draw("histE");
	hratio->GetYaxis()->SetRangeUser(-1., 2.);
	hratio->GetYaxis()->SetNdivisions(4,false);
	hratio->GetYaxis()->SetTitleFont(42);
	hratio->GetXaxis()->SetTitle("#frac{#Delta m}{m}");
	hratio->GetXaxis()->SetTitleSize(0.2);
	hratio->GetXaxis()->SetLabelSize(0.16);
	hratio->GetYaxis()->SetLabelSize(0.16);
	hratio->GetYaxis()->SetTitleSize(0.15);
	hratio->GetYaxis()->SetTitleOffset(0.65);
	hratio->GetXaxis()->SetTitleOffset(0.8);
      }
      else hratio->Draw("histEsame");
    
    }



 
  pad1->SetLogy(0);
  c1->SaveAs(TString::Format("plots/massesShapeCorrected_cat%d.png",c));
  c1->SaveAs(TString::Format("plots/massesShapeCorrected_cat%d.pdf",c));
  c1->SaveAs(TString::Format("~/www/plotsNota/massesShapeCorrected_cat%d.png",c));
  c1->SaveAs(TString::Format("~/www/plotsNota/massesShapeCorrected_cat%d.pdf",c));
 
  plot->SetAxisRange(0.0001,0.4,"Y");
  pad1->SetLogy(1);
  c1->SaveAs(TString::Format("plots/massesShapeCorrected_cat%d_LOG.png",c));
  c1->SaveAs(TString::Format("plots/massesShapeCorrected_cat%d_LOG.pdf",c));
  
  c1->SaveAs(TString::Format("~/www/plotsNota/massesShapeCorrected_cat%d_LOG.png",c));
  c1->SaveAs(TString::Format("~/www/plotsNota/massesShapeCorrected_cat%d_LOG.pdf",c));


}




void PlotSigShapeMCofficial(RooWorkspace* w, Int_t c) {

  // Variables
  RooArgSet* ntplVars1 = defineVariables();
  RooArgSet* ntplVars2 = defineVariables();
  
  RooRealVar* PhotonsMass = w->var("PhotonsMass");

  TCanvas* c1 = new TCanvas("c1", "c1", 1);
  c1->cd();
  RooPlot* plot1;
  RooPlot* plot2;
  /*  TH1F* h1 = new TH1F("h1", "h1", 60,-0.12, 0.12);
  TH1F* h2 = new TH1F("h2", "h2", 60,-0.12, 0.12);
  std::string var("PhotonsMass/150 -1");
  h1->GetXaxis()->SetTitle("#frac{#Delta m}{m} ");
  */
  
  /*  TH1F* h1 = new TH1F("h1", "h1", 40, 0, 40);
  TH1F* h2 = new TH1F("h2", "h2", 40, 0, 40 );
  std::string var("nvtx");
  h1->GetXaxis()->SetTitle("nvtx");
*/

  
  /* TH1F* h1 = new TH1F("h1", "h1", 80, 0, 400);
  TH1F* h2 = new TH1F("h2", "h2", 80, 0, 400 );
  std::string var("ph1_pt");
  h1->GetXaxis()->SetTitle("p_{T} [GeV]");*/

  /*  TH1F* h1 = new TH1F("h1", "h1", 60, 0, 1.2);
  TH1F* h2 = new TH1F("h2", "h2", 60, 0, 1.2 );
  std::string var("ph1_r9");
  h1->GetXaxis()->SetTitle("R9");*/
  
  TH1F* h1 = new TH1F("h1", "h1", 120, 100, 200);
  TH1F* h2 = new TH1F("h2", "h2", 120, 100, 200);
  std::string var("PhotonsMass");
  h1->GetXaxis()->SetTitle("m_{#gamma #gamma} [GeV]");

  TChain* sigTree1= new TChain();
  TChain* sigTree2= new TChain();
  
 
  sigTree1->Add("histograms_CMS-HGG_M150_Livia.root/ggh_m150_8TeV");
  // sigTree1->Add("histograms_CMS-HGG_HighMass_17042014.root/vbf_m150_8TeV");
  //sigTree1->Add("histograms_CMS-HGG_HighMass_17042014.root/wzh_m150_8TeV");
  //sigTree1->Add("histograms_CMS-HGG_HighMass_09042014.root/tth_m150_8TeV");
  sigTree1->SetTitle("sigTree1");
  sigTree1->SetName("sigTree1");

  sigTree2->Add("histograms_CMS-HGG_M150Moriond2013.root/ggh_m150_8TeV");
//  sigTree2->Add("histograms_CMS-HGG_Moriond_M150.root/vbf_m150_8TeV");
//sigTree2->Add("histograms_CMS-HGG_Moriond_M150.root/wzh_m150_8TeV");
//sigTree2->Add("histograms_CMS-HGG_SM_RD_M150.root/tth_m150_8TeV");
  sigTree2->SetTitle("sigTree2");
  sigTree2->SetName("sigTree2");

  std::cout<<sigTree2->GetEntries()<<std::endl;

  // common preselection cut
  //  && (ph1_r9>0.94 && ph2_r9>0.94 )
  // 1)  prime 4 cat livia
  if (c==0) sigTree1->Draw((var+">>h1").c_str(), "evweight*(PhotonsMass>0.7*150 && PhotonsMass<1.3*150)");
  /*  if (c==1) sigTree1->Draw((var+">>h1").c_str(), "evweight*(PhotonsMass>0.7*150 && PhotonsMass<1.3*150&& (abs(ph1_eta)<1.5 && abs(ph2_eta)<1.5) && (ph1_r9<0.94 || ph2_r9<0.94 )) ");
  if (c==2) sigTree1->Draw((var+">>h1").c_str(), "evweight*(PhotonsMass>0.7*150 && PhotonsMass<1.3*150&& (abs(ph1_eta)>1.5 || abs(ph2_eta)>1.5) && (ph1_r9>0.94 && ph2_r9>0.94 ))");
  if (c==3) sigTree1->Draw((var+">>h1").c_str(), "evweight*(PhotonsMass>0.7*150 && PhotonsMass<1.3*150&& (abs(ph1_eta)>1.5 || abs(ph2_eta)>1.5) && (ph1_r9<0.94 || ph2_r9<0.94 )) ");*/
  
  if (c==0) sigTree2->Draw((var+">>h2").c_str(), "evweight*(PhotonsMass>0.7*150 && PhotonsMass<1.3*150 )");
  /*  if (c==1) sigTree2->Draw((var+">>h2").c_str(), "evweight*(PhotonsMass>0.7*150 && PhotonsMass<1.3*150&& (abs(ph1_eta)<1.5 && abs(ph2_eta)<1.5) && (ph1_r9<0.94 || ph2_r9<0.94 )) ");
  if (c==2) sigTree2->Draw((var+">>h2").c_str(), "evweight*(PhotonsMass>0.7*150 && PhotonsMass<1.3*150&& (abs(ph1_eta)>1.5 || abs(ph2_eta)>1.5) && (ph1_r9>0.94 && ph2_r9>0.94 ))");
  if (c==3) sigTree2->Draw((var+">>h2").c_str(), "evweight*(PhotonsMass>0.7*150 && PhotonsMass<1.3*150&& (abs(ph1_eta)>1.5 || abs(ph2_eta)>1.5) && (ph1_r9<0.94 || ph2_r9<0.94 )) ");*/
  

  h1->SetLineColor(kBlack);
  h2->SetLineColor(kRed);
 
  

 
 
  
  h1->GetXaxis()->SetTitleFont(42);
  h1->GetXaxis()->SetTitleSize(0.05);
  h1->SetMarkerSize(0);
  h2->SetMarkerSize(0);
  h1->Sumw2();
  h2->Sumw2();
  h1->Scale(1./h1->Integral());
  h1->GetYaxis()->SetRangeUser(0.0001, 0.45);
  h1->Draw("hist");
  h2->Scale(1./h2->Integral());
  h2->Draw("histSAME");
 
  std::cout<<"Private: Mean = "<<h1->GetMean()<<" +/- "<<h1->GetMeanError()<<"  RMS: "<<h1->GetRMS()<<" +/- "<<h1->GetRMSError()<<std::endl;
  std::cout<<"Official: Mean = "<<h2->GetMean()<<" +/- "<<h2->GetMeanError()<<"  RMS: "<<h2->GetRMS()<<" +/- "<<h2->GetRMSError()<<std::endl;

  TLegend* legmc = new TLegend(0.5, 0.6, 0.87, 0.89, "", "brNDC");
  legmc->AddEntry(h1,"m_{#gamma#gamma} = 150 GeV Private","L");
  legmc->AddEntry(h2,"m_{#gamma#gamma} = 150 GeV Official","L");
  legmc->SetTextSize(0.0206044);  
  legmc->SetTextFont(42);
  legmc->SetBorderSize(0);
  legmc->SetFillStyle(0);
  legmc->Draw(); 

  TLatex *lat  = new TLatex(0.65,0.9,TString::Format("Cat: %d", c));
  
  lat->SetTextSize(0.038);
  lat->SetTextAlign(11);
  lat->SetTextFont(42);
  lat->Draw("same");
  lat->SetNDC();


 

  TPaveText* label_cms = get_labelCMS(0, "2012", true);
  TPaveText* label_sqrt = get_labelSqrt(0);
  label_cms->Draw("same");
  label_sqrt->Draw("same");

  h1->GetYaxis()->SetRangeUser(0.0001, 0.25);
  c1->SaveAs(TString::Format("plots/massesShape_Official_cat%d.png",c));
  c1->SaveAs(TString::Format("plots/massesShape_Official_cat%d.pdf",c));
  c1->SaveAs(TString::Format("~/www/plotsNota/massesShape_Official_cat%d.png",c));
  c1->SaveAs(TString::Format("~/www/plotsNota/massesShape_Official_cat%d.pdf",c));
  
  c1->SetLogy();
  
  c1->SaveAs(TString::Format("plots/massesShape_Official_cat%d_LOG.png",c));
  c1->SaveAs(TString::Format("plots/massesShape_Official_cat%d_LOG.pdf",c));

  c1->SaveAs(TString::Format("~/www/plotsNota/massesShape_Official_cat%d_LOG.png",c));
  c1->SaveAs(TString::Format("~/www/plotsNota/massesShape_Official_cat%d_LOG.pdf",c));

  std::cout<<" CAT: "<<c<< " Mean: "<< h1->GetMean()<< " Mean Official: "<<h2->GetMean()<<std::endl;
  std::cout<<"Diff %"<<(h1->GetMean()-h2->GetMean())/h2->GetMean()std::endl;
}










void PlotVtxRecoEff() {

  Int_t nmass = 7;
  Int_t masses[7] = {150, 200, 250, 300, 400, 600, 800};
  Double_t massesD[7] = {150., 200., 250., 300., 400., 600, 800.};

  TString inDir = "";

  
  TH1F* massCat0[7];
  TH1F* massCat1[7];
  TH1F* massCat2[7];
  TH1F* massCat3[7];
  TH1F* massCat4[7];

  TH1F* massCat0_vtxOK[7];
  TH1F* massCat1_vtxOK[7];
  TH1F* massCat2_vtxOK[7];
  TH1F* massCat3_vtxOK[7];
  TH1F* massCat4_vtxOK[7];


 
  Double_t effCAT0[7] ;
  Double_t effCAT1[7] ;
  Double_t effCAT2[7] ;
  Double_t effCAT3[7] ;
  Double_t effCAT4[7] ;

  
  for(int iMass=0;iMass<7;iMass++){ //per ogni massa
  //chain summing up all production modes
  TChain* sigTree1  = new TChain();

// common preselection cut
  TString maincut = TString::Format("PhotonsMass>=(%d*0.1) && PhotonsMass<=(%d*1.9)", masses[iMass], masses[iMass]);   // livia
  
  if(iMass<5){
  sigTree1->Add(TString::Format("histograms_CMS-HGG_11032014_MC.root/ggh_m%d_8TeV", masses[iMass]));
  sigTree1->Add(TString::Format("histograms_CMS-HGG_11032014_MC.root/vbf_m%d_8TeV", masses[iMass]));
  sigTree1->Add(TString::Format("histograms_CMS-HGG_11032014_MC.root/wzh_m%d_8TeV", masses[iMass]));
  sigTree1->Add(TString::Format("histograms_CMS-HGG_11032014_MC.root/tth_m%d_8TeV", masses[iMass]));


  }else {

  sigTree1->Add(TString::Format("histograms_CMS-HGG_HighMass_600800.root/ggh_m%d_8TeV", masses[iMass]));
  sigTree1->Add(TString::Format("histograms_CMS-HGG_HighMass_600800.root/vbf_m%d_8TeV", masses[iMass]));
  sigTree1->Add(TString::Format("histograms_CMS-HGG_HighMass_600800.root/wzh_m%d_8TeV", masses[iMass]));
  sigTree1->Add(TString::Format("histograms_CMS-HGG_HighMass_600800.root/tth_m%d_8TeV", masses[iMass]));

  }
  sigTree1->SetTitle("sigTree1");
  sigTree1->SetName("sigTree1");

 
  massCat0[iMass]= new TH1F(TString::Format("massCat0_m%d", masses[iMass]),TString::Format("massCat0_m%d", masses[iMass]), 120, masses[iMass]*0.3,masses[iMass]*1.7 );

  massCat1[iMass]= new TH1F(TString::Format("massCat1_m%d", masses[iMass]),TString::Format("massCat1_m%d", masses[iMass]), 120, masses[iMass]*0.3,masses[iMass]*1.7 );

  massCat2[iMass]= new TH1F(TString::Format("massCat2_m%d", masses[iMass]),TString::Format("massCat2_m%d", masses[iMass]), 120, masses[iMass]*0.3,masses[iMass]*1.7 );

  massCat3[iMass]= new TH1F(TString::Format("massCat3_m%d", masses[iMass]),TString::Format("massCat3_m%d", masses[iMass]), 120, masses[iMass]*0.3,masses[iMass]*1.7 );

  massCat4[iMass]= new TH1F(TString::Format("massCat4_m%d", masses[iMass]),TString::Format("massCat4_m%d", masses[iMass]), 120, masses[iMass]*0.3,masses[iMass]*1.7 );



  massCat0_vtxOK[iMass]= new TH1F(TString::Format("massCat0_m%d_vtxOK", masses[iMass]),TString::Format("massCat0_m%d_vtxOK", masses[iMass]), 120, masses[iMass]*0.3,masses[iMass]*1.7 );
  
  massCat1_vtxOK[iMass]= new TH1F(TString::Format("massCat1_m%d_vtxOK", masses[iMass]),TString::Format("massCat1_m%d_vtxOK", masses[iMass]), 120, masses[iMass]*0.3,masses[iMass]*1.7 );

  massCat2_vtxOK[iMass]= new TH1F(TString::Format("massCat2_m%d_vtxOK", masses[iMass]),TString::Format("massCat2_m%d_vtxOK", masses[iMass]), 120, masses[iMass]*0.3,masses[iMass]*1.7 );
  
  massCat3_vtxOK[iMass]= new TH1F(TString::Format("massCat3_m%d_vtxOK", masses[iMass]),TString::Format("massCat3_m%d_vtxOK", masses[iMass]), 120, masses[iMass]*0.3,masses[iMass]*1.7 );
 
  massCat4_vtxOK[iMass]= new TH1F(TString::Format("massCat4_m%d_vtxOK", masses[iMass]),TString::Format("massCat4_m%d_vtxOK", masses[iMass]), 120, masses[iMass]*0.3,masses[iMass]*1.7 );
 
 

  TString cutCat0 = "((abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 ))";
  TString cutCat1 = "((abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9<0.94 || ph2_r9<0.94 ))";
  TString cutCat2 = "((abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 && ph2_r9>0.94 ))";
  TString cutCat3 = "((abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 || ph2_r9>0.94 ))";
  TString cutCat4 = "";



  sigTree1->Draw(TString::Format("PhotonsMass>>massCat0_m%d", masses[iMass]),"evweight*("+ maincut+" && "+cutCat0+")");
  sigTree1->Draw(TString::Format("PhotonsMass>>massCat1_m%d", masses[iMass]),"evweight*("+ maincut+" && "+cutCat1+")");
  sigTree1->Draw(TString::Format("PhotonsMass>>massCat2_m%d", masses[iMass]),"evweight*("+ maincut+" && "+cutCat2+")");
  sigTree1->Draw(TString::Format("PhotonsMass>>massCat3_m%d", masses[iMass]),"evweight*("+ maincut+" && "+cutCat3+")");
  sigTree1->Draw(TString::Format("PhotonsMass>>massCat4_m%d", masses[iMass]),"evweight*("+ maincut+")");


  std::cout<<"----------m: "<<masses[iMass]<<" --------"<<std::endl;
  std::cout<<massCat0[iMass]->Integral()<<std::endl;
  std::cout<<massCat1[iMass]->Integral()<<std::endl;
  std::cout<<massCat2[iMass]->Integral()<<std::endl;
  std::cout<<massCat3[iMass]->Integral()<<std::endl;
  std::cout<<massCat4[iMass]->Integral()<<std::endl;

  sigTree1->Draw(TString::Format("PhotonsMass>>massCat0_m%d_vtxOK", masses[iMass]),"evweight*("+ maincut+" && abs(gv_z-vtx_z)<1 &&"+cutCat0+")");
  sigTree1->Draw(TString::Format("PhotonsMass>>massCat1_m%d_vtxOK", masses[iMass]),"evweight*("+ maincut+" && abs(gv_z-vtx_z)<1 &&"+cutCat1+")");
  sigTree1->Draw(TString::Format("PhotonsMass>>massCat2_m%d_vtxOK", masses[iMass]),"evweight*("+ maincut+" && abs(gv_z-vtx_z)<1 &&"+cutCat2+")");
  sigTree1->Draw(TString::Format("PhotonsMass>>massCat3_m%d_vtxOK", masses[iMass]),"evweight*("+ maincut+" && abs(gv_z-vtx_z)<1 &&"+cutCat3+")");
  sigTree1->Draw(TString::Format("PhotonsMass>>massCat4_m%d_vtxOK", masses[iMass]),"evweight*("+ maincut+" && abs(gv_z-vtx_z)<1)");

  std::cout<<massCat0_vtxOK[iMass]->Integral()<<std::endl;
  std::cout<<massCat1_vtxOK[iMass]->Integral()<<std::endl;
  std::cout<<massCat2_vtxOK[iMass]->Integral()<<std::endl;
  std::cout<<massCat3_vtxOK[iMass]->Integral()<<std::endl;
  std::cout<<massCat4_vtxOK[iMass]->Integral()<<std::endl;

 

  effCAT0[iMass] = massCat0_vtxOK[iMass]->Integral()/massCat0[iMass]->Integral();
  effCAT1[iMass] = massCat1_vtxOK[iMass]->Integral()/massCat1[iMass]->Integral();
  effCAT2[iMass] = massCat2_vtxOK[iMass]->Integral()/massCat2[iMass]->Integral();
  effCAT3[iMass] = massCat3_vtxOK[iMass]->Integral()/massCat3[iMass]->Integral();
  effCAT4[iMass] = massCat4_vtxOK[iMass]->Integral()/massCat4[iMass]->Integral();

  std::cout<<"------------- eff ------------ "<<std::endl;

  std::cout<<effCAT0[iMass]<<std::endl;
  std::cout<<effCAT1[iMass]<<std::endl;
  std::cout<<effCAT2[iMass]<<std::endl;
  std::cout<<effCAT3[iMass]<<std::endl;
  std::cout<<effCAT4[iMass]<<std::endl;

  }



  Double_t massesErr[7] = {0.,0.,0.,0.,0.};
  Double_t effErr[7] = {0.,0.,0.,0.,0.};
  

  TGraphErrors* effCAT0graph = new TGraphErrors(7, massesD, effCAT0, massesErr, effErr);
  TGraphErrors* effCAT1graph = new TGraphErrors(7, massesD, effCAT1, massesErr, effErr);
  TGraphErrors* effCAT2graph = new TGraphErrors(7, massesD, effCAT2, massesErr, effErr);
  TGraphErrors* effCAT3graph = new TGraphErrors(7, massesD, effCAT3, massesErr, effErr);
  TGraphErrors* effCAT4graph = new TGraphErrors(7, massesD, effCAT4, massesErr, effErr);

  effCAT0graph->SetMarkerColor(kPink-8);
  effCAT1graph->SetMarkerColor(kOrange+1);
  effCAT2graph->SetMarkerColor(kSpring-8);
  effCAT3graph->SetMarkerColor(kBlue-7);
  effCAT4graph->SetMarkerColor(kMagenta-9);
  effCAT0graph->SetLineColor(kPink-8);
  effCAT1graph->SetLineColor(kOrange+1);
  effCAT2graph->SetLineColor(kSpring-8);
  effCAT3graph->SetLineColor(kBlue-7);
  effCAT4graph->SetLineColor(kMagenta-9);

 TCanvas* c1 = new TCanvas("c1", "c1", 1);
  c1->cd();
  
  

  TMultiGraph* mg= new TMultiGraph();
  mg->Add(effCAT0graph);
  mg->Add(effCAT1graph);
  mg->Add(effCAT2graph);
  mg->Add(effCAT3graph);
  mg->Add(effCAT4graph);

  mg->Draw("APE");
  mg->GetYaxis()->SetRangeUser(0., 1.);
  mg->GetYaxis()->SetTitle("Vtx MVA Reconstruction efficiency");
  mg->GetXaxis()->SetTitle("m_{X} [GeV]");
  TLegend* legmc = new TLegend(0.6, 0.2, 0.85, 0.7, "", "brNDC");
  legmc->AddEntry(effCAT0graph, "CAT 0", "P");
  legmc->AddEntry(effCAT1graph, "CAT 1", "P");
  legmc->AddEntry(effCAT2graph, "CAT 2", "P");
  legmc->AddEntry(effCAT3graph, "CAT 3", "P");
  legmc->AddEntry(effCAT4graph, " All Categories Combined", "P");
  legmc->SetTextSize(0.0206044);  
  legmc->SetTextFont(42);
  legmc->SetBorderSize(0);
  legmc->SetFillStyle(0);
  legmc->Draw("same"); 

  TLatex *lat  = new TLatex(0.65,0.9,"Vertex Recnstruction Efficiency");
  
  lat->SetTextSize(0.038);
  lat->SetTextAlign(11);
  lat->SetTextFont(42);
  // lat->Draw("same");
  lat->SetNDC();


 

  TPaveText* label_cms = get_labelCMS(0, "2012", true);
  TPaveText* label_sqrt = get_labelSqrt(0);
  label_cms->Draw("same");
  label_sqrt->Draw("same");

  c1->SetLogy(0);
  c1->SaveAs("plots/VtxRecoEff.png");
  c1->SaveAs("plots/VtxRecoEff.pdf");
  c1->SaveAs("~/www/plotsNota/VtxRecoEff.png");
  c1->SaveAs("~/www/plotsNota/VtxRecoEff.pdf");
 }







void MakePlotMassDataMC(RooWorkspace* w, Float_t mass ) {

  Int_t ncat = NCAT;
  TString inDir = "";

  Float_t MASS(mass);

  // Luminosity:
  Float_t Lum = 19500.0;  
  RooRealVar lumi("lumi","lumi",Lum);
  w->import(lumi); 
  
  // Variables
  RooArgSet* ntplVars = defineVariables();
  RooRealVar* PhotonsMass = w->var("PhotonsMass");
   
  int iMass = abs(mass);      
 // common preselection cut
  TString mainCut = "(PhotonsMass>=130 && PhotonsMass<=1200)";   //130-2000
  RooPlot*  plotPhotonsMassDataMC[NCAT];
  
  //**********DATA***************//
  // create tree
  TFile file("histograms_CMS-HGG_24072013.root");
  TTree* dataTree = (TTree*) file.Get("Data");
   
  //**********G+jets***************//
  // create tree

  TChain* gjTree = new TChain();
  gjTree->Add("histograms_CMS-HGG_24072013.root/gjet_20_8TeV_pf");
  gjTree->Add("histograms_CMS-HGG_24072013.root/gjet_40_8TeV_pf");

  TChain* qcdTree = new TChain();
  qcdTree->Add("histograms_CMS-HGG_24072013.root/qcd_30_8TeV_pf");
  qcdTree->Add("histograms_CMS-HGG_24072013.root/qcd_40_8TeV_pf"); 
 
  TChain* jjTree = new TChain();
  jjTree->Add("histograms_CMS-HGG_24072013.root/qcd_30_8TeV_ff");
  jjTree->Add("histograms_CMS-HGG_24072013.root/qcd_40_8TeV_ff"); 
 
  std::cout<< qcdTree->GetEntries()<<std::endl;

  //**********DIPHOTJET***************//
  // create tree

  TChain* diphotjTree = new TChain();
  diphotjTree->Add("histograms_CMS-HGG_24072013.root/diphojet_8TeV");
   

  //**********DIPHOT***************//
  // create tree

  TChain* diphotTree = new TChain();
  diphotTree->Add("histograms_CMS-HGG_24072013.root/dipho_Box_25_8TeV");
  diphotTree->Add("histograms_CMS-HGG_24072013.root/dipho_Box_250_8TeV");
   

  

  TH1F* h_data[NCAT+1];
  TH1F* h_data_b[NCAT+1];
  TH1F*  h_gj[NCAT+1];
  TH1F*  h_qcd[NCAT+1];
  TH1F*  h_jj[NCAT+1];
  TH1F*  h_diphot[NCAT+1];
  TH1F*  h_diphotj[NCAT+1];
  TH1F* h_sum[NCAT+1]; 
 
  
  TCanvas* ctmp = new TCanvas("ctmp","PhotonsMass Background Categories",1);
  Int_t nbin = 87;
  Double_t  min = 130;
  Double_t  max = 1000;

   for (int c=NCAT; c<NCAT+1; ++c) {


     h_data[c]= new TH1F(TString::Format("h_data_cat%d",c), TString::Format("h_data_cat%d",c), nbin, min, max);
     h_data_b[c]= new TH1F(TString::Format("h_data_b_cat%d",c), TString::Format("h_data_b_cat%d",c), nbin, min, max);
     h_gj[c]= new TH1F(TString::Format("h_gj_cat%d",c), TString::Format("h_gj_cat%d",c), nbin, min, max);
     h_qcd[c]= new TH1F(TString::Format("h_qcd_cat%d",c), TString::Format("h_qcd_cat%d",c), nbin, min, max);
     h_jj[c]= new TH1F(TString::Format("h_jj_cat%d",c), TString::Format("h_jj_cat%d",c), nbin, min, max);
     h_diphot[c]= new TH1F(TString::Format("h_diphot_cat%d",c), TString::Format("h_diphot_cat%d",c), nbin, min, max);
     h_diphotj[c]= new TH1F(TString::Format("h_diphotj_cat%d",c), TString::Format("h_diphotj_cat%d",c), nbin, min, max);



    // 1)  prime 4 cat livia
     if (c==0){//&&(PhotonsMass<178 || PhotonsMass >402) 
   
      dataTree->Draw("PhotonsMass>>h_data_cat0", "("+mainCut+"&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 ))*1");
      dataTree->Draw("PhotonsMass>>h_data_b_cat0", "("+mainCut+"&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 ))*1");
      gjTree->Draw("PhotonsMass>>h_gj_cat0", "("+mainCut+"&&(abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 ))*evweight");
      gjTree->Draw("PhotonsMass>>h_qcd_cat0", "("+mainCut+"&&(abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 ))*evweight");
      diphotTree->Draw("PhotonsMass>>h_diphot_cat0","("+ mainCut+"&&(abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 ))*evweight");
      diphotjTree->Draw("PhotonsMass>>h_diphotj_cat0","("+ mainCut+"&&(abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 ))*evweight");
    
     }
     if (c==1){//&&(PhotonsMass<178 || PhotonsMass >402)&&

       dataTree->Draw("PhotonsMass>>h_data_cat1", "("+mainCut+"&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9<0.94 || ph2_r9<0.94 )  )*1");
       dataTree->Draw("PhotonsMass>>h_data_b_cat1", "("+mainCut+"&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9<0.94 || ph2_r9<0.94 )  )*1");
     gjTree->Draw("PhotonsMass>>h_gj_cat1", "("+mainCut+"&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9<0.94 || ph2_r9<0.94 ))*evweight ");
     qcdTree->Draw("PhotonsMass>>h_qcd_cat1", "("+mainCut+"&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9<0.94 || ph2_r9<0.94 ))*evweight ");
      diphotTree->Draw("PhotonsMass>>h_diphot_cat1", "("+mainCut+"&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9<0.94 || ph2_r9<0.94 ) )*evweight");
     diphotjTree->Draw("PhotonsMass>>h_diphotj_cat1", "("+mainCut+"&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9<0.94 || ph2_r9<0.94 ) )*evweight");

       }
     if (c==2){//&&(PhotonsMass<178 || PhotonsMass >402)

     dataTree->Draw("PhotonsMass>>h_data_cat2", "("+mainCut+"&& (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 && ph2_r9>0.94 ))*1");
     dataTree->Draw("PhotonsMass>>h_data_b_cat2", "("+mainCut+"&& (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 && ph2_r9>0.94 ))*1");
     gjTree->Draw("PhotonsMass>>h_gj_cat2", "("+mainCut+"&& (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 && ph2_r9>0.94 ))*evweight");
     qcdTree->Draw("PhotonsMass>>h_qcd_cat2", "("+mainCut+"&& (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 && ph2_r9>0.94 ))*evweight");
     diphotTree->Draw("PhotonsMass>>h_diphot_cat2", "("+mainCut+"&& (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 && ph2_r9>0.94 ))*evweight");
     diphotjTree->Draw("PhotonsMass>>h_diphotj_cat2", "("+mainCut+"&& (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 && ph2_r9>0.94 ))*evweight");
    
   }
     if (c==3){//&&(PhotonsMass<178 || PhotonsMass >402)

     dataTree->Draw("PhotonsMass>>h_data_cat3", "("+mainCut+"&& (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9<0.94 || ph2_r9<0.94 ))*1");
     dataTree->Draw("PhotonsMass>>h_data_b_cat3", "("+mainCut+"&& (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9<0.94 || ph2_r9<0.94 ))*1");
     gjTree->Draw("PhotonsMass>>h_gj_cat3", "("+mainCut+"&& (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9<0.94 || ph2_r9<0.94 ))*evweight ");
     qcdTree->Draw("PhotonsMass>>h_qcd_cat3", "("+mainCut+"&& (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9<0.94 || ph2_r9<0.94 ))*evweight ");
      diphotTree->Draw("PhotonsMass>>h_diphot_cat3","("+ mainCut+"&& (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9<0.94 || ph2_r9<0.94 ) )*evweight");
      diphotjTree->Draw("PhotonsMass>>h_diphotj_cat3","("+ mainCut+"&& (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9<0.94 || ph2_r9<0.94 ))*evweight ");
    
   }
     if (c==4){//&&(PhotonsMass<178 || PhotonsMass >402)

    dataTree->Draw("PhotonsMass>>h_data_cat4","("+ mainCut+")*1");
    dataTree->Draw("PhotonsMass>>h_data_b_cat4","("+ mainCut+")*1");
    gjTree->Draw("PhotonsMass>>h_gj_cat4", "("+mainCut+")*evweight");
    qcdTree->Draw("PhotonsMass>>h_qcd_cat4", "("+mainCut+")*evweight");
    jjTree->Draw("PhotonsMass>>h_jj_cat4", "("+mainCut+")*evweight");
    diphotTree->Draw("PhotonsMass>>h_diphot_cat4","("+ mainCut+")*evweight");
    diphotjTree->Draw("PhotonsMass>>h_diphotj_cat4","("+ mainCut+")*evweight");
    

   }

 
     //fisso la normalizzazione di gjets a quella gj+qcd xkè la shape di qcd _pf fa schifo, quindi prenidamoq uella di gj _pf
     Double_t qcdInt = h_qcd[c]->Integral();
     Double_t gJetIntFixed = qcdInt+h_gj[c]->Integral();
     h_gj[c]->Scale(gJetIntFixed/h_gj[c]->Integral());     
     h_gj[c]->SetFillColor(kAzure+8);
     h_gj[c]->Sumw2();
     
     h_diphotj[c]->SetFillColor(kTeal+9);
     h_diphotj[c]->Sumw2();
     
     //sommo i due diphot
     h_diphot[c]->Add(h_diphotj[c]);
     h_diphot[c]->SetFillColor(kSpring+7);
     h_diphot[c]->Sumw2();
     
     //provo jj per i referee PLB
     h_jj[c]->SetFillColor(kOrange+8);
     std::cout<<"/////////////////////////////"<<h_jj[c]->GetEntries()<<std::endl;
     h_jj[c]->Sumw2();
     
     h_sum[c] = (TH1F*) h_gj[c]->Clone();
     h_sum[c]->Add(h_diphot[c]);
     h_sum[c]->Add(h_jj[c]);
     h_sum[c]->Sumw2();
     h_sum[c]->SetFillColor(kBlack);
     h_sum[c]->SetFillStyle(3003);
     h_sum[c]->SetMarkerSize(0);

     //make kolmogorov test between data and MC
     Double_t CHI2ndf = h_data[c]->Chi2Test(h_sum[c], "UWPCHI2/NDF");
     TPaveText* label = new TPaveText(0.6, 0.67, 0.85, 0.7, "brNDC" );
     label->SetFillColor(kWhite);
     label->SetBorderSize(0.);
     label->SetTextSize(0.038);
     label->SetTextAlign(11);
     label->SetTextFont(42);
     // label->AddText(TString::Format("#chi^{2}/NDF: %.3f", CHI2ndf));
     
     THStack hs("hs","hs");
     hs.Add(h_jj[c]);
     hs.Add(h_gj[c]); 
     hs.Add(h_diphot[c]); 
  
     std::cout<<h_sum[c]->Integral()<<std::endl;
     std::cout<<h_diphot[c]->Integral()<<" %: "<<h_diphot[c]->Integral()/h_sum[c]->Integral()<<std::endl;
     std::cout<<h_gj[c]->Integral()<<" %: "<<h_gj[c]->Integral()/h_sum[c]->Integral()<<std::endl;
     std::cout<<h_jj[c]->Integral()<<" %: "<<h_jj[c]->Integral()/h_sum[c]->Integral()<<std::endl;
     std::cout<<h_data[c]->Integral()<<std::endl;
   
 
    ctmp->cd();
    h_data[c]->Sumw2();
    h_data_b[c]->Sumw2();
    h_sum[c]->Sumw2();

    TH1F* h1_ratio1 = (TH1F*)h_data_b[c]->Clone();
    TH1F* h1_ratio1_unblind = (TH1F*)h_data[c]->Clone();
    TH1F* h1_ratio2 = (TH1F*)h_sum[c]->Clone();

    //for(int i = 0;i < h1_ratio1->GetNbinsX();i++) std::cout<<" ratio1: "<<h1_ratio1->GetBinContent(i)<<std::endl;//if(h1_ratio1->GetBinContent(i)==0)h1_ratio1->SetBinContent(i, -1);
   
   

    //  std::cout<<"int 700-850: "<<h_data[c]->Integral(h_data[c]->FindBin(700.),h_data[c]->FindBin(950.) )<<"  int 850-1200: "<<h_data[c]->Integral(h_data[c]->FindBin(950.),h_data[c]->FindBin(1200.) )<<std::endl;



    ctmp->Clear();
    //-------pad 1-------//
    TPad * pad1 = new TPad("pad1", "pad1",0.01,0.13,0.75,1.);  
    
   
    pad1->SetRightMargin(0.1);
    
    pad1->SetLogy();
    pad1->Draw();
    pad1->cd();

    h_data_b[c]->SetMarkerSize(0.6);
    h_data_b[c]->Draw("pe");
    
    h_data_b[c]->GetYaxis()->SetTitle(TString::Format("Events /10 GeV", (max-min)/nbin));
    //h__data[c]->GetXaxis()->SetTitle("m_{#gamma #gamma} [GeV]");
    h_data_b[c]->GetYaxis()->SetRangeUser(0.1, h_data_b[c]->GetMaximum()*2);
    hs.Draw("histsame");
    h_data_b[c]->Draw("pesame");
    //  h_data[c]->Draw("pesame");
    h_sum[c]->Draw("E2same");
    // std::cout<<"--------_> "<<h_data_b[c]->Integral()<<std::endl; 
  
  

    TLegend *leg1;
    if(c!=4)leg1 = new TLegend(0.6075,0.6536441,0.8575,0.9340678, TString::Format("Category %d",c), "brNDC");
    else leg1 = new TLegend(0.5075,0.6536441,0.8075,0.9340678, TString::Format("All classes combined",c), "brNDC");
    leg1->AddEntry(h_data_b[c],"Data","PE");
    leg1->AddEntry(h_diphot[c],"#gamma + #gamma", "F");
    leg1->AddEntry(h_gj[c],"#gamma + jet","F");
    leg1->AddEntry(h_jj[c],"jet + jet", "F");   
    leg1->AddEntry(h_sum[c], "Bkg uncertainty", "F");
    
   
    
    leg1->SetTextSize(0.035);
    leg1->SetTextFont(42);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
    leg1->Draw("same");

    // label->Draw("same");
 

    TPaveText* label_cms = get_labelCMS(0, "2012", false);
    TPaveText* label_sqrt = get_labelSqrt(0);
    //label_cms->Draw("same");
    //label_sqrt->Draw("same");
    int iPos=11 ;
    CMS_lumi( pad1,false,iPos );
    pad1->SetLogy(0);
    pad1->RedrawAxis();

    ctmp->cd();

    //-------pad 2------//
    TPad * pad2 = new TPad("pad2", "pad2",0.01,0.001,0.75,0.2);
    pad2->SetGrid();
    
    //pad2->SetTopMargin(0.01);
    pad2->SetBottomMargin(0.4);
    pad2->SetRightMargin(0.1);
    pad2->Draw();
    pad2->cd();

    Double_t xmax = h1_ratio1->GetXaxis()->GetXmax();
    Double_t xmin = h1_ratio1->GetXaxis()->GetXmin();
    TLine* line = new TLine(xmin,1.,xmax,1.);
  

    h1_ratio1->SetStats(0);
    
    h1_ratio1->Divide(h_sum[c]);
    h1_ratio1_unblind->Divide(h_sum[c]);
    h1_ratio2->Divide(h_sum[c]);
    h1_ratio1->SetMarkerStyle(20);
    h1_ratio1->SetMarkerSize(0.6);
    //  h1_ratio1->GetXaxis()->SetTitle(xAxis.c_str());
    h1_ratio1->GetYaxis()->SetRangeUser(0., 2.); // cwr zoom
    h1_ratio1->GetYaxis()->SetNdivisions(2,false);
    h1_ratio1->GetYaxis()->SetTitle("Data/Bkg");
    h1_ratio1->GetYaxis()->SetTitleFont(42);
    h1_ratio1->GetXaxis()->SetTitleFont(42);
    h1_ratio1->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    
    h1_ratio1->GetXaxis()->SetTitleSize(0.2);
    h1_ratio1->GetXaxis()->SetLabelSize(0.16);
    h1_ratio1->GetXaxis()->SetLabelOffset(5.16);
    h1_ratio1->GetYaxis()->SetLabelSize(0.16);
    h1_ratio1->GetYaxis()->SetTitleSize(0.15);
    h1_ratio1->GetYaxis()->SetTitleOffset(0.45);
    h1_ratio1->GetXaxis()->SetTitleOffset(0.8);

    
    
    for(int j = 0;j<=h1_ratio1->GetNbinsX();j++){
      
      if(h_sum[c]->GetBinContent(j))  h1_ratio1->SetBinError(j,sqrt(pow(h_data[c]->GetBinError(j)/h_sum[c]->GetBinContent(j), 2)+ pow(h_data[c]->GetBinContent(j)*h_sum[c]->GetBinError(j)/(h_sum[c]->GetBinContent(j)*h_sum[c]->GetBinContent(j)),2)));
      else h1_ratio1->SetBinError(j,0.);
    }
    h1_ratio1->Draw("PEX0");
   
    for(int j = 0;j<=h1_ratio1_unblind->GetNbinsX();j++){
     
      if(h_sum[c]->GetBinContent(j))  h1_ratio1_unblind->SetBinError(j,sqrt(pow(h_data[c]->GetBinError(j)/h_sum[c]->GetBinContent(j), 2)+ pow(h_data[c]->GetBinContent(j)*h_sum[c]->GetBinError(j)/(h_sum[c]->GetBinContent(j)*h_sum[c]->GetBinContent(j)),2)));
      else h1_ratio1_unblind->SetBinError(j,0.);
    }
  
    for(int j = 0;j<=h1_ratio2->GetNbinsX();j++){
      if(h_sum[c]->GetBinContent(j)) h1_ratio2->SetBinError(j,h_sum[c]->GetBinError(j)/h_sum[c]->GetBinContent(j));
      else h1_ratio2->SetBinError(j,0.);
    }
    h1_ratio2->Draw("E2same");
    
    line->SetLineWidth(1.);
    line->Draw("same");
    
  

    ctmp->cd();
    //-------pad 3------//
    
   
    TPad * pad3 = new TPad("pad3", "pad3",0.68,0.001,1.,0.2);
    pad3->SetGrid();
    
    //pad2->SetTopMargin(0.01);
    pad3->SetBottomMargin(0.4);
    pad3->SetRightMargin(0.1);
    //pad3->Draw();
    pad3->cd();

    TH1F* h1_res = new TH1F("h1_res", "h1_res",11,-5, 5.);
   
    for(int i = 1;i < h1_ratio1_unblind->GetNbinsX();i++){
      if(h1_ratio1_unblind->GetBinError(i)==0)continue;
      //   std::cout<<((h1_ratio1_unblind->GetBinContent(i))-1)/h1_ratio1_unblind->GetBinError(i)<<std::endl;
      h1_res->SetFillColor(kBlack);
      h1_res->Fill(((h1_ratio1_unblind->GetBinContent(i))-1)/h1_ratio1_unblind->GetBinError(i));

    }
    // std::cout<<h1_res->GetEntries()<<std::endl;
   
    h1_res->GetYaxis()->SetNdivisions(4,false);
    h1_res->GetXaxis()->SetTitle("pull");
    h1_res->GetXaxis()->SetTitleSize(0.15);
    h1_res->GetXaxis()->SetTitleOffset(0.3);
    h1_res->GetXaxis()->SetLabelSize(0.1);
    h1_res->GetYaxis()->SetLabelSize(0.1);
    h1_res->Draw("hbar");
   

    Double_t mean = h1_res->GetMean();
    Double_t sigma = h1_res->GetRMS();
    Double_t meanErr = h1_res->GetMeanError();
    Double_t sigmaErr = h1_res->GetRMSError();

    //   ctmp->cd();
    TPaveText* label2 = new TPaveText(0.1528545,0.02150538,0.8593564,0.296661, "brNDC" );
    label2->SetFillColor(kWhite);
    label2->SetBorderSize(0.);
    label2->SetTextSize(0.1213922);
    label2->SetTextAlign(11);
    label2->SetTextFont(42);
    label2->AddText(TString::Format("#splitline{Mean: %.3f #pm %.3f }{RMS: %.3f #pm %.3f}", mean,meanErr,sigma, sigmaErr));


    label2->Draw("same");
  
    ctmp->SaveAs("plots/DATA_MC_MASS_"+TString::Format("cat%d.png", c));
    ctmp->SaveAs("plots/DATA_MC_MASS_"+TString::Format("cat%d.C", c));
    ctmp->SaveAs("plots/DATA_MC_MASS_"+TString::Format("cat%d.pdf", c));
    ctmp->SaveAs("~/www/plotsNota/DATA_MC_MASS_"+TString::Format("cat%d.png", c));
    ctmp->SaveAs("~/www/plotsNota/DATA_MC_MASS_"+TString::Format("cat%d.pdf", c));
  
    pad1->SetLogy(1);
    ctmp->SaveAs("plots/DATA_MC_MASS_"+TString::Format("cat%d_LOG.png", c));
    ctmp->SaveAs("plots/DATA_MC_MASS_"+TString::Format("cat%d_LOG.pdf", c));
    ctmp->SaveAs("~/www/plotsNota/DATA_MC_MASS_"+TString::Format("cat%d_LOG.png", c));
    ctmp->SaveAs("~/www/plotsNota/DATA_MC_MASS_"+TString::Format("cat%d_LOG.pdf", c));
    ctmp->SaveAs("~/www/plotsPAS/DATA_MC_MASS_"+TString::Format("cat%d_LOG.png", c));
    ctmp->SaveAs("~/www/plotsPAS/DATA_MC_MASS_"+TString::Format("cat%d_LOG.pdf", c));
    ctmp->SaveAs("~/www/plotsPAS/DATA_MC_MASS_"+TString::Format("cat%d_LOG.C", c));
    ctmp->SaveAs("~/www/plotsPAS/DATA_MC_MASS_"+TString::Format("cat%d_LOG.root", c));
  }

 
  



}







void MakePlotNVTXDataMC(RooWorkspace* w, Float_t mass ) {

  Int_t ncat = NCAT;
  TString inDir = "";

  Float_t MASS(mass);

  // Luminosity:
  Float_t Lum = 19500.0;  
  RooRealVar lumi("lumi","lumi",Lum);
  w->import(lumi); 
  
  // Variables
  RooArgSet* ntplVars = defineVariables();
  RooRealVar* nvtx = w->var("nvtx");
   
  int iMass = abs(mass);      
 // common preselection cut
  TString mainCut = "";   //130-500
 
  
  //**********DATA***************//
  // create tree
  TFile file("histograms_CMS-HGG_24072013.root");
  TTree* dataTree = (TTree*) file.Get("Data");
   
  //**********G+jets***************//
  // create tree

  TChain* gjTree = new TChain();
  gjTree->Add("histograms_CMS-HGG_24072013.root/gjet_20_8TeV_pf");
  gjTree->Add("histograms_CMS-HGG_24072013.root/gjet_40_8TeV_pf"); 

  TChain* qcdTree = new TChain();
  qcdTree->Add("histograms_CMS-HGG_24072013.root/qcd_30_8TeV_pf");
  qcdTree->Add("histograms_CMS-HGG_24072013.root/qcd_40_8TeV_pf"); 
 
  //**********DIPHOTJET***************//
  // create tree

  TChain* diphotjTree = new TChain();
  diphotjTree->Add("histograms_CMS-HGG_24072013.root/diphojet_8TeV");
   

  //**********DIPHOT***************//
  // create tree

  TChain* diphotTree = new TChain();
  diphotTree->Add("histograms_CMS-HGG_24072013.root/dipho_Box_25_8TeV");
  diphotTree->Add("histograms_CMS-HGG_24072013.root/dipho_Box_250_8TeV");
   

  

  TH1F* h_data[NCAT+1];
  TH1F*  h_gj[NCAT+1];
  TH1F*  h_qcd[NCAT+1];
  TH1F*  h_diphot[NCAT+1];
  TH1F*  h_diphotj[NCAT+1];
  TH1F* h_sum[NCAT+1]; 
 
  
  TCanvas* ctmp = new TCanvas("ctmp","nvtx Background Categories",1);
  Int_t nbin = 50;
  Double_t  min = 0;
  Double_t  max = 50;

   for (int c=0; c<NCAT+1; ++c) {


     h_data[c]= new TH1F(TString::Format("h_data_cat%d",c), TString::Format("h_data_cat%d",c), nbin, min, max);
     h_gj[c]= new TH1F(TString::Format("h_gj_cat%d",c), TString::Format("h_gj_cat%d",c), nbin, min, max);
     h_qcd[c]= new TH1F(TString::Format("h_qcd_cat%d",c), TString::Format("h_qcd_cat%d",c), nbin, min, max);
     h_diphot[c]= new TH1F(TString::Format("h_diphot_cat%d",c), TString::Format("h_diphot_cat%d",c), nbin, min, max);
     h_diphotj[c]= new TH1F(TString::Format("h_diphotj_cat%d",c), TString::Format("h_diphotj_cat%d",c), nbin, min, max);



    // 1)  prime 4 cat livia
     if (c==0){//&&(PhotonsMass<178 || PhotonsMass >402) 
   
      dataTree->Draw("nvtx>>h_data_cat0", "("+mainCut+" (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 ))*1");
      gjTree->Draw("nvtx>>h_gj_cat0", "("+mainCut+"(abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 ))*evweight");
      qcdTree->Draw("nvtx>>h_qcd_cat0", "("+mainCut+"(abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 ))*evweight");
      diphotTree->Draw("nvtx>>h_diphot_cat0","("+ mainCut+"(abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 ))*evweight");
      diphotjTree->Draw("nvtx>>h_diphotj_cat0","("+ mainCut+"(abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 ))*evweight");
    
     }
     if (c==1){//&&(nvtx<178 || nvtx >402)&&

     dataTree->Draw("nvtx>>h_data_cat1", "("+mainCut+" (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9<0.94 || ph2_r9<0.94 )  )*1");
     gjTree->Draw("nvtx>>h_gj_cat1", "("+mainCut+" (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9<0.94 || ph2_r9<0.94 ))*evweight ");
     qcdTree->Draw("nvtx>>h_qcd_cat1", "("+mainCut+" (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9<0.94 || ph2_r9<0.94 ))*evweight ");
      diphotTree->Draw("nvtx>>h_diphot_cat1", "("+mainCut+" (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9<0.94 || ph2_r9<0.94 ) )*evweight");
     diphotjTree->Draw("nvtx>>h_diphotj_cat1", "("+mainCut+" (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9<0.94 || ph2_r9<0.94 ) )*evweight");

       }
     if (c==2){//&&(nvtx<178 || nvtx >402)

     dataTree->Draw("nvtx>>h_data_cat2", "("+mainCut+" (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 && ph2_r9>0.94 ))*1");
     gjTree->Draw("nvtx>>h_gj_cat2", "("+mainCut+" (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 && ph2_r9>0.94 ))*evweight");
     qcdTree->Draw("nvtx>>h_qcd_cat2", "("+mainCut+" (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 && ph2_r9>0.94 ))*evweight");
     diphotTree->Draw("nvtx>>h_diphot_cat2", "("+mainCut+" (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 && ph2_r9>0.94 ))*evweight");
     diphotjTree->Draw("nvtx>>h_diphotj_cat2", "("+mainCut+" (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 && ph2_r9>0.94 ))*evweight");
    
   }
     if (c==3){//&&(nvtx<178 || nvtx >402)

     dataTree->Draw("nvtx>>h_data_cat3", "("+mainCut+" (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 || ph2_r9>0.94 ))*1");
     gjTree->Draw("nvtx>>h_gj_cat3", "("+mainCut+" (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 || ph2_r9>0.94 ))*evweight ");
     qcdTree->Draw("nvtx>>h_qcd_cat3", "("+mainCut+" (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 || ph2_r9>0.94 ))*evweight ");
      diphotTree->Draw("nvtx>>h_diphot_cat3","("+ mainCut+" (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 || ph2_r9>0.94 ) )*evweight");
      diphotjTree->Draw("nvtx>>h_diphotj_cat3","("+ mainCut+" (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 || ph2_r9>0.94 ))*evweight ");
    
   }
     if (c==4){//&&(nvtx<178 || nvtx >402)

    dataTree->Draw("nvtx>>h_data_cat4","("+ mainCut+"1)*1");
    gjTree->Draw("nvtx>>h_gj_cat4", "("+mainCut+"1)*evweight");
    qcdTree->Draw("nvtx>>h_qcd_cat4", "("+mainCut+"1)*evweight");
    diphotTree->Draw("nvtx>>h_diphot_cat4","("+ mainCut+"1)*evweight");
    diphotjTree->Draw("nvtx>>h_diphotj_cat4","("+ mainCut+"1)*evweight");
    

   }

 
     //fisso la normalizzazione di gjets a quella gj+qcd xkè la shape di qcd _pf fa schifo, quindi prenidamoq uella di gj _pf
     Double_t qcdInt = h_qcd[c]->Integral();
     Double_t gJetIntFixed = qcdInt+h_gj[c]->Integral();
     h_gj[c]->Scale(gJetIntFixed/h_gj[c]->Integral());   
     h_gj[c]->SetFillColor(kAzure+8);
     h_gj[c]->Sumw2();
   
     h_diphotj[c]->SetFillColor(kTeal+9);
     h_diphotj[c]->Sumw2();

   //sommo i due diphot samples
   h_diphot[c]->Add(h_diphotj[c]);
   h_diphot[c]->SetFillColor(kSpring+7);
   h_diphot[c]->Sumw2();


   
   h_sum[c] = (TH1F*) h_gj[c]->Clone();
   h_sum[c]->Add(h_diphot[c]);
   
   h_sum[c]->SetFillColor(kBlack);
   h_sum[c]->SetFillStyle(3004);
   h_sum[c]->SetMarkerSize(0);
   

   THStack hs("hs","hs");  
   
   hs.Add(h_gj[c]);
   hs.Add(h_diphot[c]); 
   
   std::cout<<h_sum[c]->Integral()<<std::endl;
   std::cout<<h_diphot[c]->Integral()<<" %: "<<h_diphot[c]->Integral()/h_sum[c]->Integral()<<std::endl;
   std::cout<<h_gj[c]->Integral()<<" %: "<<h_gj[c]->Integral()/h_sum[c]->Integral()<<std::endl;
   std::cout<<h_data[c]->Integral()<<std::endl;
   
 
    ctmp->cd();
      h_data[c]->Sumw2();
    h_sum[c]->Sumw2();
    TH1F* h1_ratio1 = (TH1F*)h_data[c]->Clone();
    TH1F* h1_ratio2 = (TH1F*)h_sum[c]->Clone();
 
    
    //-------pad 1-------//
    TPad * pad1 = new TPad("pad1", "pad1",0.01,0.13,1.,1.);  
    
    //pad1->SetTopMargin(0.1);
    //pad1->SetBottomMargin(0.01);
    pad1->SetRightMargin(0.1);
    
    pad1->SetLogy();
    pad1->Draw();
    pad1->cd();
   
   
    h_data[c]->SetMarkerSize(0.7);
    h_data[c]->Draw("pe");
    
    h_data[c]->GetYaxis()->SetTitle(TString::Format("Events/%.2f", (max-min)/nbin));
    h_data[c]->GetXaxis()->SetTitle("nvtx");
    //  h_data[c]->GetYaxis()->SetRangeUser(0., 28000);
    hs.Draw("histsame");
    h_data[c]->Draw("pesame");
    h_sum[c]->Draw("E2same");
    std::cout<<"--------_> "<<h_data[c]->Integral()<<std::endl; 
  

    TLegend *leg1;
    if(c!=4)leg1 = new TLegend(0.6075,0.7536441,0.8575,0.9340678, TString::Format("Category %d",c), "brNDC");
    else leg1 = new TLegend(0.6075,0.7536441,0.8575,0.9340678, TString::Format("All Categories Combined",c), "brNDC");
    leg1->AddEntry(h_data[c],"Data","PE");
    leg1->AddEntry(h_gj[c],"prompt + fake","F");
    leg1->AddEntry(h_diphot[c],"prompt + prompt", "F");
    leg1->AddEntry(h_sum[c], "Bkg Err", "F");
    
   
    
    leg1->SetTextSize(0.035);
    leg1->SetTextFont(42);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
    leg1->Draw("same");
 

    TPaveText* label_cms = get_labelCMS(0, "2012", false);
    TPaveText* label_sqrt = get_labelSqrt(0);
    label_cms->Draw("same");
    label_sqrt->Draw("same");

    pad1->SetLogy(0);
    pad1->RedrawAxis();

    ctmp->cd();

    //-------pad 2------//
    TPad * pad2 = new TPad("pad2", "pad2",0.01,0.001,1.,0.2);
    pad2->SetGrid();
    
    //pad2->SetTopMargin(0.01);
    pad2->SetBottomMargin(0.4);
    pad2->SetRightMargin(0.1);
    pad2->Draw();
    pad2->cd();

    Double_t xmax = h1_ratio1->GetXaxis()->GetXmax();
    Double_t xmin = h1_ratio1->GetXaxis()->GetXmin();
    TLine* line = new TLine(xmin,1.,xmax,1.);
  

    h1_ratio1->SetStats(0);
    h1_ratio1->Divide(h_sum[c]);
    h1_ratio1->SetMarkerStyle(20);
    h1_ratio1->SetMarkerSize(1.1);
    //  h1_ratio1->GetXaxis()->SetTitle(xAxis.c_str());
    h1_ratio1->GetYaxis()->SetRangeUser(0., 2.); // cwr zoom
    h1_ratio1->GetYaxis()->SetNdivisions(4,false);
    h1_ratio1->GetYaxis()->SetTitle("Data/Bkg.");
    h1_ratio1->GetYaxis()->SetTitleFont(42);
    h1_ratio1->GetXaxis()->SetTitleFont(42);
    h1_ratio1->GetXaxis()->SetTitle("nvtx");
    
    h1_ratio1->GetXaxis()->SetTitleSize(0.23);
    h1_ratio1->GetXaxis()->SetLabelSize(0.16);
    h1_ratio1->GetYaxis()->SetLabelSize(0.16);
    h1_ratio1->GetYaxis()->SetTitleSize(0.15);
    h1_ratio1->GetYaxis()->SetTitleOffset(0.25);
    h1_ratio1->GetXaxis()->SetTitleOffset(0.5);

    
    
    for(int j = 0;j<=h1_ratio1->GetNbinsX();j++){
      if(h_sum[c]->GetBinContent(j))  h1_ratio1->SetBinError(j,h_data[c]->GetBinError(j)/h_sum[c]->GetBinContent(j));
      else h1_ratio1->SetBinError(j,0.);
    }
    h1_ratio1->Draw("PEX0");
    
    h1_ratio2->Divide(h_sum[c]);
    for(int j = 0;j<=h1_ratio2->GetNbinsX();j++){
      if(h_sum[c]->GetBinContent(j)) h1_ratio2->SetBinError(j,h_sum[c]->GetBinError(j)/h_sum[c]->GetBinContent(j));
      else h1_ratio2->SetBinError(j,0.);
    }
    h1_ratio2->Draw("E2same");
    
    line->SetLineWidth(1.);
    line->Draw("same");


    ctmp->SaveAs("plots/DATA_MC_NVTX_"+TString::Format("cat%d.png", c));
    ctmp->SaveAs("plots/DATA_MC_NVTX_"+TString::Format("cat%d.pdf", c));
    ctmp->SaveAs("plots/DATA_MC_NVTX_"+TString::Format("cat%d.root", c));
  
    pad1->SetLogy(1);
    ctmp->SaveAs("plots/DATA_MC_NVTX_"+TString::Format("cat%d_LOG.png", c));
    ctmp->SaveAs("plots/DATA_MC_NVTX_"+TString::Format("cat%d_LOG.pdf", c));
    ctmp->SaveAs("plots/DATA_MC_NVTX_"+TString::Format("cat%d_LOG.root", c));
  }

 
  



}




void MakePlotRooKeys(int cat){
  TString wsDir = "BiasStudy/workspaces/";
  TFile* f = new TFile(""+wsDir+"HighMass-hgg.RooKeysPdfMCBkg_8TeV.root");
  RooWorkspace *w = f->Get("w_bias");
  
  RooAbsPdf* r2 = (RooAbsPdf*) *w->pdf(TString::Format("BkgMCKeyPdf_bw2_cat%",cat));
  /*RooKeysPdf* r2 = (RooKeysPdf*) *w->pdf(TString::Format("BkgMCKeyPdf_bw3_cat%",cat));
  RooKeysPdf* r2 = (RooKeysPdfx*) *w->pdf(TString::Format("BkgMCKeyPdf_bw3_cat%",cat));*/
  RooDataSet* mc =(RooDataSet*) *w->data(TString::Format("BkgMCWeight_cat%d",cat));

  TCanvas* ctmp = new TCanvas("ctmp","PhotonsMass Background RooKeysPdf",0,0,500,500);
  RooPlot* plotPhotonsMassBkgMC;
  Int_t nBinsMass(320);
  plotPhotonsMassBkgMC = w->var("PhotonsMass")->frame(130, 1000,nBinsMass);

  r2->plotOn(plotPhotonsMassBkgMC,"L", LineColor(kBlack), LineWidth(2));
  /* r3->plotOn(plotPhotonsMassBkgMC,"L", LineColor(8), LineWidth(2));
  r4->plotOn(plotPhotonsMassBkgMC,"L", LineColor(kOrange), LineWidth(2));*/
  mc->plotOn(plotPhotonsMassBkgMC);
  plotPhotonsMassBkgMC->SetAxisRange(0.001,plotPhotonsMassBkgMC->GetMaximum()*30.,"Y");
    

  TLegend *leg1 = new TLegend(0.2175,0.8236441,0.6575,0.9240678, TString::Format("Category %d",cat), "brNDC");
  leg1->AddEntry(plotPhotonsMassBkgMC->getObject(0),"Bkg MC","LPE");
  TLegend *leg2 = new TLegend(0.4375,0.7236441,0.85,0.9240678, TString::Format("RooKeysPdf",cat), "brNDC");
  
  leg2->AddEntry(plotPhotonsMassBkgMC->getObject(1),"Bw x 2","L");
  leg2->AddEntry(plotPhotonsMassBkgMC->getObject(2),"Bw x 3","L");
  leg2->AddEntry(plotPhotonsMassBkgMC->getObject(3),"Bw x 4","L");
  
  leg1->SetTextSize(0.035);
  leg1->SetTextFont(42);
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  
  leg2->SetTextSize(0.035);
  leg2->SetTextFont(42);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  
  TPaveText* label_cms = get_labelCMS(0, "2012", true);
  TPaveText* label_sqrt = get_labelSqrt(0);
  
  ctmp->cd();
  
  plotPhotonsMassBkgMC->Draw();  
  leg1->Draw("same");
  leg2->Draw("same");
  
  label_cms->Draw("same");
  label_sqrt->Draw("same");

  
  ctmp->SetLogy(1);
  
  ctmp->SaveAs("plots/BKG_MC_rookeyspdf_"+TString::Format("cat%d_LOG.png", cat));
  ctmp->SaveAs("plots/BKG_MC_rookeyspdf_"+TString::Format("cat%d_LOG.root", cat));
}




