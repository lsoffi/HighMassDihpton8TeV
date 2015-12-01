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
RooFitResult*  BkgModelFitDiJetFunc(RooWorkspace*, Bool_t, Float_t, bool);
RooFitResult*  BkgModelFitExpolFunc(RooWorkspace*, Bool_t, Float_t, bool);
RooHistFunc* getRooHistFunc(int cat, RooRealVar* var);
void SetConstantParams(const RooArgSet* params);



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
  RooRealVar* ph1_eta = new RooRealVar("ph1_eta", "eta(g1)",-10,10,"");
  RooRealVar* ph2_eta = new RooRealVar("ph2_eta", "eta(g2)",-10,10,"");
  RooRealVar* ph1_r9 = new RooRealVar("ph1_r9", "R9(g1)",-10,10,"");
  RooRealVar* ph2_r9 = new RooRealVar("ph2_r9", "R9(g2)",-10,10,"");
  RooRealVar* evweight = new RooRealVar("evweight","Reweightings",0,1000,"");
  RooRealVar* nvtx = new RooRealVar("nvtx", "nvtx", 0, 50, "");
  //  RooRealVar* btagCategory = new RooRealVar("btagCategory","event category",0.9,2.1,"") ;
  //  RooRealVar* newweight = new RooRealVar("newweight","Reweightings",0., 1000,"");
  RooArgSet* ntplVars = new RooArgSet(*PhotonsMass,*PhotonsMassTrue, *ph1_eta, *ph2_eta, *ph1_r9, *ph2_r9, *evweight, *nvtx);
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



TTree *dataset2tree(RooDataSet *dataset, RooArgSet* args, Double_t newW){

  
  RooArgList argList(*args);
  argList.Print();
  Double_t variables[50];
  Long64_t nEntries= dataset->numEntries();
  
  TTree *tree = new TTree("tree","tree");
  tree->SetDirectory(0);
  TIterator *it1=NULL; 
  it1 = argList.createIterator();
  int index1=0;
  for(RooRealVar *var = (RooRealVar *) it1->Next(); var!=NULL;
      var = (RooRealVar *) it1->Next(),index1++){
    TString name(var->GetName());
    name.ReplaceAll("-","_");
    tree->Branch(name, &(variables[index1]), name+"/D");
  }

  Double_t newweight;
  tree->Branch("newweight", &newweight, "newweight/D");
  tree->Print();
 
  for(Long64_t jentry=0; jentry<nEntries; jentry++){
   
    TIterator *it=NULL; 
    RooArgList argList1(*(dataset->get(jentry)));
    it = argList1.createIterator();
    int index=0;
    for(RooRealVar *var = (RooRealVar *) it->Next(); var!=NULL;
	var = (RooRealVar *) it->Next(), index++){
      variables[index]=var->getVal();
      
    }

    newweight = dataset->weight()*newW;
   
    delete it;
    tree->Fill();
  }
  tree->ResetBranchAddresses();
  return tree;
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
  TString fileBkgName("HighMass-hgg.inputbkg");
  
  TString card_name("HighMass-hgg_models_Bkg_8TeV_test.rs");
  HLFactory hlf("HLFactory", card_name, false);
  RooWorkspace* w = hlf.GetWs();
 
  RooFitResult* fitresults;



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

 

  cout << endl; cout << "Now AddBkgData" << endl;
  AddBkgData(w, mass);         
  

  cout << endl; cout << "Now BkgModelFit" << endl;    
  bool blind = false; 
  bool dobandsHere= false;
  MakeRooKeysPDFMCBkg( w, mass, true);
 
  return;
}


// Data dataset
void AddBkgData(RooWorkspace* w, Float_t mass) {

  // initializations
  Int_t ncat = NCAT;

  Float_t minMassFit, maxMassFit;
  minMassFit = MINmass;
  maxMassFit = MAXmass;

  // retrieve the data tree; no common preselection cut applied yet; 
  TString inDir = "";
  // TFile dataFile("histograms_CMS-HGG_19032013.root");
  TFile dataFile("histograms_CMS-HGG_26022014_DATA.root");//mass range fino a 1000
  TTree* dataTree = (TTree*) dataFile.Get("Data");

  // Variables
  RooArgSet* ntplVars = defineVariables();

  // common preselection cut
  TString mainCut = TString::Format("PhotonsMass>=(%.1f) && PhotonsMass<=(%.1f)", minMassFit, maxMassFit);   // livia
  
  // Create dataset
  RooDataSet Data("Data","dataset",dataTree,*ntplVars,mainCut,"evweight");
  cout << endl;
  cout << "Data, everything: " << endl;
  Data.Print("v");
  cout << "---- nX:  " << Data.sumEntries() << endl;
  cout << endl;

  // split into NCAT categories
  RooDataSet* dataToFit[NCAT];  
  for (int c=0; c<4; ++c) {

  // 1)  prime 4 cat livia
    if (c==0) dataToFit[c] = (RooDataSet*) Data.reduce(*w->var("PhotonsMass"),mainCut+TString::Format("&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 )"));
    if (c==1) dataToFit[c] = (RooDataSet*) Data.reduce(*w->var("PhotonsMass"),mainCut+TString::Format("&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9<0.94 || ph2_r9<0.94 ) "));
    if (c==2) dataToFit[c] = (RooDataSet*) Data.reduce(*w->var("PhotonsMass"),mainCut+TString::Format("&& (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 && ph2_r9>0.94 )"));
    if (c==3) dataToFit[c] = (RooDataSet*) Data.reduce(*w->var("PhotonsMass"),mainCut+TString::Format("&& (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9<0.94 || ph2_r9<0.94 ) "));   


    cout << endl; cout << "for category = " << c << endl;
    dataToFit[c]->Print("v");
    cout << "---- nX:  " << dataToFit[c]->sumEntries() << endl;

     w->import(*dataToFit[c],Rename(TString::Format("Data_cat%d",c)));
  }

  cout << "data, no split" << endl;
  // Create full data set without categorization
  RooDataSet* data = (RooDataSet*) Data.reduce(*w->var("PhotonsMass"),mainCut);
  w->import(*data, Rename("Data"));
  cout << endl;
 
  data->Print("v");
  cout << "---- nX:  " << data->sumEntries() << endl; 
}








RooFitResult* BkgModelFitExpPARFunc(RooWorkspace* w, Bool_t dobands, Float_t mass,Int_t c,  bool blind) {

  Int_t ncat = NCAT;

  RooDataSet* data;
 
  RooFitResult* fitresult;


  RooPlot* plotPhotonsMassBkg;

  Float_t minMassFit, maxMassFit;
 
    minMassFit = MINmass;
    maxMassFit = MAXmass;
  
  // Fit data with background pdf for data limit
  RooRealVar* PhotonsMass = w->var("PhotonsMass");  
  PhotonsMass->setUnit("GeV");
 
  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
  
  
    data = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    PhotonsMass->setRange("bkg range", MINmass, MAXmass);
    
    // fit con expol 
    RooFormulaVar *p1mod= new RooFormulaVar(TString::Format("par1ExpPAR_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_ExpPAR1_cat%d",c)));
    RooFormulaVar *p2mod= new RooFormulaVar(TString::Format("par2ExpPAR_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_ExpPAR2_cat%d",c)));
    
    RooAbsPdf* PhotonsMassBkg = new RooGenericPdf(TString::Format("PhotonsMassBkg_cat%d",c), "exp(-@1*@0)*pow(@0, @2)", RooArgList(*PhotonsMass, *p1mod, *p2mod));
    
    fitresult = PhotonsMassBkg->fitTo(*data, Range(minMassFit,maxMassFit),RooFit::FitOptions("MHTER"), SumW2Error(kTRUE), Save(kTRUE));   
   
    std::cout<<TString::Format("******************************** Background Fit results mass %f cat %d ***********************************", mass, c)<<std::endl;
    fitresult->Print("V");
    w->import(*PhotonsMassBkg);
 

    //************************************************
    // Plot PhotonsMass background fit results per categories 
    TCanvas* ctmp = new TCanvas("ctmp","PhotonsMass Background ",0,0,500,500);
    Int_t nBinsMass(60);
    plotPhotonsMassBkg = PhotonsMass->frame(minMassFit, maxMassFit,nBinsMass);
  
    data->plotOn(plotPhotonsMassBkg,RooFit::Invisible());    
   
    PhotonsMassBkg->plotOn(plotPhotonsMassBkg,LineColor(kBlue),Range(minMassFit,maxMassFit),NormRange("bkg range")); 
   
    double chi2 = plotPhotonsMassBkg->chiSquare(3);
    Int_t ndof = nBinsMass-2;
    std::cout<<"------> "<< ndof<<std::endl;
    double prob = TMath::Prob(chi2*ndof, ndof);
    std::cout<<prob<<std::endl;

    
 

    blind=false;
    if( blind ) {
     
      RooDataSet* data_down = (RooDataSet*) data.reduce(*w->var("PhotonsMass"),"PhotonsMass < 173.5");
      RooDataSet* data_up = (RooDataSet*) data.reduce(*w->var("PhotonsMass"),"PhotonsMass >850");

      data_up->plotOn(plotPhotonsMassBkg);    
      data_down->plotOn(plotPhotonsMassBkg); 


   
    } else {
      data->plotOn(plotPhotonsMassBkg);    
      } 
       
   
    plotPhotonsMassBkg->GetXaxis()->SetTitle("m_{#gamma #gamma} [GeV]");
    plotPhotonsMassBkg->SetAxisRange(0.001,plotPhotonsMassBkg->GetMaximum()*1.5,"Y");
    plotPhotonsMassBkg->Draw();  

    TLegend *legdata = new TLegend(0.3790323,0.7775424,0.6290323,0.9279661, TString::Format("Category %d",c), "brNDC");
    legdata->AddEntry(plotPhotonsMassBkg->getObject(2),"Data","LPE");
    legdata->AddEntry(plotPhotonsMassBkg->getObject(1),"Parametric Model: ExpPAR","L");
  
    legdata->SetTextSize(0.035);
    legdata->SetTextFont(42);
    legdata->SetBorderSize(0);
    legdata->SetFillStyle(0);
    legdata->Draw("same");

    TPaveText* label_cms = get_labelCMS(0, "2012", false);
    TPaveText* label_sqrt = get_labelSqrt(0);
    label_cms->Draw("same");
    label_sqrt->Draw("same");

    //write down the chi2 of the fit on the
      
    TPaveText* label_chi2 = new TPaveText(0.5524194,0.6419492,0.796371,0.7690678,"brNDC");
    label_chi2->SetFillColor(kWhite);
    label_chi2->SetTextSize(0.035);
    label_chi2->SetTextFont(42);
    label_chi2->SetTextAlign(31); // align right
    label_chi2->AddText(TString::Format("Fit chi square/dof = %.3f", chi2));
    label_chi2->AddText(TString::Format("Chi square Prob = %.3f", prob));
    label_chi2->Draw("same");
    
 
    int massI(mass);
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPPAR_M%d.png",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPPAR_M%d.pdf",c,massI));

    ctmp->SetLogy();
    plotPhotonsMassBkg->SetAxisRange(1.3,plotPhotonsMassBkg->GetMaximum()*1.5,"Y");
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPPAR_LOG_M%d.png",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPPAR_LOG_M%d.pdf",c,massI));


  RooFitResult* r;

  return r;
}







void SetConstantParams(const RooArgSet* params) {

  cout << endl; cout << "Entering SetConstantParams" << endl;
  TIterator* iter(params->createIterator());
  for (TObject *a = iter->Next(); a != 0; a = iter->Next()) {
    RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
    if (rrv) { rrv->setConstant(true); std::cout << " " << rrv->GetName(); }
  }  
}







// Signal Data Set
void MakeRooKeysPDFMCBkg(RooWorkspace* w, Float_t mass, Bool_t isMirror) {


  TString wsDir = "BiasStudy/workspaces/"+filePOSTfix;
  
  RooWorkspace *wBias = new RooWorkspace("w_bias","w_bias");  


  Int_t ncat = NCAT;
  TString inDir = "";

  Float_t MASS(mass);

  // Luminosity:
  Float_t Lum = 19500.0;  
  RooRealVar lumi("lumi","lumi",Lum);
  wBias->import(lumi); 
  
  // Variables
  RooArgSet* ntplVars = defineVariables();
  ntplVars->Print("V");
  RooArgSet* ntplVars_newweight = defineVariables_newWeight();
  int iMass = abs(mass);  
  RooRealVar* PhotonsMass = w->var("PhotonsMass");  
  PhotonsMass->setRange(MINmass, MAXmass);  
  TFile sigFile1("histograms_CMS-HGG_24072013.root");   //ggh prod mode tree livia
  
  // common preselection cut
  TString mainCut = "PhotonsMass>=100 && PhotonsMass<=1000";   // livia
  
  //get sumEntries of QCD 
  TChain* qcdTree=new TChain();
  qcdTree->Add("histograms_CMS-HGG_24072013.root/qcd_30_8TeV_pf");
  qcdTree->Add("histograms_CMS-HGG_24072013.root/qcd_40_8TeV_pf");

 
  RooDataSet qcdMCWeighted("qcdMCWeighted","MC qcd weighted",qcdTree,*ntplVars,mainCut,"evweight");
  cout << endl;
  cout << "qcdMCWeighted" << endl;
  qcdMCWeighted.Print("v");
  Double_t qcdInt_ = qcdMCWeighted.sumEntries();
  cout << "---- nX: qcd Int " << qcdInt_ << endl; 


  //get sumEntries of QCD 
  TChain* gjTree=new TChain();
  gjTree->Add("histograms_CMS-HGG_24072013.root/gjet_20_8TeV_pf");
  gjTree->Add("histograms_CMS-HGG_24072013.root/gjet_40_8TeV_pf");

 

  RooDataSet gjMCWeighted("gjMCWeighted","MC gj weighted",gjTree,*ntplVars,mainCut,"evweight");
  cout << endl;
  cout << "gjMCWeighted" << endl;
  gjMCWeighted.Print("v");
  Double_t gjInt_ = gjMCWeighted.sumEntries();
  cout << "---- nX: gj Int " << gjInt_ << endl; 



 
  //chain summing up all production modes
  TChain* sigTree1  = new TChain();
  sigTree1->Add("histograms_CMS-HGG_24072013.root/gjet_20_8TeV_pf");
  sigTree1->Add("histograms_CMS-HGG_24072013.root/gjet_40_8TeV_pf");
  sigTree1->Add("histograms_CMS-HGG_24072013.root/diphojet_8TeV");
  sigTree1->Add("histograms_CMS-HGG_24072013.root/dipho_Box_25_8TeV");
  //sigTree1->Add("histograms_CMS-HGG_24072013.root/dipho_Box_250_8TeV");
  

  sigTree1->SetTitle("sigTree1");
  sigTree1->SetName("sigTree1");
  
 


  
  // Create signal dataset composed with different productions, the weight is already applied in our ntuples
  RooDataSet BkgMCWeighted("BkgMCWeighted","MC BKG weighted",sigTree1,*ntplVars,mainCut,"evweight");
  cout << endl;
  cout << "BkgMCWeighted" << endl;
  BkgMCWeighted.Print("v");
  cout << "---- nX:  " << BkgMCWeighted.sumEntries() << endl; 
  
  // apply a common preselection cut; split in categories
  cout << endl;
 
  RooDataSet* data[NCAT];
  RooDataSet* qcdMC[NCAT];
  RooDataSet* gjMC[NCAT];
  TTree* gjAndQcdTree[NCAT];
  RooDataSet* QcdToGjMC[NCAT];

  RooDataSet* BkgMC[NCAT];
  TTree* BkgMCcopyTree[NCAT];
  RooDataSet* BkgMCcopy[NCAT];
  RooKeysPdf* BkgMCKeyPdf[NCAT];
  RooKeysPdf* BkgMCKeyPdf_bw2[NCAT];
  RooKeysPdf* BkgMCKeyPdf_bw3[NCAT];
  RooKeysPdf* BkgMCKeyPdf_bw4[NCAT];
  RooKeysPdf* BkgMCKeyPdf_bw2_noMirr[NCAT];


  RooDerivative* BkgMCKeyPdf_D1[NCAT];
  RooDerivative* BkgMCKeyPdf_bw2_D1[NCAT];
  RooDerivative* BkgMCKeyPdf_bw3_D1[NCAT];
  RooDerivative* BkgMCKeyPdf_bw4_D1[NCAT];
  RooDerivative* BkgMCKeyPdf_bw2_D1_noMirr[NCAT];

  RooDerivative* BkgMCKeyPdf_D2[NCAT];
  RooDerivative* BkgMCKeyPdf_bw2_D2[NCAT];
  RooDerivative* BkgMCKeyPdf_bw3_D2[NCAT];
  RooDerivative* BkgMCKeyPdf_bw4_D2[NCAT];
  RooDerivative* BkgMCKeyPdf_bw2_D2_noMirr[NCAT];

  RooPlot* plotPhotonsMassBkgMC[NCAT];
  RooPlot* plotPhotonsMassBkgMC_D1[NCAT];
  RooPlot* plotPhotonsMassBkgMC_D2[NCAT];
  TCanvas* ctmp = new TCanvas("ctmp","PhotonsMass Background Categories",0,0,500,500);
  Int_t nBinsMass(200);
  Double_t  minMassFit = MINmass;
  Double_t  maxMassFit = MAXmass;

  w->Print("V");
  //  RooArgSet* argset_ = new RooArgSet(*w->var("PhotonsMass"), *w->var("evweight"));

  for (int c=2; c<3; ++c) {


    data[c] = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));

    // 1)  prime 4 cat livia
  
   
    //reduce QCD dataset
    if (c==0) qcdMC[c] = (RooDataSet*) qcdMCWeighted.reduce(*ntplVars, mainCut+TString::Format("&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 )"));
    if (c==1) qcdMC[c] = (RooDataSet*) qcdMCWeighted.reduce(*ntplVars,mainCut+TString::Format("&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9<0.94 || ph2_r9<0.94 ) "));
    if (c==2) qcdMC[c] = (RooDataSet*) qcdMCWeighted.reduce(*ntplVars,mainCut+TString::Format("&& (abs(ph1_eta)>1.4442 || abs(ph2_eta)>1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 )"));
    if (c==3) qcdMC[c] = (RooDataSet*) qcdMCWeighted.reduce(*ntplVars,mainCut+TString::Format("&& (abs(ph1_eta)>1.4442 || abs(ph2_eta)>1.4442) && (ph1_r9<0.94 || ph2_r9<0.94 ) "));

 
    //reduce gj dataset
    if (c==0) gjMC[c] = (RooDataSet*) gjMCWeighted.reduce(*ntplVars, mainCut+TString::Format("&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 )"));
    if (c==1) gjMC[c] = (RooDataSet*) gjMCWeighted.reduce(*ntplVars,mainCut+TString::Format("&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9<0.94 || ph2_r9<0.94 ) "));
    if (c==2) gjMC[c] = (RooDataSet*) gjMCWeighted.reduce(*ntplVars,mainCut+TString::Format("&& (abs(ph1_eta)>1.4442 || abs(ph2_eta)>1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 )"));
    if (c==3) gjMC[c] = (RooDataSet*) gjMCWeighted.reduce(*ntplVars,mainCut+TString::Format("&& (abs(ph1_eta)>1.4442 || abs(ph2_eta)>1.4442) && (ph1_r9<0.94 || ph2_r9<0.94 ) "));
 
    gjMC[c]->Print();



    Float_t qcdInt = qcdMC[c]->sumEntries();
    Float_t gjInt = gjMC[c]->sumEntries();
    Float_t gjEntries = gjMC[c]->numEntries();
    std::cout<<"qcd: "<<qcdInt<<" gj: "<<gjInt<<" gjEntries: "<<gjEntries<<std::endl;
    Float_t qcdWeight = qcdInt/gjInt;
 
    
    gjAndQcdTree[c] = (TTree*)dataset2tree(gjMC[c], ntplVars, qcdWeight);
    
   
    QcdToGjMC[c] = new RooDataSet(gjMC[c]->GetName(),gjMC[c]->GetTitle(),gjAndQcdTree[c],*ntplVars_newweight, mainCut, "newweight"  );  
    QcdToGjMC[c]->Print("");
  
    


    if (c==0) BkgMC[c] = (RooDataSet*) BkgMCWeighted.reduce(*ntplVars,mainCut+TString::Format("&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 )"));
    if (c==1) BkgMC[c] = (RooDataSet*) BkgMCWeighted.reduce(*ntplVars,mainCut+TString::Format("&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9<0.94 || ph2_r9<0.94 ) "));
    if (c==2) BkgMC[c] = (RooDataSet*) BkgMCWeighted.reduce(*ntplVars,mainCut+TString::Format("&& (abs(ph1_eta)>1.4442 || abs(ph2_eta)>1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 )"));
    if (c==3) BkgMC[c] = (RooDataSet*) BkgMCWeighted.reduce(*ntplVars,mainCut+TString::Format("&& (abs(ph1_eta)>1.4442 || abs(ph2_eta)>1.4442) && (ph1_r9<0.94 || ph2_r9<0.94 ) "));


    BkgMC[c]->Print();   
    BkgMCcopyTree[c] = (TTree*) dataset2tree(BkgMC[c], ntplVars, 1.);

    BkgMCcopy[c] =  new RooDataSet(BkgMC[c]->GetName(),BkgMC[c]->GetTitle(),BkgMCcopyTree[c] ,*ntplVars_newweight, mainCut, "newweight");     
    BkgMCcopy[c]->append(*QcdToGjMC[c]);
    BkgMCcopy[c]->Print();


   

    wBias->import(*BkgMCcopy[c],Rename(TString::Format("BkgMCWeight_cat%d",c)));
    
    cout << "cat " << c << ", BkgMC[c]: " << endl;
    BkgMCcopy[c]->Print("v");
    cout << "---- for category " << c << ", nX for [c]:  " << BkgMCcopy[c]->sumEntries() << endl; 
    cout << endl;


    isMirror = true;
    if(isMirror){
    
    //create the rookeyspdf
      BkgMCKeyPdf[c] = new RooKeysPdf(TString::Format("BkgMCKeyPdf_cat%d",c),TString::Format("BkgMCKeyPdf_cat%d",c),*PhotonsMass,*BkgMCcopy[c],RooKeysPdf::MirrorBoth );
      BkgMCKeyPdf_bw2[c] = new RooKeysPdf(TString::Format("BkgMCKeyPdf_bw2_cat%d",c),TString::Format("BkgMCKeyPdf_bw2_cat%d",c),*PhotonsMass,*BkgMCcopy[c],RooKeysPdf::MirrorBoth,2 );
      BkgMCKeyPdf_bw3[c] = new RooKeysPdf(TString::Format("BkgMCKeyPdf_bw3_cat%d",c),TString::Format("BkgMCKeyPdf_bw3_cat%d",c),*PhotonsMass,*BkgMCcopy[c],RooKeysPdf::MirrorBoth,3 );
      //  BkgMCKeyPdf_bw4[c] = new RooKeysPdf(TString::Format("BkgMCKeyPdf_bw4_cat%d",c),TString::Format("BkgMCKeyPdf_bw4_cat%d",c),*PhotonsMass,*BkgMCcopy[c],RooKeysPdf::MirrorBoth,4 );
      
    }else{

      //create the rookeyspdf
      BkgMCKeyPdf[c] = new RooKeysPdf(TString::Format("BkgMCKeyPdf_cat%d",c),TString::Format("BkgMCKeyPdf_cat%d",c),*PhotonsMass,*BkgMCcopy[c],RooKeysPdf::NoMirror );
      BkgMCKeyPdf_bw2[c] = new RooKeysPdf(TString::Format("BkgMCKeyPdf_bw2_cat%d",c),TString::Format("BkgMCKeyPdf_bw2_cat%d",c),*PhotonsMass,*BkgMCcopy[c],RooKeysPdf::NoMirror,2 );
      BkgMCKeyPdf_bw3[c] = new RooKeysPdf(TString::Format("BkgMCKeyPdf_bw3_cat%d",c),TString::Format("BkgMCKeyPdf_bw3_cat%d",c),*PhotonsMass,*BkgMCcopy[c],RooKeysPdf::NoMirror,3 );
      //  BkgMCKeyPdf_bw4[c] = new RooKeysPdf(TString::Format("BkgMCKeyPdf_bw4_cat%d",c),TString::Format("BkgMCKeyPdf_bw4_cat%d",c),*PhotonsMass,*BkgMCcopy[c],RooKeysPdf::NoMirror,4 );
    

    }


    wBias->import(*BkgMCKeyPdf_bw2[c]);
    //   wBias->import(*BkgMCKeyPdf_bw[c]);
    //  wBias->import(*BkgMCKeyPdf_bw3[c]);
    //wBias->import(*BkgMCKeyPdf_bw4[c]);
    
    plotPhotonsMassBkgMC[c] = PhotonsMass->frame(minMassFit, maxMassFit,nBinsMass);
    // data[c]->plotOn(plotPhotonsMassBkgMC[c]);
    //BkgMC[c]->plotOn(plotPhotonsMassBkgMC[c],"PE", MarkerColor(kRed), LineColor(kRed), MarkerSize(1.));
    BkgMCcopy[c]->plotOn(plotPhotonsMassBkgMC[c],"PE", MarkerColor(kRed), LineColor(kRed), MarkerSize(0.8));
    BkgMCKeyPdf[c]->plotOn(plotPhotonsMassBkgMC[c],"L", LineColor(kBlue), LineWidth(2));
    BkgMCKeyPdf_bw2[c]->plotOn(plotPhotonsMassBkgMC[c],"L", LineColor(kBlack), LineWidth(2));
    BkgMCKeyPdf_bw3[c]->plotOn(plotPhotonsMassBkgMC[c],"L", LineColor(8), LineWidth(2));
    //  BkgMCKeyPdf_bw4[c]->plotOn(plotPhotonsMassBkgMC[c],"L", LineColor(kOrange), LineWidth(2));

    plotPhotonsMassBkgMC[c]->SetAxisRange(0.001,plotPhotonsMassBkgMC[c]->GetMaximum()*30.,"Y");
    
    /*
    //create the rookeyspdf_D1

    BkgMCKeyPdf_D1[c] =   (RooDerivative*) BkgMCKeyPdf[c]->derivative(*PhotonsMass, 1, 0.01);
    BkgMCKeyPdf_bw2_D1[c] =  (RooDerivative*)   BkgMCKeyPdf_bw2[c]->derivative(*PhotonsMass, 1, 0.01);
    BkgMCKeyPdf_bw3_D1[c] =   (RooDerivative*) BkgMCKeyPdf_bw3[c]->derivative(*PhotonsMass, 1, 0.01);
    BkgMCKeyPdf_bw4_D1[c] =   (RooDerivative*) BkgMCKeyPdf_bw4[c]->derivative(*PhotonsMass, 1, 0.01);


    //create the rookeyspdf_D2
    BkgMCKeyPdf_D2[c] =   (RooDerivative*) BkgMCKeyPdf[c]->derivative(*PhotonsMass, 2, 0.01);
    BkgMCKeyPdf_bw2_D2[c] = (RooDerivative*)  BkgMCKeyPdf_bw2[c]->derivative(*PhotonsMass, 2,0.01);
    BkgMCKeyPdf_bw3_D2[c] = (RooDerivative*)   BkgMCKeyPdf_bw3[c]->derivative(*PhotonsMass, 2,0.01);
    BkgMCKeyPdf_bw4_D2[c] =  (RooDerivative*)  BkgMCKeyPdf_bw4[c]->derivative(*PhotonsMass, 2,0.01);

    
    plotPhotonsMassBkgMC_D1[c] = PhotonsMass->frame(minMassFit, maxMassFit);
    BkgMCKeyPdf_D1[c]->plotOn(plotPhotonsMassBkgMC_D1[c],"L", LineColor(kBlue), LineWidth(2));
    BkgMCKeyPdf_bw2_D1[c]->plotOn(plotPhotonsMassBkgMC_D1[c],"L", LineColor(kBlack), LineWidth(2));
    BkgMCKeyPdf_bw3_D1[c]->plotOn(plotPhotonsMassBkgMC_D1[c],"L", LineColor(8), LineWidth(2));
    BkgMCKeyPdf_bw4_D1[c]->plotOn(plotPhotonsMassBkgMC_D1[c],"L", LineColor(kOrange), LineWidth(2));

    plotPhotonsMassBkgMC_D1[c]->SetAxisRange(plotPhotonsMassBkgMC_D1[c]->GetMinimum()*1.3,0.15,"Y");   
    
   
    plotPhotonsMassBkgMC_D2[c] = PhotonsMass->frame(minMassFit, maxMassFit);
    BkgMCKeyPdf_D2[c]->plotOn(plotPhotonsMassBkgMC_D2[c],"", LineColor(kBlue), LineWidth(2));
    BkgMCKeyPdf_bw2_D2[c]->plotOn(plotPhotonsMassBkgMC_D2[c],"", LineColor(kBlack), LineWidth(2));
    BkgMCKeyPdf_bw3_D2[c]->plotOn(plotPhotonsMassBkgMC_D2[c],"", LineColor(8), LineWidth(2));
    BkgMCKeyPdf_bw4_D2[c]->plotOn(plotPhotonsMassBkgMC_D2[c],"", LineColor(kOrange), LineWidth(2));

    plotPhotonsMassBkgMC_D2[c]->SetAxisRange(plotPhotonsMassBkgMC_D2[c]->GetMinimum()*1.3,0.04,"Y");
    
    */
 
    TLegend *leg1 = new TLegend(0.2175,0.8236441,0.6575,0.9240678, TString::Format("Category %d",c), "brNDC");
    leg1->AddEntry(plotPhotonsMassBkgMC[c]->getObject(0),"Bkg MC","LPE");

    TLegend *leg2 = new TLegend(0.4375,0.7236441,0.85,0.9240678, TString::Format("RooKeysPdf",c), "brNDC");

    TLegend *leg3 = new TLegend(0.2175,0.8236441,0.6575,0.9240678, TString::Format("Category %d",c), "brNDC");
    if(isMirror){
      leg2->AddEntry(plotPhotonsMassBkgMC[c]->getObject(1),"Default Bw","L");
      leg2->AddEntry(plotPhotonsMassBkgMC[c]->getObject(2),"Bw x 2","L");
      leg2->AddEntry(plotPhotonsMassBkgMC[c]->getObject(3),"Bw x 3","L");
      /* leg2->AddEntry(plotPhotonsMassBkgMC[c]->getObject(4),"Bw x 4","L");
      leg2->AddEntry(plotPhotonsMassBkgMC[c]->getObject(4),"Bw x 4","L");
      leg2->AddEntry(plotPhotonsMassBkgMC[c]->getObject(4),"Bw x 4","L");*/
    }else{
      leg2->AddEntry(plotPhotonsMassBkgMC[c]->getObject(1),"Default bw NoMirr","L");
      leg2->AddEntry(plotPhotonsMassBkgMC[c]->getObject(2),"Bw x 2 NoMirr","L");
      leg2->AddEntry(plotPhotonsMassBkgMC[c]->getObject(3),"Bw x 3 NoMirr","L");
      //leg2->AddEntry(plotPhotonsMassBkgMC[c]->getObject(4),"Bw x 4 NoMirr","L");
    }
    
    leg1->SetTextSize(0.035);
    leg1->SetTextFont(42);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
   
    leg2->SetTextSize(0.035);
    leg2->SetTextFont(42);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);

    leg3->SetTextSize(0.035);
    leg3->SetTextFont(42);
    leg3->SetBorderSize(0);
    leg3->SetFillStyle(0);
   
   
    TPaveText* label_cms = get_labelCMS(0, "2012", true);
    TPaveText* label_sqrt = get_labelSqrt(0);

  
    ctmp->cd();
    plotPhotonsMassBkgMC[c]->GetXaxis()->SetTitle("m_{#gamma #gamma} [GeV]");
    plotPhotonsMassBkgMC[c]->Draw();  
    leg1->Draw("same");
    leg2->Draw("same");
    
    label_cms->Draw("same");
    label_sqrt->Draw("same");

   
    ctmp->SetLogy(1);
    if(isMirror){
    ctmp->SaveAs("plots/BKG_MC_rookeyspdf_"+TString::Format("cat%d_LOG.png", c));
    ctmp->SaveAs("plots/BKG_MC_rookeyspdf_"+TString::Format("cat%d_LOG.pdf", c));
    ctmp->SaveAs("plots/BKG_MC_rookeyspdf_"+TString::Format("cat%d_LOG.root", c));
    }else{
      ctmp->SaveAs("plots/BKG_MC_rookeyspdf_NoMirr_"+TString::Format("cat%d_LOG.png", c));
      ctmp->SaveAs("plots/BKG_MC_rookeyspdf_NoMirr_"+TString::Format("cat%d_LOG.pdf", c));
      ctmp->SaveAs("plots/BKG_MC_rookeyspdf_NoMirr_"+TString::Format("cat%d_LOG.root", c));
  
    }
       ctmp->SetLogy(0);
       /*  //plot D1
       plotPhotonsMassBkgMC_D1[c]->GetYaxis()->SetTitle("First Derivative ");
       plotPhotonsMassBkgMC_D1[c]->GetXaxis()->SetTitle("m_{#gamma #gamma} [GeV]");
       plotPhotonsMassBkgMC_D1[c]->Draw();  
       leg2->SetHeader("RooKeysPDF - 1st Derivative");
       leg2->Draw("same");
       leg3->Draw("same");
       
       label_cms->Draw("same");
       label_sqrt->Draw("same");
       if(isMirror){
	 ctmp->SaveAs("plots/BKG_MC_rookeyspdf_D1_"+TString::Format("cat%d.png", c));
	 ctmp->SaveAs("plots/BKG_MC_rookeyspdf_D1_"+TString::Format("cat%d.pdf", c));
	 ctmp->SaveAs("plots/BKG_MC_rookeyspdf_D1_"+TString::Format("cat%d.root", c));
       }else{
	 ctmp->SaveAs("plots/BKG_MC_rookeyspdf_D1_NoMirr_"+TString::Format("cat%d_LOG.png", c));
	 ctmp->SaveAs("plots/BKG_MC_rookeyspdf_D1_NoMirr_"+TString::Format("cat%d_LOG.pdf", c));
	 ctmp->SaveAs("plots/BKG_MC_rookeyspdf_D1_NoMirr_"+TString::Format("cat%d_LOG.root", c));
	 
       }  
       ctmp->SetLogy(0);
       ctmp->Clear();

    //plot D2
    plotPhotonsMassBkgMC_D2[c]->GetYaxis()->SetTitle("Second Derivative");
    plotPhotonsMassBkgMC_D2[c]->GetXaxis()->SetTitle("m_{#gamma #gamma} [GeV]");
    plotPhotonsMassBkgMC_D2[c]->Draw();  
    leg2->SetHeader("RooKeysPDF - 2nd Derivative");
    leg2->Draw("same");
    leg3->Draw("same");
	
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    if(isMirror){
    ctmp->SaveAs("plots/BKG_MC_rookeyspdf_D2_"+TString::Format("cat%d.png", c));
    ctmp->SaveAs("plots/BKG_MC_rookeyspdf_D2_"+TString::Format("cat%d.pdf", c));
    ctmp->SaveAs("plots/BKG_MC_rookeyspdf_D2_"+TString::Format("cat%d.root", c));
    }else{
      ctmp->SaveAs("plots/BKG_MC_rookeyspdf_D2_NoMirr_"+TString::Format("cat%d_LOG.png", c));
      ctmp->SaveAs("plots/BKG_MC_rookeyspdf_D2_NoMirr_"+TString::Format("cat%d_LOG.pdf", c));
      ctmp->SaveAs("plots/BKG_MC_rookeyspdf_D2_NoMirr_"+TString::Format("cat%d_LOG.root", c));
      
      } */
  }

  
    
  std::cout << "done with importing MC background datasets and RooKeysPdfs" << std::endl;
  

  TString filename(wsDir+"HighMass-hgg.RooKeysPdfMCBkg_8TeV_cat2.root");
  wBias->writeToFile(filename);
  cout << "Write background RooKeys workspace in: " << filename << " file" << endl;

  return;

}



RooHistFunc* getRooHistFunc(int cat, RooRealVar* var){


  double mass[5] = {150., 200., 250., 300., 400.};
  double c0[5] = {57.2743,73.98, 86.1266,99.099,115.707};
  double c1[5] = {69.246,76.9358,75.4159,72.8432,69.0506 };
  double c2[5] = {27.6973,33.0972,34.1901,37.6118,36.0121 };
  double c3[5] = {54.2063,60.2946,59.8451,60.5352,53.4554 };
  double all[5] = {187.921, 217.376, 226.456,237.32,241.692};

  
  TH1F* h_all = new TH1F("h_all", "h_all", 15, 130, 875);
  for(int i=0;i<5;i++){
    std::cout<<cat<<std::endl;
    std::cout<<mass[i]<<"  "<<h_all->FindBin(mass[i])<<"   "<<c0[i]/19620.<<std::endl;
    if(cat==0) h_all->SetBinContent(h_all->FindBin(mass[i]),c0[i]/19620.);
    if(cat==1) h_all->SetBinContent(h_all->FindBin(mass[i]),c1[i]/19620.);
    if(cat==2) h_all->SetBinContent(h_all->FindBin(mass[i]),c2[i]/19620.);
    if(cat==3) h_all->SetBinContent(h_all->FindBin(mass[i]),c3[i]/19620.);
    if(cat==4) h_all->SetBinContent(h_all->FindBin(mass[i]),all[i]/19620.);
   
  }
  std::cout<<"--------------------------->"<<cat<<std::endl;

  if(cat==0) h_all->SetBinContent(h_all->FindBin(350),(c0[3]+c0[4])/2/19620.);
  if(cat==1) h_all->SetBinContent(h_all->FindBin(350),(c1[3]+c1[4])/2/19620.);
  if(cat==2) h_all->SetBinContent(h_all->FindBin(350),(c2[3]+c2[4])/2/19620.);
  if(cat==3) h_all->SetBinContent(h_all->FindBin(350),(c3[3]+c3[4])/2/19620.);
  if(cat==4) h_all->SetBinContent(h_all->FindBin(350),(all[3]+all[4])/2/19620.);
 
  double highmass[9]={450., 500., 550., 600., 650., 700., 750., 800., 850.};
  for(int i = 0; i<9;i++){
  if(cat==0) h_all->SetBinContent(h_all->FindBin(highmass[i]),(c0[4])/19620.);
  if(cat==1) h_all->SetBinContent(h_all->FindBin(highmass[i]),(c1[4])/19620.);
  if(cat==2) h_all->SetBinContent(h_all->FindBin(highmass[i]),(c2[4])/19620.);
  if(cat==3) h_all->SetBinContent(h_all->FindBin(highmass[i]),(c3[4])/19620.);
  if(cat==4) h_all->SetBinContent(h_all->FindBin(highmass[i]),(all[4])/19620.);
  }
  
  TCanvas* c = new TCanvas("c", "c", 1);
  c->cd();
  //  TPaveText* label_cms = get_labelCMS(0, false);
  //TPaveText* label_sqrt = get_labelSqrt(0);  
  h_all->GetYaxis()->SetRangeUser(0.0001, 0.1);
  h_all->GetXaxis()->SetTitle("m_{H} [GeV]");
  h_all->GetYaxis()->SetTitle("Signal Yield/1pb");
  h_all->Draw();

  RooDataHist* rooData_all = new RooDataHist("rooData_all", "roData_all",*var,h_all);
  RooHistFunc* rooFunc_all = new RooHistFunc(TString::Format("PhotonsMassSig_cat%d_norm",cat),TString::Format("PhotonsMassSig_cat%d_norm",cat), *var,*rooData_all, 3);
  RooPlot* plot = var->frame();
  rooData_all->plotOn(plot);
  rooFunc_all->plotOn(plot);
  plot->Draw("same");
  c->SetLogy();

  c->SaveAs("plots/signalYield.png");
  c->SaveAs("plots/signalYield.pdf");

  return rooFunc_all;

}






RooHistFunc* getRooHistFuncFitMIN(int cat, RooRealVar* var){


  double mass[8] = {150., 250, 350, 450, 550, 650, 750, 850};
  double c0[8] = {130., 200., 250., 300., 300., 450., 450., 550. };
  double c1[8] = {130., 200., 250., 300., 300., 450., 450., 550. };
  double c2[8] = {130., 200., 250., 300., 300., 450., 450., 550. };
  double c3[8] = {130., 200., 250., 300., 300., 450., 450., 550. };
 
  
  TH1F* h_all = new TH1F("h_all", "h_all", 8, 100, 900);
  for(int i=0;i<8;i++){
    std::cout<<cat<<std::endl;
    std::cout<<mass[i]<<"  "<<h_all->FindBin(mass[i])<<"   "<<c0[i]<<std::endl;
    if(cat==0) h_all->SetBinContent(h_all->FindBin(mass[i]),c0[i]);
    if(cat==1) h_all->SetBinContent(h_all->FindBin(mass[i]),c1[i]);
    if(cat==2) h_all->SetBinContent(h_all->FindBin(mass[i]),c2[i]);
    if(cat==3) h_all->SetBinContent(h_all->FindBin(mass[i]),c3[i]);
      
  }
  TCanvas* c = new TCanvas("c", "c", 1);
  c->cd();
  TPaveText* label_cms = get_labelCMS(0,"2012", false);
  TPaveText* label_sqrt = get_labelSqrt(0);  

  h_all->GetXaxis()->SetTitle("m_{H} [GeV]");
  h_all->GetYaxis()->SetTitle("Fit Range: Minimum");
  h_all->Draw("");
 
  RooDataHist* rooData_all = new RooDataHist("rooData_all", "roData_all",*var,h_all);
  RooHistFunc* rooFunc_all = new RooHistFunc(TString::Format("PhotonsMassBkg_cat%d_min",cat),TString::Format("PhotonsMassBkg_cat%d_min",cat), *var,*rooData_all, 3);
  RooPlot* plot = var->frame();
  rooData_all->plotOn(plot);
  rooFunc_all->plotOn(plot);
  plot->Draw("same");
  // c->SetLogy();

  c->SaveAs("plots/.png");
  c->SaveAs("plots/Bkg_fitMinimum.pdf");

  return rooFunc_all;

}



TF1* getFuncFitMIN(int cat, RooRealVar* var){


  double mass[8] = {150., 250, 350, 450, 550, 650, 750, 850};
  double c0[8] = {130., 200., 250., 300., 300., 450., 450., 550. };
  double c1[8] = {130., 200., 250., 300., 300., 450., 450., 550. };
  double c2[8] = {130., 200., 250., 300., 300., 450., 450., 550. };
  double c3[8] = {130., 200., 250., 300., 300., 450., 450., 550. };
 
  TF1* f = new TF1("f", "[0]+[1]*x", 150, 850);
  f->SetParameter(0,1000 );
  f->SetParameter(1, -0.5);
  f->SetParameter(2, -1000.);
  
  
  TH1F* h_all = new TH1F("h_all", "h_all", 8, 100, 900);
  for(int i=0;i<8;i++){
    if(cat==0) h_all->SetBinContent(h_all->FindBin(mass[i]),c0[i]);
    if(cat==1) h_all->SetBinContent(h_all->FindBin(mass[i]),c1[i]);
    if(cat==2) h_all->SetBinContent(h_all->FindBin(mass[i]),c2[i]);
    if(cat==3) h_all->SetBinContent(h_all->FindBin(mass[i]),c3[i]);
      
  }
  TCanvas* c = new TCanvas("c", "c", 1);
  c->cd();
  TPaveText* label_cms = get_labelCMS(0,"2012", false);
  TPaveText* label_sqrt = get_labelSqrt(0);  

  h_all->GetXaxis()->SetTitle("m_{H} [GeV]");
  h_all->GetYaxis()->SetTitle("Fit Range: Minimum");
  h_all->Draw("");
  h_all->Fit("f");
  f->Draw("Lsame");
  for(int i=0;i<8;i++){
    std::cout<<cat<<std::endl;
    std::cout<<mass[i]<<"  "<<h_all->FindBin(mass[i])<<"   "<<c0[i]<< "    "<<f->Eval(mass[i])<<std::endl;
       
  }
 
  return f;

}






RooHistFunc* getRooHistFuncFitMAX(int cat, RooRealVar* var){

  
  double mass[8] = {150,250, 350, 450, 550, 650, 750, 850};
  double c0[8] = {230, 550, 650, 700, 800, 800, 1000, 1000};
  double c1[8] = {230, 550, 650, 700, 800, 800, 1000, 1000};
  double c2[8] = {230, 550, 650, 700, 800, 800, 1000, 1000};
  double c3[8] = {230, 550, 650, 700, 800, 800, 1000, 1000};


  TH1F* h_all = new TH1F("h_all", "h_all", 8, 100, 900);
  for(int i=0;i<8;i++){
     std::cout<<cat<<std::endl;
     std::cout<<mass[i]<<"  "<<h_all->FindBin(mass[i])<<"   "<<c0[i]<<std::endl;
     if(cat==0) h_all->SetBinContent(h_all->FindBin(mass[i]),c0[i]);
     if(cat==1) h_all->SetBinContent(h_all->FindBin(mass[i]),c1[i]);
     if(cat==2) h_all->SetBinContent(h_all->FindBin(mass[i]),c2[i]);
     if(cat==3) h_all->SetBinContent(h_all->FindBin(mass[i]),c3[i]);
      
  }
  TCanvas* c = new TCanvas("c", "c", 1);
  c->cd();
  TPaveText* label_cms = get_labelCMS(0,"2012", false);
  TPaveText* label_sqrt = get_labelSqrt(0);  

  h_all->GetXaxis()->SetTitle("m_{H} [GeV]");
  h_all->GetYaxis()->SetTitle("Fit Range: Maximum");
  h_all->Draw();
  

  RooDataHist* rooData_all = new RooDataHist("rooData_all", "roData_all",*var,h_all);
  RooHistFunc* rooFunc_all = new RooHistFunc(TString::Format("PhotonsMassBkg_cat%d_max",cat),TString::Format("PhotonsMassBkg_cat%d_max",cat), *var,*rooData_all, 7);

  RooPlot* plot = var->frame();
  rooData_all->plotOn(plot);
  rooFunc_all->plotOn(plot);
  plot->Draw("same");
  c->SetLogy();

  c->SaveAs("plots/.png");
  c->SaveAs("plots/Bkg_fitMaximum.pdf");

  return rooFunc_all;

}



TF1* getFuncFitMAX(int cat, RooRealVar* var){

  
  double mass[8] = {150,250, 350, 450, 550, 650, 750, 850};
  double c0[8] = {230, 550, 650, 700, 800, 800, 1000, 1000};
  double c1[8] = {230, 550, 650, 700, 800, 800, 1000, 1000};
  double c2[8] = {230, 550, 650, 700, 800, 800, 1000, 1000};
  double c3[8] = {230, 550, 650, 700, 800, 800, 1000, 1000};

  TF1* f = new TF1("f", "[0]+[2]*pow(x,[1])", 150, 850);
  f->SetParameter(0,1000 );
  f->SetParameter(1, -0.5);
  f->SetParameter(2, -1000.);
  

  TH1F* h_all = new TH1F("h_all", "h_all", 8, 100, 900);
  for(int i=0;i<8;i++){
    if(cat==0) h_all->SetBinContent(h_all->FindBin(mass[i]),c0[i]);
    if(cat==1) h_all->SetBinContent(h_all->FindBin(mass[i]),c1[i]);
    if(cat==2) h_all->SetBinContent(h_all->FindBin(mass[i]),c2[i]);
    if(cat==3) h_all->SetBinContent(h_all->FindBin(mass[i]),c3[i]);
      
  }
  TCanvas* c = new TCanvas("c", "c", 1);
  c->cd();
  TPaveText* label_cms = get_labelCMS(0,"2012", false);
  TPaveText* label_sqrt = get_labelSqrt(0);  

  h_all->GetXaxis()->SetTitle("m_{H} [GeV]");
  h_all->GetYaxis()->SetTitle("Fit Range: Maximum");
  h_all->Draw();
  h_all->Fit("f");
  f->Draw("Lsame");
  for(int i=0;i<8;i++){
    std::cout<<cat<<std::endl;
    std::cout<<mass[i]<<"  "<<h_all->FindBin(mass[i])<<"   "<<c0[i]<< "    "<<f->Eval(mass[i])<<std::endl;
       
  }
  

 
  return f;

}



makeRooHistPlot(RooHistFunc* rooFitMin,RooHistFunc* rooFitMax ,TF1* fmin, TF1* fmax, RooRealVar* var){

  TCanvas* cfit = new TCanvas("cfit", "cfit", 1);
  cfit->cd();
  
  RooPlot* p = var->frame();
  rooFitMin->plotOn(p, LineColor(kAzure+9));
  rooFitMax->plotOn(p, LineColor(kViolet+9));
  p->Draw();
  fmax->SetLineColor(kViolet+9);
  fmax->SetLineStyle(kDashed);
  fmin->SetLineColor(kAzure+9);
  fmin->SetLineStyle(kDashed);
  fmax->Draw("same");
  fmin->Draw("same");
  p->GetYaxis()->SetTitle("m_{#gamma #gamma} [GeV]");
  p->GetXaxis()->SetTitle("m_{H} [GeV]");

  TLegend* leg = new TLegend(0.2, 0.65, 0.45, 0.89, "", "brNDC");
  leg->SetTextSize(0.0206044);  
  leg->SetTextFont(42);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(p->getObject(1), "Fit Maximum", "L");
  leg->AddEntry(p->getObject(0), "Fit Minimum", "L");
  leg->Draw("same");
  

  TPaveText* label_cms = get_labelCMS(0, "2012", false);
  TPaveText* label_sqrt = get_labelSqrt(0);
  label_cms->Draw("same");
  label_sqrt->Draw("same");
  cfit->SaveAs("plots/bkg_fitRange.png");
  cfit->SaveAs("plots/bkg_fitRange.pdf");
  cfit->SaveAs("~/www/BkgFit/bkg_fitRange.pdf");
  cfit->SaveAs("~/www/BkgFit/bkg_fitRange.png");
 



}
