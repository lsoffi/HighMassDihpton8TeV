
using namespace RooFit;
using namespace RooStats;
//#include "RooCruijff.h"
//#include "HZZ4LRooPdfs.h"
//#include "RooCPSHighMassGGHNoInterf.h"

static const Int_t NCAT = 1;  // chiara
Int_t MINmass;
Int_t MAXmass;
Float_t MINmassT;
Float_t MAXmassT;
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



TPaveText* get_labelCMS( int legendQuadrant = 0 , bool sim=false) {

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
    x1 = 0.105;
    y1 = 0.953;//953
    x2 = 0.6;
    y2 = 0.975;//975

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
    y1 = 0.953;//953
    x2 = 0.94;
    y2 = 0.975;//975
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
  // RooRealVar* PhotonsMassTrue  = new RooRealVar("PhotonsMassTrue", "M(gg)",MINmassT,MAXmassT,"GeV");
  RooRealVar* ph1_eta = new RooRealVar("ph1_eta", "eta(g1)",-10,10,"");
  RooRealVar* ph2_eta = new RooRealVar("ph2_eta", "eta(g2)",-10,10,"");
  RooRealVar* ph1_r9 = new RooRealVar("ph1_r9", "R9(g1)",-10,10,"");
  RooRealVar* ph2_r9 = new RooRealVar("ph2_r9", "R9(g2)",-10,10,"");
  RooRealVar* evweight = new RooRealVar("evweight","weightings",0,1000,"");
  RooRealVar* nvtx = new RooRealVar("nvtx", "nvtx", 0, 50, "");
  RooArgSet* ntplVars = new RooArgSet(*PhotonsMass,/**PhotonsMassTrue*,*/*ph1_eta, *ph2_eta, *ph1_r9, *ph2_r9, *evweight);
  
  return ntplVars;
}

RooArgSet* defineVariablesM250() {

  // define variables of the input ntuple //livia
  RooRealVar* PhotonsMass  = new RooRealVar("PhotonsMass", "M(gg)",130, 1000,"GeV");
  RooRealVar* PhotonsMassTrue  = new RooRealVar("PhotonsMassTrue", "M(gg)",MINmassT,MAXmassT,"GeV");
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



void runfits(const Float_t mass=150, Bool_t dobands = false, Float_t width=0.1, std::string model) {

  //******************************************************************//
  //  Running mode  corresponds to the following cases
  //         - full run set:
  //         - create signal and background data sets 
  //         - make and fit signal and background  models 
  //         - write signal and background workspaces in root files
  //         - write data card
  //*******************************************************************//

 
  TString fileBaName("HighMass-hgg");    
  TString fileBkgName("HighMass-hgg.inputbkg");
  
  TString card_name("HighMass-hgg_models_Bkg_8TeV_test.rs");
  HLFactory hlf("HLFactory", card_name, false);
  RooWorkspace* w = hlf.GetWs();
 
  RooFitResult* fitresults;


  // Luminosity:
  Float_t Lum = 19500.0;  
  RooRealVar lumi("lumi","lumi",Lum);
  w->import(lumi); 
  
   //compute roohistfunc with min and max for the fir range 
  RooRealVar* var = new RooRealVar("var","var", 100, 900 );

  RooHistFunc* rooFitMin = getRooHistFuncFitMIN(0, var);
  RooHistFunc* rooFitMax = getRooHistFuncFitMAX(0, var);
  TF1* fFitMin = getFuncFitMIN(0, var);
  TF1* fFitMax = getFuncFitMAX(0, var);
  
  
 

  //  w->Print("v");
  
  cout << endl; cout << "Now AddSigData" << endl;

  double widths[4] = {0.1,2.,5,10};//, 2., 5., 10.
  double masses[7] = {150., 200., 250., 300., 400., 600., 800} ; //

 
 
  double p0[NCAT][4][7];
  double p0e[NCAT][4][7];
  double p1[NCAT][4][7];
  double p1e[NCAT][4][7];
  double p2[NCAT][4][7];
  double p2e[NCAT][4][7];

  TFile* file = new TFile("BW_corrections.root","RECREATE");
  file->cd();
    for(int v=0; v<1;v++){

      std::cout<<"  WIDTH: ---------------------------------------------------------> "<<widths[v]<<std::endl;
      for (int m=0; m<1;m++){
	std::cout<<"  MASS: ---------------------------------------------------------> "<<masses[m]<<std::endl;
	/*		
	float min;
	float max;
	float wind ;
	if(width<1)wind=0.05;
	else wind = masses[m]*widths[v]/100;
	min =masses[m]-4.5*wind;
	max = masses[m]+4.5*wind;
	if(min<130)min=130.;
	if(max>1000)max=1000.;
	MINmassT=min;
	MAXmassT=max;
	MINmass=min;
	MAXmass=max;*/	
	std::cout<<fFitMin->Eval(masses[m])<<"    "<<fFitMax->Eval(masses[m])<<std::endl;
	// makeRooHistPlot(rooFitMin,rooFitMax,fFitMin,fFitMax,var);
  
	var->setVal(masses[m]);
	double newmin = fFitMin->Eval(masses[m]);
	double newmax = fFitMax->Eval(masses[m]);
	
	if(masses[m]==150){
	  newmin=130;
	  newmax=230;
	}
	
	MINmassT=newmin;
	MAXmassT=newmax;
	MINmass=newmin;
	MAXmass=newmax;
	
	AddSigData(w, masses[m],widths[v], model);
	std::cout<<MINmassT<<"  "<<MAXmassT<<std::endl;

	for(int c=0; c<NCAT; c++){
	  std::cout<<"  CAT: ---------------------------------------------------------> "<<c<<std::endl;
	  double res[6];
	  //	  SigModelTheoFcnFitCPS(w);  
	  SigModelTheoFcnFit(w,c,masses[m], widths[v],model,res);  

	p0[c][v][m]= res[0] ;
	p0e[c][v][m]= res[1] ;
	p1[c][v][m]= res[2] ;
	p1e[c][v][m]= res[3] ;
	p2[c][v][m]= res[4] ;
	p2e[c][v][m]= res[5] ;
	 
      }

    }
  }

  
    /*
    for(int c=0; c<NCAT; c++){
      TH2D* hp0 = new TH2D(TString::Format("hp0_cat%d",c),TString::Format("hp0_cat%d",c),6, 130., 850., 6, 0.,10);
      TH2D* hp1 = new TH2D(TString::Format("hp1_cat%d",c),TString::Format("hp1_cat%d",c),6, 130., 850., 6, 0.,10);
      TH2D* hp2 = new TH2D(TString::Format("hp2_cat%d",c),TString::Format("hp2_cat%d",c),6, 130., 850., 6, 0.,10 );
       
      int i = 0;
      for(int v=0; v<4;v++){
	std::endl;std::cout<<" --------- "<<std::endl;
	std::endl;std::cout<<" ---- p0 ----- "<<std::endl;
	for (int m=0; m<7;m++) std::cout<<p0[c][v][m]<<"   ";
	std::cout<<" "<<std::endl;
	std::cout<<" ---- p1 ----- "<<std::endl;
	for (int m=0; m<7;m++) std::cout<<p1[c][v][m]<<"   ";
	std::cout<<" "<<std::endl;
	std::cout<<" ---- p2 ----- "<<std::endl;
	for (int m=0; m<7;m++) std::cout<<p2[c][v][m]<<"   ";
	std::cout<<" "<<std::endl;

      }
      //m200-400-800-150->1-4-6-0
      double vec_p0[16] = {p0[c][0][1],p0[c][1][1],p0[c][2][1],p0[c][3][1], p0[c][0][4],p0[c][1][4],p0[c][2][4],p0[c][3][4],p0[c][0][6],p0[c][1][6],p0[c][2][6],p0[c][3][6],p0[c][0][0],p0[c][1][0],p0[c][2][0],p0[c][3][0] }; 
      double vec_p1[16] = {p1[c][0][1],p1[c][1][1],p1[c][2][1],p1[c][3][1], p1[c][0][4],p1[c][1][4],p1[c][2][4],p1[c][3][4],p1[c][0][6],p1[c][1][6],p1[c][2][6],p1[c][3][6],p1[c][0][0],p1[c][1][0],p1[c][2][0],p1[c][3][0] };  


      //   TGraph2D* g0 = makeTG2D(0, c,vec_p0);
      //TGraph2D* g1 = makeTG2D(1, c,vec_p1);
      //g0->Write();
      //g1->Write();
    }
    file->Write();
    file->Close();*/
  return;
}


TGraph2D* makeTG2D(int p, int c, double* vec){

 
  TH2D* hp0 = new TH2D(TString::Format("hp%d_cat%d",p,c),TString::Format("hp%d_cat%d",p,c),7, 140., 860., 7, -1.,10.5);
   
      TGraph2D* g0 = new TGraph2D(hp0);

      /*  g0->SetPoint(12, 150.,0,1	   );
      g0->SetPoint(13, 150.,2,-0.0716112 );
      g0->SetPoint(14, 150.,5,-0.0196842 );
      g0->SetPoint(15, 150.,10,-0.0124951 );
      g0->SetPoint(0, 200., 0,1	   );
      g0->SetPoint(1, 200., 2,-0.0129597 );
      g0->SetPoint(2, 200., 5,-0.0092817  );
      g0->SetPoint(3, 200., 10,-0.0077557);
      g0->SetPoint(4, 400., 0,1	   );
      g0->SetPoint(5, 400., 2,-0.0297586  );
      g0->SetPoint(6, 400., 5,-0.0061914  );
      g0->SetPoint(7, 400., 10,-0.0047460  );
      g0->SetPoint(8, 850., 0,1	   );
      g0->SetPoint(9, 850., 2,0.0170733   );
      g0->SetPoint(10, 850.,5,-0.041108  );
      g0->SetPoint(11, 850.,10,-0.0286136 );
   */




      //    for(int i = 0; i<16; i++) std::cout<<vec[i]<<std::endl;
      g0->SetPoint(12, 150., 0, vec[12]);
      g0->SetPoint(13, 150., 2, vec[13] );
      g0->SetPoint(14, 150., 5,vec[14] );
      g0->SetPoint(15, 150., 10.1, vec[15] );
      g0->SetPoint(0, 200., 0,vec[0]);
      g0->SetPoint(1, 200., 2, vec[1]);
      g0->SetPoint(2, 200., 5,  vec[2]);
      g0->SetPoint(3, 200., 10.1,vec[3]);
      g0->SetPoint(4, 400., 0, vec[4]);
      g0->SetPoint(5, 400., 2,vec[5]);
      g0->SetPoint(6, 400., 5,vec[6]);
      g0->SetPoint(7, 400., 10.1, vec[7]);
      g0->SetPoint(8, 850.1., 0, vec[8]);
      g0->SetPoint(9, 850.1, 2, vec[9]);
      g0->SetPoint(10, 850.1, 5,vec[10]);
      g0->SetPoint(11, 850.1, 10.1, vec[11]);
    
      g0->SetTitle(TString::Format("p%d_cat%d",p,c));
      g0->SetName(TString::Format("p%d_cat%d",p,c));
      TCanvas* ca = new TCanvas("c", "c",600);
      ca->cd();
      TPad* pad = new TPad("pad","pad", 0.01, 0.03, 0.9, 0.9.);
      pad->Draw();
      pad->cd();
      TPaveText* label_cms = get_labelCMS(0, false);
      TPaveText* label_sqrt = get_labelSqrt(0);  
      
      hp0->SetStats(0);
      //  hp0->SetTitle("");
     
      hp0->GetXaxis()->SetTitle("m_{H} [GeV]");
      hp0->GetYaxis()->SetTitle("#Gamma_{H} [GeV]");
      hp0->Draw("axis");
      g0->Draw("COLZsame");
      label_cms->Draw("same");
      label_sqrt->Draw("same");
     
      g0->GetZaxis()->SetTitleSize(0.04);
      g0->GetZaxis()->SetTitleOffset(1.6);
      g0->GetZaxis()->SetLabelOffset(-0.1);
      std::cout<<g0->Interpolate(300, 5)<<std::endl;
      pad->RedrawAxis();
      ca->RedrawAxis();
      ca->SaveAs(TString::Format("/afs/cern.ch/user/s/soffi/www/plotsNota/p%d_cat%d.png",p,c));
      ca->SaveAs(TString::Format("/afs/cern.ch/user/s/soffi/www/plotsNota/p%d_cat%d.pdf",p,c));

      return g0;

}

// Signal Data Set
void AddSigData(RooWorkspace* w,  Float_t mass, Float_t width, std::string model) {

  Int_t ncat = NCAT;
  TString inDir = "";

  Float_t MASS(mass);

  
  // Variables
  RooArgSet* ntplVars = defineVariables();
  ntplVars->add(*w->var("PhotonsMassTrue"));
  (*w->var("PhotonsMassTrue"))->setRange(MINmassT, MAXmassT);
  // (*w->var("PhotonsMassTrue"))->setBinning();
  (*w->var("PhotonsMassTrue")).Print("V");
  int iMass(mass);      
  int iWidth(width);      

 
  //chain summing up all production modes
  TChain* sigTree1  = new TChain();

  
   if(width<1)  sigTree1->Add(TString::Format("histograms_CMS-HGG_08052014_MC.root/ggh_m%d_8TeV",iMass));
   else  sigTree1->Add(TString::Format("histograms_CMS-HGG_06062014_MC_w%d.root/ggh_m%d_8TeV",  iWidth, iMass));
   sigTree1->SetTitle("sigTree1");
   sigTree1->SetName("sigTree1");


  // common preselection cut
  //  TString mainCut = TString::Format("PhotonsMassTrue>=%f && PhotonsMassTrue<=%f",MINmassT, MAXmassT);   // livia 
  TString mainCut = TString::Format("PhotonsMassTrue>%f && PhotonsMassTrue<%f ",MINmassT, MAXmassT );   // livia
  
  
  // Create signal dataset composed with different productions, the weight is already applied in our ntuples
  RooDataSet sigWeighted("sigWeighted","dataset",sigTree1,*ntplVars,mainCut,"evweight");
  cout << endl;
  cout << "sigWeighted" << endl;
  sigWeighted.Print("v");
  cout << "---- nX:  " << sigWeighted.sumEntries() << endl; 
  
  // apply a common preselection cut; split in categories
  cout << endl;
  RooDataSet* signal;
 
  for(int c = 0; c <NCAT; c++){
    // 1)  prime 4 cat livia
    if (c==0) signal = (RooDataSet*) sigWeighted.reduce(*w->var("PhotonsMassTrue"),mainCut+TString::Format("&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 )"));
    if (c==1) signal = (RooDataSet*) sigWeighted.reduce(*w->var("PhotonsMassTrue"),mainCut+TString::Format("&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9<0.94 || ph2_r9<0.94 ) "));
    if (c==2) signal = (RooDataSet*) sigWeighted.reduce(*w->var("PhotonsMassTrue"),mainCut+TString::Format("&& (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 && ph2_r9>0.94 )"));
    if (c==3) signal = (RooDataSet*) sigWeighted.reduce(*w->var("PhotonsMassTrue"),mainCut+TString::Format("&& (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9<0.94 || ph2_r9<0.94 ) "));
    int imass(mass);
    int iwidth(width);
    w->import(*signal,Rename(TString::Format("SigWeight",c, imass, iwidth)));
    
    cout << "cat " << c << ", signal: " << endl;
    signal->Print("v");
    cout << "---- for category " << c << ", nX for signal[c]:  " << signal->sumEntries() << endl; 
    cout << endl;
  }

  // Create full weighted signal data set without categorization
  RooDataSet* signalAll = (RooDataSet*) sigWeighted.reduce(*w->var("PhotonsMassTrue"),mainCut);
  w->import(*signalAll, Rename("SigWeight_cat4"));
  cout << "now signalAll" << endl;
  signalAll->Print("v");
  cout << "---- nX for signalAll:  " << signalAll->sumEntries() << endl; 
  cout << endl;


  // w->Print("V");
}




void SigModelTheoFcnFit(RooWorkspace* w, int c, double mass, double width,std::string model, double* p) {


  
    
  RooDataSet* signal;
  RooRealVar* PhotonsMassTrue = w->var("PhotonsMassTrue"); 
  
  TCanvas c1("c1", "c1", 1);
 

 
  float min;
  float max;
  
  min =MINmassT;
  max = MAXmassT;
  
  int nbin=60;
   
    int imass(mass);
    int iwidth(width); 
    std::cout<<"---------------------> "<<TString::Format("SigWeight",c, imass, iwidth)<<std::endl;
    RooDataSet* signal;
    signal= (RooDataSet*) w->data(TString::Format("SigWeight",c, imass, iwidth));

    
   //BW
  
    RooRealVar sigmaBW_var(TString::Format("sigmaBW_var_cat%d_m%d_w%d",c, imass, iwidth), TString::Format("sigmaBW_var_cat%dm%d_w%d",c, imass, iwidth), 0., 0.2*mass);
    RooRealVar meanBW_var(TString::Format("meanBW_var_cat%d_m%d_w%d",c, imass, iwidth), TString::Format("meanBW_var_cat%d_m%d_w%d",c, imass, iwidth),0.7*mass, 1.3*mass);
    std::cout<<" width:--------> "<<width<<std::endl;
    if(width==10.)sigmaBW_var.setVal(0.1*mass);
    if(width==5.)sigmaBW_var.setVal(0.05*mass);
    if(width==2.)sigmaBW_var.setVal(0.02*mass);
    if(width<1.)sigmaBW_var.setVal(0.1);
    meanBW_var.setVal(mass);
    sigmaBW_var.setConstant();
    meanBW_var.setConstant();
    w->import(sigmaBW_var);
    w->import(meanBW_var);
    
    RooFormulaVar meanBW(TString::Format("meanBW_cat%d_m%d_w%d",c, imass, iwidth),"","@0",*w->var(TString::Format("meanBW_var_cat%d_m%d_w%d",c, imass, iwidth)));    
    RooFormulaVar sigmaBW (TString::Format("sigmaBW_cat%d_m%d_w%d",c, imass, iwidth),"","@0",*w->var(TString::Format("sigmaBW_var_cat%d_m%d_w%d",c, imass, iwidth))); 
    RooBreitWigner BW(TString::Format("BW_cat%d_m%d_w%d",c,imass, iwidth),TString::Format("BW_cat%d_m%d_w%d",c,imass, iwidth), *PhotonsMassTrue, meanBW, sigmaBW);
   
    
    
    //evaluate disagreement
    TH1F* h_datar = new TH1F(TString::Format("h_datar_cat%d_m%d_w%d",c, imass, iwidth),TString::Format("h_datar%d_m%d_w%d",c, imass, iwidth), nbin,MINmassT, MAXmassT);
    TH1F* h_pdfr = new TH1F(TString::Format("h_pdfr%d_m%d_w%d",c, imass, iwidth),TString::Format("h_pdf%d_m%d_w%d",c, imass, iwidth), nbin,MINmassT, MAXmassT);
    h_datar = (TH1F*) signal->createHistogram(TString::Format("h_datar_cat%d_m%d_w%d",c, imass, iwidth),*PhotonsMassTrue);
    h_pdfr =  (TH1F*) BW->createHistogram(TString::Format("h_pdfr_cat%d_m%d_w%d",c, imass, iwidth),*PhotonsMassTrue);
    h_pdfr->Sumw2();
    h_pdfr ->Scale(h_datar->Integral()/h_pdfr->Integral());
    h_datar->Sumw2();
    h_pdfr->Sumw2();
     
    std::cout<<"---------------------> "<<h_datar->GetNbinsX()<<std::endl;
    for(int j = 1; j<h_datar->GetNbinsX(); j++ ){
      double datacont= h_datar->GetBinContent(j);
      double dataerr= h_datar->GetBinError(j);
      h_datar->SetBinError(j, dataerr/h_pdfr->GetBinContent(j));
      h_datar->SetBinContent(j, datacont/h_pdfr->GetBinContent(j));     
    }

    h_datar->GetXaxis()->SetRangeUser(min, max);

    double f_min =  mass*(1-4*width/100);
    if(f_min<MINmassT) f_min = MINmassT;
    double f_max =  mass*(1-4*width/100);
    if(f_max>MAXmassT) f_max = MAXmassT;


    TF1* f = new TF1(TString::Format("f_cat%d",c), "exp([0]*(x-[1]))", mass*(1-4*width/100),mass*(1+4*width/100));
   //   f->SetParLimits(1, 0., 100);
    f->SetParLimits(0, -4.,0);
    f->FixParameter(1, mass);
   
     h_datar->Fit(f, "BMR");
     RooDataHist hist_datar(TString::Format("hist_datar_cat%d",c),TString::Format("hist_datar_cat%d",c), *PhotonsMassTrue, h_datar, 1.);
    c1.cd();
    c1.Clear();
    TPad pad1("pad1", "pad1",0.01256281,0.23945,0.5741206,1.);  
    pad1.Draw();    
    pad1.SetLogy();  
    pad1.cd();   
    
    RooPlot* plotG = (RooPlot*)PhotonsMassTrue->frame( Range(148, 152),Bins(nbin));
    signal->plotOn(plotG, "L");
    BW.plotOn(plotG, LineColor(kGreen), LineStyle(kDashed));
 
    if(width<1){
      f->SetParameter(0,0);
    }

    
    p[0] = f->GetParameter(0);
    p[1] = f->GetParError(0);
    p[2] = f->GetParameter(0);
    p[3] = f->GetParError(0);
    p[4] = f->GetParameter(0);
    p[5] = f->GetParError(0);


   

    RooRealVar v1(TString::Format("v1_cat%d_m%d_w%d",c,imass,iwidth),TString::Format("v1_cat%d_m%d_w%d",c,imass,iwidth),f->GetParameter(0));    
    w->import(v1);
   RooFormulaVar f1(TString::Format("f1%d_m%d_w%d",c,imass,iwidth), TString::Format("f1%d_m%d_w%d",c,imass,iwidth), "@0",*w->var(TString::Format("v1_cat%d_m%d_w%d",c,imass,iwidth)));
 

    RooRealVar mh("mh", "mh",mass);
    mh.setConstant(),
    RooGenericPdf p0f("p0f",  "exp(@1*(@0-@2))", RooArgList(*PhotonsMassTrue, f1,mh));
   
       
    if(width>1) RooProdPdf prod(TString::Format("prod_cat%d_m%d_w%d",c,imass,iwidth),TString::Format("prod_cat%d_m%d_w%d",c,imass,iwidth), RooArgList(p0f, BW));

    std::cout<<"flag"<<std::endl;
    if(width>1)prod.plotOn(plotG, LineColor(kBlue));

    plotG->GetYaxis()->SetRangeUser(0.01,plotG->GetMaximum()*10 );
    plotG->GetXaxis()->SetTitle("m_{#gamma #gamma}^{true} [GeV]");
    plotG->GetXaxis()->SetTitleFont(42);
    plotG->GetXaxis()->SetTitleSize(0.05);
    TPaveText* label_cms = get_labelCMS(0,  true);
    TPaveText* label_sqrt = get_labelSqrt(0);


    TLatex *lat  = new TLatex(0.6,0.8,TString::Format("Cat: %d", c));  
    lat->SetTextSize(0.038);
    lat->SetTextAlign(11);
    lat->SetTextFont(42); 
    lat->SetNDC();
    TLegend* legmc = new TLegend(0.6, 0.55, 0.9, 0.8, "", "brNDC");
    legmc->SetTextSize(0.0206044);  
    legmc->SetTextFont(42);
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->SetHeader(("Model: "+model).c_str());
    int imass(mass)
    legmc->AddEntry(plotG->getObject(0),TString::Format("m_{#gamma#gamma} = %d GeV",imass),"LPE");    
    legmc->AddEntry(plotG->getObject(1),"BW ","L");
    if(width>1)legmc->AddEntry(plotG->getObject(2),"BW #times Exp","L");   
  
    
    if(width>-1)plotG->Draw();    
    plotG->SetTitle("");
    
    lat->Draw("same");
    label_cms->Draw("same");
    label_sqrt->Draw("same"); 
    legmc->Draw("same");
    plotG->Draw("same");
    c1.cd();
    
    if(width>1.){
      // c1.cd();
      //-------pad 3------//DIFF 0.01256281,0.264945,0.5741206,1
      TPad pad3("pad3", "pad3",0.01256281,0.02510989,0.5741206,0.22989);
      pad3.SetGrid();
      pad3.Draw();
      
      pad3.cd();
 

      h_datar->GetYaxis()->SetRangeUser(-2.,h_datar->GetMaximum()*1.5 );
      h_datar->GetYaxis()->SetNdivisions(505);
      h_datar->SetMarkerSize(0.4);
      //   h_datar->GetXaxis()->SetTitle("m_{#gamma #gamma}^{true}");
      h_datar->GetXaxis()->SetTitleSize(0.15);
      h_datar->GetXaxis()->SetLabelSize(0.1);
      h_datar->GetYaxis()->SetLabelSize(0.1);
      h_datar->GetYaxis()->SetTitleSize(0.13);
      h_datar->GetYaxis()->SetTitle("DATA/FIT");
      h_datar->GetXaxis()->SetTitle("");
r      h_datar->GetYaxis()->SetTitleOffset(0.45);
      h_datar->GetXaxis()->SetTitleOffset(0.8);
      
      TH1F* h_dataCopyr =(TH1F*) h_datar->Clone();
      for (int i = 0; i< h_dataCopyr->GetNbinsX()+1;i++) if (h_dataCopyr->GetBinContent(i)==0) h_dataCopyr->SetBinContent(i, 1);
      h_datar->SetMarkerSize(0);
      h_datar->Draw("P");
      h_dataCopyr->Divide(h_dataCopyr);
      h_dataCopyr->SetLineColor(kAzure-2);
      h_dataCopyr->SetFillStyle(3002);
      h_dataCopyr->SetMarkerSize(0.);
      std::cout<<"flag1"<<std::endl;   
      
      h_dataCopyr->Draw("HISTSAME");
      h_datar->SetStats(kFALSE);
      h_datar->Draw("PSAME");
    
      std::cout<<"flag2"<<std::endl;  

      c1.cd();
    }
    pad1->Draw();
    
    c1.SaveAs(TString::Format("plots/Prod_BW_cat%d_LOG_M%d_w%d.png",c,imass, iwidth)); 
    c1.SaveAs(TString::Format("~/www/plotsNota/Prod_BW_cat%d_LOG_M%d_w%d.png",c,imass, iwidth)); 
    c1.SaveAs(TString::Format("~/www/plotsNota/Prod_BW_cat%d_LOG_M%d_w%d.pdf",c,imass, iwidth)); 
    c1.SaveAs(TString::Format("plots/Prod_BW_cat%d_LOG_M%d_w%d.pdf",c,imass, iwidth)); 
    std::cout<<"flag3"<<std::endl;   
    
    std::cout<<"flag4"<<std::endl;
    
}





void SigModelBWFit(double mass) {

  
  
  TCanvas c1("c1", "c1", 1);
 
  float min;
  float max;
  
  min =0;
  max = 1000;
  
  int nbin=100;
  int c = 0;

  double width=0.1*mass;
    int imass(mass);
    int iwidth(width); 
    std::cout<<"---------------------> "<<TString::Format("SigWeight_cat%d_m%d_w%d",c, imass, iwidth)<<std::endl;
      if(mass==200){
    minF=190;
    maxF=210;
    min=180;
    max=220;
    }else if(mass==400){
      minF=330;
      maxF=470;
      min=250;
      max=550;
    }else if(mass==600){
      minF=200;
      maxF=1000;
      min=200;
      max=1000;
    }

      RooDataSet* signal;
      RooRealVar* PhotonsMassTrue = new RooRealVar("PhotonsMassTrue", "PhotonsMassTrue", min, max); 
  
    //evaluate disagreement
      //  TFile* f = new TFile("ggH_HToGG_M_400_TuneZ2star_8TeV_pythia6_tauola_cff_py_GEN_RECO_PU_V1_soffi-ggH_HToGG_M_400_TuneZ2star_8TeV_pythia6_tauola_cff_py_GEN_RECO_PU_V1-3c38936f8952bfd1e4113f6445a6ad95_2_1_pGW.root", "READ");
    //  TFile* f = new TFile("ggH_HToGG_M_200_TuneZ2star_8TeV_pythia6_tauola_cff_py_GEN_RECO_PU_V1_soffi-ggH_HToGG_M_200_TuneZ2star_8TeV_pythia6_tauola_cff_py_GEN_RECO_PU_V1-b15f6399671ac21b9c6c7c8e6161203a_10_1_oRg.root", "READ");
      //TFile* f = new TFile("root:://eoscms//eos/cms/store/group/phys_higgs/soffi/MC/HIG-Summer12PLHE-00021_genvalid_M400.root", "READ");
      TFile* f = TFile::Open("root:://eoscms//eos/cms/store/group/phys_higgs/soffi/MC/HIG-Summer12PLHE-00021_genvalid_M400.root");
     // TTree* tree = f->Get("event");   
     TTree* tree = f->Get("Events");   
    TH1F* h = new TH1F("h", "h", 100, min, max);
    // tree->Draw("gp_p4->M()>>h", "gp_pdgid==25");
    tree->Draw("recoGenParticles_genParticles__VALIDATION.obj.mass()>>h", "recoGenParticles_genParticles__VALIDATION.obj.pdgId()==25");
   
    RooDataHist hist("hist", "hist", *PhotonsMassTrue, h, 1.);    
   //BW
  
    RooRealVar sigmaBW_var(TString::Format("sigmaBW_var_cat%d_m%d_w%d",c, imass, iwidth), TString::Format("sigmaBW_var_cat%dm%d_w%d",c, imass, iwidth), 0., mass);//0., 0.5*mass
    RooRealVar meanBW_var(TString::Format("meanBW_var_cat%d_m%d_w%d",c, imass, iwidth), TString::Format("meanBW_var_cat%d_m%d_w%d",c, imass, iwidth),0., 1000.);//
    std::cout<<" width:--------> "<<width<<std::endl;
    
    RooFormulaVar meanBW(TString::Format("meanBW_cat%d_m%d_w%d",c, imass, iwidth),"","@0",meanBW_var);    
    RooFormulaVar sigmaBW (TString::Format("sigmaBW_cat%d_m%d_w%d",c, imass, iwidth),"","@0",sigmaBW_var); 
    RooBreitWigner BW(TString::Format("BW_cat%d_m%d_w%d",c,imass, iwidth),TString::Format("BW_cat%d_m%d_w%d",c,imass, iwidth), *PhotonsMassTrue, meanBW, sigmaBW); 
  
    BW.fitTo(hist);

    TH1F* h_pdf = new TH1F("h_pdf", "h_pdf", 100, min, max);
    h_pdf->Sumw2();
 
    h->Sumw2();
    h_pdf->Sumw2();
    h_pdf = (TH1F*) BW.createHistogram("PhotonsMassTrue", 100);
    h_pdf ->Scale(h->Integral()/h_pdf->Integral());
    for(int j = 1; j<h->GetNbinsX(); j++ ){
      double datacont= h->GetBinContent(j);
      double dataerr= h->GetBinError(j);
    
      h->SetBinError(j, dataerr/h_pdf->GetBinContent(j));
      h->SetBinContent(j, datacont/h_pdf->GetBinContent(j));     
    }
 
    h->GetXaxis()->SetRangeUser(min, max);
    //   TF1* func = new TF1(TString::Format("f_cat%d",c), "exp([0]*(x-[1]))",minF, maxF);
    TF1* func = new TF1(TString::Format("f_cat%d",c), "[0]+[1]*x",minF, maxF); //+[2]*x*x+[3]*x*x*x
    // func->SetParLimits(0, -4.,0);
    // func->FixParameter(1, mass);   
    // h->Fit(func, "BMR");
    c1.cd();

    for(int i = 0; i<4;i++)  std::cout<<func->GetParameter(i)<<std::endl;
    
    RooRealVar v1(TString::Format("v1_cat%d_m%d_w%d",c,imass,iwidth),TString::Format("v1_cat%d_m%d_w%d",c,imass,iwidth),func->GetParameter(0));    
    RooRealVar v2(TString::Format("v2_cat%d_m%d_w%d",c,imass,iwidth),TString::Format("v2_cat%d_m%d_w%d",c,imass,iwidth),func->GetParameter(1));    
    RooRealVar v3(TString::Format("v3_cat%d_m%d_w%d",c,imass,iwidth),TString::Format("v3_cat%d_m%d_w%d",c,imass,iwidth),func->GetParameter(2));    
    RooRealVar v4(TString::Format("v4_cat%d_m%d_w%d",c,imass,iwidth),TString::Format("v4_cat%d_m%d_w%d",c,imass,iwidth),func->GetParameter(3));    
    RooFormulaVar f1(TString::Format("f1%d_m%d_w%d",c,imass,iwidth), TString::Format("f1%d_m%d_w%d",c,imass,iwidth), "@0",v1);    
    RooFormulaVar f2(TString::Format("f2%d_m%d_w%d",c,imass,iwidth), TString::Format("f2%d_m%d_w%d",c,imass,iwidth), "@0",v2);    
    RooFormulaVar f3(TString::Format("f3%d_m%d_w%d",c,imass,iwidth), TString::Format("f3%d_m%d_w%d",c,imass,iwidth), "@0",v3);     
    RooFormulaVar f4(TString::Format("f4%d_m%d_w%d",c,imass,iwidth), TString::Format("f4%d_m%d_w%d",c,imass,iwidth), "@0",v4);    
    RooRealVar mh("mh", "mh",mass);
  
    //   mh.setConstant(),
    // RooGenericPdf p0f("p0f",  "exp(@1*(@0-@2))", RooArgList(*PhotonsMassTrue, f1,mh));
    RooGenericPdf p0f("p0f",  "@@1+@0*@2", RooArgList(*PhotonsMassTrue, f1,f2));//+@3*@0*@0+@4*@0*@0*@0
    // RooProdPdf prod(TString::Format("prod_cat%d_m%d_w%d",c,imass,iwidth),TString::Format("prod_cat%d_m%d_w%d",c,imass,iwidth), RooArgList( BW, p0f));
    
    c1.Clear();
    TPad pad1("pad1", "pad1",0.01256281,0.23945,0.91206,1.);  
    pad1.Draw();    
    pad1.SetLogy();  
    pad1.cd();
    RooPlot* plotG = (RooPlot*)PhotonsMassTrue->frame( Range(min, max),Bins(100)); 
    hist.plotOn(plotG);   
    BW.plotOn(plotG, LineColor(kGreen), LineStyle(kDashed));
    // prod.plotOn(plotG, LineColor(kBlue));
    plotG->GetYaxis()->SetRangeUser(0.05,plotG->GetMaximum()*10. );
    plotG->GetXaxis()->SetTitle("m_{#gamma #gamma}^{true} [GeV]");
    plotG->GetXaxis()->SetTitleFont(42);
    plotG->GetXaxis()->SetTitleSize(0.05);
    TPaveText* label_cms = get_labelCMS(0,  true);
    TPaveText* label_sqrt = get_labelSqrt(0);
    plotG->Draw();
    
    TLegend* legmc = new TLegend(0.6, 0.6, 0.9, 0.8, "", "brNDC");
    legmc->SetTextSize(0.0406044);  
    legmc->SetTextFont(42);
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    std::string model("GGH");
    legmc->SetHeader(("Model: "+model).c_str());
    int imass(mass)
    legmc->AddEntry(plotG->getObject(0),TString::Format("m_{#gamma#gamma} = %d GeV",imass),"L");    
    legmc->AddEntry(plotG->getObject(1)," BW ","L");    
    
    label_cms->Draw("same");
    label_sqrt->Draw("same"); 
    legmc->Draw("same");
    plotG->Draw("same");
    c1.cd();
    
    
    TPad pad3("pad3", "pad3",0.01256281,0.012510989,0.91206,0.22989);
    pad3.SetGrid();
    pad3.Draw();
    
    pad3.cd();
    
    
    h->GetYaxis()->SetRangeUser(-1.,3 );
    h->GetYaxis()->SetNdivisions(505);
    h->SetMarkerSize(0.4);
    h->GetXaxis()->SetTitleSize(0.15);
    h->GetXaxis()->SetLabelSize(0.1);
    h->GetYaxis()->SetLabelSize(0.1);
    h->GetYaxis()->SetTitleSize(0.13);
    h->GetYaxis()->SetTitle("DATA/FIT");
    h->GetXaxis()->SetTitle("");
    h->GetYaxis()->SetTitleOffset(0.45);
    h->GetXaxis()->SetTitleOffset(0.8);
    for (int i = 1; i< h->GetNbinsX();i++)std::cout<<h->GetBinContent(i)<<std::endl;
    TH1F* h_dataCopyr =(TH1F*) h->Clone();
    for (int i = 0; i< h_dataCopyr->GetNbinsX()+1;i++) if (h_dataCopyr->GetBinContent(i)==0) h_dataCopyr->SetBinContent(i, 1);
    h->SetMarkerSize(0);
    h->Draw("P");
    h->GetYaxis()->SetRangeUser(-1.,3 );
    h_dataCopyr->Divide(h_dataCopyr);
    h_dataCopyr->SetLineColor(kAzure-2);
    h_dataCopyr->SetFillStyle(3002);
    h_dataCopyr->SetMarkerSize(0.);
    std::cout<<"flag1"<<std::endl;   
      
    h_dataCopyr->Draw("HISTSAME");
    h->SetStats(kFALSE);
    h->Draw("PSAME");
    
    std::cout<<"flag2"<<std::endl;  
    
    
    c1.SaveAs("plots/BW_test.png"); 
    c1.SaveAs(TString::Format("~/www/plotsNota/BW_test_M%d.png",imass)); 
    c1.SaveAs("plots/BW_test.pdf"); 
    std::cout<<"flag3"<<std::endl;   
    
    std::cout<<"flag4"<<std::endl;
    
}





void SigModelBWFitHIST(int mass, bool ispythia6, bool ispythia8) {

  
  
  TCanvas c1("c1", "c1", 1);
 
  float min;
  float max;
  
  min =0;
  max = 1000;
  
  int nbin=100;
  int c = 0;

  double width=0.1*mass;
    int imass(mass);
    int iwidth(width); 
    std::cout<<"---------------------> "<<TString::Format("SigWeight_cat%d_m%d_w%d",c, imass, iwidth)<<std::endl;
      if(mass==200){
    minF=190;
    maxF=210;
    min=180;
    max=220;
    }else if(mass==400){
      minF=250;
      maxF=600;
      min=250;
      max=600;
    }else if(mass==600){
	minF=400;//200
	maxF=800;//900
      min=400;
      max=800;
    }
      else if(mass==800){
      minF=450;
      maxF=950;
      min=450;
      max=950;
    }

      TFile* f;
      if(ispythia6&&mass==600) f = TFile::Open("/afs/cern.ch/work/s/soffi/CMSSW536/src/ggH_HToGG_M_600_Width_60p0_TuneZ2star_8TeV_pythia6_tauola_cff_py_GEN_SIM_.root"); 
      if(ispythia6&&mass==800) f = TFile::Open("/afs/cern.ch/work/s/soffi/CMSSW536/src/ggH_HToGG_M_800_Width_80p0_TuneZ2star_8TeV_pythia6_tauola_cff_py_GEN.root"); 
       if(ispythia6&&mass==400) f = TFile::Open("/afs/cern.ch/work/s/soffi/CMSSW536/src/ggH_HToGG_M_400_Width_40p0_TuneZ2star_8TeV_pythia6_tauola_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_PU.root");

       if(ispythia8&&mass==600) f = TFile::Open("/afs/cern.ch/work/s/soffi/CMSSW_5_3_21/src/test_ggH_HToGG_pythia8_M600_G60.root"); 
       if(ispythia8&&mass==800) f = TFile::Open("/afs/cern.ch/work/s/soffi/CMSSW_5_3_21/src/test_ggH_HToGG_pythia8_M800_G80.root"); 
       if(ispythia8&&mass==400) f = TFile::Open("/afs/cern.ch/work/s/soffi/CMSSW_5_3_21/src/test_ggH_HToGG_pythia8_M400_G40.root");

       if(!ispythia6&&!ispythia8) f = TFile::Open(TString::Format("root://eoscms//eos/cms/store/group/phys_higgs/soffi/MC/HIG-Summer12PLHE-00021_genvalid_M%d.root",mass));
      TTree* tree = f->Get("Events");   
      TH1F* h = new TH1F("h", "h", 100, min, max);
      TH1F* h2 = new TH1F("h2", "h2", 900/h->GetBinWidth(3),100, 1000);
      std::string var;
      if(ispythia6)var="recoGenParticles_genParticles__AODSIMoneStep";
      if(ispythia8)var="recoGenParticles_genParticles__SIM";
      if(!ispythia8&&!ispythia6) var="recoGenParticles_genParticles__VALIDATION";
      
      tree->Draw((var+".obj.mass()>>h").c_str(), (var+".obj.pdgId()==25").c_str());
      tree->Draw((var+".obj.mass()>>h2").c_str(), (var+".obj.pdgId()==25").c_str());

         
      TF1 *bw = new TF1("bw"," ([0])/ ( TMath::Power( TMath::Power(x,2) - TMath::Power([1],2) , 2 ) + TMath::Power([1],2)*TMath::Power([2],2) )",min, max);
      bw->SetParameter(0, 1e14);
    //bw->FixParameter(1, 800);//722 600  400
      if(mass==400){
	bw->SetParLimits(1, 370,430 );//  700, 900  -- 570,630 ---  370,430
	bw->SetParLimits(2,10,60);// 200,400 ---80,160---   10,40
	if(ispythia6)	bw->FixParameter(2,40.);
      }else if(mass==600){
	bw->SetParLimits(1, 500,680 );//  700, 900  -- 570,630 ---  370,430
	bw->SetParLimits(2,20,2 so00);// 200,400 ---80,160---   10,40
	if(ispythia6)	bw->FixParameter(2,60.);
	  }else if(mass==800){
	bw->SetParLimits(1, 700,900 );//  700, 900  -- 570,630 ---  370,430
	bw->SetParLimits(2,40,400);// 200,400 ---80,160---   10,40
	if(ispythia6)	bw->FixParameter(2,80.);   
      }
      TH1F* h_bw = new TH1F("h_bw", "h_bw",100, min, max);
      
      h->Fit("bw","BR");
      h_bw = (TH1F*)bw->GetHistogram();
   
      /* double intg = h_bw->Integral(h_bw->FindBin(100), h_bw->FindBin(1000));
    double intg2p = h_bw->Integral(h_bw->FindBin(400), h_bw->FindBin(950));
    

    double ratio = intg2p/intg;
    std::cout<<"------> "<<intg<<"   "<<intg2p<<"   "<<TString::Format("%.3f",ratio)<<std::endl;
    */
    
    TPaveText* label_cms = get_labelCMS(0,  true);
    TPaveText* label_sqrt = get_labelSqrt(0);
    c1.Clear();
    TPad pad1("pad1", "pad1",0.01256281,0.23945,0.91206,1.);  
    pad1.Draw();    
    //pad1.SetLogy();  
    pad1.cd();
    
    TLegend* legmc = new TLegend(0.6, 0.6, 0.9, 0.8, "", "brNDC");
    legmc->SetTextSize(0.0406044);  
    legmc->SetTextFont(42);
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    std::string model("GGH");
    legmc->SetHeader(("Model: "+model).c_str());
    int imass(mass)
    legmc->AddEntry(h,TString::Format("m_{#gamma#gamma} = %d GeV",imass),"L");    
    legmc->AddEntry(bw," BW ","L");    
    
    bw->SetLineColor(kBlue);
    h2->SetStats(kFALSE);
    h2->GetXaxis()->SetTitle("m_{#gamma #gamma}^{true} [GeV]");
    h2->GetXaxis()->SetTitleFont(42);
    h2->GetXaxis()->SetTitleSize(0.05);
    h2->GetXaxis()->SetRangeUser(100, 1000);
    h2->Draw("hist");
    //h->Draw("histsame");
    bw->Draw("same");
    //  h_bw->Draw("same");


   
    label_cms->Draw("same");
    label_sqrt->Draw("same"); 
    legmc->Draw("same");    
    c1.cd();

    TPad pad3("pad3", "pad3",0.01256281,0.012510989,0.91206,0.22989);
    pad3.SetGrid();
    pad3.Draw();    
    pad3.cd();
    

    
    TH1F* hcopy = new TH1F("hcopy", "hcopy",100,100, 1000);
    for(int i = 1; i<h->GetNbinsX(); i++){
      hcopy->SetBinContent(hcopy->FindBin(h->GetBinCenter(i)), 1-((h->GetBinContent(i)-h_bw->GetBinContent(i))/h->Integral()));
      std::cout<<hcopy->GetBinContent(i)<<std::endl;
    }
    //    hcopy->Divide(h_bw);
    // TF1* func = new TF1("func", "[0]+[1]*(x-600)",min, max); //+[2]*x*x+[3]*x*x*x
    //TF1* func = new TF1("func", "[0]+[1]*(x-400)+[2]*pow(x-400, 2)+[3]*pow(x-400, 3)+[4]*pow(x-400, 4)+[5]*pow(x-400,5)+[6]*pow(x-400,6)+[7]*pow(x-400,7)",min, max); //+[2]*x*x+[3]*x*x*x
    //func->FixParameter(0, 1);
    //   hcopy->Fit("func", "BR");
    hcopy->SetStats(kFALSE);
 
    hcopy->GetYaxis()->SetRangeUser(0.995,1.005 );
    hcopy->GetYaxis()->SetNdivisions(505);
    hcopy->SetMarkerSize(0.4);
    hcopy->GetXaxis()->SetTitleSize(0.15);
    hcopy->GetXaxis()->SetLabelSize(0.1);
    hcopy->GetYaxis()->SetLabelSize(0.1);
    hcopy->GetYaxis()->SetTitleSize(0.13);
    hcopy->GetYaxis()->SetTitle("1- #frac{DATA_{i}-FIT_{i}}{DATA_{TOT}}");
    hcopy->GetXaxis()->SetTitle("");
    hcopy->GetYaxis()->SetTitleOffset(0.45);
    hcopy->GetXaxis()->SetTitleOffset(0.8);
    hcopy->Draw("hist");
    // func->Draw("same");
    c1.cd();

   
    c1.cd();
    c1.SetLogy();

    if(ispythia6)  c1.SaveAs(TString::Format("~/www/plotsNota/BW_test_M%d_HIST_pythia6.png",imass)); 
    if(ispythia8)  c1.SaveAs(TString::Format("~/www/plotsNota/BW_test_M%d_HIST_pythia8.png",imass)); 
    if(!ispythia6&&!ispythia8)  c1.SaveAs(TString::Format("~/www/plotsNota/BW_test_M%d_HIST_poweghMinlo.png",imass)); 
 
    bw->Delete();
}








void SigModelINTERF(int mass, int width) {

  
  
  TCanvas c1("c1", "c1", 1);
 
  float min;
  float max;
  
  min =0;
  max = 1000;
  
  int nbin=100;
  int c = 0;

  //  double width=0.1*mass;
    int imass(mass);
    int iwidth(width); 
    std::cout<<"---------------------> "<<TString::Format("SigWeight_cat%d_m%d_w%d",c, imass, iwidth)<<std::endl;
     if(mass==125){
    minF=124.99;
    maxF=125.01;
    min=124.99;
    max=125.01;
    }else if(mass==200){
    minF=100;
    maxF=600;
    min=100;
    max=600;
    }else if(mass==400){
      minF=200;
      maxF=600;
      min=200;
      max=600;
    }else if(mass==600){
	minF=400;//200
	maxF=800;//900
      min=400;
      max=800;
    }
      else if(mass==800){
      minF=400;
      maxF=1000;
      min=400;
      max=1000;
    }

      TFile* f_sig;
      TFile* f_int;
      f_sig = TFile::Open(TString::Format("/afs/cern.ch/work/s/soffi/CMSSW_5_3_21/src/dump/SoffiDumper/output_SIG_M%d_G%d.root",imass,iwidth)); 
      f_int = TFile::Open(TString::Format("/afs/cern.ch/work/s/soffi/CMSSW_5_3_21/src/dump/SoffiDumper/output_INT_M%d_G%d.root",imass,iwidth)); 
     
      TTree* tree_sig = f_sig->Get("soffidumper/myTree");   
      TTree* tree_int = f_int->Get("soffidumper/myTree");   
      TH1F* h_sig = new TH1F("h_sig", "h_sig", 100, min, max);
      TH1F* h_int = new TH1F("h_int", "h_int", 100, min, max);
      h_sig->Sumw2();
      h_int->Sumw2();
      
      tree_sig->Draw("mX>>h_sig", TString::Format("evweight*(mX>%f&&mX<%f)",min,max));
      tree_int->Draw("mX>>h_int", TString::Format("evweight*(mX>%f&&mX<%f)",min,max));


               
      TF1 *bw = new TF1("bw"," ([0])/ ( TMath::Power( TMath::Power(x,2) - TMath::Power([1],2) , 2 ) + TMath::Power([1],2)*TMath::Power([2],2) )",min, max);
      bw->SetParameter(0, h_sig->Integral());
    //bw->FixParameter(1, 800);//722 600  400
      if(mass==125){
	bw->FixParameter(2,0.004);
	bw->FixParameter(1,125);
      }else if(mass==200){
	bw->SetParLimits(1, 170,230 );//  700, 900  -- 570,630 ---  370,430
	//	bw->SetParLimits(2,10,60);// 200,400 ---80,160---   10,40
	bw->FixParameter(2,20.);
	bw->FixParameter(1,200.);
      }else if(mass==400){
	bw->SetParLimits(1, 385,415 );//  700, 900  -- 570,630 ---  370,430
	//	bw->SetParLimits(2,10,60);// 200,400 ---80,160---   10,40
	bw->FixParameter(2,40.);
	bw->FixParameter(1,400.);
      }else if(mass==600){
	bw->SetParLimits(1, 500,680 );//  700, 900  -- 570,630 ---  370,430
	//	bw->SetParLimits(2,20,2 so00);// 200,400 ---80,160---   10,40
	bw->FixParameter(2,60.);
	  }else if(mass==800){
	bw->SetParLimits(1, 700,900 );//  700, 900  -- 570,630 ---  370,430
	//	bw->SetParLimits(2,40,400);// 200,400 ---80,160---   10,40
       	bw->FixParameter(2,80.);   
	bw->FixParameter(1,800.);   
      }
  
     
    
      double scale;
      if(mass==125) scale=1;
      if(mass==200) scale=0.00354;
      if(mass==400) scale=0.00152;
      if(mass==800) scale=0.0000293;
 
      double sig_intg = h_sig->Integral();
      double int_intg = h_int->Integral();

      //  h_sig->Scale(scale*scale);
      h_int->Scale(scale);
    

      TH1F* h_sigint = h_sig->Clone(); 
      h_sigint->Sumw2();
      h_sigint->Add(h_int);
  

      std::cout<<h_sigint->Integral()<<std::endl;
      //bw->FixParameter(0,h_sigint->Integral());
      h_sigint->Scale(1./ h_sigint->Integral());
      TH1F* h_b = new TH1F("h_b", "h_b",h_sigint->GetNbinsX(), min, max); 
      h_b->Sumw2();    
      h_sigint->Fit("bw","BR");
      h_b = (TH1F*)bw->GetHistogram();
  
    
      TPaveText* label_cms = get_labelCMS(0,  true);
      TPaveText* label_sqrt = get_labelSqrt(0);
      c1.Clear();
      TPad pad1("pad1", "pad1",0.01256281,0.23945,0.91206,1.);  
      pad1.Draw();    
      //pad1.SetLogy();  
      pad1.cd();
    
      TLegend* legmc;
      if(mass==800)legmc = new TLegend(0.25, 0.6, 0.55, 0.8, "", "brNDC");
      else legmc = new TLegend(0.55, 0.6, 0.85, 0.8, "", "brNDC");
      legmc->SetTextSize(0.0406044);  
      legmc->SetTextFont(42);
      legmc->SetBorderSize(0);
      legmc->SetFillStyle(0);
      std::string model("GGH");
      legmc->SetHeader(("Model: "+model).c_str());
      int imass(mass);
      // legmc->AddEntry(h_sig,TString::Format("#splitline{m^{true}_{#gamma#gamma} = %d [GeV]}{  }",imass),"L");    
      legmc->AddEntry(bw,TString::Format("#splitline{BW M %d [GeV]}{#Gamma %d [GeV]} ",imass, iwidth),"L");    
    
      // bw->SetLineColor(kBlue);
      h_sigint->SetStats(kFALSE);
      h_sigint->GetXaxis()->SetTitle("m_{#gamma #gamma}^{true} [GeV]");
      h_sigint->GetXaxis()->SetTitleFont(42);
      h_sigint->GetXaxis()->SetTitleSize(0.05);
      h_sigint->GetXaxis()->SetRangeUser(100, 1000);
     
      h_b->Draw("L");
      // bw->Draw("Lsame");
      h_sigint->Draw("histsame");
         
      label_cms->Draw("same");
      label_sqrt->Draw("same"); 
      legmc->Draw("same");    
      c1.cd();
      
      TPad pad3("pad3", "pad3",0.01256281,0.072510989,0.91206,0.27989);
      pad3.SetGrid();
      pad3.Draw();    
      pad3.cd();
    
      
     
      h_b->Scale(1./ h_b->Integral());
      TH1F* hcopy=h_sigint->Clone();
      // hcopy->Scale(1./ hcopy->Integral());
      std::cout<<h_sigint->GetNbinsX()<<" "<<h_b->GetNbinsX()<<std::endl;
      hcopy->Add(h_b, -1);
      //hcopy->Scale(1./hcopy->Integral());
     
    hcopy->SetStats(kFALSE);
 
    // hcopy->GetYaxis()->SetRangeUser(0.95,1.05 );
    hcopy->GetYaxis()->SetRangeUser(-0.02, 0.02 );
    hcopy->GetYaxis()->SetNdivisions(505);
    hcopy->SetMarkerSize(0.4);
    hcopy->GetXaxis()->SetTitleSize(0.15);
    hcopy->GetXaxis()->SetLabelSize(0.1);
    hcopy->GetYaxis()->SetLabelSize(0.1);
    hcopy->GetYaxis()->SetTitleSize(0.13);
    hcopy->GetYaxis()->SetTitle(" #frac{DATA_{i}-FIT_{i}}{DATA_{TOT}}");
    hcopy->GetXaxis()->SetTitle("m_{#gamma #gamma}^{true} [GeV]");
    hcopy->GetYaxis()->SetTitleOffset(0.45);
    hcopy->GetXaxis()->SetTitleOffset(0.8);
    hcopy->Draw("hist");

    TF1* func = new TF1("func", "[0]+[1]*x+[2]*pow(x, 2)+[3]*pow(x, 3)+[4]*pow(x, 4)+[5]*pow(x,5)+[6]*pow(x,6)+[7]*pow(x,7)",min, max); //+[2]*x*x+[3]*x*x*x   +[4]*pow(x-400, 4)+[5]*pow(x-400,5)+[6]*pow(x-400,6)+[7]*pow(x-400,7)
    //func->FixParameter(0, 1);
    // hcopy->Fit("func", "BR");


    // func->Draw("same");

    
    c1.cd();

   
    c1.cd();
    c1.SetLogy();
    c1.SaveAs(TString::Format("~/www/plotsNota/SIGINT_test_M%d_G%d_HIST.png",imass, width)); 
   

}






 

void SigModelTheoFcnFitCPS() {

  /* RooDataSet* signal;

    
  signal= (RooDataSet*) w->data(TString::Format("SigWeight",c, 800,5));

 */
  RooRealVar* PhotonsMassTrue = new RooRealVar("PhotonsMassTrue", "PhotonsMassTrue", 0, 1500);//w->var("PhotonsMassTrue"); 
  
  TCanvas c1("c1", "c1", 1);
 
  int nbin=100;
  int c = 0;
  std::string model("GGH");
  double mass=800;
  double width=20;
    int imass(mass);
    int iwidth(width); 
   
   //BW1 
    RooRealVar sigmaBW1_var(TString::Format("sigmaBW1_var_cat%d_m%d_w%d",c, imass, iwidth), TString::Format("sigmaBW1_var_cat%dm%d_w%d",c, imass, iwidth), 0., 0.2*mass);
    RooRealVar meanBW_var(TString::Format("meanBW_var_cat%d_m%d_w%d",c, imass, iwidth), TString::Format("meanBW_var_cat%d_m%d_w%d",c, imass, iwidth),0.7*mass, 1.3*mass);
    sigmaBW1_var.setVal(0.1*mass);
    meanBW_var.setVal(mass);
    sigmaBW1_var.setConstant();
    RooFormulaVar sigmaBW1 (TString::Format("sigmaBW1_cat%d_m%d_w%d",c, imass, iwidth),"","@0",sigmaBW1_var); 
    RooFormulaVar meanBW (TString::Format("meanBW_cat%d_m%d_w%d",c, imass, iwidth),"","@0",meanBW_var); 
    RooBreitWigner BW1(TString::Format("BW1_cat%d_m%d_w%d",c,imass, iwidth),TString::Format("BW1_cat%d_m%d_w%d",c,imass, iwidth), *PhotonsMassTrue, meanBW, sigmaBW1);

    //BW2 
    RooRealVar sigmaBW2_var(TString::Format("sigmaBW2_var_cat%d_m%d_w%d",c, imass, iwidth), TString::Format("sigmaBW2_var_cat%dm%d_w%d",c, imass, iwidth), 0., 0.2*mass);
    sigmaBW2_var.setVal(0.05*mass);
    sigmaBW2_var.setConstant();
    RooFormulaVar sigmaBW2 (TString::Format("sigmaBW2_cat%d_m%d_w%d",c, imass, iwidth),"","@0",sigmaBW2_var); 
    RooBreitWigner BW2(TString::Format("BW2_cat%d_m%d_w%d",c,imass, iwidth),TString::Format("BW2_cat%d_m%d_w%d",c,imass, iwidth), *PhotonsMassTrue, meanBW, sigmaBW2);
   
    //BW3 
    RooRealVar sigmaBW3_var(TString::Format("sigmaBW3_var_cat%d_m%d_w%d",c, imass, iwidth), TString::Format("sigmaBW3_var_cat%dm%d_w%d",c, imass, iwidth), 0., 0.2*mass);
    sigmaBW3_var.setVal(103);
    sigmaBW3_var.setConstant();
    RooFormulaVar sigmaBW3 (TString::Format("sigmaBW3_cat%d_m%d_w%d",c, imass, iwidth),"","@0",sigmaBW3_var); 
    RooBreitWigner BW3(TString::Format("BW3_cat%d_m%d_w%d",c,imass, iwidth),TString::Format("BW3_cat%d_m%d_w%d",c,imass, iwidth), *PhotonsMassTrue, meanBW, sigmaBW3);
 

    TFile* f = new TFile("/tmp/soffi/HIG-Summer12PLHE-00021_genvalid_M800_P15.root", "READ");
    TTree* tree = f->Get("Events");   
    TH1F* h = new TH1F("h", "h", 100,0,1500);
    tree->Draw("recoGenParticles_genParticles__VALIDATION.obj.mass()>>h", "recoGenParticles_genParticles__VALIDATION.obj.pdgId()==25");
    RooDataHist hist("hist", "hist", *PhotonsMassTrue, h, 1.);  

  
    
    //CPS
    RooRealVar cp1_var(TString::Format("cp1_var_cat%d_m%d_w%d",c, imass, iwidth), TString::Format("cp1_var_cat%d_m%d_w%d",c, imass, iwidth),0.,100.);
    cp1_var.setVal(13.3);
    cp1_var.setConstant();
    RooRealVar cp2_var(TString::Format("cp2_var_cat%d_m%d_w%d",c, imass, iwidth), TString::Format("cp2_var_cat%d_m%d_w%d",c, imass, iwidth),0.,100.);
    cp2_var.setVal(6.66);
    cp2_var.setConstant();
    RooRealVar cp3_var(TString::Format("cp3_var_cat%d_m%d_w%d",c, imass, iwidth), TString::Format("cp3_var_cat%d_m%d_w%d",c, imass, iwidth),0.,100.);
    cp3_var.setVal(1.);
    cp3_var.setConstant();
    RooFormulaVar cp1(TString::Format("cp1_cat%d_m%d_w%d",c, imass, iwidth),"","@0",cp1_var);
    RooFormulaVar cp2(TString::Format("cp2_cat%d_m%d_w%d",c, imass, iwidth),"","@0",cp2_var);
    RooFormulaVar cp3(TString::Format("cp3_cat%d_m%d_w%d",c, imass, iwidth),"","@0",cp3_var);
    RooRealVar br_var(TString::Format("br_var_cat%d_m%d_w%d",c, imass, iwidth), TString::Format("br_var_cat%d_m%d_w%d",c, imass, iwidth),0.);
    br_var.setVal(0.);
    br_var.setConstant();
    RooFormulaVar br(TString::Format("br_cat%d_m%d_w%d",c, imass, iwidth),"","@0",br_var);
    RooRealVar is8_var(TString::Format("is8_var_cat%d_m%d_w%d",c, imass, iwidth), TString::Format("is8_var_cat%d_m%d_w%d",c, imass, iwidth),1.);
    is8_var.setVal(1.);
    is8_var.setConstant();
    Bool_t is8 =  true;
  //  RooFormulaVar is8(TString::Format("is8_cat%d_m%d_w%d",c, imass, iwidth),"","@0",is8_var);
   
    RooCPSHighMassGGHNoInterf cps1(TString::Format("CPS1_cat%d_m%d_w%d",c,imass, iwidth),TString::Format("CPS1_cat%d_m%d_w%d",c,imass, iwidth),*PhotonsMassTrue, meanBW, cp1,br,is8);
    RooCPSHighMassGGHNoInterf cps2(TString::Format("CPS2_cat%d_m%d_w%d",c,imass, iwidth),TString::Format("CPS2_cat%d_m%d_w%d",c,imass, iwidth),*PhotonsMassTrue, meanBW, cp2,br,is8);
    RooCPSHighMassGGHNoInterf cps3(TString::Format("CPS3_cat%d_m%d_w%d",c,imass, iwidth),TString::Format("CPS3_cat%d_m%d_w%d",c,imass, iwidth),*PhotonsMassTrue, meanBW, cp3,br,is8);

    //cps3.fitTo(hist);


    RooPlot* plotG = PhotonsMassTrue->frame();
    //hist->plotOn(plotG, RooFit::Invisible());
    hist->plotOn(plotG);
    cps3.plotOn(plotG, LineColor(kBlue));
    
    //cps.GetYaxis()->SetRangeUser(0.01,plotG->GetMaximum()*10 );
    plotG->GetXaxis()->SetTitle("m_{#gamma #gamma}^{true} [GeV]");
    plotG->GetXaxis()->SetTitleFont(42);
    plotG->GetXaxis()->SetTitleSize(0.05);
    TPaveText* label_cms = get_labelCMS(0,  true);
    TPaveText* label_sqrt = get_labelSqrt(0);

    TF1 *bw = new TF1("bw"," ([0])/ ( TMath::Power( TMath::Power(x,2) - TMath::Power([1],2) , 2 ) + TMath::Power([1],2)*TMath::Power([2],2) )",0, 1500);
    bw->FixParameter(0, 100000);
    bw->FixParameter(1, 600);//722 600  400
    // bw->SetParLimits(1,  700, 900 );//  700, 900  -- 570,630 ---  370,430
    bw->FixParameter(2,60);//283 120 30
    //bw->SetParLimits(2, 200,400);// 200,400 ---80,160---   10,40
    bw->SetLineColor(kRed);
    TH1F* h_bw = new TH1F("h_bw", "h_bw",300, min, max );
    // h->Fit("bw","BR");
    h_bw = (TH1F*)bw->GetHistogram();
    h_bw->Scale(h->Integral()/h_bw->Integral());

    TLegend* legmc = new TLegend(0.2, 0.55, 0.45, 0.8, "", "brNDC");
    legmc->SetTextSize(0.0206044);  
    legmc->SetTextFont(42);
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->SetHeader(("Model: "+model).c_str());
    int imass(mass);
    legmc->AddEntry(plotG->getObject(0),"POWHEG 15","L");    
    legmc->AddEntry(plotG->getObject(1),"CPS (M 800, c' 1)","L");
    //legmc->AddEntry(h_bw,"BW(M 600, W 60","L");
     
    plotG->Draw();    
    plotG->SetTitle("");
    label_cms->Draw("same");
    label_sqrt->Draw("same"); 
    legmc->Draw("same");

    
    
    
    
    //h_bw->Draw("histsame");
    std::cout<<    h_bw->Integral()<<std::endl;
    c1.cd();
    c1.SetLogy();
  
    c1.SaveAs(TString::Format("plots/CPS_BW_cat%d_LOG_M%d_w%d.png",c,imass, iwidth)); 
    c1.SaveAs(TString::Format("~/www/plotsNota/CPS_BW_cat%d_M%d_w%d.png",c,imass, iwidth)); 
    c1.SaveAs(TString::Format("plots/CPS_BW_cat%d_LOG_M%d_w%d.pdf",c,imass, iwidth)); 
   
}







TGraph2D* makeDELTA(){

  double m[10] = {150., 200.,300.,400, 500, 600, 700, 800, 900, 1000};
  double cp[24] = {0.1, 0.2, 0.4, 0.6, 0.8, 1., 2, 5, 7, 10, 12, 15,17,  20,22, 25, 27, 30, 32, 35, 37, 40, 42, 45};
  TH2D* h0 = new TH2D("h0", "h0",9, m, 23,cp);
  
  double a_delta[7][5];
  
  a_delta[0][0] = 0.015081577003;
  a_delta[0][1] = 0.0576512292027;
 a_delta[0][2] = 0.131081297994;
 a_delta[0][3] = 0.243380963802;
 a_delta[0][4] = 0.39500105381;
 a_delta[1][0] = 0.0972077324986;
 a_delta[1][1] = 0.434791147709;
 a_delta[1][2] = 0.861800074577;
 a_delta[1][3] = 1.58613061905;
 a_delta[1][4] = 2.41115784645;
 a_delta[2][0] = 0.242783889174;
 a_delta[2][1] = 0.827510118484;
 a_delta[2][2] = 1.78883588314;
 a_delta[2][3] = 3.17804527283;
 a_delta[2][4] = 4.9168419838;
 a_delta[3][0] = 0.368439763784;
 a_delta[3][1] = 1.45783913136;
 a_delta[3][2] = 3.27468585968;
 a_delta[3][3] = 5.93577003479;
 a_delta[3][4] = 9.49385356903;
 a_delta[4][0] = 0.679009974003;
 a_delta[4][1] = 2.96808123589;
 a_delta[4][2] = 6.60964488983;
 a_delta[4][3] = 11.7470684052;
 a_delta[4][4] = 18.1083869934;
 a_delta[5][0] = 1.51261091232;
 a_delta[5][1] = 5.85586261749;
 a_delta[5][2] = 12.3102693558;
 a_delta[5][3] = 21.393032074;
 a_delta[5][4] = 31.7420043945;
 a_delta[6][0] = 1.7427560091;
 a_delta[6][1] = 7.33418560028;
 a_delta[6][2] = 15.4683561325;
 a_delta[6][3] = 26.5165843964;
 a_delta[6][4] = 37.535987854;
 
 double mOld[10] = {400, 500, 600, 700, 800, 900, 1000};
 double cpOld[24] = { 0.2, 0.4, 0.6, 0.8, 1.};

 double deltaM6[6] = {a_delta[6][0], a_delta[6][1], a_delta[6][2], a_delta[6][3], a_delta[6][4]};
 for(int j=0;j<5;j++)deltaM6[j]/=a_delta[6][4];
 TGraph* gr6 = new TGraph(5, cpOld, deltaM6);
 double deltaM5[6] = {a_delta[5][0], a_delta[5][1], a_delta[5][2], a_delta[5][3], a_delta[5][4]};
 for(int j=0;j<5;j++)deltaM5[j]/=a_delta[5][4];
 TGraph* gr5 = new TGraph(5, cpOld, deltaM5);
 double deltaM4[6] = {a_delta[4][0], a_delta[4][1], a_delta[4][2], a_delta[4][3], a_delta[4][4]};
 for(int j=0;j<5;j++)deltaM4[j]/=a_delta[4][4];
 TGraph* gr4 = new TGraph(5, cpOld, deltaM4);
 double deltaM3[6] = {a_delta[3][0], a_delta[3][1], a_delta[3][2], a_delta[3][3], a_delta[3][4]};
 for(int j=0;j<5;j++)deltaM3[j]/=a_delta[3][4];
 TGraph* gr3 = new TGraph(5, cpOld, deltaM3);
 double deltaM2[5] = {a_delta[2][0], a_delta[2][1], a_delta[2][2], a_delta[2][3], a_delta[2][4]};
 for(int j=0;j<5;j++)deltaM2[j]/=a_delta[2][4];
 TGraph* gr2 = new TGraph(5, cpOld, deltaM2);
 double deltaM1[5] = {a_delta[1][0], a_delta[1][1], a_delta[1][2], a_delta[1][3], a_delta[1][4]};
 for(int j=0;j<5;j++)deltaM1[j]/=a_delta[1][4];
 TGraph* gr1 = new TGraph(5, cpOld, deltaM1);
 double deltaM0[5] = {a_delta[0][0], a_delta[0][1], a_delta[0][2], a_delta[0][3], a_delta[0][4]};
 for(int j=0;j<5;j++)deltaM0[j]/=a_delta[0][4];
 TGraph* gr0 = new TGraph(5, cpOld, deltaM0);

 double deltacp1[7] = {a_delta[0][4], a_delta[1][4], a_delta[2][4], a_delta[3][4], a_delta[4][4], a_delta[5][4], a_delta[6][4]};
 TGraph* grcp1 = new TGraph(7, mOld, deltacp1);

    TLegend* legmc = new TLegend(0.2, 0.65, 0.55, 0.8, "", "brNDC");
    legmc->SetTextSize(0.0206044);  
    legmc->SetTextFont(42);
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    std::string model("GGH");
    legmc->SetHeader(("Model: "+model).c_str());

    legmc->AddEntry(gr0,"M 400","LP");  
    legmc->AddEntry(gr1,"M 500","LP");  
    legmc->AddEntry(gr2,"M 600","LP");  
    legmc->AddEntry(gr3,"M 700","LP");  
    legmc->AddEntry(gr4,"M 800","LP");  
    legmc->AddEntry(gr5,"M 900","LP");  
    legmc->AddEntry(gr6,"M 1000","LP"); 
    legmc->SetHeader("a_delta[i][j]/a_delta[i][4] -> (a_delta @ M XX && C' YY) / (a_delta @ C' = 1)"); 

 
 for(int i=0;i<7;i++) for(int j=0;j<5;j++) h0->SetBinContent(h0->FindBin(mOld[i],cpOld[j]), a_delta[i][j]);
      
 //   TGraph2D* g0 = new TGraph2D(h0);
    
      gr0->SetTitle("delta");
      gr0->SetName("delat");
      TCanvas* ca = new TCanvas("c", "c",600);
      ca->cd();
      TPad* pad = new TPad("pad","pad", 0.01, 0.03, 0.9, 0.9.);
      pad->Draw();
      pad->cd();
      TPaveText* label_cms = get_labelCMS(0, false);
      TPaveText* label_sqrt = get_labelSqrt(0);  
      
      h0->SetStats(0);
      //  hp0->SetTitle("");
     
      gr0->GetXaxis()->SetTitle("m_{H} [GeV]");
      gr0->GetYaxis()->SetTitle("C'");
      //  h0->Draw("axis");

      gr0->SetLineColor(15);
      gr1->SetLineColor(1);
      gr2->SetLineColor(2);
      gr3->SetLineColor(3);
      gr4->SetLineColor(4);
      gr5->SetLineColor(5);
      gr6->SetLineColor(6);
      gr0->SetMarkerColor(15);
      gr1->SetMarkerColor(1);
      gr2->SetMarkerColor(2);
      gr3->SetMarkerColor(3);
      gr4->SetMarkerColor(4);
      gr5->SetMarkerColor(5);
      gr6->SetMarkerColor(6);
      gr0->Draw("APL");
      TF1* f4 = new TF1("f4", "[0]+[1]*x+[2]*x*x");
      gr4->Fit(f4);
      gr1->Draw("LPsame");
      gr2->Draw("LPsame");
      gr3->Draw("LPsame");
      gr4->Draw("LPsame");
      gr5->Draw("LPsame");
      gr6->Draw("LPsame");
         


      legmc->Draw("same");
      label_cms->Draw("same");
      label_sqrt->Draw("same");
     
      ca->SaveAs(TString::Format("/afs/cern.ch/user/s/soffi/www/plotsNota/delta.png"));
      ca->SaveAs(TString::Format("/afs/cern.ch/user/s/soffi/www/plotsNota/delta.pdf"));

      grcp1->GetXaxis()->SetTitle("m_{H} [GeV]");
      grcp1->GetYaxis()->SetTitle("C'");
      legmc->Clear();
      legmc->SetHeader("a_delta[i][4] -> a_delta @ C' = 1");
      TF1* fcp1 = new TF1("fcp1", "[0]+[1]*x+[2]*x*x");
      grcp1->Fit(fcp1);
      grcp1->Draw("APE");
      label_cms->Draw("same");
      label_sqrt->Draw("same");
      legmc->Draw("same");
     
      ca->SaveAs(TString::Format("/afs/cern.ch/user/s/soffi/www/plotsNota/deltaCP1.png"));
      ca->SaveAs(TString::Format("/afs/cern.ch/user/s/soffi/www/plotsNota/deltaCP1.pdf"));

      for(int i=0;i<7;i++) for(int j=0;j<5;j++)std::cout<<"M : "<<mOld[i]<<" C': "<<cpOld[j]<< " a_delta: "<<a_delta[i][j]<< "   FIT: "<< f4->Eval(cpOld[j])<< "   "<<fcp1->Eval(mOld[i])<<"   "<<f4->Eval(cpOld[j])*fcp1->Eval(mOld[i])<<std::endl;


}



TGraph2D* makeWIDTH(){

  double a_width[7][5];

    a_width[0][0] = 26.5976848602;
    a_width[0][1] = 26.5976848602;
    a_width[0][2] = 26.7212715149;
    a_width[0][3] = 26.5424861908;
    a_width[0][4] = 26.4101047516;
    a_width[1][0] = 58.7029762268;
    a_width[1][1] = 58.7029762268;
    a_width[1][2] = 58.7458114624;
    a_width[1][3] = 58.1389541626;
    a_width[1][4] = 58.612991333;
    a_width[2][0] = 103.93271637;
    a_width[2][1] = 103.93271637;
    a_width[2][2] = 103.803878784;
    a_width[2][3] = 103.800476074;
    a_width[2][4] = 103.74822998;
    a_width[3][0] = 162.97102356;
    a_width[3][1] = 162.97102356;
    a_width[3][2] = 162.841094971;
    a_width[3][3] = 162.802963257;
    a_width[3][4] = 162.688156128;
    a_width[4][0] = 235.570587158;
    a_width[4][1] = 235.570587158;
    a_width[4][2] = 235.355682373;
    a_width[4][3] = 235.014694214;
    a_width[4][4] = 235.070922852;
    a_width[5][0] = 320.552429199;
    a_width[5][1] = 320.552429199;
    a_width[5][2] = 320.081390381;
    a_width[5][3] = 319.721862793;
    a_width[5][4] = 319.82611084;
    a_width[6][0] = 416.118835449;
    a_width[6][1] = 416.118835449;
    a_width[6][2] = 430.821716309;
    a_width[6][3] = 416.076812744;
    a_width[6][4] = 415.9637146;

 
    

    
    
    double mOld[10] = {0.,400, 500, 600, 700, 800, 900, 1000};
    double width[8];
    width[0] = 0;
    for(int i=1;i<8;i++){
      float sum=0;
      
      for(int j=0;j<5;j++){
	sum+=a_width[i-1][j]; 
      }
      width[i] = sum/5;
      std::cout<<"m[i]: "<<mOld[i]<<" w[i]: "<<width[i]<<std::endl;
    }
    
    TGraph* gr = new TGraph(8, mOld, width);

    TLegend* legmc = new TLegend(0.2, 0.65, 0.55, 0.8, "", "brNDC");
    legmc->SetTextSize(0.0206044);  
    legmc->SetTextFont(42);
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    std::string model("GGH");
    legmc->SetHeader(("Model: "+model).c_str());

    legmc->AddEntry(gr,"width","LP");  
    
//legmc->SetHeader("a_delta[i][j]/a_delta[i][4] -> (a_delta @ M XX && C' YY) / (a_delta @ C' = 1)"); 

      gr->SetTitle("delta");
      gr->SetName("delat");
      TCanvas* ca = new TCanvas("c", "c",600);
      ca->cd();
      TPad* pad = new TPad("pad","pad", 0.01, 0.03, 0.9, 0.9.);
      pad->Draw();
      pad->cd();
      TPaveText* label_cms = get_labelCMS(0, false);
      TPaveText* label_sqrt = get_labelSqrt(0);  
     
      gr->GetXaxis()->SetTitle("m_{H} [GeV]");
//gr->GetYaxis()->SetTitle("C'");
      //  h0->Draw("axis");

      gr->Draw("APL");
      TF1* f = new TF1("f", "[0]+[1]*x+[2]*x*x+[3]*x*x*x", 400, 1000);
      TF1* flow = new TF1("f", "[0]+[1]*x", 0, 400);
      gr->Fit(f, "R");     
      legmc->Draw("same");
      std::cout<<f->Eval(400.)<<std::endl;
      gr->Fit(flow, "R");
      std::cout<<flow->Eval(400.)<<std::endl;
      f->Draw("same");
      flow->Draw("same");
      label_cms->Draw("same");
      label_sqrt->Draw("same");
     
      ca->SaveAs(TString::Format("/afs/cern.ch/user/s/soffi/www/plotsNota/width.png"));
      ca->SaveAs(TString::Format("/afs/cern.ch/user/s/soffi/www/plotsNota/width.pdf"));

}



// Fit signal with model gauss pdfs
void SigModelFitConvBW(RooWorkspace* w, Float_t mass, Double_t width, std::string model) {

  Int_t ncat = NCAT;
  RooDataSet* sigToFit[NCAT];
  
  
  Float_t MASS(mass);  
  Float_t minMassFit(mass*0.8);
  Float_t maxMassFit(mass*1.2); 
  RooRealVar* PhotonsMass = w->var("PhotonsMass"); 
  Double_t scaleSyst;
  Double_t smearSyst;

  std::cout<<"----------------------------------------------------------------------------------------"<<std::endl;
  //  PhotonsMass->setRange("sigrange",minMassFit-20,maxMassFit+20); 
  TPaveText* label_cms = get_labelCMS(0, "2012", true);
  TPaveText* label_sqrt = get_labelSqrt(0);
  TFile* f = new TFile("sigShapeCorrections.root", "READ");
  // Fit Signal 
  for (int c=0; c<ncat; ++c) {
    cout << "---------- Category = " << c << endl;
    //  sigToFit[c] = (RooDataSet*) w->data(TString::Format("SigWeight_cat%d",c));
    //  w->import(*sigToFit[c]);

    //introduce systs 
    if(c==0 || c==1)scaleSyst = 0.005;
    if(c==2 || c==3)scaleSyst = 0.007;
    if(c==0)smearSyst = 0.005;
    if(c==1)smearSyst = 0.0058;
    if(c==2 || c==3)smearSyst = 0.01;

    //get sigma from TF1:   
    TF1* fcn = (TF1*)f->Get(TString::Format("f%d",c));
    Float_t massF = (Float_t) (*w->var("MH")).getVal();
    Float_t sigmaCorr = fcn->Eval(massF);
    if(massF==150)sigmaCorr=1;
    //std::cout<<"Mass: "<<massF<<" corr: "<<sigmaCorr<<std::endl;
    RooRealVar rooSigmaCorr (TString::Format("rooSigmaCorr_cat%d",c), TString::Format("rooSigmaCorr_cat%d",c), sigmaCorr, "");
    rooSigmaCorr.setConstant();
    w->import(rooSigmaCorr);


     ( *w->var(TString::Format("mShift_cat%d",c))).setConstant();
     ( *w->var(TString::Format("mSmear_cat%d",c))).setConstant();
    //cb
    RooFormulaVar CBpos_mean_draw(TString::Format("CBpos_mean_draw_cat%d",c),"","@0+@1",RooArgList(*w->var(TString::Format("ReducedMass_sig_mean_cat%d",c)),*w->var("MH")));
    RooFormulaVar CBpos_mean(TString::Format("CBpos_mean_cat%d",c),"",TString::Format("@0+%f*@1", scaleSyst),RooArgList(*w->var(TString::Format("ReducedMass_sig_mean_cat%d",c)),*w->var(TString::Format("mShift_cat%d",c))));
    //RooFormulaVar CBpos_mean(TString::Format("CBpos_mean_cat%d",c),"",TString::Format("@0"),*w->var(TString::Format("ReducedMass_sig_mean_cat%d",c)));   
    RooFormulaVar CBpos_sigma(TString::Format("CBpos_sigma_cat%d",c),"",TString::Format("(sqrt(@0*@0*@3*@3+%f*%f*@2)*@1)",smearSyst,smearSyst),RooArgList(*w->var(TString::Format("ReducedMass_sig_sigma_cat%d",c)),*w->var("MH"),*w->var(TString::Format("mSmear_cat%d",c)),*w->var(TString::Format("rooSigmaCorr_cat%d",c)) ) );

    //    std::cout<<"-------------------> SIGMA: "<<CBpos_sigma->getVal()<<"    MASS: "<<(*w->var("MH")).getVal()<<std::endl;


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
  
    RooRealVar sigmaBW_var(TString::Format("sigmaBW_var_cat%d",c), TString::Format("sigmaBW_var_cat%d",c), 30., 50.);
    RooRealVar meanBW_var(TString::Format("meanBW_var_cat%d",c), TString::Format("meanBW_var_cat%d",c),350., 450.);
    std::cout<<" width:--------> "<<width<<std::endl;
    sigmaBW_var.setVal(40);
    meanBW_var.setVal(400);
    //sigmaBW_var.setConstant();
    //  meanBW_var.setConstant();
    w->import(sigmaBW_var);
    w->import(meanBW_var);
    
    RooFormulaVar meanBW(TString::Format("meanBW_cat%d",c),"","@0",*w->var(TString::Format("meanBW_var_cat%d",c)));    
    RooFormulaVar sigmaBW (TString::Format("sigmaBW_cat%d",c),"","@0",*w->var(TString::Format("sigmaBW_var_cat%d",c))); 
    //else if(width==2) sigmaBW = new RooFormulaVar(TString::Format("sigmaBW_cat%d",c),"","@0*0.02",*w->var("MH"));   
    //else if(width==5) sigmaBW = new RooFormulaVar(TString::Format("sigmaBW_cat%d",c),"","@0*0.05",*w->var("MH"));   
    //else if(width==7) sigmaBW = new RooFormulaVar(TString::Format("sigmaBW_cat%d",c),"","@0*0.07",*w->var("MH"));   
    //   else if(width==10) sigmaBW = new RooFormulaVar(TString::Format("sigmaBW_cat%d",c),"","@0*0.10",*w->var("MH"));   
    RooBreitWigner SigModelBW(TString::Format("SigModelBW_cat%d",c),TString::Format("SigModelBW_cat%d",c), *PhotonsMass, meanBW, sigmaBW);

  
    //CONV 
    RooFFTConvPdf  ConvolutedRes_CB(TString::Format("PhotonsMassSig_cat%d",c),TString::Format("PhotonsMassSig_cat%d",c), *PhotonsMass,SigModelBW, ResAddPdf);
    //RooFitResult* fitresults_CB = (RooFitResult* ) ConvolutedRes_CB.fitTo(*sigToFit[c],RooFit::Range("sigrange"), RooFit::Save(kTRUE));
    // std::cout<<TString::Format("******************************** Signal Fit results CB mass %f cat %d***********************************", mass, c)<<std::endl;
    // fitresults_CB->Print("V");
    w->import(ConvolutedRes_CB);
    // std::cout<<".............> "<<c<<std::endl;
    
    RooHistFunc* rooFunc_norm = getRooHistFunc(c,w->var("MH"), model );
    w->import(*rooFunc_norm);

    std::cout<<"SIG NORM ----->"<<rooFunc_norm->getVal(*w->var("MH"));
    // w->Print("V");

    if(width ==10 ){ //if i want to plot the fit
      sigToFit[c] = (RooDataSet*) w->data(TString::Format("SigWeight_cat%d",c));

      RooFitResult* fitresults_CB = (RooFitResult* ) ConvolutedRes_CB.fitTo(*sigToFit[c],SumW2Error(kTRUE), RooFit::Save(kTRUE));
      fitresults_CB->Print("V");
      w->import(ConvolutedRes_CB);
      RooPlot* plotOnlyResPdf = w->var("PhotonsMass")->frame(Range(minMassFit-20,maxMassFit+20),Bins(100));
      sigToFit[c]->plotOn(plotOnlyResPdf, LineColor(kRed), LineStyle(kDashed));
      double max = plotOnlyResPdf->GetMaximum();
      
      // Plot to verify everything is ok
      RooPlot* plotPhotonsMassAll = w->var("PhotonsMass")->frame(Range(minMassFit-20,maxMassFit+20),Bins(60));
      sigToFit[c]->plotOn(plotPhotonsMassAll);
      SigModelBW.plotOn(plotPhotonsMassAll, LineColor(kGreen), LineStyle(kDashed));
      //  ResAddPdf_draw.plotOn(plotPhotonsMassAll, LineColor(kRed), LineStyle(kDashed));
      ConvolutedRes_CB.plotOn(plotPhotonsMassAll, LineColor(kBlue));

      
      TCanvas* c1 = new TCanvas("c1","PhotonsMass",0,0,800,800);
      c1->cd(1);
      
      plotPhotonsMassAll->Draw();  
      plotPhotonsMassAll->GetYaxis()->SetRangeUser(0.01, max*1.6);
      plotPhotonsMassAll->GetXaxis()->SetRangeUser(210, 290);

      TLatex *lat  = new TLatex(0.55,0.9,TString::Format("Cat: %d", c));  
      lat->SetTextSize(0.038);
      lat->SetTextAlign(11);
      lat->SetTextFont(42); 
      lat->SetNDC();

      TLegend *legmc = new TLegend(0.55, 0.6, 0.87, 0.88, ("Model: "+model).c_str(), "brNDC");
      legmc->AddEntry(plotPhotonsMassAll->getObject(0),"Simulation","LPE");
      legmc->AddEntry(plotPhotonsMassAll->getObject(1),"BW","L");
      //legmc->AddEntry(plotPhotonsMassAll->getObject(2)," CB + CB ","L");
      legmc->AddEntry(plotPhotonsMassAll->getObject(2),"BW #otimes Resolution","L");
      legmc->SetTextSize(0.0206044);
      legmc->SetTextFont(42);
      legmc->SetBorderSize(0);
      legmc->SetFillStyle(0);
      legmc->Draw();
      lat->Draw("same");
      label_cms->Draw("same");
      label_sqrt->Draw("same");
      

      
      int massI(mass);
      c1->SetLogy();
     
      plotPhotonsMassAll->GetXaxis()->SetTitle("m_{#gamma #gamma}[GeV]");
     
      c1->SetLogy(0);
      c1->SaveAs("plots/prelimSignalCBCconvBW"+TString::Format(("_M%d_cat%d_"+model+"_FixW.png").c_str(),massI, c));
      c1->SaveAs("plots/prelimSignalCBCconvBW"+TString::Format(("_M%d_cat%d_"+model+"_FixW.pdf").c_str(),massI, c));
      c1->SaveAs("~/www/plotsNota/prelimSignalCBCconvBW"+TString::Format(("_M%d_cat%d_"+model+"_FixW.png").c_str(),massI, c));
      c1->SaveAs("~/www/plotsNota/prelimSignalCBCconvBW"+TString::Format(("_M%d_cat%d_"+model+"_FixW.pdf").c_str(),massI, c));
      
      c1->SetLogy();
      plotPhotonsMassAll->GetYaxis()->SetRangeUser(0.01,max*10. );
      plotPhotonsMassAll->GetXaxis()->SetRangeUser(210, 290);
      c1->SaveAs("plots/prelimSignalCBCconvBW"+TString::Format(("_M%d_cat%d_LOG_"+model+"_FixW.png").c_str(),massI,c));
      c1->SaveAs("plots/prelimSignalCBCconvBW"+TString::Format(("_M%d_cat%d_LOG_"+model+"_FixW.pdf").c_str(),massI,c));
      c1->SaveAs("~/www/plotsNota/prelimSignalCBCconvBW"+TString::Format(("_M%d_cat%d_LOG_"+model+"_FixW.png").c_str(),massI,c));
      c1->SaveAs("~/www/plotsNota/prelimSignalCBCconvBW"+TString::Format(("_M%d_cat%d_LOG_"+model+"_FixW.pdf").c_str(),massI,c));
   


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
									  *w->var(TString::Format("meanBW_var_cat%d",c)),
									  *w->var(TString::Format("sigmaBW_var_cat%d",c))));
    
        
    SetConstantParams(w->set(TString::Format("ConvolutedPdfParam_cat%d",c)));
    
    //w->Print("V");
    
  }

}


void SigModelFitCBC(RooWorkspace* w, Float_t mass) {

  Int_t ncat = NCAT;
  RooDataSet* sigToFit[NCAT];
  
  
  Float_t MASS(mass);  
  Float_t minMassFit(mass*0.8);
  Float_t maxMassFit(mass*1.2); 
  
  if(mass==150.) minMassFit = MINmass;
  std::cout<<"----------------------------------------------------------------------------------------"<<std::endl;
  //PhotonsMass->setRange("sigrange",minMassFit-20,maxMassFit+20); 
  TPaveText* label_cms = get_labelCMS(0, "2012", true);
  TPaveText* label_sqrt = get_labelSqrt(0);
 
  // Fit Signal 
  for (int c=0; c<ncat; ++c) {
    cout << "---------- Category = " << c << endl;
    sigToFit[c] = (RooDataSet*) w->data(TString::Format("SigWeight_cat%d",c));
   
    RooRealVar* PhotonsMass = w->var("PhotonsMass"); 
    //cb
   
    RooFormulaVar CBpos_mean_draw(TString::Format("CBpos_mean_draw_cat%d",c),"","@0+250",*w->var(TString::Format("ReducedMass_sig_mean_cat%d",c)));
    RooFormulaVar CBpos_sigma(TString::Format("CBpos_sigma_cat%d",c),"","@0*250",*w->var(TString::Format("ReducedMass_sig_sigma_cat%d",c) ));
    RooFormulaVar CBpos_alphaCB(TString::Format("CBpos_alphaCB_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_alphaCBpos_cat%d",c)) );
    RooFormulaVar CBneg_alphaCB(TString::Format("CBneg_alphaCB_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_alphaCBneg_cat%d",c)) );
    RooFormulaVar CBpos_n(TString::Format("CBpos_n_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_Npos_cat%d",c)) );
    RooFormulaVar CBneg_n(TString::Format("CBneg_n_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_Nneg_cat%d",c)) );
    RooFormulaVar CBpos_frac(TString::Format("CBpos_frac_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_frac_cat%d",c)) );
    
    RooCBShape ResCBpos_draw(TString::Format("ResCBpos_draw_cat%d",c),TString::Format("ResCBpos_draw_cat%d",c) , *PhotonsMass, CBpos_mean_draw, CBpos_sigma, CBpos_alphaCB, CBpos_n) ;
    RooCBShape ResCBneg_draw(TString::Format("ResCBneg_draw_cat%d",c),TString::Format("ResCBneg_draw_cat%d",c) , *PhotonsMass, CBpos_mean_draw, CBpos_sigma, CBneg_alphaCB, CBneg_n) ;
   
    //PhotonsMass->setBins(40000, "cache");  
   

    //add CB neg + Gauss
    RooFormulaVar Gauss_frac(TString::Format("Gauss_frac_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_fracGauss_cat%d",c)));    
    RooFormulaVar Gauss_sigma(TString::Format("Gauss_sigma_cat%d",c),"","@0*250",*w->var(TString::Format("ReducedMass_sig_sigmaGauss_cat%d",c)));
    RooGaussian ResGauss_draw(TString::Format("ResGauss_draw_cat%d",c),TString::Format("ResGauss_draw_cat%d",c),*PhotonsMass, CBpos_mean_draw, Gauss_sigma );
    RooAddPdf ResAddGaussPdf_draw(TString::Format("PhotonsMassSig_cat%d",c),TString::Format("PhotonsMassSig_cat%d",c), RooArgList(ResCBneg_draw, ResGauss_draw), Gauss_frac);
    // w->import(ResAddGaussPdf_draw);
    
    //CBC
    RooFormulaVar CBC_mean(TString::Format("CBC_mean_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_sig_mean_cat%d",c)) );
    RooFormulaVar CBC_sigma(TString::Format("CBC_sigma_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_sig_sigma_cat%d",c)) );
    RooFormulaVar CBC_alphaC(TString::Format("CBC_alphaC_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_sig_alphaC_cat%d",c)) );
    RooFormulaVar CBC_alphaCB(TString::Format("CBC_alphaCB_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_sig_alphaCB_cat%d",c)) );
    RooFormulaVar CBC_n(TString::Format("CBC_n_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_sig_n_cat%d",c)) );

    
    RooCBCrujffPdf ResCBCPdf_draw(TString::Format("PhotonsMassSig_cat%d",c),TString::Format("PhotonsMassSig_cat%d",c) , *PhotonsMass, CBC_mean, CBC_sigma, CBC_alphaC, CBC_alphaCB, CBC_n) ; 
    w->import(ResCBCPdf_draw);


    double width=0.1;
    if(width < 2.){ //if i want to plot the fit
      //RooFitResult* fitresults_Gauss = (RooFitResult* ) ResAddGaussPdf_draw.fitTo(*sigToFit[c], RooFit::Save(kTRUE));
      RooFitResult* fitresults_Gauss = (RooFitResult* ) ResCBCPdf_draw.fitTo(*sigToFit[c],Range(minMassFit,maxMassFit),SumW2Error(kTRUE), RooFit::Save(kTRUE));
      std::cout<<TString::Format("******************************** Signal Fit results Gauss mass %f cat %d***********************************", mass, c)<<std::endl;
      fitresults_Gauss->Print("V");
      
      RooPlot* plotOnlyResPdf = w->var("PhotonsMass")->frame(Range(minMassFit-20,maxMassFit+20),Bins(100));
      sigToFit[c]->plotOn(plotOnlyResPdf, LineColor(kRed), LineStyle(kDashed));
      double max = plotOnlyResPdf->GetMaximum();
      
      // Plot to verify everything is ok
      RooPlot* plotPhotonsMassAll = w->var("PhotonsMass")->frame(Range(minMassFit-20,maxMassFit+20),Bins(100));
      sigToFit[c]->plotOn(plotPhotonsMassAll);

      // ResAddGaussPdf_draw.plotOn(plotPhotonsMassAll, LineColor(kRed), LineStyle(kDashed), NormRange("sigrange"));
      //ResAddGaussPdf_draw.plotOn(plotPhotonsMassAll,Components(TString::Format("ResGauss_draw_cat%d",c)), LineColor(kOrange), LineStyle(kDashed));
      //ResAddGaussPdf_draw.plotOn(plotPhotonsMassAll,Components(TString::Format("ResCBneg_draw_cat%d",c)), LineColor(kViolet), LineStyle(kDashed));
      ResCBCPdf_draw.plotOn(plotPhotonsMassAll, LineColor(kRed), LineStyle(kDashed));

      
      TCanvas* c1 = new TCanvas("c1","PhotonsMass",0,0,800,800);
      c1->cd(1);
      
      plotPhotonsMassAll->Draw();  
      
      TLegend *legmc = new TLegend(0.5491457,0.75,0.801457,0.9340659, TString::Format("Category %d",c), "brNDC");
      legmc->AddEntry(plotPhotonsMassAll->getObject(0),"Simulation","LPE");
      //legmc->AddEntry(plotPhotonsMassAll->getObject(2),"CB + Gauss","L");
      //legmc->AddEntry(plotPhotonsMassAll->getObject(3),"Gauss","L");
      //legmc->AddEntry(plotPhotonsMassAll->getObject(4),"CB","L");
      
      legmc->SetTextSize(0.0206044);
      legmc->SetTextFont(42);
      legmc->SetBorderSize(0);
      legmc->SetFillStyle(0);
      legmc->Draw();
      
      label_cms->Draw("same");
      label_sqrt->Draw("same");
      

      
      int massI(mass);
      c1->SetLogy();
      plotPhotonsMassAll->GetYaxis()->SetRangeUser(0.0000000001,max*10. );
      c1->SaveAs("plots/prelimSignalCBGaussCconvBW"+TString::Format("_M%d_cat%d_LOG.png",massI,c));
      c1->SaveAs("plots/prelimSignalCBGaussCconvBW"+TString::Format("_M%d_cat%d_LOG.root",massI,c));
      
      /*  plotPhotonsMassAll->GetYaxis()->SetRangeUser(0.0001, max*1.1 );
      c1->SetLogy(0);
      c1->SaveAs("plots/prelimSignalCBGaussCconvBW"+TString::Format("_M%d_cat%d.png",massI, c));
      c1->SaveAs("plots/prelimSignalCBGaussCconvBW"+TString::Format("_M%d_cat%d.root",massI, c));
      
        c1->SaveAs("plots/prelimSignalCBCconvBW"+TString::Format("_M%d_cat%d.png",massI, c));
	   c1->SaveAs("plots/prelimSignalCBCconvBW"+TString::Format("_M%d_cat%d.root",massI, c));
	   
	   c1->SetLogy();
	   c1->SaveAs("plots/prelimSignalCBCconvBW"+TString::Format("_M%d_cat%d_LOG.png",massI,c));
	   c1->SaveAs("plots/prelimSignalCBCconvBW"+TString::Format("_M%d_cat%d_LOG.root",massI,c));
      */
      
    }
    
    // IMPORTANT: fix all pdf parameters to constant
    
    /*  w->defineSet(TString::Format("ConvolutedPdfGaussParam_cat%d",c),RooArgSet( *w->var(TString::Format("ReducedMass_sig_sigma_cat%d",c)), 
									       *w->var(TString::Format("ReducedMass_sig_alphaCBneg_cat%d",c)),
									       *w->var(TString::Format("ReducedMass_sig_Nneg_cat%d",c)),	   
									       *w->var(TString::Format("ReducedMass_sig_fracGauss_cat%d",c)),  
									       *w->var(TString::Format("ReducedMass_sig_sigmaGauss_cat%d",c)),  
									       *w->var(TString::Format("ReducedMass_sig_mean_cat%d",c)),
									       *w->var(TString::Format("sigmaBW_var_cat%d",c))));
    
    SetConstantParams(w->set(TString::Format("ConvolutedPdfGaussParam_cat%d",c)));
    */

    w->defineSet(TString::Format("CBCParam_cat%d",c),RooArgSet(  *w->var(TString::Format("PhotonsMass_sig_alphaCB_cat%d",c)),
										*w->var(TString::Format("PhotonsMass_sig_n_cat%d",c)),	   
										*w->var(TString::Format("PhotonsMass_sig_alphaC_cat%d",c)),  
										*w->var(TString::Format("PhotonsMass_sig_mean_cat%d",c)),  
										*w->var(TString::Format("PhotonsMass_sig_sigma_cat%d",c))));  
										
    
    SetConstantParams(w->set(TString::Format("CBCParam_cat%d",c)));
    //w->Print("V");
    
  }

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
    /*   TCanvas* ctmp = new TCanvas("ctmp","PhotonsMass Background ",0,0,500,500);
    Int_t nBinsMass(60);
    plotPhotonsMassBkg = PhotonsMass->frame(minMassFit, maxMassFit,nBinsMass);
  
    data->plotOn(plotPhotonsMassBkg,RooFit::Invisible());    
   
    PhotonsMassBkg->plotOn(plotPhotonsMassBkg,LineColor(kBlue),Range(minMassFit,maxMassFit),NormRange("bkg range")); 
   
    double chi2 = plotPhotonsMassBkg->chiSquare(2);
    Int_t ndof = nBinsMass-2;
    std::cout<<"------> "<< ndof<<std::endl;
    double prob = TMath::Prob(chi2*ndof, ndof);
    std::cout<<prob<<std::endl;


    blind=false;
    if( blind ) {
     
      RooDataSet* data_down = (RooDataSet*) data.reduce(*w->var("PhotonsMass"),"PhotonsMass < 327.666 ");
      RooDataSet* data_up = (RooDataSet*) data.reduce(*w->var("PhotonsMass"),"PhotonsMass >650");

      data_up->plotOn(plotPhotonsMassBkg);    
      data_down->plotOn(plotPhotonsMassBkg); 


   
    } else {
      data->plotOn(plotPhotonsMassBkg);    
      } 
      TH1F* h_data = new TH1F("h_data","h_data", 60,minMassFit, maxMassFit);
    TH1F* h_pdf = new TH1F("h_pdf","h_pdf", 60,minMassFit, maxMassFit);
    h_data = (TH1F*) data->createHistogram("PhotonsMass", 60, 0, 0);
    h_pdf =  (TH1F*) PhotonsMassBkg->createHistogram("PhotonsMass",60);
    h_pdf ->Scale(h_data->Integral()/h_pdf->Integral());

    //-------pad 1-------//
    TPad * pad1 = new TPad("pad1", "pad1",0.01,0.13,0.75,1.);  
    
    
    pad1->SetRightMargin(0.1);
    
    pad1->SetLogy();
    pad1->Draw();
    pad1->cd();
    
    pad1->Clear();
    
   
    plotPhotonsMassBkg->GetXaxis()->SetTitle("m_{#gamma #gamma} [GeV]");
    plotPhotonsMassBkg->SetAxisRange(0.001,plotPhotonsMassBkg->GetMaximum()*1.5,"Y");
    plotPhotonsMassBkg->Draw();  

    TLegend *legdata = new TLegend(0.5334677,0.680339,0.8245968,0.8958475, TString::Format("Category %d",c), "brNDC");
    if(mass<950.)legdata = new TLegend(0.2334677,0.300339,0.5645968,0.4958475, TString::Format("Category %d",c), "brNDC");
    legdata->AddEntry(plotPhotonsMassBkg->getObject(2),"Data","LPE");
    legdata->AddEntry(plotPhotonsMassBkg->getObject(1),"Fit Model","L");
  
  

    dobands=false;*/
    //********************************************************************************//
    /*    if (dobands) {

      RooAbsPdf *cpdf; cpdf = PhotonsMassBkg;
      TGraphAsymmErrors onesigma;
      TGraphAsymmErrors twosigma;
      
      RooRealVar *nlim = new RooRealVar(TString::Format("nlim%d",c),"",0.0,0.0,10.0);
      nlim->removeRange();
      
      RooCurve *nomcurve = dynamic_cast<RooCurve*>(plotPhotonsMassBkg->getObject(1));
      
      double el1;
      double eh1;
      double el2;
      double eh2;
  
      int j = 0;
      for (int i=1; i<(plotPhotonsMassBkg->GetXaxis()->GetNbins()+1); ++i) {
	double lowedge = plotPhotonsMassBkg->GetXaxis()->GetBinLowEdge(i);
	double upedge  = plotPhotonsMassBkg->GetXaxis()->GetBinUpEdge(i);
	double center  = plotPhotonsMassBkg->GetXaxis()->GetBinCenter(i);
	
	double nombkg = nomcurve->interpolate(center);


	nlim->setVal(nombkg);
	PhotonsMass->setRange("errRange",lowedge,upedge);
	RooAbsPdf *epdf = 0;
	epdf = new RooExtendPdf("epdf","",*cpdf,*nlim,"errRange");
	
	RooAbsReal *nll = epdf->createNLL(*(data),Extended());
	RooMinimizer minim(*nll);
	minim.setStrategy(0);
	double clone = 1.0 - 2.0*RooStats::SignificanceToPValue(1.0);
	double cltwo = 1.0 - 2.0*RooStats::SignificanceToPValue(2.0);
	
	minim.migrad();
	minim.minos(*nlim);
	printf("errlo = %5f, errhi = %5f\n",nlim->getErrorLo(),nlim->getErrorHi());


	el1 = nlim->getErrorLo();
	eh1= nlim->getErrorHi();
	//	std::cout<<"-----------------------------------------------------------------> "<< minim.minos(*nlim)<<std::endl;
	minim.setErrorLevel(0.5*pow(ROOT::Math::normal_quantile(1-0.5*(1-cltwo),1.0), 2)); // the 0.5 is because qmu is -2*NLL
	minim.migrad();
	minim.minos(*nlim);
	el2 = nlim->getErrorLo();
	eh2= nlim->getErrorHi();


	delete nll;
	delete epdf;
	if( el1 != 0. && eh1 != 0. && el2 != 0. && eh2 != 0. &&  el1 != 1. && eh1 != 1.  && el2 != 1.  && eh2 != 1. ) {
	onesigma.SetPoint(j,center,nombkg);
	twosigma.SetPoint(j,center,nombkg);
	onesigma.SetPointError(j,0.,0.,-el1,eh1);
	twosigma.SetPointError(j,0.,0.,-el2,eh2);
	j++;
	}

	}

      PhotonsMass->setRange("errRange",minMassFit,maxMassFit);
      
      twosigma.SetLineColor(kGreen);
      twosigma.SetFillColor(kGreen);
      twosigma.SetMarkerColor(kGreen);
      twosigma.Draw("C3 SAME");
  
      onesigma.SetLineColor(kYellow);
      onesigma.SetFillColor(kYellow);
      onesigma.SetMarkerColor(kYellow);
      onesigma.Draw("C3 SAME");
   
      legdata->AddEntry(&onesigma, "#pm 1 #sigma", "F" );
      legdata->AddEntry(&twosigma, "#pm 2 #sigma","F" );
      plotPhotonsMassBkg->Draw("SAME"); 
     
      }
    */
    /* legdata->SetTextSize(0.035);
    legdata->SetTextFont(42);
    legdata->SetBorderSize(0);
    legdata->SetFillStyle(0);
    legdata->Draw("same");


    TPaveText* label_cms = get_labelCMS(0, "2012", false);
    TPaveText* label_sqrt = get_labelSqrt(0);
    label_cms->Draw("same");
    label_sqrt->Draw("same");

    //write down the chi2 of the fit on the
      
    TPaveText* label_chi2= new TPaveText(0.5744355,0.5750847,0.80583871,0.65822034,"brNDC");
    if(mass<950.)label_chi2 = new TPaveText(0.4344355,0.1550847,0.5983871,0.2622034,"brNDC");
    label_chi2->SetFillColor(kWhite);
    label_chi2->SetTextSize(0.035);
    label_chi2->SetTextFont(42);
    label_chi2->SetTextAlign(31); // align right
    label_chi2->AddText(TString::Format("Fit chi square/dof = %.3f", chi2));
    label_chi2->AddText(TString::Format("Chi square Prob = %.3f", prob));
    label_chi2->Draw("same");
    

    ctmp->cd();
    //-------pad 2------//
    TPad * pad2 = new TPad("pad2", "pad2",0.01,0.001,0.75,0.2);
    pad2->SetGrid();
    
    //pad2->SetTopMargin(0.01);
    pad2->SetBottomMargin(0.4);
    pad2->SetRightMargin(0.1);
    pad2->Draw();
    pad2->cd();
    h_data->Sumw2();
    h_pdf->Sumw2();
    TH1F* h_data1 = h_data->Clone();
    h_data->Add(h_pdf, -1);
    for(int i = 0; i< h_data->GetNbinsX();i++) h_data->SetBinError(i,h_data1->GetBinError(i));
    h_data->GetYaxis()->SetRangeUser(-50., 50.);
    h_data->GetYaxis()->SetNdivisions(505);
    h_data->SetMarkerSize(0.4);
    h_data->GetXaxis()->SetTitle("m_{#gamma #gamma}");
    h_data->GetXaxis()->SetTitleSize(0.2);
    h_data->GetXaxis()->SetLabelSize(0.16);
    h_data->GetYaxis()->SetLabelSize(0.16);
    h_data->GetYaxis()->SetTitleSize(0.15);
    h_data->GetYaxis()->SetTitle("");
    h_data->GetYaxis()->SetTitleOffset(0.65);
    h_data->GetXaxis()->SetTitleOffset(0.8);
    //    for(int i=0; i<h_data->GetNbinsX();i++) h_data->SetBinError(i, h_data->GetBinError(i)/h_pdf->GetBinContent(i));
    h_data->Draw("pe");
    TH1F* h_dataCopy = h_data->Clone();
    for (int i = 0; i< h_dataCopy->GetNbinsX();i++) if (h_dataCopy->GetBinContent(i)==0) h_dataCopy->SetBinContent(i, 1);
    h_dataCopy->Add(h_dataCopy, -1);
    h_dataCopy->SetLineColor(kRed);
    h_dataCopy->Draw("histsame");
 
    int massI(mass);
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPPAR_M%d.png",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPPAR_M%d.pdf",c,massI));

    ctmp->SetLogy();
    plotPhotonsMassBkg->SetAxisRange(1.3,plotPhotonsMassBkg->GetMaximum()*1.5,"Y");
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPPAR_LOG_M%d.png",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPPAR_LOG_M%d.pdf",c,massI));
    ctmp->SaveAs("~/www/plotsNota/prelimBkg"+TString::Format("_cat%d_EXPPAR_LOG_M%d.png",c,massI));
    ctmp->SaveAs("~/www/plotsNota/prelimBkg"+TString::Format("_cat%d_EXPPAR_LOG_M%d.pdf",c,massI));
    */
    
  RooFitResult* r;

  return r;
}







RooFitResult* BkgModelFitDiJetFunc(RooWorkspace* w, Bool_t dobands, Float_t mass, bool blind) {

  Int_t ncat = NCAT;
  std::cout<<"isBlind: "<<blind<<std::endl;
  // retrieve pdfs and datasets from workspace to fit with pdf models
  RooDataSet* data[NCAT];
 
  RooFitResult* fitresult[NCAT];;
  RooPlot* plotPhotonsMassBkg[NCAT];

  // dobands and dosignal
  RooDataSet* signal[NCAT];

  RooAbsPdf* PhotonsMassSig[NCAT];
  
  Float_t minMassFit, maxMassFit;
  minMassFit = MINmass;
  maxMassFit = MAXmass;
   
  // Fit data with background pdf for data limit
  RooRealVar* PhotonsMass = w->var("PhotonsMass");  
  PhotonsMass->setUnit("GeV");
 
  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
  
  for (int c = 0; c < ncat; ++c) {
    data[c] = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    
  
    // fit a la dijets
    RooFormulaVar *p0mod = new RooFormulaVar(TString::Format("par0DiJet_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_norm_cat%d",c)));
    RooFormulaVar *p1mod = new RooFormulaVar(TString::Format("par1DiJet_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_slope1_cat%d",c)));
    RooFormulaVar *p2mod = new RooFormulaVar(TString::Format("par2DiJet_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_slope2_cat%d",c)));
    RooFormulaVar *p3mod = new RooFormulaVar(TString::Format("par3DiJet_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_slope3_cat%d",c)));
    PhotonsMass->setRange("bkg range", MINmass,MAXmass);
    RooFormulaVar *x     = new RooFormulaVar(TString::Format("xDiJet_cat%d",c),"","@0/8000.",*w->var("PhotonsMass"));

   
    RooAbsPdf* PhotonsMassBkgTmp0 = new RooGenericPdf(TString::Format("PhotonsMassBkg_cat%d",c), "pow(1-@0, @2)/pow(@0, @1+@3*log(@0))", RooArgList(*x, *p1mod, *p2mod,*p3mod));
   

    fitresult[c] = PhotonsMassBkgTmp0->fitTo(*data[c], Range(minMassFit,maxMassFit),RooFit::FitOptions("MHTER"), SumW2Error(kTRUE), Save(kTRUE));
    w->import(*PhotonsMassBkgTmp0);
   
    std::cout<<TString::Format("******************************** Background Fit results mass %f cat %d***********************************", mass, c)<<std::endl;
    fitresult[c]->Print("V");
   

    //************************************************
    // Plot PhotonsMass background fit results per categories 
    TCanvas* ctmp = new TCanvas("ctmp","PhotonsMass Background Categories",0,0,500,500);
    Int_t nBinsMass(60);
    plotPhotonsMassBkg[c] = PhotonsMass->frame(minMassFit, maxMassFit,nBinsMass);
  
    data[c]->plotOn(plotPhotonsMassBkg[c],RooFit::Invisible());    
   
    PhotonsMassBkgTmp0->plotOn(plotPhotonsMassBkg[c],LineColor(kBlue),Range(minMassFit,maxMassFit),NormRange("bkg range")); 
    double chi2 = plotPhotonsMassBkg[c]->chiSquare(3);
    Int_t ndof = nBinsMass-3;
    std::cout<<"------> "<< ndof<<std::endl;
    double prob = TMath::Prob(chi2*ndof, ndof);
    std::cout<<prob<<std::endl;
    blind = true;
    if( blind ) {
     
      RooDataSet* data_down = (RooDataSet*) data[c].reduce(*w->var("PhotonsMass"),"PhotonsMass < 178.");
      RooDataSet* data_up = (RooDataSet*) data[c].reduce(*w->var("PhotonsMass"),"PhotonsMass >402");

      data_up->plotOn(plotPhotonsMassBkg[c]);    
      data_down->plotOn(plotPhotonsMassBkg[c]); 


   
    } else {
      data[c]->plotOn(plotPhotonsMassBkg[c]);    
      } 
       
   
    plotPhotonsMassBkg[c]->GetXaxis()->SetTitle("m_{#gamma #gamma}[GeV]");
    plotPhotonsMassBkg[c]->SetAxisRange(0.001,plotPhotonsMassBkg[c]->GetMaximum()*1.5,"Y");
    plotPhotonsMassBkg[c]->Draw();  

    TLegend *legdata = new TLegend(0.3790323,0.7775424,0.6290323,0.9279661, TString::Format("Category %d",c), "brNDC");
    legdata->AddEntry(plotPhotonsMassBkg[c]->getObject(2),"Data","LPE");
    legdata->AddEntry(plotPhotonsMassBkg[c]->getObject(1),"Parametric Model: DiJet","L");
    legdata->SetTextSize(0.035);
    legdata->SetTextFont(42);
    // legdata->SetTextAlign(31);
    legdata->SetBorderSize(0);
    legdata->SetFillStyle(0);
    legdata->Draw("same");

    TPaveText* label_cms = get_labelCMS(0, "2012", false);
    TPaveText* label_sqrt = get_labelSqrt(0);
    label_cms->Draw("same");
    label_sqrt->Draw("same");

    //write down the chi2 of the fit on the
 
    TPaveText* label_chi2 = new TPaveText(0.5524194,0.6419492,0.796371,0.7690678, "brNDC");
    label_chi2->SetFillColor(kWhite);
    label_chi2->SetTextSize(0.035);
    label_chi2->SetTextFont(42);
    label_chi2->SetTextAlign(31); // align right
    label_chi2->AddText(TString::Format("Fit chi square/dof = %.3f", chi2));
    label_chi2->AddText(TString::Format("Chi square Prob = %.3f", prob));
    label_chi2->Draw("same");

    

    //********************************************************************************//
    if (dobands) {

      RooAbsPdf *cpdf; cpdf = PhotonsMassBkgTmp0;
      TGraphAsymmErrors *onesigma = new TGraphAsymmErrors();
      TGraphAsymmErrors *twosigma = new TGraphAsymmErrors();
      
      RooRealVar *nlim = new RooRealVar(TString::Format("nlim%d",c),"",0.0,0.0,10.0);
      nlim->removeRange();
      
      RooCurve *nomcurve = dynamic_cast<RooCurve*>(plotPhotonsMassBkg[c]->getObject(1));
      
      for (int i=1; i<(plotPhotonsMassBkg[c]->GetXaxis()->GetNbins()+1); ++i) {
	double lowedge = plotPhotonsMassBkg[c]->GetXaxis()->GetBinLowEdge(i);
	double upedge  = plotPhotonsMassBkg[c]->GetXaxis()->GetBinUpEdge(i);
	double center  = plotPhotonsMassBkg[c]->GetXaxis()->GetBinCenter(i);
	
	double nombkg = nomcurve->interpolate(center);
	nlim->setVal(nombkg);
	PhotonsMass->setRange("errRange",lowedge,upedge);
	RooAbsPdf *epdf = 0;
	epdf = new RooExtendPdf("epdf","",*cpdf,*nlim,"errRange");
	
	RooAbsReal *nll = epdf->createNLL(*(data[c]),Extended());
	RooMinimizer minim(*nll);
	minim.setStrategy(0);
	double clone = 1.0 - 2.0*RooStats::SignificanceToPValue(1.0);
	double cltwo = 1.0 - 2.0*RooStats::SignificanceToPValue(2.0);
	
	minim.migrad();
	minim.minos(*nlim);
	// printf("errlo = %5f, errhi = %5f\n",nlim->getErrorLo(),nlim->getErrorHi());
	
	onesigma->SetPoint(i-1,center,nombkg);
	onesigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
	
	minim.setErrorLevel(0.5*pow(ROOT::Math::normal_quantile(1-0.5*(1-cltwo),1.0), 2)); // the 0.5 is because qmu is -2*NLL
	// eventually if cl = 0.95 this is the usual 1.92!      
	
	minim.migrad();
	minim.minos(*nlim);
	
	twosigma->SetPoint(i-1,center,nombkg);
	twosigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
	
	delete nll;
	delete epdf;
      }

      PhotonsMass->setRange("errRange",minMassFit,maxMassFit);
      
      twosigma->SetLineColor(kGreen);
      twosigma->SetFillColor(kGreen);
      twosigma->SetMarkerColor(kGreen);
      twosigma->Draw("L3 SAME");
      
      onesigma->SetLineColor(kYellow);
      onesigma->SetFillColor(kYellow);
      onesigma->SetMarkerColor(kYellow);
      onesigma->Draw("L3 SAME");
      
      plotPhotonsMassBkg[c]->Draw("SAME"); 
    }

    int massI(mass);
    ctmp->SaveAs("plots/prelimBkg"+TString::Format("_cat%d_DIJET_M%d.png",c,massI));
    ctmp->SaveAs("plots/prelimBkg"+TString::Format("_cat%d_DIJET_M%d.pdf",c,massI));
    ctmp->SaveAs("plots/prelimBkg"+TString::Format("_cat%d_DIJET_M%d.root",c,massI));

    ctmp->SetLogy();
    plotPhotonsMassBkg[c]->SetAxisRange(1.3,plotPhotonsMassBkg[c]->GetMaximum()*1.5,"Y");
    ctmp->SaveAs("plots/prelimBkg"+TString::Format("_cat%d_DIJET_LOG_M%d.png",c,massI));
    ctmp->SaveAs("plots/prelimBkg"+TString::Format("_cat%d_DIJET_LOG_M%d.pdf",c,massI));
    ctmp->SaveAs("plots/prelimBkg"+TString::Format("_cat%d_DIJET_LOG_M%d.root",c,massI));

  }



  return fitresult;
}




RooFitResult* BkgModelFitDiJetEXPFunc(RooWorkspace* w, Bool_t dobands, Float_t mass,Int_t c, bool blind) {

  Int_t ncat = NCAT;
  std::cout<<"isBlind: "<<blind<<std::endl;
  // retrieve pdfs and datasets from workspace to fit with pdf models
  RooDataSet* data;
 
  RooFitResult* fitresult;;
  
  RooPlot* plotPhotonsMassBkg;

  // dobands and dosignal
  RooDataSet* signal;

  RooAbsPdf* PhotonsMassSig;
  
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
    
  
    // fit a la dijets
    
    RooFormulaVar *p1mod = new RooFormulaVar(TString::Format("par1DiJetEXP_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_slope1_3_cat%d",c)));
    RooFormulaVar *p2mod = new RooFormulaVar(TString::Format("par2DiJetEXP_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_slope2_3_cat%d",c)));
    RooFormulaVar *p3mod = new RooFormulaVar(TString::Format("par3DiJetEXP_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_slope3_3_cat%d",c)));
    RooFormulaVar *exp1 = new RooFormulaVar(TString::Format("expDiJetEXP_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_exp1DiJetEXP_cat%d",c)));
    RooFormulaVar *pFrac1 = new RooFormulaVar(TString::Format("fracDiJetEXP_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_fracDiJetEXP_cat%d",c)));
   
   
    RooGenericPdf* PhotonsMassBkgTmp0DiJet = new RooGenericPdf(TString::Format("PhotonsMassBkg_DIJETE_truth_cat%d",c), "pow(1-@0/8000., @2)/pow(@0/8000., @1+@3*log(@0/8000.))", RooArgList(*w->var("PhotonsMass"), *p1mod, *p2mod,*p3mod));
  
    RooExponential* PhotonsMassBkgTmp0Exp = new RooExponential(TString::Format("PhotonsMassBkg_EXP_truth_cat%d",c),"", *w->var("PhotonsMass"),  *exp1);
    
   
    RooAddPdf* PhotonsMassBkgTmpAdd = new RooAddPdf(TString::Format("PhotonsMassBkg_cat%d",c),TString::Format("PhotonsMassBkg_cat%d",c) , RooArgList(*PhotonsMassBkgTmp0DiJet, *PhotonsMassBkgTmp0Exp), RooArgList(*pFrac1));
    
    fitresult = PhotonsMassBkgTmpAdd->fitTo(*data,RooFit::FitOptions("MHTR"), Save(kTRUE));//RooFit::FitOptions("MHTER"), Range(minMassFit,maxMassFit),    
    w->import(*PhotonsMassBkgTmpAdd);
  
    std::cout<<TString::Format("******************************** Background DiJetEXP Fit results mass %f cat %d ***********************************", mass, c)<<std::endl;
    fitresult->Print("V");
   

    //************************************************
    // Plot PhotonsMass background fit results per categories 
    TCanvas* ctmp = new TCanvas("ctmp","PhotonsMass Background Categories",0,0,500,500);
    Int_t nBinsMass(60);
    plotPhotonsMassBkg = PhotonsMass->frame(minMassFit, maxMassFit,nBinsMass);
  
    data->plotOn(plotPhotonsMassBkg,RooFit::Invisible());    
   
    PhotonsMassBkgTmpAdd->plotOn(plotPhotonsMassBkg,LineColor(kBlue)); 
    PhotonsMassBkgTmpAdd->plotOn(plotPhotonsMassBkg,Components(TString::Format("PhotonsMassBkg_EXP_truth_cat%d",c)),LineColor(kViolet),LineStyle(kDashed)); 
    PhotonsMassBkgTmpAdd->plotOn(plotPhotonsMassBkg,Components(TString::Format("PhotonsMassBkg_DIJETE_truth_cat%d",c)),LineColor(kOrange),LineStyle(kDashed));   

    double chi2 = plotPhotonsMassBkg->chiSquare(3);
    Int_t ndof = nBinsMass-5;
    std::cout<<"------> "<< ndof<<std::endl;
    double prob = TMath::Prob(chi2*ndof, ndof);
    std::cout<<prob<<std::endl;
   
    if( blind ) {
 
      RooDataSet* data_down = (RooDataSet*) data.reduce(*w->var("PhotonsMass"),"PhotonsMass < 173.5");
      RooDataSet* data_up = (RooDataSet*) data.reduce(*w->var("PhotonsMass")," PhotonsMass >850");
      TH1F* h_up= new TH1F("h_up", "h_up",nBinsMass, 130, 1000);
      h_up->Sumw2();
      data_up->fillHistogram(h_up, RooArgList(*PhotonsMass));
      TH1F* h_down= new TH1F("h_down", "h_down",nBinsMass, 130, 1000);
      h_down->Sumw2();
      data_down->fillHistogram(h_down, RooArgList(*PhotonsMass));
   
    } else {
      data->plotOn(plotPhotonsMassBkg);    
      } 
       
   
    plotPhotonsMassBkg->GetXaxis()->SetTitle("m_{#gamma #gamma} [GeV]");
    plotPhotonsMassBkg->SetAxisRange(0.1,10000,"Y");
    plotPhotonsMassBkg->Draw();  
   if( blind ) {
       h_up->Draw("sameP");
       h_down->Draw("sameP");
     }
  
    TLegend *legdata = new TLegend(0.3790323,0.7775424,0.6290323,0.9279661, TString::Format("Category %d",c), "brNDC");
    legdata->AddEntry(plotPhotonsMassBkg->getObject(4),"Data","LPE");
    legdata->AddEntry(plotPhotonsMassBkg->getObject(1),"Parametric Model: DiJetEXP","L");  
    legdata->AddEntry(plotPhotonsMassBkg->getObject(2),"Parametric Model: EXP","L");
    legdata->AddEntry(plotPhotonsMassBkg->getObject(3),"Parametric Model: DiJet","L");
  
    legdata->SetTextSize(0.035);
    legdata->SetTextFont(42);
    // legdata->SetTextAlign(31);
    legdata->SetBorderSize(0);
    legdata->SetFillStyle(0);
    legdata->Draw("same");

    TPaveText* label_cms = get_labelCMS(0, "2012", false);
    TPaveText* label_sqrt = get_labelSqrt(0);
    label_cms->Draw("same");
    label_sqrt->Draw("same");

    //write down the chi2 of the fit on the
 
    TPaveText* label_chi2 = new TPaveText(0.5524194,0.6419492,0.796371,0.7690678, "brNDC");
    label_chi2->SetFillColor(kWhite);
    label_chi2->SetTextSize(0.035);
    label_chi2->SetTextFont(42);
    label_chi2->SetTextAlign(31); // align right
    label_chi2->AddText(TString::Format("Fit chi square/dof = %.3f", chi2));
    label_chi2->AddText(TString::Format("Chi square Prob = %.3f", prob));
    label_chi2->Draw("same");

    

    //********************************************************************************//
    if (dobands) {

      RooAbsPdf *cpdf; cpdf = PhotonsMassBkgTmp0;
      TGraphAsymmErrors *onesigma = new TGraphAsymmErrors();
      TGraphAsymmErrors *twosigma = new TGraphAsymmErrors();
      
      RooRealVar *nlim = new RooRealVar(TString::Format("nlim%d",c),"",0.0,0.0,10.0);
      nlim->removeRange();
      
      RooCurve *nomcurve = dynamic_cast<RooCurve*>(plotPhotonsMassBkg[c]->getObject(1));
      
      for (int i=1; i<(plotPhotonsMassBkg->GetXaxis()->GetNbins()+1); ++i) {
	double lowedge = plotPhotonsMassBkg->GetXaxis()->GetBinLowEdge(i);
	double upedge  = plotPhotonsMassBkg->GetXaxis()->GetBinUpEdge(i);
	double center  = plotPhotonsMassBkg->GetXaxis()->GetBinCenter(i);
	
	double nombkg = nomcurve->interpolate(center);
	nlim->setVal(nombkg);
	PhotonsMass->setRange("errRange",lowedge,upedge);
	RooAbsPdf *epdf = 0;
	epdf = new RooExtendPdf("epdf","",*cpdf,*nlim,"errRange");
	
	RooAbsReal *nll = epdf->createNLL(*(data),Extended());
	RooMinimizer minim(*nll);
	minim.setStrategy(0);
	double clone = 1.0 - 2.0*RooStats::SignificanceToPValue(1.0);
	double cltwo = 1.0 - 2.0*RooStats::SignificanceToPValue(2.0);
	
	minim.migrad();
	minim.minos(*nlim);
	// printf("errlo = %5f, errhi = %5f\n",nlim->getErrorLo(),nlim->getErrorHi());
	
	onesigma->SetPoint(i-1,center,nombkg);
	onesigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
	
	minim.setErrorLevel(0.5*pow(ROOT::Math::normal_quantile(1-0.5*(1-cltwo),1.0), 2)); // the 0.5 is because qmu is -2*NLL
	// eventually if cl = 0.95 this is the usual 1.92!      
	
	minim.migrad();
	minim.minos(*nlim);
	
	twosigma->SetPoint(i-1,center,nombkg);
	twosigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
	
	delete nll;
	delete epdf;
      }

      PhotonsMass->setRange("errRange",minMassFit,maxMassFit);
      
      twosigma->SetLineColor(kGreen);
      twosigma->SetFillColor(kGreen);
      twosigma->SetMarkerColor(kGreen);
      twosigma->Draw("L3 SAME");
      
      onesigma->SetLineColor(kYellow);
      onesigma->SetFillColor(kYellow);
      onesigma->SetMarkerColor(kYellow);
      onesigma->Draw("L3 SAME");
      
      plotPhotonsMassBkg->Draw("SAME"); 
    }

    int massI(mass);
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETEXP_M%d.png",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETEXP_M%d.pdf",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETEXP_M%d.root",c,massI));

    ctmp->SetLogy();
    //  plotPhotonsMassBkg->SetAxisRange(1.3,plotPhotonsMassBkg[c]->GetMaximum()*1.5,"Y");
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETEXP_LOG_M%d.png",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETEXP_LOG_M%d.pdf",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETEXP_LOG_M%d.root",c,massI));

 

  return fitresult;
}






RooFitResult* BkgModelFitDiJetEXPOLFunc(RooWorkspace* w, Bool_t dobands, Float_t mass, Int_t c, bool blind) {

  Int_t ncat = NCAT;
 
  // retrieve pdfs and datasets from workspace to fit with pdf models
  RooDataSet* data;
 
  RooFitResult* fitresult;;

  RooPlot* plotPhotonsMassBkg;

  // dobands and dosignal
  RooDataSet* signal;

  RooAbsPdf* PhotonsMassSig;
  
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
    
  
    // fit a la dijets
    
    RooFormulaVar *p1mod = new RooFormulaVar(TString::Format("par1DiJetEXPOL_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_slope1_4_cat%d",c)));
    RooFormulaVar *p2mod = new RooFormulaVar(TString::Format("par2DiJetEXPOL_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_slope2_4_cat%d",c)));
    RooFormulaVar *p3mod = new RooFormulaVar(TString::Format("par3DiJetEXPOL_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_slope3_4_cat%d",c)));
    RooFormulaVar *expol1 = new RooFormulaVar(TString::Format("expol1DiJetEXPOL_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_expol1_4_cat%d",c)));
    RooFormulaVar *expol2 = new RooFormulaVar(TString::Format("expol2DiJetEXPOL_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_expol2_4_cat%d",c)));
    RooFormulaVar *pFrac1 = new RooFormulaVar(TString::Format("fracDiJetEXPOL_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_fracDiJetEXPOL_cat%d",c)));
   
   
    RooGenericPdf* PhotonsMassBkgTmp0DiJet = new RooGenericPdf(TString::Format("PhotonsMassBkg_DiJetEx_cat%d",c), "pow(1-@0/8000., @2)/pow(@0/8000., @1+@3*log(@0/8000.))", RooArgList(*w->var("PhotonsMass"), *p1mod, *p2mod,*p3mod));
    RooGenericPdf* PhotonsMassBkgTmp0Expol = new RooGenericPdf(TString::Format("PhotonsMassBkg_ExpolDiJ_cat%d",c), "exp(-@0/(@1+@2*@0))", RooArgList(*w->var("PhotonsMass"), *expol1, *expol2));
    
   
    RooAddPdf* PhotonsMassBkgTmpAdd = new RooAddPdf(TString::Format("PhotonsMassBkg_cat%d",c),TString::Format("PhotonsMassBkg_cat%d",c) , RooArgList(*PhotonsMassBkgTmp0DiJet, *PhotonsMassBkgTmp0Expol), RooArgList(*pFrac1));
    
    fitresult = PhotonsMassBkgTmpAdd->fitTo(*data,RooFit::FitOptions("MHTR"), Save(kTRUE));//RooFit::FitOptions("MHTER"), Range(minMassFit,maxMassFit),    
    w->import(*PhotonsMassBkgTmpAdd);
  
    std::cout<<TString::Format("******************************** Background DiJetEXPOL Fit results mass %f cat %d ***********************************", mass, c)<<std::endl;
    fitresult->Print("V");
   

    //************************************************
    // Plot PhotonsMass background fit results per categories 
    TCanvas* ctmp = new TCanvas("ctmp","PhotonsMass Background Categories",0,0,500,500);
    Int_t nBinsMass(60);
    plotPhotonsMassBkg = PhotonsMass->frame(minMassFit, maxMassFit,nBinsMass);
  
    data->plotOn(plotPhotonsMassBkg,RooFit::Invisible());    
   
    PhotonsMassBkgTmpAdd->plotOn(plotPhotonsMassBkg,LineColor(kBlue)); 
    PhotonsMassBkgTmpAdd->plotOn(plotPhotonsMassBkg,Components(TString::Format("PhotonsMassBkg_ExpolDiJ_cat%d",c)),LineColor(kViolet),LineStyle(kDashed)); 
    PhotonsMassBkgTmpAdd->plotOn(plotPhotonsMassBkg,Components(TString::Format("PhotonsMassBkg_DiJetEx_cat%d",c)),LineColor(kOrange),LineStyle(kDashed));   

    double chi2 = plotPhotonsMassBkg->chiSquare(3);
    Int_t ndof = nBinsMass-6;
    std::cout<<"------> "<< ndof<<std::endl;
    double prob = TMath::Prob(chi2*ndof, ndof);
    std::cout<<prob<<std::endl;
    blind = false;
    if( blind ) {
     
      RooDataSet* data_down = (RooDataSet*) data.reduce(*w->var("PhotonsMass"),"PhotonsMass < 178.");
      RooDataSet* data_up = (RooDataSet*) data.reduce(*w->var("PhotonsMass"),"PhotonsMass >402");

      data_up->plotOn(plotPhotonsMassBkg);    
      data_down->plotOn(plotPhotonsMassBkg); 


   
    } else {
      data->plotOn(plotPhotonsMassBkg);    
      } 
       
   
    plotPhotonsMassBkg->GetXaxis()->SetTitle("m_{#gamma #gamma} [GeV]");
    plotPhotonsMassBkg->SetAxisRange(0.001,plotPhotonsMassBkg->GetMaximum()*1.5,"Y");
    plotPhotonsMassBkg->Draw();  

    TLegend *legdata = new TLegend(0.3790323,0.7775424,0.6290323,0.9279661, TString::Format("Category %d",c), "brNDC");
    legdata->AddEntry(plotPhotonsMassBkg->getObject(4),"Data","LPE");
    legdata->AddEntry(plotPhotonsMassBkg->getObject(1),"Parametric Model: DiJetEXPOL","L");  
    legdata->AddEntry(plotPhotonsMassBkg->getObject(2),"Parametric Model: Expol","L");
    legdata->AddEntry(plotPhotonsMassBkg->getObject(3),"Parametric Model: DiJet","L");
  
    legdata->SetTextSize(0.035);
    legdata->SetTextFont(42);
    // legdata->SetTextAlign(31);
    legdata->SetBorderSize(0);
    legdata->SetFillStyle(0);
    legdata->Draw("same");

    TPaveText* label_cms = get_labelCMS(0, "2012", false);
    TPaveText* label_sqrt = get_labelSqrt(0);
    label_cms->Draw("same");
    label_sqrt->Draw("same");

    //write down the chi2 of the fit on the
 
    TPaveText* label_chi2 = new TPaveText(0.5524194,0.6419492,0.796371,0.7690678, "brNDC");
    label_chi2->SetFillColor(kWhite);
    label_chi2->SetTextSize(0.035);
    label_chi2->SetTextFont(42);
    label_chi2->SetTextAlign(31); // align right
    label_chi2->AddText(TString::Format("Fit chi square/dof = %.3f", chi2));
    label_chi2->AddText(TString::Format("Chi square Prob = %.3f", prob));
    label_chi2->Draw("same");

    

    //********************************************************************************//
    if (dobands) {

      RooAbsPdf *cpdf; cpdf = PhotonsMassBkgTmp0;
      TGraphAsymmErrors *onesigma = new TGraphAsymmErrors();
      TGraphAsymmErrors *twosigma = new TGraphAsymmErrors();
      
      RooRealVar *nlim = new RooRealVar(TString::Format("nlim%d",c),"",0.0,0.0,10.0);
      nlim->removeRange();
      
      RooCurve *nomcurve = dynamic_cast<RooCurve*>(plotPhotonsMassBkg->getObject(1));
      
      for (int i=1; i<(plotPhotonsMassBkg->GetXaxis()->GetNbins()+1); ++i) {
	double lowedge = plotPhotonsMassBkg->GetXaxis()->GetBinLowEdge(i);
	double upedge  = plotPhotonsMassBkg->GetXaxis()->GetBinUpEdge(i);
	double center  = plotPhotonsMassBkg->GetXaxis()->GetBinCenter(i);
	
	double nombkg = nomcurve->interpolate(center);
	nlim->setVal(nombkg);
	PhotonsMass->setRange("errRange",lowedge,upedge);
	RooAbsPdf *epdf = 0;
	epdf = new RooExtendPdf("epdf","",*cpdf,*nlim,"errRange");
	
	RooAbsReal *nll = epdf->createNLL(*(data),Extended());
	RooMinimizer minim(*nll);
	minim.setStrategy(0);
	double clone = 1.0 - 2.0*RooStats::SignificanceToPValue(1.0);
	double cltwo = 1.0 - 2.0*RooStats::SignificanceToPValue(2.0);
	
	minim.migrad();
	minim.minos(*nlim);
	// printf("errlo = %5f, errhi = %5f\n",nlim->getErrorLo(),nlim->getErrorHi());
	
	onesigma->SetPoint(i-1,center,nombkg);
	onesigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
	
	minim.setErrorLevel(0.5*pow(ROOT::Math::normal_quantile(1-0.5*(1-cltwo),1.0), 2)); // the 0.5 is because qmu is -2*NLL
	// eventually if cl = 0.95 this is the usual 1.92!      
	
	minim.migrad();
	minim.minos(*nlim);
	
	twosigma->SetPoint(i-1,center,nombkg);
	twosigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
	
	delete nll;
	delete epdf;
      }

      PhotonsMass->setRange("errRange",minMassFit,maxMassFit);
      
      twosigma->SetLineColor(kGreen);
      twosigma->SetFillColor(kGreen);
      twosigma->SetMarkerColor(kGreen);
      twosigma->Draw("L3 SAME");
      
      onesigma->SetLineColor(kYellow);
      onesigma->SetFillColor(kYellow);
      onesigma->SetMarkerColor(kYellow);
      onesigma->Draw("L3 SAME");
      
      plotPhotonsMassBkg->Draw("SAME"); 
    }

    int massI(mass);
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETEXPOL_M%d.png",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETEXPOL_M%d.pdf",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETEXPOL_M%d.root",c,massI));

    ctmp->SetLogy();
    plotPhotonsMassBkg->SetAxisRange(1.3,plotPhotonsMassBkg->GetMaximum()*1.5,"Y");
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETEXPOL_LOG_M%d.png",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETEXPOL_LOG_M%d.pdf",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETEXPOL_LOG_M%d.root",c,massI));

  



  return fitresult;
}





RooFitResult* BkgModelFitExpolFunc(RooWorkspace* w, Bool_t dobands, Float_t mass, bool blind) {

  Int_t ncat = NCAT;
  std::cout<<"isBlind: "<<blind<<std::endl;
  // retrieve pdfs and datasets from workspace to fit with pdf models
  RooDataSet* data[NCAT];
 
  RooFitResult* fitresult[NCAT];;
  RooPlot* plotPhotonsMassBkg[NCAT];

  // dobands and dosignal
  RooDataSet* signal[NCAT];

  RooAbsPdf* PhotonsMassSig[NCAT];
  
  Float_t minMassFit, maxMassFit;
 
    minMassFit = MINmass;
    maxMassFit = MAXmass;
  
  // Fit data with background pdf for data limit
  RooRealVar* PhotonsMass = w->var("PhotonsMass");  
  PhotonsMass->setUnit("GeV");
 
  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
  
  for (int c = 0; c < ncat; ++c) {
    data[c] = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    
  
    // fit con expo pol
   
    RooFormulaVar *p1mod = new RooFormulaVar(TString::Format("par1Expol_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_expol1_cat%d",c)));
    RooFormulaVar *p2mod = new RooFormulaVar(TString::Format("par2Expol_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_expol2_cat%d",c)));
    RooFormulaVar *p3mod = new RooFormulaVar(TString::Format("par3Expol_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_expol3_cat%d",c)));
 
    PhotonsMass->setRange("bkg range", MINmass,MAXmass);
    RooFormulaVar *x     = new RooFormulaVar(TString::Format("xExpol_cat%d",c),"","@0",*w->var("PhotonsMass"));

   
    RooAbsPdf* PhotonsMassBkgTmp0 = new RooGenericPdf(TString::Format("PhotonsMassBkg_cat%d",c), "exp(-@0*@0/(@1+@2*@0+@3*@0*@0))", RooArgList(*x, *p1mod, *p2mod,*p3mod));
   

    fitresult[c] = PhotonsMassBkgTmp0->fitTo(*data[c], Range(minMassFit,maxMassFit),RooFit::FitOptions("MHTER"), SumW2Error(kTRUE), Save(kTRUE));
    w->import(*PhotonsMassBkgTmp0);
   
    std::cout<<TString::Format("******************************** Background Fit results mass %f cat %d***********************************", mass, c)<<std::endl;
    fitresult[c]->Print("V");
   

    //************************************************
    // Plot PhotonsMass background fit results per categories 
    TCanvas* ctmp = new TCanvas("ctmp","PhotonsMass Background Categories",0,0,500,500);
    Int_t nBinsMass(60);
    plotPhotonsMassBkg[c] = PhotonsMass->frame(minMassFit, maxMassFit,nBinsMass);
  
    data[c]->plotOn(plotPhotonsMassBkg[c],RooFit::Invisible());    
   
    PhotonsMassBkgTmp0->plotOn(plotPhotonsMassBkg[c],LineColor(kBlue),Range(minMassFit,maxMassFit),NormRange("bkg range")); 
    double chi2 = plotPhotonsMassBkg[c]->chiSquare(3);
    Int_t ndof = nBinsMass-2;
    std::cout<<"------> "<< ndof<<std::endl;
    double prob = TMath::Prob(chi2*ndof, ndof);
    std::cout<<prob<<std::endl;
    blind =false;
    if( blind ) {
     
      RooDataSet* data_down = (RooDataSet*) data[c].reduce(*w->var("PhotonsMass"),"PhotonsMass < 178.");
      RooDataSet* data_up = (RooDataSet*) data[c].reduce(*w->var("PhotonsMass"),"PhotonsMass >402");

      data_up->plotOn(plotPhotonsMassBkg[c]);    
      data_down->plotOn(plotPhotonsMassBkg[c]); 


   
    } else {
      data[c]->plotOn(plotPhotonsMassBkg[c]);    
      } 
       
  
    plotPhotonsMassBkg[c]->GetXaxis()->SetTitle("m_{#gamma #gamma}[GeV]");
    plotPhotonsMassBkg[c]->SetAxisRange(0.001,plotPhotonsMassBkg[c]->GetMaximum()*1.5,"Y");
    plotPhotonsMassBkg[c]->Draw();  

    TLegend *legdata = new TLegend(0.3790323,0.7775424,0.6290323,0.9279661, TString::Format("Category %d",c), "brNDC");
    legdata->AddEntry(plotPhotonsMassBkg[c]->getObject(2),"Data","LPE");
    legdata->AddEntry(plotPhotonsMassBkg[c]->getObject(1),"Parametric Model: Expol","L");
    legdata->SetTextSize(0.035);
    legdata->SetTextFont(42);
    // legdata->SetTextAlign(31);
    legdata->SetBorderSize(0);
    legdata->SetFillStyle(0);
    legdata->Draw("same");

    TPaveText* label_cms = get_labelCMS(0, "2012", false);
    TPaveText* label_sqrt = get_labelSqrt(0);
    label_cms->Draw("same");
    label_sqrt->Draw("same");

    //write down the chi2 of the fit on the
      
    TPaveText* label_chi2 = new TPaveText(0.5524194,0.6419492,0.796371,0.7690678, "brNDC");
    label_chi2->SetFillColor(kWhite);
    label_chi2->SetTextSize(0.035);
    label_chi2->SetTextFont(42);
    label_chi2->SetTextAlign(31); // align right
    label_chi2->AddText(TString::Format("Fit chi square/dof = %.3f", chi2));
    label_chi2->AddText(TString::Format("Chi square Prob = %.3f", prob));
    label_chi2->Draw("same");

    
    dobands = false;
    //********************************************************************************//
    if (dobands) {

      RooAbsPdf *cpdf; cpdf = PhotonsMassBkgTmp0;
      TGraphAsymmErrors *onesigma = new TGraphAsymmErrors();
      TGraphAsymmErrors *twosigma = new TGraphAsymmErrors();
      
      RooRealVar *nlim = new RooRealVar(TString::Format("nlim%d",c),"",0.0,0.0,10.0);
      nlim->removeRange();
      
      RooCurve *nomcurve = dynamic_cast<RooCurve*>(plotPhotonsMassBkg[c]->getObject(1));
      
      for (int i=1; i<(plotPhotonsMassBkg[c]->GetXaxis()->GetNbins()+1); ++i) {
	double lowedge = plotPhotonsMassBkg[c]->GetXaxis()->GetBinLowEdge(i);
	double upedge  = plotPhotonsMassBkg[c]->GetXaxis()->GetBinUpEdge(i);
	double center  = plotPhotonsMassBkg[c]->GetXaxis()->GetBinCenter(i);
	
	double nombkg = nomcurve->interpolate(center);
	nlim->setVal(nombkg);
	PhotonsMass->setRange("errRange",lowedge,upedge);
	RooAbsPdf *epdf = 0;
	epdf = new RooExtendPdf("epdf","",*cpdf,*nlim,"errRange");
	
	RooAbsReal *nll = epdf->createNLL(*(data[c]),Extended());
	RooMinimizer minim(*nll);
	minim.setStrategy(0);
	double clone = 1.0 - 2.0*RooStats::SignificanceToPValue(1.0);
	double cltwo = 1.0 - 2.0*RooStats::SignificanceToPValue(2.0);
	
	minim.migrad();
	minim.minos(*nlim);
	// printf("errlo = %5f, errhi = %5f\n",nlim->getErrorLo(),nlim->getErrorHi());
	
	onesigma->SetPoint(i-1,center,nombkg);
	onesigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
	
	minim.setErrorLevel(0.5*pow(ROOT::Math::normal_quantile(1-0.5*(1-cltwo),1.0), 2)); // the 0.5 is because qmu is -2*NLL
	// eventually if cl = 0.95 this is the usual 1.92!      
	
	minim.migrad();
	minim.minos(*nlim);
	
	twosigma->SetPoint(i-1,center,nombkg);
	twosigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
	
	delete nll;
	delete epdf;
      }

      PhotonsMass->setRange("errRange",minMassFit,maxMassFit);
      
      twosigma->SetLineColor(kGreen);
      twosigma->SetFillColor(kGreen);
      twosigma->SetMarkerColor(kGreen);
      twosigma->Draw("L3 SAME");
      
      onesigma->SetLineColor(kYellow);
      onesigma->SetFillColor(kYellow);
      onesigma->SetMarkerColor(kYellow);
      onesigma->Draw("L3 SAME");
      
      plotPhotonsMassBkg[c]->Draw("SAME"); 
    }
    int massI(mass);
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPOL_M%d.png",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPOL_M%d.pdf",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPOL_M%d.root",c,massI));

    ctmp->SetLogy();
    plotPhotonsMassBkg[c]->SetAxisRange(1.3,plotPhotonsMassBkg[c]->GetMaximum()*1.5,"Y");
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPOL_LOG_M%d.png",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPOL_LOG_M%d.pdf",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPOL_LOG_M%d.root",c,massI));

  }



  return fitresult;
}




void SetConstantParams(const RooArgSet* params) {

  cout << endl; cout << "Entering SetConstantParams" << endl;
  TIterator* iter(params->createIterator());
  for (TObject *a = iter->Next(); a != 0; a = iter->Next()) {
    RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
    if (rrv) { rrv->setConstant(true); std::cout << " " << rrv->GetName(); }
  }  
}

// Write signal pdfs and datasets into the workspace 
void MakeSigWS(RooWorkspace* w, const char* fileBaseName, Float_t width, std::string model){
  
  TString wsDir = "workspaces/"+filePOSTfix;
  Int_t ncat = NCAT;

  RooWorkspace *wAll = new RooWorkspace("w_all","w_all");  

  //********************************//
  // Retrieve P.D.F.s
   //w->Print("V");
  for (int c=0; c<ncat; ++c) {
    std::cout<<"flag"<<std::endl;
      wAll->import(*w->pdf("PhotonsMassSig"+TString::Format("_cat%d",c)));//*w->pdf("PhotonsMassSigCBCExt"+TString::Format("_cat%d",c))
    
      wAll->import(*w->data(TString::Format("SigWeight_cat%d",c)));
      wAll->import(*w->function("PhotonsMassSig"+TString::Format("_cat%d_norm",c)));
                                                 
  }
  std::cout << "done with importing signal pdfs" << std::endl;

  // (2) Systematics on energy scale and resolution // chiara: per ora tutte le sistematiche non hanno senso
  // wAll->factory("CMS_hgg_sig_m0_absShift[1,1.0,1.0]");
  // wAll->factory("CMS_hgg_sig_m0_absShift_cat0[1,1.0,1.0]");
  // wAll->factory("CMS_hgg_sig_m0_absShift_cat1[1,1.0,1.0]");
  // wAll->factory("prod::CMS_hgg_sig_m0_cat0(massggnewvtx_sig_m0_cat0, CMS_hgg_sig_m0_absShift)");
  // wAll->factory("prod::CMS_hgg_sig_m0_cat1(massggnewvtx_sig_m0_cat1, CMS_hgg_sig_m0_absShift)");

  // (3) Systematics on resolution: create new sigmas
  // wAll->factory("CMS_hgg_sig_sigmaScale[1,1.0,1.0]");
  // wAll->factory("CMS_hgg_sig_sigmaScale_cat0[1,1.0,1.0]");
  // wAll->factory("CMS_hgg_sig_sigmaScale_cat1[1,1.0,1.0]");
  // wAll->factory("prod::CMS_hgg_sig_sigma_cat0(massggnewvtx_sig_sigma0_cat0, CMS_hgg_sig_sigmaScale)");
  // wAll->factory("prod::CMS_hgg_sig_gsigma_cat0(massggnewvtx_sig_sigma1_cat0, CMS_hgg_sig_sigmaScale)");
  // wAll->factory("prod::CMS_hgg_sig_sigma_cat1(massggnewvtx_sig_sigma0_cat1, CMS_hgg_sig_sigmaScale)");
  // wAll->factory("prod::CMS_hgg_sig_gsigma_cat1(massggnewvtx_sig_sigma1_cat1, CMS_hgg_sig_sigmaScale)")

  TString filename(wsDir+TString(fileBaseName)+TString::Format(("_m%.2f_w%.2f.inputsig_"+model+"_FixW.root").c_str(),w->var("MH")->getVal(),width));
  wAll->writeToFile(filename);
  cout << "Write signal workspace in: " << filename << " file" << endl;
  
  return;
}

// Write background pdfs and datasets into the workspace 
void MakeBkgWS(RooWorkspace* w, const char* fileBaseName, double mass) {

  TString wsDir = "workspaces/"+filePOSTfix;
  Int_t ncat = NCAT;  

  RooWorkspace *wAll = new RooWorkspace("w_all","w_all");
  
  //********************************//
  // Retrieve the datasets and PDFs
  RooDataSet* data[NCAT];
 
  for (int c=0; c<ncat; ++c) {
  
    data[c] = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    ((RooRealVar*) data[c]->get()->find("PhotonsMass"))->setBins(320) ;
 
    RooDataHist* dataBinned = data[c]->binnedClone();
 
    wAll->import(*w->pdf(TString::Format("PhotonsMassBkg_cat%d",c)));
 
    wAll->import(*dataBinned, Rename(TString::Format("data_obs_cat%d",c)));
 
    wAll->import(*w->data(TString::Format("Data_cat%d",c)), Rename(TString::Format("data_unbinned_obs_cat%d",c)));
 
  }
  std::cout << "done with importing background pdfs" << std::endl;
  

  TString filename;
  filename = (wsDir+TString(fileBaseName)+TString::Format("_m%.2f.root",w->var("MH")->getVal()));


  wAll->writeToFile(filename);
  cout << "Write background workspace in: " << filename << " file" << endl;

  std::cout << std::endl; 
  std::cout << "observation:" << std::endl;
  for (int c=0; c<ncat; ++c) {
    std::cout << "cat " << c << ", " << wAll->data(TString::Format("data_obs_cat%d",c))->sumEntries() << endl;
    wAll->data(TString::Format("data_obs_cat%d",c))->Print();
  }
  std::cout << std::endl;
  
  return;
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

// preparing datacards
void MakeDataCard_1Channel(RooWorkspace* w, const char* fileBaseName, const char* fileBkgName, Float_t width,int iChan, std::string model) {

  TString cardDir = "/afs/cern.ch/work/s/soffi/CMSSW611-Analysis/src/h2gglobe/ChiaraFitLimits/datacardWithAllSyst/"+filePOSTfix;
  Int_t ncat = NCAT;
  TString wsDir   = "/afs/cern.ch/work/s/soffi/CMSSW611-Analysis/src/h2gglobe/ChiaraFitLimits/workspaces/"+filePOSTfix;

  // **********************
  // Retrieve the datasets
  cout << "Start retrieving dataset" << endl;
  
  RooDataSet* data[9];
  RooDataSet* signal[9];
  for (int c=0; c<ncat; ++c) {
    data[c]   = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    signal[c] = (RooDataSet*) w->data(TString::Format("SigWeight_cat%d",c));
  }

  RooRealVar*  lumi = w->var("lumi");
  /*
  // *****************************
  // Print Expected event yields
  cout << "======== Expected Events Number =====================" << endl;  
  cout << ".........Measured Data for L = " << lumi->getVal() << " pb-1 ............................" << endl;  
  cout << "#Events data: " <<  w->data("Data")->sumEntries()  << endl;
  for (int c = 0; c < ncat; ++c) {
    cout << TString::Format("#Events data cat%d:   ",c) << data[c]->sumEntries()  << endl;
  }
  cout << ".........Expected Signal for L = " << lumi->getVal() << " pb-1 ............................" << endl;  
  cout << "#Events Signal:      " << w->data("SigWeight")->sumEntries()  << endl;
  Float_t siglikeErr[6];
  for (int c = 0; c < ncat; ++c) {
    cout << TString::Format("#Events Signal cat%d: ",c) << signal[c]->sumEntries() << endl;
    siglikeErr[c]=0.6*signal[c]->sumEntries();
  }
  cout << "====================================================" << endl;  

 */
  // *****************************
  // Printdata Data Card int file
  TString filename(cardDir+TString(fileBaseName)+"_"+"8TeV"+Form(("_m%.2f_w%.2f_channel%d_"+model+"_FixW.txt").c_str(),w->var("MH")->getVal(),width,iChan));
  ofstream outFile(filename);

  outFile << "#CMS-HGG HighMass DataCard for Unbinned Limit Setting, " << lumi->getVal() <<  " pb-1 " << endl;
  outFile << "#Run with: combine -d datacardName.txt -U -m *mass* -H ProfileLikelihood -M MarkovChainMC --rMin=0 --rMax=20.0  -b 3000 -i 50000 --optimizeSim=1 --tries 30" << endl;
  outFile << "# Lumi =  " << lumi->getVal() << " pb-1" << endl;
  outFile << "imax *" << endl;
  outFile << "jmax *" << endl;
  outFile << "kmax *" << endl;
  outFile << "---------------" << endl;

  outFile << "shapes data_obs * " << wsDir+TString(fileBkgName)+TString::Format("_m%.2f.root",w->var("MH")->getVal()) << " w_all:data_obs_$CHANNEL" << endl;
  outFile << "shapes sig * "      << wsDir+TString(fileBaseName)+"_8TeV"+TString::Format(("_m%.2f_w%.2f.inputsig_"+model+"_FixW.root").c_str(),w->var("MH")->getVal(),width) << " w_all:PhotonsMassSig_$CHANNEL" << endl;
  outFile << "shapes bkg * "      << wsDir+TString(fileBkgName)+TString::Format("_m%.2f.root",w->var("MH")->getVal()) << " w_all:PhotonsMassBkg_$CHANNEL" << endl;

  outFile << "---------------" << endl;
  outFile << Form("bin          cat%d", iChan) << endl;
  outFile <<  "observation   "  <<  Form("%.10lg",data[iChan]->sumEntries()) << endl;
  outFile << "------------------------------" << endl;
  outFile << "bin                 " << Form("cat%d      cat%d", iChan, iChan) << endl;
  outFile << "process                 sig      bkg" << endl;
  outFile << "process                   0        1" << endl;
  // if(signalScaler==1.)
  // signalScaler=1./signal[2]->sumEntries()*20;
  outFile << "rate                   " 
    //	  << signal[iChan]->sumEntries()*signalScaler << " " << data[iChan]->sumEntries() << endl;
    // << 1 << " " << data[iChan]->sumEntries() << endl;
	  <<19620. << "  "<<data[iChan]->sumEntries() << endl;
  outFile << "--------------------------------" << endl;
  outFile << "# signal scaled by " << signalScaler << endl;

 
  outFile << "lumi_8TeV     lnN     1.026000  - " << endl;
  outFile << "eff_trig     lnN     1.010000  - " << endl;
  outFile << "global_syst     lnN     1.050000  - " << endl;
  if(iChan==0){    
    outFile << "id_eff_eb     lnN     1.005000  - " << endl;    
    outFile << "id_eff_ee     lnN     1.000000  - " << endl;    
    outFile << "r9Eff   lnN   1.0145/0.9915   - " << endl;
    outFile << "vtxEff   lnN   0.996/1.008   - " << endl; 
  }else if(iChan==1){    
    outFile << "id_eff_eb     lnN     1.005000  - " << endl;    
    outFile << "id_eff_ee     lnN     1.000000  - " << endl;    
    outFile << "r9Eff   lnN   0.985/1.0085   - " << endl;
    outFile << "vtxEff   lnN   0.998/1.005   - " << endl; 
  }else if(iChan==2){    
    outFile << "id_eff_eb     lnN     1.00000 - " << endl;    
    outFile << "id_eff_ee     lnN     1.02600  - " << endl;    
    outFile << "r9Eff   lnN   1.0193/0.964   - " << endl;
    outFile << "vtxEff   lnN   0.996/1.007   - " << endl; 
  }else if(iChan==3){    
    outFile << "id_eff_eb     lnN     1.000000  - " << endl;    
    outFile << "id_eff_ee     lnN     1.026000  - " << endl;    
    outFile << "r9Eff   lnN   0.981/1.0357   - " << endl;
    outFile << "vtxEff   lnN   0.998/1.003   - " << endl; 
  }
    
  outFile << Form("mShift_cat%d    param   0 1 ",iChan) << endl;
  outFile << Form("mSmear_cat%d    param   0 1 ",iChan) << endl;



  // outFile << "CMS_VV_eff_g         lnN  0.8/1.20      - # Signal Efficiency" << endl;
  // outFile << "# Parametric shape uncertainties, entered by hand." << endl;
  // outFile << Form("CMS_hgg_sig_m0_absShift    param   1   0.0125   # displacement of the mean w.r.t. nominal in EB*EX category, good R9",iChan) << endl;
  // outFile << Form("CMS_hgg_sig_sigmaScale     param   1   0.1   # multiplicative correction to sigmas in EB*EX category, good R9",iChan) << endl;
  // outFile << Form("rooHistFunc_cat%d_norm       ",iChan) << endl;
  // outFile << Form("CMS_hgg_bkg_8TeV_slope2_cat%d         flatParam  # Mean and absolute uncertainty on background slope",iChan) << endl;
  // outFile << Form("CMS_hgg_bkg_8TeV_slope3_cat%d         flatParam  # Mean and absolute uncertainty on background slope",iChan) << endl;
  // if (iChan != 2 )  outFile << Form("CMS_hgg_bkg_8TeV_slope1_cat%d         flatParam  # Mean and absolute uncertainty on background slope",iChan) << endl;
  outFile.close();

  cout << "Write data card in: " << filename << " file" << endl;

  return;
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
      if((mass>1000)&&(mass<2000))
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
  // sigTree1->Add("histograms_CMS-HGG_24072013.root/dipho_Box_250_8TeV");
  

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

  //  RooArgSet* argset_ = new RooArgSet(*w->var("PhotonsMass"), *w->var("evweight"));

  for (int c=1; c<2; ++c) {

    // 0) chiara: 1cat only
    // signal[c] =  (RooDataSet*) sigWeighted.reduce(*w->var("massggnewvtx"),mainCut);   

    data[c] = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));

    // 1)  prime 4 cat livia
  
   
    //reduce QCD dataset
    if (c==0) qcdMC[c] = (RooDataSet*) qcdMCWeighted.reduce(*ntplVars, mainCut+TString::Format("&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 )"));
    if (c==1) qcdMC[c] = (RooDataSet*) qcdMCWeighted.reduce(*ntplVars,mainCut+TString::Format("&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9<0.94 || ph2_r9<0.94 ) "));
    if (c==2) qcdMC[c] = (RooDataSet*) qcdMCWeighted.reduce(*ntplVars,mainCut+TString::Format("&& (abs(ph1_eta)>1.4442 || abs(ph2_eta)>1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 )"));
    if (c==3) qcdMC[c] = (RooDataSet*) qcdMCWeighted.reduce(*ntplVars,mainCut+TString::Format("&& (abs(ph1_eta)>1.4442 || abs(ph2_eta)>1.4442) && (ph1_r9>0.94 || ph2_r9>0.94 ) "));

 
    //reduce gj dataset
    if (c==0) gjMC[c] = (RooDataSet*) gjMCWeighted.reduce(*ntplVars, mainCut+TString::Format("&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 )"));
    if (c==1) gjMC[c] = (RooDataSet*) gjMCWeighted.reduce(*ntplVars,mainCut+TString::Format("&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9<0.94 || ph2_r9<0.94 ) "));
    if (c==2) gjMC[c] = (RooDataSet*) gjMCWeighted.reduce(*ntplVars,mainCut+TString::Format("&& (abs(ph1_eta)>1.4442 || abs(ph2_eta)>1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 )"));
    if (c==3) gjMC[c] = (RooDataSet*) gjMCWeighted.reduce(*ntplVars,mainCut+TString::Format("&& (abs(ph1_eta)>1.4442 || abs(ph2_eta)>1.4442) && (ph1_r9>0.94 || ph2_r9>0.94 ) "));
 
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
    if (c==3) BkgMC[c] = (RooDataSet*) BkgMCWeighted.reduce(*ntplVars,mainCut+TString::Format("&& (abs(ph1_eta)>1.4442 || abs(ph2_eta)>1.4442) && (ph1_r9>0.94 || ph2_r9>0.94 ) "));


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
      //   BkgMCKeyPdf[c] = new RooKeysPdf(TString::Format("BkgMCKeyPdf_cat%d",c),TString::Format("BkgMCKeyPdf_cat%d",c),*PhotonsMass,*BkgMCcopy[c],RooKeysPdf::MirrorBoth );
      BkgMCKeyPdf_bw2[c] = new RooKeysPdf(TString::Format("BkgMCKeyPdf_bw2_cat%d",c),TString::Format("BkgMCKeyPdf_bw2_cat%d",c),*PhotonsMass,*BkgMCcopy[c],RooKeysPdf::MirrorBoth,2 );
      BkgMCKeyPdf_bw3[c] = new RooKeysPdf(TString::Format("BkgMCKeyPdf_bw3_cat%d",c),TString::Format("BkgMCKeyPdf_bw3_cat%d",c),*PhotonsMass,*BkgMCcopy[c],RooKeysPdf::MirrorBoth,3 );
      BkgMCKeyPdf_bw4[c] = new RooKeysPdf(TString::Format("BkgMCKeyPdf_bw4_cat%d",c),TString::Format("BkgMCKeyPdf_bw4_cat%d",c),*PhotonsMass,*BkgMCcopy[c],RooKeysPdf::MirrorBoth,4 );
      
    }else{

      //create the rookeyspdf
      BkgMCKeyPdf[c] = new RooKeysPdf(TString::Format("BkgMCKeyPdf_cat%d",c),TString::Format("BkgMCKeyPdf_cat%d",c),*PhotonsMass,*BkgMCcopy[c],RooKeysPdf::NoMirror );
      BkgMCKeyPdf_bw2[c] = new RooKeysPdf(TString::Format("BkgMCKeyPdf_bw2_cat%d",c),TString::Format("BkgMCKeyPdf_bw2_cat%d",c),*PhotonsMass,*BkgMCcopy[c],RooKeysPdf::NoMirror,2 );
      BkgMCKeyPdf_bw3[c] = new RooKeysPdf(TString::Format("BkgMCKeyPdf_bw3_cat%d",c),TString::Format("BkgMCKeyPdf_bw3_cat%d",c),*PhotonsMass,*BkgMCcopy[c],RooKeysPdf::NoMirror,3 );
      BkgMCKeyPdf_bw4[c] = new RooKeysPdf(TString::Format("BkgMCKeyPdf_bw4_cat%d",c),TString::Format("BkgMCKeyPdf_bw4_cat%d",c),*PhotonsMass,*BkgMCcopy[c],RooKeysPdf::NoMirror,4 );
    

    }


    wBias->import(*BkgMCKeyPdf_bw2[c]);
    //   wBias->import(*BkgMCKeyPdf_bw[c]);
    wBias->import(*BkgMCKeyPdf_bw3[c]);
    wBias->import(*BkgMCKeyPdf_bw4[c]);
    
    plotPhotonsMassBkgMC[c] = PhotonsMass->frame(minMassFit, maxMassFit,nBinsMass);
    data[c]->plotOn(plotPhotonsMassBkgMC[c]);
    //BkgMC[c]->plotOn(plotPhotonsMassBkgMC[c],"PE", MarkerColor(kRed), LineColor(kRed), MarkerSize(1.));
    //BkgMCcopy[c]->plotOn(plotPhotonsMassBkgMC[c],"PE", MarkerColor(kGreen), LineColor(kGreen), MarkerSize(0.5));
    //   BkgMCKeyPdf[c]->plotOn(plotPhotonsMassBkgMC[c],"L", LineColor(kBlue), LineWidth(2));
    BkgMCKeyPdf_bw2[c]->plotOn(plotPhotonsMassBkgMC[c],"L", LineColor(kBlack), LineWidth(2));
    BkgMCKeyPdf_bw3[c]->plotOn(plotPhotonsMassBkgMC[c],"L", LineColor(8), LineWidth(2));
    BkgMCKeyPdf_bw4[c]->plotOn(plotPhotonsMassBkgMC[c],"L", LineColor(kOrange), LineWidth(2));

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
      leg2->AddEntry(plotPhotonsMassBkgMC[c]->getObject(4),"Bw x 4 NoMirr","L");
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
  

  TString filename(wsDir+"HighMass-hgg.RooKeysPdfMCBkg_8TeV_cat123.root");
  wBias->writeToFile(filename);
  cout << "Write background RooKeys workspace in: " << filename << " file" << endl;

  return;

}



RooHistFunc* getRooHistFunc(int cat, RooRealVar* var, std::string model){


  double mass[7] = {150., 200., 250., 300., 400., 600, 800};

  double c0[7];
  double c1[7];
  double c2[7];
  double c3[7];
  double all[7];
  //changed 10/04
  /*  double c0g[5] = {0.00300125 , 0.00377226 , 0.00446634 , 0.0051454 , 0.00592569};
  double c1g[5] = {0.00416221 , 0.00503098 , 0.00494583 , 0.00483801 , 0.00465511 };
  double c2g[5] = {0.0015748 , 0.00188656 , 0.00215004 , 0.00246588 , 0.00260045 };
  double c3g[5] = {0.00209722 , 0.00255168 , 0.0025782 , 0.00251383 , 0.00218195 };
  double allg[5] = {0.0108564 , 0.0132374 , 0.0141401 , 0.014937 , 0.0153788 };  
  double c0h[5] = {0.00302 , 0.003772 , 0.0044645 , 0.005133 , 0.00595};
  double c1h[5] = {0.00414 , 0.004974 , 0.0049136 , 0.0048243 , 0.00463 };
  double c2h[5] = {0.00158 , 0.001885 , 0.0021557 , 0.0024568 , 0.0026 };
  double c3h[5] = {0.002092 , 0.002533 , 0.00256584 , 0.0024928 , 0.002173 };
  double allh[5] = {0.010853 , 0.01316 , 0.0140993 , 0.014881 , 0.0153686};*/

  double c0g[7] = {0.00281259889626 ,0.00360405499841 ,0.0042512557194 ,0.00480861314285 ,0.00555793387016 ,0.00711359319501 ,0.00823700633031};
  double c1g[7] = { 0.00384737104156 ,0.00409835498596 ,0.00396547370455 ,0.003904278683 ,0.00368952356883 ,0.00312694674604 ,0.00264070618382   };
  double c2g[7] = {0.00147564281183 ,0.00190380886793 ,0.00200106257462 ,0.00237975631151 ,0.00252602251139 ,0.00252633118101 ,0.00243976040669 };
  double c3g[7] = {0.00179971742008 ,0.001998562755 ,0.00209254398016 ,0.00195718051131 ,0.00166752834056 ,0.00124744484509 ,0.000929175591328};
  double allg[7] = { 0.00993533016973 ,0.0116047816073 ,0.0123103359787 ,0.0130498286487 ,0.0134410082909 ,0.0140143159671 ,0.0142466485122};  

  double c0v[7] = {0.00295597,   0.00351098,   0.00404037,   0.00454414,   0.00547479,0.00702859,   0.00817238};
  double c1v[7] = {0.00402382,   0.00410415,   0.00415866,   0.00418736,   0.00416732,0.00381745,   0.00305454   };
  double c2v[7] = {0.00144884,   0.00165231,   0.00183088,   0.00198455,   0.0022172,0.00238368 ,  0.00215175 };
  double c3v[7] = {0.00205624,   0.00204125,   0.00201627,   0.00198131,   0.00188142,0.0015618 ,  0.00108238 };
  double allv[7] = {0.0104849,   0.0113087,   0.0120462,   0.0126974.   0.0137407,0.0147915,   0.0144611   };
  
  for(int i = 0; i<7; i++){
    if (model=="GGH"){
      c0[i] = c0g[i];
      c1[i] = c1g[i];
      c2[i] = c2g[i];
      c3[i] = c3g[i];
      all[i] = allg[i];

    }else if(model=="VBF"){
      c0[i] = c0v[i];
      c1[i] = c1v[i];
      c2[i] = c2v[i];
      c3[i] = c3v[i];
      all[i] = allv[i];


    }
  }
  
  TH1F* h_all = new TH1F("h_all", "h_all", 15, 130., 875.);
  //m150-400
  for(int i=0;i<5;i++){
    std::cout<<"i: "<<i<<" mass: "<<mass[i]<<std::endl;
    if(cat==0) h_all->SetBinContent(h_all->FindBin(mass[i]),c0[i]);
    if(cat==1) h_all->SetBinContent(h_all->FindBin(mass[i]),c1[i]);
    if(cat==2) h_all->SetBinContent(h_all->FindBin(mass[i]),c2[i]);
    if(cat==3) h_all->SetBinContent(h_all->FindBin(mass[i]),c3[i]);
    if(cat==4) h_all->SetBinContent(h_all->FindBin(mass[i]),all[i]);
   
  }
  std::cout<<"--------------------------->"<<cat<<std::endl;
  //m350
  if(cat==0) h_all->SetBinContent(h_all->FindBin(350),(c0[3]+c0[4])/2);
  if(cat==1) h_all->SetBinContent(h_all->FindBin(350),(c1[3]+c1[4])/2);
  if(cat==2) h_all->SetBinContent(h_all->FindBin(350),(c2[3]+c2[4])/2);
  if(cat==3) h_all->SetBinContent(h_all->FindBin(350),(c3[3]+c3[4])/2);
  if(cat==4) h_all->SetBinContent(h_all->FindBin(350),(all[3]+all[4])/2);

  //m450-550 
  double r1mass[3]={450, 500, 550};
  for(int i = 0; i<3;i++){
  if(cat==0) h_all->SetBinContent(h_all->FindBin(r1mass[i]),(c0[4]+c0[5])/2);
  if(cat==1) h_all->SetBinContent(h_all->FindBin(r1mass[i]),(c1[4]+c1[5])/2);
  if(cat==2) h_all->SetBinContent(h_all->FindBin(r1mass[i]),(c2[4]+c2[5])/2);
  if(cat==3) h_all->SetBinContent(h_all->FindBin(r1mass[i]),(c3[4]+c3[5])/2);
  if(cat==4) h_all->SetBinContent(h_all->FindBin(r1mass[i]),(all[4]+all[5])/2);
  }

  //m600 e 850
  double highmass[2]={ 600., 850.};
  for(int i = 0; i<2;i++){
  if(cat==0) h_all->SetBinContent(h_all->FindBin(highmass[i]),c0[i+5]);
  if(cat==1) h_all->SetBinContent(h_all->FindBin(highmass[i]),c1[i+5]);
  if(cat==2) h_all->SetBinContent(h_all->FindBin(highmass[i]),c2[i+5]);
  if(cat==3) h_all->SetBinContent(h_all->FindBin(highmass[i]),c3[i+5]);
  if(cat==4) h_all->SetBinContent(h_all->FindBin(highmass[i]),all[i+5]);
  }
  

  //m650-800 
  double r2mass[4]={650, 700, 750, 800};
  for(int i = 0; i<4;i++){
  if(cat==0) h_all->SetBinContent(h_all->FindBin(r2mass[i]),(c0[5]+c0[6])/2);
  if(cat==1) h_all->SetBinContent(h_all->FindBin(r2mass[i]),(c1[5]+c1[6])/2);
  if(cat==2) h_all->SetBinContent(h_all->FindBin(r2mass[i]),(c2[5]+c2[6])/2);
  if(cat==3) h_all->SetBinContent(h_all->FindBin(r2mass[i]),(c3[5]+c3[6])/2);
  if(cat==4) h_all->SetBinContent(h_all->FindBin(r2mass[i]),(all[5]+all[6])/2);
  }

  for(int i=0;i<16;i++){
    std::cout<<cat<<std::endl;
    // std::cout<<mass[i]<<"  "<<h_all->FindBin(mass[i])<<"   "<<c0[i]<<std::endl;
    std::cout<<i<<"  "<<h_all->GetBinCenter(i)<<"  "<<h_all->GetBinContent(i)<<std::endl;

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
  /*  RooPlot* plot = var->frame();
  rooData_all->plotOn(plot);
  rooFunc_all->plotOn(plot);
  plot->Draw("same");
  c->SetLogy();

  c->SaveAs("plots/signalYield.png");
  c->SaveAs("plots/signalYield.pdf");*/
  
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
  TPaveText* label_cms = get_labelCMS(0, false);
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
  TPaveText* label_cms = get_labelCMS(0, false);
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
  TPaveText* label_cms = get_labelCMS(0,false);
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
  TPaveText* label_cms = get_labelCMS(0, false);
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
