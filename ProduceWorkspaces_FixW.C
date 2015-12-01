#include "HighMass-HggFitter_mgg_FixW.cc"
//

void ProduceWorkspaces(int mass, Double_t width, std::string model){
  // gSystem->SetIncludePath("-I/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.32.00-cms5/include");
  /* gSystem->SetIncludePath("-I/afs/cern.ch/cms/slc5_amd64_gcc472/lcg/roofit/5.34.04-cms2/include");
  gSystem->Load("libRooFit");
  //
gROOT->ProcessLine(".L RooCruijff.cxx+");
  gROOT->ProcessLine(".x RooCBCBPdf.cxx+");
  gROOT->ProcessLine(".x RooCPSHighMassGGHNoInterf.cxx+");*/
  //gSystem->SetIncludePath("-I/afs/cern.ch/cms/slc5_amd64_gcc472/lcg/roofit/5.34.04-cms2/include");
  //  gROOT->ProcessLine(".L RooCPSHighMassGGHNoInterf.cxx++");
  cout << "mass = " << mass << endl; 
  // runfits(mass, true);
  runfits(mass, true, width, model);
  
 
}
