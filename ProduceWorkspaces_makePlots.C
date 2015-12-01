#include "HighMass-HggFitter_mgg_makePlots.cc"
//

void ProduceWorkspaces_makePlots(int mass){
  cout << "Now plot Data vs MC bkg" << endl;
  gSystem->SetIncludePath("-I/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.32.00-cms5/include");
  cout << "Now plot Data vs MC bkg" << endl;
  gSystem->Load("libRooFit");
  cout << "Now plot Data vs MC bkg" << endl;
  // gROOT->ProcessLine(".x RooCBCrujffPdf.cxx+");
  cout << "Now plot Data vs MC bkg" << endl;
  cout << "mass = " << mass << endl; 
  // runfits(mass, true);
  runfits(mass, true);
  
 
}
