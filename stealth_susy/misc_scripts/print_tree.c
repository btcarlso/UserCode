{
  gSystem->Load("libSusyEvent.so");
  gROOT->LoadMacro("GenTreeViewerRA3.cc+");

  TFile* source = TFile::Open("dcap:///pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/w1JetsToLNu_v1/susyEvents_145_1_FoQ.root");
  TTree* susyTree = (TTree*)source->Get("susyTree");
  susy::Event* susyEvent = new susy::Event();

  susyTree->SetBranchAddress("susyEvent", &susyEvent);

  float minPt = 3.;

  for(int i = 0; i < 3; ++i){
    susyTree->GetEntry(i);
    viewGenTreeRA3(*susyEvent, minPt);
  }
}
