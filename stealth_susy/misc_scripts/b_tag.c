void b_tag(){

  vector<bool> *jet_btag=0;
  TBranch *b_jet_btag;
  Int_t nJets=-9; 

  TH1F *nJ=new TH1F("nJ","n-jets",10,0.5,10.5);
  TH1F *nb_jet=new TH1F("nb-jet","b-jet", 11,-0.5,10.5);
  
  TTree *tree = (TTree*)_file0->Get("tree"); 
  tree->SetBranchAddress("jets_n",&nJets); 
  tree->SetBranchAddress("jet_btag",&jet_btag, &b_jet_btag); 
  int entries=tree->GetEntries();
  for(int i=0; i<entries;i++){
    tree->GetEntry(i);
    //cout << "nJets: " << nJets; 
    int nB=0; 
    for(int ib=0; ib<jet_btag->size();ib++){
      bool bT=jet_btag->at(ib);
      if(bT==1) nB++;
    }
    //cout << " btags: " << nB << endl; 
    nJ->Fill(nJets);
    nb_jet->Fill(nB);
  }
  TCanvas *C = new TCanvas();
  C->Divide(1,2);
  C->cd(1);
  nJ->Draw();
  C->cd(2);
  nb_jet->Draw();
}
