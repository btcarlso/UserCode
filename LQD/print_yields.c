print_yields(){
  float stmin300[]={475,325,575};
  float stmin900[]={1425,1425,1525}; 

  TFile *fbkg = new TFile("output_file_allMC.root","READ"); 
  TFile *fsig = new TFile("output_file_RPV_LQD221_M300.root","READ"); 

  for (int nJ=5; nJ<=7; nJ++){
    TH1F *hbkg=(TH1F*)fbkg->Get(Form("h_nJets%d_stmin",nJ)); 
    TH1F *hsig=(TH1F*)fsig->Get(Form("h_nJets%d_stmin",nJ)); 
    int ibin=hbkg->FindBin(stmin300[nJ-5]);
    cout << "nJ: " << nJ << " "; 
    cout << "bkg: " << hbkg->GetBinContent(ibin) << "  "; 
    cout << "sig: " << hsig->GetBinContent(ibin) << endl; 


    }
}
