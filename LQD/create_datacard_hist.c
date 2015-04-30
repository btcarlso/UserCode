void create_datacard_hist(int M){
  TFile *fbkg =new TFile("output_file_allMC.root","READ"); 
  TFile *fsig = new TFile(Form("output_file_RPV_LQD221_M%d.root",M),"READ"); 
  TFile *fdata = new TFile("output_file_singleMu.root","READ");

  TFile *out = new TFile("out.root","RECREATE"); 
  
  for(int nJ=5; nJ<=7; nJ++){
    out->mkdir(Form("njets_%d",nJ)); 
    out->cd(Form("njets_%d",nJ)); 
    TH1F *hbkg=(TH1F*)fbkg->Get(Form("h_nJets%d_st",nJ)); 
    hbkg->SetName("bkg");
    
    TH1F *hsig=(TH1F*)fsig->Get(Form("h_nJets%d_st",nJ));
    hsig->SetName("sig");

    TH1F *hdata=(TH1F*)fdata->Get(Form("h_nJets%d_st",nJ));
    hdata->SetName("data_obs");

    cout << hdata->Integral() << endl; 
    cout << hsig->Integral() << " " << hbkg->Integral() << endl; 
    
    hbkg->Write();
    hsig->Write();
    hdata->Write(); 

  }
}
