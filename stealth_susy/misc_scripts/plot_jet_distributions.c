{
  gROOT->SetBatch();
  TFile *input = new TFile("output_file_wJets.root","READ"); 
  TFile *output = new TFile("output_jets.root","RECREATE"); 
  for(int nJ=2; nJ<=5; nJ++){
    TCanvas *C = new TCanvas(Form("nJets_%d",nJ)); 
    C->Divide(2,2); 
    for(int j=1; j<=4; j++){
      C->cd(j); 
      TH1F *h=(TH1F*)input->Get(Form("jet%d_pt_nJet%d",j,nJ)); 
      h->DrawCopy(); 

    }
    output->cd();
    C->Write();
  }



}
