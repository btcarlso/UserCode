void print_table(){
  TFile *fbkg = new TFile("output_file_ttFullLept.root","READ"); 
  TFile *fdata = new TFile("output_file_singleMu.root","READ"); 

  TH1F *hdata=(TH1F*)fdata->Get("h_misQID"); 
  TH1F *hmc=(TH1F*)fbkg->Get("h_misQID");

  hdata->LabelsDeflate();
  for(int i=1; i<=hdata->GetNbinsX(); i++){
    cout << hdata->GetXaxis()->GetBinLabel(i) << " " << hdata->GetBinContent(i) << endl; 
  }
  for(int i=1; i<=hdata->GetNbinsX(); i++){
    cout << hmc->GetXaxis()->GetBinLabel(i) << " " <<hmc->GetBinContent(i) << endl;
  }

}
