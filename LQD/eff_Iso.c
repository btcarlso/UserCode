void eff_Iso(){
  TFile *fdata =new TFile("output_file_singleEle.root","READ");
  TFile *fmc =new TFile("output_file_test.root","READ");

  TH1F *pass_data=(TH1F*)fdata->Get("nJets_Iso_pass_1btag");
  TH1F *total_data=(TH1F*)fdata->Get("nJets_Iso_total_1btag");

  TH1F *pass=(TH1F*)fmc->Get("nJets_Iso_pass_1btag"); 
  TH1F *total=(TH1F*)fmc->Get("nJets_Iso_total_1btag");

  for(int nJ=2; nJ<=7; nJ++){
    double effData=pass_data->GetBinContent(nJ)/total_data->GetBinContent(nJ); 
    double eff=pass->GetBinContent(nJ)/total->GetBinContent(nJ); 

    cout << nJ << " jets eff: " << eff << " +/- " << TMath::Sqrt(eff*(1-eff)/total->GetBinContent(nJ)) << " data  " << effData <<  " +/- " << TMath::Sqrt(effData*(1-effData)/total_data->GetBinContent(nJ)) << endl; 
    
  }

}
