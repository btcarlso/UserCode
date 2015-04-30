void plot_pt(){

  TH2F *wJets=(TH2F*)_file0->Get("jetpt_nJets"); 
  TH2F *ttbar=(TH2F*)_file1->Get("jetpt_nJets");
  TH2F *allMC=(TH2F*)_file2->Get("jetpt_nJets");
  TH2F *data=(TH2F*)_file3->Get("jetpt_nJets");

  TProfile *wJets_px=(TProfile*)wJets->ProfileX("wJets_px");
  TProfile *ttbar_px=(TProfile*)ttbar->ProfileX("ttbar_px");
  TProfile *allMC_px=(TProfile*)allMC->ProfileX("allMC_px");
  TProfile *data_px=(TProfile*)data->ProfileX("data_px");

  wJets_px->SetLineColor(kRed);
  ttbar_px->SetLineColor(kBlue);
  allMC_px->SetLineColor(kBlack);
  
  allMC_px->Draw("histo");
  wJets_px->Draw("histo same");
  ttbar_px->Draw("histo same");
  data_px->Draw("E1 same");

}
