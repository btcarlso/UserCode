void k_factor(){

  TFile *wJets_HT=new TFile("output_file_trigger_wJetsHTbinnedLO_1mu.root","READ"); //HT binned LO
  TFile *wJets_jet = new TFile("output_file_trigger_wJetsLO_1mu.root","READ"); //Jet binned LO
  
  TCanvas *C = new TCanvas(); 

  for(int nJ=2; nJ<=6; nJ++){
    TH1F *ht=(TH1F*)wJets_HT->Get(Form("ht_nJets%d",nJ)); 
    TH1F *ht_b=(TH1F*)wJets_jet->Get(Form("ht_nJets%d",nJ)); 

    ht_b->Divide(ht); 
    ht_b->GetYaxis()->SetTitle("Jet binned / HT binned sample");
    ht_b->SetAxisRange(300,3000); 
    ht_b->SetLineColor(nJ); 
    ht_b->SetMarkerColor(nJ); 

    if(nJ==2){
      ht_b->Draw(); 
      //      ht->Draw("same");
    }


  }


}
