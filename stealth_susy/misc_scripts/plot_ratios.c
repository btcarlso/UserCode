void plot_ratios(){

  TFile *wJets_st= new TFile("output_file_wJets.root","READ"); 
  TFile *wJets_genst = new TFile("output_file_wJetsGen.root","READ"); 

  TGraphAsymmErrors *r43G=(TGraphAsymmErrors*)wJets_genst->Get("Ratio4-3"); 
  TGraphAsymmErrors *r53G=(TGraphAsymmErrors*)wJets_genst->Get("Ratio5-3");
  TGraphAsymmErrors *r63G=(TGraphAsymmErrors*)wJets_genst->Get("Ratio6-3");

  TGraphAsymmErrors *r43=(TGraphAsymmErrors*)wJets_st->Get("Ratio4-3");
  TGraphAsymmErrors *r53=(TGraphAsymmErrors*)wJets_st->Get("Ratio5-3");
  TGraphAsymmErrors *r63=(TGraphAsymmErrors*)wJets_st->Get("Ratio6-3");


  r43->SetMarkerStyle(8); 
  r43->SetMarkerSize(0.5); 
  r43->SetMarkerColor(kRed); 

  r53->SetMarkerStyle(8);
  r53->SetMarkerSize(0.5); 
  r53->SetMarkerColor(kRed);

  r63->SetMarkerStyle(8);
  r63->SetMarkerSize(0.5);
  r63->SetMarkerColor(kRed);


  TCanvas *C = new TCanvas();
  C->Divide(2,2); 
  C->cd(1); 
  //  r43->Draw("ap");
  r43G->Draw("ap");
  C->cd(2);
  // r53->Draw("ap"); 
  r53G->Draw("ap"); 
  C->cd(3);
  //r63->Draw("ap"); 
  r63G->Draw("ap"); 

}
