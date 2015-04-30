void opt_met(){

  TFile *Fsig=new TFile("/uscms_data/d3/btcarlso/TEST_DIR/CMSSW_5_3_8_patch3/src/SusyAnalysis/SusyNtuplizer/macro/output_file_stealth_300_200.root","READ"); 
  TFile *Fbkg=new TFile("/uscms_data/d3/btcarlso/TEST_DIR/CMSSW_5_3_8_patch3/src/SusyAnalysis/SusyNtuplizer/macro/output_file_allMC.root","READ");

  TH1F *metSig=(TH1F*)Fsig->FindObjectAny("met_1Mu_1El_0btag_SF"); 
  TH1F *metBkg=(TH1F*)Fbkg->FindObjectAny("met_1Mu_1El_0btag_SF");

  TH1F *s_bkg=(TH1F*)metSig->Clone("s_bkg"); 
  TH1F *s_bkg_metmax=(TH1F*)metSig->Clone("s_bkg_metmax");

  s_bkg->GetYaxis()->SetTitle("S/#sqrt{B}"); 

  int Nx=metSig->GetNbinsX();


  for(int i=1; i<=Nx; i++){
    double met=metSig->GetBinCenter(i); 
  
    double Smax=metSig->Integral(1,i); 
    double Bmax=metBkg->Integral(1,i); 

    double S=metSig->Integral(i,Nx); 
    double B=metBkg->Integral(i,Nx); 

    if(B>0)s_bkg->SetBinContent(i,S/TMath::Sqrt(B)); 
    else s_bkg->SetBinContent(i,0); 

    if(Bmax>0)s_bkg_metmax->SetBinContent(i,Smax/TMath::Sqrt(Bmax));
    else s_bkg_metmax->SetBinContent(i,0);


  }

  TCanvas *c = new TCanvas();
  c->cd(); 

  s_bkg->SetAxisRange(0,300); 

  s_bkg_metmax->SetLineColor(kRed);
  s_bkg_metmax->SetLineWidth(2); 
  s_bkg_metmax->SetLineStyle(kDashed); 

  TLegend *Leg = new TLegend(0.55,0.35,0.75,0.55);
  Leg->SetFillColor(0); 
  Leg->SetLineColor(0);
  Leg->SetBorderSize(0); 

  Leg->AddEntry(s_bkg,"#slash{E}_{T} > #slash{E}_{T}^{min}","l");
  Leg->AddEntry(s_bkg_metmax,"#slash{E}_{T} < #slash{E}_{T}^{max}","l");

  s_bkg->SetMaximum(25); 
  s_bkg->Draw("histo");
  s_bkg_metmax->Draw("histo same"); 
  TLatex txt; 
  txt.SetNDC(kTRUE); 
  txt.DrawLatex(0.20,0.8,"0 b-tag, #mu,e, M_{#tilde{q}}=300 GeV");

  Leg->DrawClone("same"); 
  c->Print("met_optimization.pdf"); 
}
