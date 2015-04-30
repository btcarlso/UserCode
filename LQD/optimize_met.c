
double getZBI(double sig, double bkg, double bkge){
  double obs=sig+bkg; //obs = signal + bkg
  // obs = expected number of _observed_ events (i.e. signal+background)
  // back = expected number of background events
  // backe = uncertainty on the expected number of background events
  if(bkg<=0) {
    //cout << "bkg must be greater than 0." << endl;
    return 0;
  }
  if(obs==0) return 0;

  double tau = bkg/(bkge*bkge);
  double n_off = tau*bkg;
  double P_BI = TMath::BetaIncomplete(1./(1.+tau), obs, n_off+1);
  double ZBI=sqrt(2.0)*TMath::ErfcInverse(2*P_BI);

  return ZBI;
}



void optimize_met(){
  TFile *fbkg = new TFile("output_file_allMC.root","READ"); 
  TFile *fsig = new TFile("output_file_RPV_LQD221_M600.root","READ"); 
  

  //TH1F *hbkg = (TH1F*)fbkg->Get("h_metmax"); //("h_dimuomassmin");
  //TH1F *hsig = (TH1F*)fsig->Get("h_metmax");//("h_dimuomassmin"); 
  TH1F *hbkg = (TH1F*)fbkg->Get("h_dimuomassmin");
  TH1F *hsig = (TH1F*)fsig->Get("h_dimuomassmin");

  TH1F *h_sb = (TH1F*)hbkg->Clone("h_ssb");
  h_sb->GetXaxis()->SetTitle("#slash{E_{T}}^{max}"); 
  // h_sb->GetXaxis()->SetTitle("M(ll)^{min}");
  h_sb->GetYaxis()->SetTitle("Z_{BI}"); 
  h_sb->SetMinimum(0);

  for(int i=1; i<=hbkg->GetNbinsX(); i++){
    double bkg = hbkg->GetBinContent(i); 
    double sig = hsig->GetBinContent(i); 
    double sens =0; 
    // if(bkg>0)sens=sig/TMath::Sqrt(bkg); 
    sens = getZBI(sig,bkg,bkg*0.5); 
    h_sb->SetBinContent(i,0);
    h_sb->SetBinContent(i,sens); 
    
  }
  h_sb->Draw("histo");
  int Iopt=h_sb->GetMaximumBin(); 
  cout << "optimum Met(max): " << h_sb->GetBinCenter(Iopt) << endl; 
  cout << "ZBi at optimum: " << h_sb->GetBinContent(Iopt) << endl; 
  cout << "used: " << h_sb->GetBinContent(h_sb->FindBin(100)) << endl; 

  cout << "bkg optimized: "<< hbkg->GetBinContent(Iopt) <<  " sig: " << hsig->GetBinContent(Iopt) << endl; 
  cout << "bkg 100: " << hbkg->GetBinContent(h_sb->FindBin(130)) << " sig: " << hsig->GetBinContent(h_sb->FindBin(130)) << endl; 


}
