void optZ(int nJ){

  cout << "nJets: " << nJ << endl; 
  TFile *fMC=new TFile("/uscms_data/d3/btcarlso/TEST_DIR/CMSSW_5_3_8_patch3/src/SusyAnalysis/SusyNtuplizer/macro/output_file_allMC.root","READ");
  TFile *fsig=new TFile("/uscms_data/d3/btcarlso/TEST_DIR/CMSSW_5_3_8_patch3/src/SusyAnalysis/SusyNtuplizer/macro/output_file_stealth_300_200.root","READ"); 

  TH1F *mMC=(TH1F*)fMC->Get(Form("h_mumu_nJet%d_0p5GeV",nJ)); 
  TH1F *msig=(TH1F*)fsig->Get(Form("h_mumu_nJet%d_0p5GeV",nJ));

  double sig=0;
  double bkg=0; 
  double BW=msig->GetBinWidth(1);
  double bkgTot=0;

  int N=30;
  TGraph *bkg_frac=new TGraph(N);
  TGraph *signf = new TGraph(N); 

  for(int i=1; i<=mMC->GetNbinsX();i++)bkgTot+=mMC->GetBinContent(i);
  cout << "BW: " << BW << endl; 
  int iZ=msig->FindBin(91); 
  for(int ii=0; ii<N; ii++){
    for(int i=iZ-ii; i<=iZ+ii;i++){
      //cout << msig->GetBinCenter(i) << endl; 
      sig+=msig->GetBinContent(i);
      bkg+=mMC->GetBinContent(i); 
    }
    //cout << Form("M: %.1f-%.1f ",mMC->GetBinCenter(iZ-ii)-BW/2,mMC->GetBinCenter(iZ+ii)+BW/2) << bkg/bkgTot << " " <<sig/(sig+bkg)<< endl; 
    bkg_frac->SetPoint(ii,BW*ii,bkg/(sig+bkg));
    signf->SetPoint(ii,BW*ii,sig/(sig+bkg));
  }

  bkg_frac->Draw("ap");
  signf->Draw("p same");
  bkg_frac->RemovePoint(0);

  int iMax=TMath::LocMax(bkg_frac->GetN(),bkg_frac->GetY()); 
  cout << "iMax: " << iMax << endl; 
  cout << "maxX: " << bkg_frac->GetX()[iMax] << endl; 
}
