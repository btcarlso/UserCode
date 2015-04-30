
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

double sRTB(double sig, double bkg){
  if(bkg<=0)bkg=0.0000001; 
  double sys=0.5*bkg; 
  sys=0; 
  return sig/TMath::Sqrt(bkg+sys*sys);
}
TFile *out; 
void sensitivity(){
  out=new TFile("optimization_out.root", "RECREATE");
  out->cd();

  for(int M=300; M<=1000; M+=100){
    if(M!=900) continue;
    sensitivity_bin(5,M); 
    sensitivity_bin(6,M);
    sensitivity_bin(7,M);
  }

}

void tokenize(TString text, double *x){
  if(text=="") return; 
  TObjArray *a = (TObjArray*)text.Tokenize(">,"); 
  int N=6; 
  int ii=0; 

  for(int i=0; i<N; i++){
    TObjString *tmp = (TObjString*)a->At(i); 
    TString tmp_string = tmp->GetString(); 
    if(i==1 || i==3 || i==5){
      x[ii]=tmp_string.Atof();
	ii++;
    }
  }
}

void sensitivity_bin(int nJ, int M){
  double sys=0.5; 

  TFile *fbkg=new TFile("output_file_allMC.root","READ"); 
  TFile *fsig=new TFile(Form("output_file_RPV_LQD221_M%d.root",M),"READ");
  
  TString name=Form("h_optimization_nJetsg%d_g1b",nJ); 
  if(nJ<7)name=Form("h_optimization_nJets%d_g1b",nJ); 

  TH1F *hbkg=(TH1F*)fbkg->Get(name); 
  TH1F *hsig=(TH1F*)fsig->Get(name); 


  int N=hbkg->GetNbinsX(); 

  std::vector<double> zbi; 
  std::vector<double> sig; 
  std::vector<double> bkg; 

  std::vector<double> sig_2; 
  std::vector<double> bkg_2; 
  std::vector<double> zbi_2; 

  std::vector<double> zbi_3; 
  std::vector<double> sig_3; 
  std::vector<double> bkg_3; 

  double sqrtb[5000];
  
  for(int i=1; i<4000; i++){
    TString BL=hbkg->GetXaxis()->GetBinLabel(i);
    TString BL2=hsig->GetXaxis()->GetBinLabel(i); 
    if(BL!=BL2) cout << "ISSUE:!!" << endl; 

  }
  TH2F *hMll_st = new TH2F(Form("hMll_st_M%d_nJ%d",M,nJ),";S_{T}^{min} (GeV); M_{ll}^{min} (GeV)",18,250,2050,17,55-12.5,450+12.5); 
  TH2F *hMmuj_st = new TH2F(Form("hMmuj_st_M%d_nJ%d",M,nJ),";S_{T}^{min} (GeV);M_{#muj}^{min} (GeV)",18,250,2050,11,-50,1050);
 
  for(int i=1; i<hbkg->GetNbinsX(); i++){
    TString BL=hbkg->GetXaxis()->GetBinLabel(i);
    
    double bins[3]={0,0,0}; 
    tokenize(BL,bins); 

    double sig_=hsig->GetBinContent(i); 
    double bkg_=hbkg->GetBinContent(i); 
    double zbi_=getZBI(sig_,bkg_,bkg_*sys);
    double srtb_=sRTB(sig_,bkg_); 
   
    

    if(zbi_>0){
      if(bins[0]>0 && bins[1]>0 && bins[2]==0){
	//cout << hMll_st->FindBin(bins[0],bins[1]) << endl; 
	hMll_st->Fill(bins[0],bins[1],zbi_); 
      }
       if(bins[0]>0 && bins[1]==105){
	//cout << hMll_st->FindBin(bins[0],bins[1]) << endl; 
	hMmuj_st->Fill(bins[0],bins[2],zbi_); 
      }

       zbi.push_back(zbi_); 
       sig.push_back(sig_); 
       bkg.push_back(bkg_); 

      if(BL.Contains("M>105") && BL.Contains("Mmuj>0") ){
	//cout << "3: " << BL << " " << sig << " " << bkg << " " << zbi_ << endl; 
	zbi_3.push_back(zbi_); //ZBI optimized for Mmuj>0 and M>100
	sig_3.push_back(sig_); 
	bkg_3.push_back(bkg_); 
      }
      if(BL.Contains("M>105")){
	//cout << "2: " << BL << " " << sig << " " << bkg << " " << zbi_ << endl; 
	zbi_2.push_back(zbi_); //ZBI optimized only for M>100
	sig_2.push_back(sig_);
	bkg_2.push_back(bkg_); 
      }
    }
    sqrtb[i]=srtb_;

  }
  //int maxI=TMath::LocMax(N,zbi); 
  //double zbiMax=TMath::MaxElement(N,zbi); 
  //  TString optBinTitle=hbkg->GetXaxis()->GetBinLabel(maxI); 

  double *zbi_a = &zbi[0];
  double *zbi_2a = &zbi_2[0]; 
  double *zbi_3a = &zbi_3[0]; 
  int N1=zbi.size(); 
  int N2=zbi_2.size(); 
  int N3=zbi_3.size(); 

  int maxI=TMath::LocMax(N1,zbi_a);
  double zbiMax=TMath::MaxElement(N1,zbi_a);

  int maxI_2=TMath::LocMax(N2,zbi_2a); 
  double zbiMax_2=TMath::MaxElement(N2, zbi_2a); 
  TString optBinTitle_2=hbkg->GetXaxis()->GetBinLabel(maxI_2); 

  int maxI_3=TMath::LocMax(N3,zbi_3a); 
  double zbiMax_3=TMath::MaxElement(N3, zbi_3a); 
  TString optBinTitle_3=hbkg->GetXaxis()->GetBinLabel(maxI_3); 

  //cout << Form("%d & %d & %.1f & %.1f & %.1f \\\\", M, nJ, zbiMax, zbiMax_2, zbiMax_3 ) << endl; 

  cout << Form("%.2f \t %.2f \t",sig_2.at(maxI_2),bkg_2.at(maxI_2)); 


  int maxISQRTB=TMath::LocMax(N,sqrtb); 
  double sqrtbMax=TMath::MaxElement(N,sqrtb); 
  /*
  TString stopt="St>1400"; 
  TString Mopt="M>125"; 
  TString Mmujopt="Mmuj>0"; 

  TH1F *hMuJ=new TH1F(Form("hMuJ_M%d_nJ%d",M,nJ),";M(#mu,j) (GeV);Z_{BI}",11,-50,1050); 

  float MuJ=0; 
  for(int i=1; i<hbkg->GetNbinsX(); i++){
    TString binL=hbkg->GetXaxis()->GetBinLabel(i);
    double sig_=hsig->GetBinContent(i);
    double bkg_=hbkg->GetBinContent(i);
    double zbi_=getZBI(sig,bkg,bkg*sys);

    if(binL.Contains(stopt)==0) continue; 
    if(binL.Contains(Mopt)==0) continue;
    hMuJ->Fill(MuJ,zbi_); 
    MuJ+=100.0; 
    //cout << binL << " " << zbi_ << endl; 
  }
  //hMuJ->Draw("histo"); 

  TH1F *hSt=new TH1F(Form("hSt_M%d_nJ%d",M,nJ),";S_{T}^{min} (GeV);Z_{BI}",18,250,2050);
  float stTmp=300;
  for(int i=1; i<hbkg->GetNbinsX(); i++){
    TString binL=hbkg->GetXaxis()->GetBinLabel(i);
    double sig=hsig->GetBinContent(i);
    double bkg=hbkg->GetBinContent(i);
    double zbi_=getZBI(sig,bkg,bkg*sys);

    if(binL.Contains(Mopt)==0) continue;
    if(binL.Contains(Mmujopt)==0) continue;
    if(zbi_>0)hSt->Fill(stTmp,zbi_);
    stTmp+=100.0;
  }
  //cout << "ST optimized, Mll>125, M: " << M << " StMax: " << hSt->GetBinCenter(hSt->GetMaximumBin()) << " " << hSt->GetMaximum() << endl;  
 

  TH1F *hdM = new TH1F(Form("hdM_M%d_nJ%d",M,nJ),";M(#mu#mu) (GeV); Z_{BI}",15,50-12.5,400+12.5); 

  float dM=50; 
  for(int i=1; i<hbkg->GetNbinsX(); i++){
    TString binL=hbkg->GetXaxis()->GetBinLabel(i);
    double sig=hsig->GetBinContent(i);
    double bkg=hbkg->GetBinContent(i);
    double zbi_=getZBI(sig,bkg,bkg*sys);

    if(binL.Contains(stopt)==0) continue;
    if(binL.Contains(Mmujopt)==0) continue;
    if(zbi_>0)hdM->Fill(dM,zbi_);
    //cout << dM << " " << zbi_  << endl; 
    dM+=25.0;

    //cout << binL << " " << zbi_ << endl;
  }
  */
  //cout << "--------------------" << endl; 

  //
  out->cd();
  hMll_st->Write();
  hMmuj_st->Write();

  gROOT->SetBatch();

  TString C1Name=hMll_st->GetName(); 
  C1Name="c_"+C1Name;
  TString C2Name=hMmuj_st->GetName(); 
  C2Name="c_"+C2Name; 

  TCanvas *c1 = new TCanvas(C1Name,"",800,600); 
  TCanvas *c2 =new TCanvas(C2Name,"",800,600);

  c1->cd();
  hMll_st->Draw("COLZ");
  c2->cd(); 
  hMmuj_st->Draw("COLZ"); 

  c1->Print(C1Name+".pdf"); 
  c2->Print(C2Name+".pdf"); 
  c1->Write();
  c2->Write();
  //hSt->Write();
  //hMuJ->Write();
  //hdM->Write();

}
