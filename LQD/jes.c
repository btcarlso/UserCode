void jes_bin(int nJ, float stmin, int M, double &eff_unc, double &unc){
  TFile *f = new TFile(Form("output_file_RPV_LQD221_M%d.root",M),"READ"); 
  TFile *fP = new TFile(Form("output_file_RPV_LQD221_M%d_jecP.root",M),"READ");
  TFile *fM = new TFile(Form("output_file_RPV_LQD221_M%d_jecM.root",M),"READ");
  TH1F *hevents=(TH1F*)f->Get("Nevents"); 
  double Ngen=hevents->GetBinContent(1); 
  if(nJ==-1){
    int nJJ_=5; 
    TH1F *ho=(TH1F*)f->Get(Form("h_nJets%d_st_uW",nJJ_));
    TH1F *hP=(TH1F*)fP->Get(Form("h_nJets%d_st_uW",nJJ_)); 
    TH1F *hM=(TH1F*)fM->Get(Form("h_nJets%d_st_uW",nJJ_));
    for(int nJ_=6; nJ_<=7; nJ_++){
      TH1F *ho_=(TH1F*)f->Get(Form("h_nJets%d_st_uW",nJ_));
      TH1F *hP_=(TH1F*)fP->Get(Form("h_nJets%d_st_uW",nJ_));
      TH1F *hM_=(TH1F*)fM->Get(Form("h_nJets%d_st_uW",nJ_));
      ho->Add(ho_);
      hP->Add(hP_);
      hM->Add(hM_); 
    }

  }
  else {
    nJJ_=nJ;
    TH1F *ho=(TH1F*)f->Get(Form("h_nJets%d_st_uW",nJJ_));
    TH1F *hP=(TH1F*)fP->Get(Form("h_nJets%d_st_uW",nJJ_));
    TH1F *hM=(TH1F*)fM->Get(Form("h_nJets%d_st_uW",nJJ_));

  }


  int i1=ho->FindBin(stmin); 
  int Nx=ho->GetNbinsX(); 

  double eff=ho->Integral(i1,Nx)/Ngen; 
  eff_unc=TMath::Sqrt(eff*(1-eff)/Ngen); 

  double uncP=TMath::Abs(1-hP->Integral(i1,Nx)/ho->Integral(i1,Nx) );
  double uncM=TMath::Abs(1-hM->Integral(i1,Nx)/ho->Integral(i1,Nx) ); 

  double effP=hP->Integral(i1,Nx)/Ngen;
    double effM=hM->Integral(i1,Nx)/Ngen;

  unc=TMath::Max(uncP,uncM);

  cout << "nJ: " << nJ <<  " stmin  " << stmin << " eff " << eff <<  " eff+: " << effP << " effM " << effM << " eff stat unc: " << eff_unc << endl; 
  //cout << "nJ: " << nJ << " stmin: " << stmin << " lept Pt " << unc << " stat " << eff_unc <<  endl; 
}


void jes(){
  int M=400; 
  //double stmin[]={325,475,325}; //300
  double stmin[]={525,575,625}; //400
  //double stmin[]={1425,1375,1375};  //900
  double jesunc=0; 
  double eff_unc=0; 

  jes_bin(-1,stmin[0],M,eff_unc,jesunc); 
  //jes_bin(6,stmin[1],M,eff_unc,jesunc);
  //jes_bin(7,stmin[2],M,eff_unc,jesunc);

}
