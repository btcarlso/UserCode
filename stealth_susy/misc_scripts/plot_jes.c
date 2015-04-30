int category_bin(int nJ, int ist){
  int nJetmax=7;
  int iBin=-1;
  if(ist==0) iBin=nJ;
  if(ist==1) iBin=nJ+nJetmax;
  if(ist==2) iBin=nJ+2*nJetmax;
  return iBin;
}

void compute_JES(TString sample, int ist){
  cout << sample << " ist: " << endl; 
  TString dir="/uscms_data/d3/btcarlso/TEST_DIR/CMSSW_5_3_8_patch3/src/SusyAnalysis/SusyNtuplizer/macro/"; 
  TFile *F = new TFile(dir+sample+".root","READ"); 
  TFile *FM = new TFile(dir+sample+"_jecM.root","READ"); 
  TFile *FP = new TFile(dir+sample+"_jecP.root","READ"); 

  TH1F *Ngen=(TH1F*)F->FindObjectAny("h_cutflow_table"); 
  TH1F *NgenP=(TH1F*)FP->FindObjectAny("h_cutflow_table");
  TH1F *NgenM=(TH1F*)FM->FindObjectAny("h_cutflow_table");

  TH1F *h=(TH1F*)F->FindObjectAny("Acceptance_uW"); //EventCategories_1Mu_1El_0btag_SF"); //"Acceptance_uW" 
  TH1F *hM=(TH1F*)FM->FindObjectAny("Acceptance_uW"); //"EventCategories_1Mu_1El_0btag_SF"); 
  TH1F *hP=(TH1F*)FP->FindObjectAny("Acceptance_uW"); //EventCategories_1Mu_1El_0btag_SF");

  int i1=category_bin(4,ist); 
  int i2=category_bin(7,ist); 

  TCanvas *c = new TCanvas(); 
  c->cd(); 

  h->SetAxisRange(i1,i2); 
  h->GetYaxis()->SetTitle("A");
  h->SetStats(kFALSE); 

  //  h->SetMaximum(0.001);
  h->Draw("histo");
  hP->SetLineColor(kRed); 
  hM->SetLineColor(kGreen); 
  hP->DrawClone("histo same"); 
  hM->DrawClone("histo same"); 
  TCanvas *c2=new TCanvas(); 
  c2->cd(); 

  for(int i=i1; i<=i2; i++){
    double a=h->GetBinContent(i); 
    double aP=hP->GetBinContent(i); 
    double aM=h->GetBinContent(i); 
    //    cout << 100*TMath::Sqrt(a*(1-a)/Ngen)/a << " " << 100*TMath::Sqrt(aP*(1-aP)/NgenP)/aP << " " << 100*TMath::Sqrt(aM*(1-aM)/NgenM)/aM << endl; 
  }

  hP->Divide(h); 
  hM->Divide(h); 
  hP->SetMaximum(1.1);
  hP->SetMinimum(0.9); 
  hP->GetYaxis()->SetTitle("A(JES^{#pm})/A"); 
  hP->GetXaxis()->SetTitle("n_{J}"); 
  hP->SetAxisRange(i1,i2);
  hP->Draw("histo"); 
  hM->Draw("histo same"); 
  c->Print(sample+"_jes_variation.pdf"); 
  c2->Print(sample+"_jes_varation_ratio.pdf"); 
  
  for(int i=i1; i<=i2; i++){
    //cout << "JES: " << i << " "<< TMath::Max(fabs(1-hP->GetBinContent(i)),fabs(1-hM->GetBinContent(i))) << endl; 
    cout << 100*TMath::Max(fabs(1-hP->GetBinContent(i)),fabs(1-hM->GetBinContent(i))) << " "; 

    //hP->GetBinContent(i) << " " << hM->GetBinContent(i) << endl; 
  }
  cout << endl; 
}

void plot_jes(){
  gROOT->SetBatch(); 
  compute_JES("output_file_stealth_300_200",0); 
  compute_JES("output_file_stealth_400_200",1);
  compute_JES("output_file_stealth_500_300",1);
  compute_JES("output_file_stealth_600_300",2);
  compute_JES("output_file_stealth_700_400",2);
  compute_JES("output_file_stealth_800_400",2);
  compute_JES("output_file_stealth_900_500",2);
}
