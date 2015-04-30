void mass_cnt(){
  TH1F *h=(TH1F*)_file0->Get("h_mumu_nJet2_0btag");
  TH1F *hs=(TH1F*)_file1->Get("h_mumu_nJet2_0btag");
  TH1F *hsb=(TH1F*)h->Clone("hsb");
  hsb->Add(hs);

  double Z=0;
  double A1=0;
  double A2=0;

  double tot=0;

  double Zsig=0;
  double A1sig=0;
  double A2sig=0;

  double totsig=0;
  TCanvas *C = new TCanvas();
  C->Divide(1,2);
  C->cd(1);
  h->Draw("histo");
  hs->SetLineColor(kRed);
  hs->Draw("histo same");
  C->cd(2); 
  hsb->Divide(h);
  hsb->Draw("histo");


  for(int i=1;i<=h->GetNbinsX();i++){
    double m=h->GetBinCenter(i);
    if(m>50 && m<=85) A1+=h->GetBinContent(i);
    if(m>85 && m<=105)Z+=h->GetBinContent(i);
    if(m>105 && m<=250)A2+=h->GetBinContent(i);
    tot+=h->GetBinContent(i);
  }

  for(int i=1;i<=h->GetNbinsX();i++){
    double m=h->GetBinCenter(i);
    if(m>50 && m<=85) A1sig+=hs->GetBinContent(i);
    if(m>85 && m<=105)Zsig+=hs->GetBinContent(i);
    if(m>105 && m<=250)A2sig+=hs->GetBinContent(i);
    totsig+=hs->GetBinContent(i);
  }

  cout << "A1: " << A1/tot << " Z " << Z/tot << " A2: " << A2/tot << endl; 
  cout << "A1: " << A1sig/totsig << " Z " << Zsig/totsig << " A2: "<< A2sig/totsig << endl; 
  
  A1=A1+A1sig;
  A2=A2+A2sig;
  Z=Z+Zsig;
  tot=tot+totsig;

  cout << "A1: " << A1/tot << " Z " << Z/tot << " A2: " << A2/tot << endl;
}
