void read_file(){
  ifstream in;
  //in.open("CopperPlateSim.csv");
  in.open("CopperPlateTemperaturewithHeatersPlaced.csv"); 
  Float_t x,y,z;
  
  Int_t nlines = 0;
  int Np=0; 
  TFile *f = new TFile("histogramFile.root","RECREATE");
  //  TNtuple *ntuple = new TNtuple("ntuple","data from ascii file","x:y:z");
  TGraph2D *gr=new TGraph2D();
  gr->SetName("OutputGraph"); 
  char str1[1]; 
  char str2[1]; 
  string tmp; 
  while (!in.eof()) {
    x=0;
    y=0;
    z=0; 
    getline(in,tmp); 
    if(tmp==",,,") {
      cout << tmp << " break: " << endl;
      break; 
    }
    if(Np > 151820){
      cout << "more than L151821 points: break " << endl; 
      break; 
      }
    Np++;
    if(Np%30!=0)continue;
    TString tmp_=tmp; 
    TObjArray *xyz= (TObjArray*)tmp_.Tokenize(","); 
    for(int i=0; i<4; i++){
      TObjString *tmpString = (TObjString*)xyz->At(i); 
      TString tmpi=tmpString->GetString();
      if(i==0)x=tmpi.Atof(); 
      if(i==1)y=tmpi.Atof();
      //if(i==2)z=tmpi.Atof();
      if(i==3)z=tmpi.Atof(); 
    }
    /* cout << "nlines: " << nlines << endl; 
    cout << tmp << endl; 
    printf("x=%8f, y=%8f, z=%8f\n",x,y,z);*/
    if(nlines < 5) cout << tmp << endl; 
    if (nlines < 5) printf("x=%8f, y=%8f, z=%8f\n",x,y,z);
    //    h->Fill(x,y,TMath::Abs(z)); 
    gr->SetPoint(nlines,x,y,z); 
    
    nlines++;
  }
  cout <<" Sucessfully read file: " << endl; 
  cout << "GetN:" << gr->GetN() << endl; 
  /*
  gr->Draw("COLZ"); 

  //cout << "Eval(0.27,1.29): " << gr->Interpolate(0.27,1.29) << endl; 
  double x_[]={1.298702,1.13411,0.9692640.805942,0.639826,0.47498,0.394462,0.689356,0.860298,1.174242,1.25476,1.338834,1.475486};
  double y_[]={-0.27432,-0.22987,-0.185674,-0.14224,-0.097282,-0.053086,-0.024638,-0.030734,-0.076708,-0.081026,-0.109474,-0.125222,0.088646};

  double To=gr->Interpolate(0.27,1.29); 
  for(int i=0; i<13; i++){
    cout << Form("%.2f,%.2f",y_[i],x_[i]) <<  gr->Interpolate(-y_[i],x_[i])-To << endl; 
    //cout <<  gr->Interpolate(-y_[i],x_[i])-To << ","; 
  }

  for(int i=0; i<13; i++){
    cout <<  gr->Interpolate(-y_[i],x_[i])-To << ",";                                                                                                    
  }


  */

  double x_[]={48.1314,45.0613,42.7596,36.8930,33.8229,33.5212,25.6546,22.5845,20.2828,50.6854,49.2437,54.6154,39.4470,38.0053,43.3770,28.2026,26.7668,32.1386,16.9702,15.5284,20.9002,55.7305,57.8876,56.5073,45.9338,44.4921,49.8639,34.6954,33.2537,38.6255,23.4570,22.0153,27.3871,52.4207,50.9790,56.3507,48.1823,39.7406,45.1123,29.9439,28.5022,33.8739,18.7050,19.5650,22.8135,47.6691,46.2274,51.5992,36.4307,34.9890,40.3608,25.1918,26.0519,29.3003,54.1560,52.7143,58.0860,42.9176,41.4759,46.8476,31.6787,32.5387,35.7872,49.4044,47.9627,53.3345,38.1655,39.0256,42.2740,54.4496,56.6067,55.8913,44.6524,45.5125,48.7609,51.1392,51.9993,55.2478,52.2447,38.2656};
  double y_[]={11.2134,11.9317,9.7768,8.2078,8.9260,6.7711,5.2021,5.9204,3.7655,13.4163,8.0459,3.4826,10.4107,5.0403,6.4769,7.4050,2.0346,3.4713,4.3994,-0.9710,0.4656,6.3045,8.5430,11.4971,8.6693,3.2988,4.7355,5.6360,0.2932,1.7298,2.6580,-2.7125,-1.2758,6.9278,1.5574,2.9941,3.9222,-1.4482,-0.0116,0.9165,-4.4539,-3.0172,-2.0910,-5.0783,-5.3600,2.1808,-3.1896,-1.7530,0.8249,-6.1953,-4.7587,-3.8324,-6.8198,-7.1014,0.4393,-4.9311,-3.4940,-2.5663,-7.9367,-6.5001,-5.5739,-8.5612,-8.8429,-4.3077,-9.6781,-8.2415,-7.3153,-10.3026,-10.5843,-11.4196,-9.1199,-6.0492,-9.0567,-12.0440,-12.3257,-10.9781,-13.7855,-14.0671,0.0813,-0.1811}; 
    double T_exp[]={-22.90,-22.68,-23.11,-23.46,-23.27,-23.35,-23.91,-23.83,-24.21,-22.80,-22.46,-22.45,-23.14,-23.26,-22.84,-23.68,-24.00,-23.51,-24.08,-23.98,-23.98,-22.18,-22.26,-23.09,-23.11,-23.33,-22.98,-23.50,-23.11,-23.04,-24.08,-23.80,-23.51,-22.44,-22.12,-22.50,-23.39,-22.18,-22.07,-23.39,-23.50,-23.26,-23.88,-23.54,-23.42,-22.79,-22.93,-21.89,-23.10,-23.06,-22.63,-23.69,-23.22,-23.14,-21.87,-22.96,-22.73,-22.91,-22.71,-22.06,-23.21,-22.70,-22.78,-22.94,-22.40,-21.59,-22.18,-22.23,-22.37,-21.70,-21.94,-22.29,-22.66,-21.92,-21.75,-22.47,-21.42,-21.70,-22.00,-22.71}; 
    /* 
  double T_exp[]={19.90,19.80,19.90,19.80,19.80,20.00,19.90,19.90,19.90,19.70,20.00,19.90,19.90,19.90,19.90,20.00,20.00,20.10,19.90,19.90,19.90,19.90,19.90,20.00,19.90,19.90,20.00,20.00,20.10,20.00,20.00,20.00,20.00,20.00,20.00,20.10,19.90,20.00,20.10,20.00,20.00,20.10,20.10,20.00,20.10,19.80,20.00,20.10,19.80,19.90,19.90,20.00,19.90,20.00,20.00,19.70,20.00,20.00,20.00,20.00,20.10,20.10,20.10,20.10,20.00,20.10,20.00,20.10,20.00,20.10,19.90,20.10,20.00,19.90,20.10,19.90,20.00,20.00,19.90,20.10};//RT
  
  double T_exp[]={-22.88,-22.66,-23.08,-23.42,-23.24,-23.31,-23.85,-23.78,-24.13,-22.78,-22.39,-22.36,-23.11,-23.18,-22.75,-23.64,-23.97,-23.44,-24.01,-23.92,-23.92,-22.03,-22.14,-23.03,-23.06,-23.12,-22.70,-23.44,-23.08,-22.98,-24.02,-23.76,-23.52,-22.34,-22.02,-22.36,-23.30,-22.03,-22.03,-23.39,-23.46,-23.29,-23.83,-23.50,-23.39,-22.68,-23.06,-21.90,-23.10,-23.07,-22.80,-23.66,-23.21,-23.12,-21.79,-23.18,-22.79,-22.98,-22.85,-22.27,-23.21,-22.72,-22.81,-23.12,-22.76,-21.97,-22.18,-22.33,-22.55,-22.37,-22.40,-22.54,-22.88,-22.22,-22.18,-22.96,-22.01,-22.39,-21.94,-22.71}; //2nd run all heat on ?
    double T_exp[]={-24.11,-24.11,-24.24,-24.47,-24.48,-24.27,-24.65,-24.70,-24.80,-24.00,-24.16,-24.07,-24.26,-24.62,-24.44,-24.56,-24.85,-24.64,-24.65,-24.63,-24.89,-23.70,-23.95,-23.95,-24.39,-24.42,-24.25,-24.60,-24.62,-24.58,-24.84,-24.74,-24.77,-24.10,-24.23,-23.92,-24.52,-24.05,-24.26,-24.65,-24.70,-24.52,-24.64,-24.67,-24.48,-24.46,-24.56,-24.14,-24.69,-24.50,-24.42,-24.73,-24.68,-24.44,-24.07,-24.61,-24.21,-24.50,-24.38,-24.30,-24.48,-24.36,-24.30,-24.44,-24.18,-23.98,-23.57,-24.16,-24.21,-23.46,-23.98,-23.98,-24.30,-24.14,-23.78,-24.12,-23.79,-23.57,-24.31,-24.48};
   */
    int Io=75;//normally 75
  double xRef=y_[Io]*2.54/100; 
  double yRef=x_[Io]*2.54/100; 
  double TRef=gr->Interpolate(xRef,yRef); 
  double TRefData=T_exp[Io]; 
  double TDataCenter = T_exp[38]; 
  int Ndata=sizeof(x_)/sizeof(double); 

  TH2F *hData = new TH2F("data","Data: T-To #circC;X(m);Y(m);",25,-.5,.5,50,0.15,1.65); 
  TH2F *hSim = new TH2F("hSim","Sim: T-To #circC;X(m);Y(m)",25,-.5,.5,50,0.15,1.65); 
  TH2F *hRatio = new TH2F("hRatio","(T-To)^{Sim}-(T-To)^{Data};X(m);Y(m)",25,-.5,.5,50,0.15,1.65); 
  TH2F *hDataCenter = new TH2F("data","Data: T-To(center) #circC;X(m);Y(m);",25,-.5,.5,50,0.15,1.65);

  for(int i=1; i<=hData->GetNbinsX(); i++){
    for(int iy=1; iy<=hData->GetNbinsY(); iy++){
      hData->SetBinContent(i,iy,-100);
      hSim->SetBinContent(i,iy,-100);
      hRatio->SetBinContent(i,iy,-100);
      hDataCenter->SetBinContent(i,iy,-100); 
  }
  }

  double chi=0; 

  for(int i=0; i<Ndata; i++){
    double xo=y_[i]*2.54/100; 
    double yo=x_[i]*2.54/100; 
    double Tsim=gr->Interpolate(xo,yo)-TRef; 
    cout << "<x,y>" << xo << yo << " T: " <<Tsim << " " << T_exp[i]-TRefData << endl; 
    int iBin=hData->FindBin(xo,yo); 
    hData->SetBinContent(iBin,T_exp[i]-TRefData); 
    hSim->SetBinContent(iBin,Tsim); 
    if(Tsim < 10)hRatio->SetBinContent(iBin,Tsim-(T_exp[i]-TRefData)); 
    hDataCenter->SetBinContent(iBin,T_exp[i]-TDataCenter); 
    chi+=TMath::Power(Tsim-(T_exp[i]-TRefData),2)/(0.5*0.5); 

  }

  cout << "chi2: " << chi << " Ndof: " << Ndata << endl; 
  gStyle->SetPalette(1); 
  gStyle->SetCanvasColor(10);
  gStyle->SetHistFillColor(10);
  gStyle->SetPaintTextFormat("4.2f"); 
  TCanvas *C1 = new TCanvas();
  TCanvas *C2 = new TCanvas();
  TCanvas *C3 = new TCanvas();
  TCanvas *C4 = new TCanvas(); 
  C1->cd();
  hData->SetMinimum(-2);
  hData->SetStats(kFALSE); 
  hData->DrawCopy("COLZ0 text1");
  C2->cd();
  hSim->SetMinimum(-2);
  hSim->SetStats(kFALSE); 
  hSim->DrawCopy("COLZ0 text1"); 
  C3->cd();
  hRatio->SetMinimum(-2);
  hRatio->SetStats(kFALSE);
  hRatio->DrawCopy("COLZ text1"); 
  C4->cd();
  hDataCenter->SetMinimum(-2.2);
  hDataCenter->SetStats(kFALSE); 
  hDataCenter->DrawCopy("COLZ text1"); 
  printf(" found %d data points\n",nlines);
  in.close();
  f->cd();
  gr->Write(); 
  hData->Write(); 
  hSim->Write();
  hRatio->Write(); 
  f->Close();
  C1->Print("Data.pdf");
  C2->Print("Sim.pdf");
  C3->Print("Difference.pdf");
  C4->Print("DeltaData.pdf");
}
