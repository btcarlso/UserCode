void compute_limit(){

  float  Nbkg[]={1415,730,384.8,205.3,121.6,68.1,44.7,28,18.6,9.32,6.53,3.88,1.47,0.83,0.383,0.383,0.383,0.383,0.383}; 
  float  NbkgSys[]={45,16,9.3,5.5,4.8,2.7,2,1.3,1.3,0.87,0.85,0.67,0.43,0.29,0.031,0.031,0.031,0.031,0.031}; 
  float  MLQ[]={300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000,1050,1100,1150,1200}; 

  int N=sizeof(MLQ)/sizeof(float);
  cout << "N :" << N << endl; 

  for(int mass=300; mass<=1000; mass+=100){

    TFile *f=new TFile(Form("output_file_RPV_LQD221_M%d.root",mass),"READ"); 
  TH1F *h=(TH1F*)f->Get("h_mumujj_LQ"); 
  int IM=1; 

  float fom[19]; 
  float tmp=0; 
  for(int im=0; im<N; im++){
    float sig=h->GetBinContent(IM);
    if(MLQ[im]<1000)IM++;
    //cout << MLQ[im] << " sig: " << sig << " sig/sqrt(B): " << sig/TMath::Sqrt(Nbkg[im]) << endl; 
     fom[im]=sig/TMath::Sqrt(Nbkg[im]); 

  }
  int iOpt=TMath::LocMax(N,fom); 
  cout << "MLQ: " << MLQ[iOpt] << endl; 
    cout << mass << " " <<  h->GetBinContent(iOpt+1) <<  " " << Nbkg[iOpt] << " " << 1+NbkgSys[iOpt]/Nbkg[iOpt] << endl; 
  delete h; 
  delete f; 

  }

}
