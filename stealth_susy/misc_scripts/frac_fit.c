#include<TFractionFitter.h>
#include<TH1F.h>


void frac_fit(){
  double Ndata []  = {315239.0, 40524.1, 8625.2};
  double NdataE[]  = {   579.8,   202.8,   95.4};
  double NW []     = {258150.0,    0.38,   14.8};
  double NWE[]     = {   502.8,    0.38,    3.1};
  double NZ []     = { 19473.9, 37494.8,  181.4};
  double NZE[]     = {    61.7,    75.4,    5.2};
  double Ntt []    = { 33156.3,  1637.8, 8214.0};
  double NttE[]    = {    32.9,     8.2,   18.3};
  double Nother [] = { 20934.8,   606.9,  470.8};
  double NotherE[] = {   494.4,     5.3,   10.7};

  TH1F *W     = new TH1F("W"    ,"W"    , 3, 0.5, 3.5);
  TH1F *Z     = new TH1F("Z"    ,"Z"    , 3, 0.5, 3.5);
  TH1F *tt    = new TH1F("tt"   ,"tt"   , 3, 0.5, 3.5);
  TH1F *other = new TH1F("other","other", 3, 0.5, 3.5);
  TH1F *data  = new TH1F("data" ,"data" , 3, 0.5, 3.5);
	
  TH1F *weightW  = new TH1F("W_weight" , "W_weight" , 3, 0.5, 3.5);	
  TH1F *weightZ  = new TH1F("Z_weight" , "Z_weight" , 3, 0.5, 3.5);	
  TH1F *weighttt = new TH1F("tt_weight", "tt_weight", 3, 0.5, 3.5);	

  for(int i=0; i<3; i++){
    double Nw_prime     = TMath::Power(NW[i]/NWE[i]        , 2); 
    double NZ_prime     = TMath::Power(NZ[i]/NZE[i]        , 2); 
    double Ntt_prime    = TMath::Power(Ntt[i]/NttE[i]      , 2); 
    //double Nother_prime = TMath::Power(Nother[i]/NotherE[i], 2); 
    W ->SetBinContent(i+1, Nw_prime );
    Z ->SetBinContent(i+1, NZ_prime );
    tt->SetBinContent(i+1, Ntt_prime);

    weightW -> SetBinContent(i+1, NW[i] /W ->GetBinContent(i+1));   
    weightZ -> SetBinContent(i+1, NZ[i] /Z ->GetBinContent(i+1));   
    weighttt-> SetBinContent(i+1, Ntt[i]/tt->GetBinContent(i+1));   
    //    other   -> SetBinContent(i+1, Nother[i] );
    data    -> SetBinContent(i+1, Ndata[i]);
  }

  W ->SetLineColor(kRed  );
  Z ->SetLineColor(kBlue );
  tt->SetLineColor(kGreen);

	for(int i=1; i<=3; i++){
	//	cout << "data: " << data->GetBinContent(i) << " W " << W->GetBinContent(i) << " Z " << Z->GetBinContent(i) << " tt " << tt->GetBinContent(i) << endl; 	
	//	cout << "data weight=1, W weight "  << weightW->GetBinContent(i) << " Z weight: " << weightZ->GetBinContent(i) << " weight tt " << weighttt->GetBinContent(i) << endl; 
	}
	
  TObjArray *mc= new TObjArray(4);
  mc->Add(W);
  mc->Add(Z);
  mc->Add(tt);
  //  mc->Add(other);
  TFractionFitter *fit = new TFractionFitter(data,mc);
  //fit->GetFitter()->FixParameter(3);
	fit->SetWeight(0,weightW);
	fit->SetWeight(1,weightZ);
	fit->SetWeight(2,weighttt);
	//fit->Constrain(0,0,2);
	//fit->Constrain(1,0,2);
	//fit->Constrain(2,0,2);

  fit->Fit();
  double R[3]; 

  for(int i=0; i<3; i++){
    double x=0;
    double xE=0;
    fit->GetResult(i,x,xE);
    double N=0; 
    if (i==0) N = NW[i] ; 
    if (i==1) N = NZ[i] ;
    if (i==2) N = Ntt[i];
	  cout << "x " << x << endl; 
	  cout << " initial fraction: " << Ndata[i]/N << endl; 
    R[i]=x*(Ndata[i]/N);
    cout << "R: " << R[i] << " +/- " << (xE/x)*R[i] << endl; 

  }
}
