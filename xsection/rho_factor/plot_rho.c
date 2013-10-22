

#include "plot_rho.h"




double MC_shape(double *x, double *par){

 // return par[0]*rho_ptMC->Interpolate(x[0]);
	return 0;

}
bool excl=false; 
double SB_func(double *x, double *par){
if(excl && x[0]>9.0 && x[0]<10.6){
TF1::RejectPoint();
return 0;
}

return par[0]+par[1]*x[0]; 

}

void CreateCanvas(TString Name,TString Title, int x, int y ){
	TCanvas *createC = new TCanvas(Name, Title,x,y);
	CName[Name]=createC;
}


void eff(TH1F *P, TH1F *T, float m1,float m2, float &Np, float &NpE, float &Nt, float &NtE, bool constant, int can_){
	if(binning=="dR")CreateCanvas(Form("Sideband_Subtraction_%d",can_),"Sideband subtraction",600,700); 
	TH1F *passed =(TH1F*)P->Clone("passed");
	TH1F *total = (TH1F*)T->Clone("total"); 
	
	TF1 *SBfunc = new TF1("SBfunc",SB_func,8.5,11.5,2);
	TF1 *SBfuncP = new TF1("SBfuncP",SB_func,8.5,11.5,2);
	TF1 *SBfuncM = new TF1("SBfuncM",SB_func,8.5,11.5,2);
	excl=true;

	if(binning=="pt") CName["mass_slices"]->cd(can_);
	if(binning=="dR") CName["mass_DR"]->cd(can_);
	
	TH1F *passP=(TH1F*)passed->Clone("pass_sysP"); 
	TH1F *passM=(TH1F*)passed->Clone("pass_sysM"); 

	TH1F *totalP=(TH1F*)passed->Clone("total_sysP"); 
	TH1F *totalM=(TH1F*)passed->Clone("total_sysM"); 

	if(constant)SBfunc->FixParameter(1,0);
	passed->Fit("SBfunc","RLQ0");
	total->Fit("SBfunc","RLQ0");
	excl=false;  
		
	if(can_>0){
		passed->DrawCopy(); 
		SBfunc->SetLineColor(kBlue);
		SBfunc->DrawClone("same"); 
		total->DrawCopy("same"); 
		if(binning=="dR"){
			CName[Form("Sideband_Subtraction_%d",can_)]->cd();
			passed->DrawCopy();
			SBfunc->DrawClone("same");
			total->DrawCopy("same");
		}
	}
	passed->Add(SBfunc,-1); 
	if(binning=="pt") CName["mass_slices_sub"]->cd(can_);
	if(binning=="dR") CName["mass_DR_sub"]->cd(can_);
	if(can_>0)passed->DrawCopy(); 
	
	SBfuncP->SetParameters(SBfunc->GetParameter(0)+SBfunc->GetParError(0),SBfunc->GetParameter(1)+SBfunc->GetParError(1)); 
	SBfuncM->SetParameters(SBfunc->GetParameter(0)-SBfunc->GetParError(0),SBfunc->GetParameter(1)-SBfunc->GetParError(1)); 
	
	passP->Add(SBfuncP,-1);
	passM->Add(SBfuncM,-1); 
	
	excl=false; 
	total->Add(SBfunc,-1);
	if(can_>0)total->DrawCopy("same");
	
	SBfuncP->SetParameters(SBfunc->GetParameter(0)+SBfunc->GetParError(0),SBfunc->GetParameter(1)+SBfunc->GetParError(1)); 
	SBfuncM->SetParameters(SBfunc->GetParameter(0)-SBfunc->GetParError(0),SBfunc->GetParameter(1)-SBfunc->GetParError(1)); 
	
	totalP->Add(SBfuncP,-1);
	totalM->Add(SBfuncM,-1); 
	
	Np=passed->Integral(passed->FindBin(m1),passed->FindBin(m2));
	Nt=total->Integral(total->FindBin(m1),total->FindBin(m2));
	
	NpE=(TMath::Abs(passP->Integral(passP->FindBin(m1),passP->FindBin(m2))-Np)+TMath::Abs(passM->Integral(passP->FindBin(m1),passP->FindBin(m2))-Np))/2;
	NtE=(TMath::Abs(totalP->Integral(total->FindBin(m1),total->FindBin(m2))-Np)+TMath::Abs(totalM->Integral(total->FindBin(m1),total->FindBin(m2))-Np))/2;
	
	cout << "Total: " << Nt << " " << NtE << " Passed: " << Np << " " << NpE << endl; 
	
	delete passed; 
	delete total; 
	delete passP; 
	delete passM;
	delete totalP;
	delete totalM; 
	delete SBfunc;
	delete SBfuncP;
	delete SBfuncM;
}



void plot_rho(){
	gROOT->SetBatch();
  double pt_bins[]={10,20,30,40,50,60,70,100};
  //double dR_bins[]={0,0.25,0.3,0.35,0.4,0.5,0.8,1,1.5};
	double dR_bins[]={0,0.2,0.3,0.4,0.8,1.0,1.5};


  TFile *input = new TFile("output_rho.root","READ");
  TH1F *tot=(TH1F*)input->FindObjectAny("hpt_den");
  TH1F *num=(TH1F*)input->FindObjectAny("hpt_num");
  TH1F *num_uW=(TH1F*)input->FindObjectAny("hpt_num_uW");

  TH1F *mu1=(TH1F*)input->FindObjectAny("hpt_lead");
  TH1F *mu2=(TH1F*)input->FindObjectAny("hpt_trail");

  TH2F *m_pt=(TH2F*)input->FindObjectAny("mass_pt");
  TH2F *m_ptnum = (TH2F*)input->FindObjectAny("mass_pt_num");
  TH2F *m_ptnum_uW = (TH2F*)input->FindObjectAny("mass_pt_num_uW");

	
	TH2F *m_dR=(TH2F*)input->FindObjectAny("mass_y0_DeltaRPtE_denBin");
	TH2F *m_dRnum = (TH2F*)input->FindObjectAny("mass_y0_DeltaRPtE_numBin");
	TH2F *m_dRnum_uW = (TH2F*)input->FindObjectAny("mass_y0_DeltaRPtE_numBin_uW");
	

//  TEfficiency *rho=(TEfficiency*)input->FindObjectAny("TEff_rho_pt");
//  rho->SetStatisticOption(4);
  mu1->SetStats(kFALSE);
  mu2->SetStats(kFALSE);
  mu2->SetLineColor(kRed);

  //  for(int i=1; i<=num->GetNbinsX();i++)tot->SetBinError(i,0);
  num->Divide(tot);
  TGraphAsymmErrors *rho_graph = new TGraphAsymmErrors(num); 
/*
  for(int i=1; i<=num_uW->GetNbinsX();i++){
    double diff=tot->GetBinContent(i)-num_uW->GetBinContent(i);
    //double err=1/TMath::Sqrt(tot->GetBinContent(i))+1/TMath::Sqrt(diff);
    double eff=rho->GetEfficiency(rho->FindFixBin(num->GetBinCenter(i)));
    double errL=rho->GetEfficiencyErrorLow(rho->FindFixBin(num->GetBinCenter(i)))/eff;
    double errH=rho->GetEfficiencyErrorUp(rho->FindFixBin(num->GetBinCenter(i)))/eff;
    
    cout << "eff: " << eff << " err: "<< errL << endl; 
    rho_graph->SetPointEYhigh(i-1,errH*num->GetBinContent(i));
    rho_graph->SetPointEYlow(i-1,errL*num->GetBinContent(i));
  }
*/
 

  TF1 *MC_shapefunc = new TF1("MC_shapefunc",MC_shape,20,100,1);
  //rho_graph->Fit(MC_shapefunc,"R0");
   cout << "Chi2: " << MC_shapefunc->GetChisquare() << " NDOF: " << MC_shapefunc->GetNDF() << endl; 

  TCanvas * mucan = new TCanvas();
  gPad->SetLogy();
  mu2->SetTitle("");
  mu2->Draw();
  mu1->Draw("Same");
  mucan->Print("Muon_distribution.pdf");
  cout << "Print slices: " << endl; 
  CreateCanvas("mass_slices","binned mass plots",800,600); 
  CreateCanvas("mass_slices_sub","binned mass plots bkg sub",800,600); 
	
  CreateCanvas("mass_DR","binned Mass in dR", 800,600);
  CreateCanvas("mass_DR_sub","binned mass in dR subtracted",800,600); 
	
	
  CName["mass_DR"]->Divide(2,3);
  CName["mass_DR_sub"]->Divide(2,3); 
  CName["mass_slices"]->Divide(3,3);
  CName["mass_slices_sub"]->Divide(3,3);
  
  TH1F *pt=(TH1F*)m_pt->ProjectionY();
  TEfficiency *rho_bkgsub = new TEfficiency("rho_bkgsub","#rho p_{T};p_{T}(#mu#mu) [GeV]; #rho",7,pt_bins);

	TH1F *tmp = new TH1F("tmp","tmp",7,pt_bins);
	TGraphAsymmErrors *rho_pt = new TGraphAsymmErrors(tmp); 
	rho_pt->SetName("rho_pt"); 
	
	TF1 *SBfunc = new TF1("SBfunc",SB_func,8.5,11.5,2);
	int can_=1; 
 	
	
for(int i=1;i<=pt->GetNbinsX(); i++){
	float T=0, P=0, TE=0, PE=0;	
	float Tw=0, Pw=0, TEw=0, PEw=0;	
	bool constant=false; 
	binning = "pt";

	
	cout << "Bin: " << i << endl; 

	float pT[]={pt->GetBinLowEdge(i),pt->GetBinCenter(i)+pt->GetBinWidth(i)/2 };   
	  
	excl=true;

	TH1F *_px=(TH1F*)m_pt->ProjectionX(Form("_px%d",i),i,i+1); 
	TH1F *_pxN=(TH1F*)m_ptnum->ProjectionX(Form("_pxN%d",i),i,i+1); 
	TH1F *_pxN_uW=(TH1F*)m_ptnum_uW->ProjectionX(Form("_pxN%d_uW",i),i,i+1); 


  	_pxN->SetLineColor(kRed);

	cout << "pT bin: " << pT[0] << "-" << pT[1] << endl; 
	_px->SetTitle(Form("%.0f<p_{T}<%.0f",pT[0],pT[1]));
	
	
	if(pt->GetBinCenter(i)>70)constant =true; 
	eff(_pxN_uW, _px,9.26,9.66,P,PE,T,TE,constant,0); 
	rho_bkgsub->SetTotalEvents(i,T); 
	rho_bkgsub->SetPassedEvents(i,P);
	
	eff(_pxN, _px,9.26,9.66,Pw,PEw,Tw,TEw,constant,can_); 
	
	float rho=Pw/Tw; 
	
	rho_pt->SetPoint(i,rho_pt->GetX()[i],rho); 
	float statH=rho_bkgsub->GetEfficiencyErrorUp(i)/rho_bkgsub->GetEfficiency(i);
	float statL=rho_bkgsub->GetEfficiencyErrorLow(i)/rho_bkgsub->GetEfficiency(i);
	
	float totH=rho_pt->GetY()[i]*TMath::Sqrt(TMath::Power(PEw/Pw,2)+TMath::Power(TEw/Tw,2)+statH*statH); 
	float totL=rho_pt->GetY()[i]*TMath::Sqrt(TMath::Power(PEw/Pw,2)+TMath::Power(TEw/Tw,2)+statL*statL);
	
	rho_pt->SetPointEYhigh(i,statH); 
	rho_pt->SetPointEYlow(i,statL); 
	
	cout << "rho: " << rho << " + " << totH << " - " << totL << endl;

	can_++;
}
	
	cout << "Binned in delta R: " << endl; 
	can_=1;
	TH1F *dR=(TH1F*)m_dR->ProjectionY();
	int Nbins=sizeof(dR_bins)/sizeof(double)-1;
	cout << "Nbins: " << Nbins << endl; 
	
	TEfficiency *rho_dR_uW = new TEfficiency("rho_dR_uW","#rho dR;; #rho",Nbins,dR_bins);
	
	TH1F *tmpR = new TH1F("tmpR","tmp",Nbins,dR_bins);
	TGraphAsymmErrors *rho_dR = new TGraphAsymmErrors(tmpR); 
	rho_dR->SetName("rho_dRPtE");
	
	TProfile *emu=(TProfile*)input->FindObjectAny("pr_eff_1S_y0"); 
	
	for(int i=2;i<=dR->GetNbinsX(); i++){
		float T=0, P=0, TE=0, PE=0;	
		float Tw=0, Pw=0, TEw=0, PEw=0;	

		bool constant=false; 
		
		binning="dR";
		
		cout << "Bin: " << i << endl; 
		
		float dR_[]={dR->GetBinLowEdge(i),dR->GetBinCenter(i)+dR->GetBinWidth(i)/2 };   
		
		excl=true;
		
		TH1F *_px=(TH1F*)m_dR->ProjectionX(Form("_dR%d",i),i,i+1); 
		TH1F *_pxN=(TH1F*)m_dRnum->ProjectionX(Form("_dRN%d",i),i,i+1); 
		TH1F *_pxN_uW=(TH1F*)m_dRnum_uW->ProjectionX(Form("_dRN%d_uW",i),i,i+1); 

		
		_pxN->SetLineColor(kRed);
		
		cout << "dR bin: " << dR_[0] << "-" << dR_[1] << endl; 
		_px->SetTitle(Form("%.2f<#DeltaR_{p_{T}}^{elliptic}<%.2f",dR_[0],dR_[1]));
		_pxN->SetTitle(Form("%.2f<#DeltaR_{p_{T}}^{elliptic}<%.2f",dR_[0],dR_[1]));
		
		eff(_pxN_uW, _px,9.26,9.66,P,PE,T,TE,constant,0); 
		rho_dR_uW->SetTotalEvents(i,T); 
		rho_dR_uW->SetPassedEvents(i,P);
		eff(_pxN, _px,9.26,9.66,Pw,PEw,Tw,TEw,constant,can_); 

		//float rho=Pw/Tw; 
		float E=emu->GetBinContent(i);
		cout << "emu: " << E << endl;
		
		float rho=(P/T)/E;
		
		if(i>0){
			rho_dR->SetPoint(i,rho_dR->GetX()[i],rho); 
			float statH=rho_dR_uW->GetEfficiencyErrorUp(i)/rho_dR_uW->GetEfficiency(i);
			float statL=rho_dR_uW->GetEfficiencyErrorLow(i)/rho_dR_uW->GetEfficiency(i);
			
			float totH=rho_dR->GetY()[i]*TMath::Sqrt(TMath::Power(PEw/Pw,2)+TMath::Power(TEw/Tw,2)+statH*statH); 
			float totL=rho_dR->GetY()[i]*TMath::Sqrt(TMath::Power(PEw/Pw,2)+TMath::Power(TEw/Tw,2)+statL*statL);
			
			rho_dR->SetPointEYhigh(i,statH); 
			rho_dR->SetPointEYlow(i,statL); 
		}
		
		
		can_++;
	}	
	TFile *out=new TFile("rho_data.root","RECREATE");
	out->cd();
	
	CName["mass_DR"]->Write();
	CName["mass_DR_sub"]->Write();
	
	rho_dR->Write();
	rho_pt->Write();

	
 rho_bkgsub->SetLineColor(kRed);
 TCanvas * c = new TCanvas();
	c->Divide(1,2);
	c->cd(1);
/*
  rho_ptMC->SetMinimum(0);
  rho_ptMC->GetXaxis()->SetTitle("p_{T} [GeV]");
  rho_ptMC->GetYaxis()->SetTitle("#rho");
  rho_ptMC->SetStats(kFALSE);
  rho_ptMC->SetTitle("");
  //cout << MC_shapefunc->GetChisquare() << " " << MC_shapefunc->GetNDF() << endl; 
  rho_ptMC->Draw("E1");
*/	
 // rho_graph->Draw("P same");
  rho_pt->Draw("P same");	
	
  TLegend *L = new TLegend(0.2,0.2,0.6,0.6);
  // L->AddEntry(rho_ptMC,"#rho MC","LE");
   L->AddEntry(rho_pt,"#rho bkg sub","LE");
   L->Draw("same");
	
	c->cd(2);
	
	rho_dR->Draw("ap");
	
	for (std::map<TString,TCanvas*>::iterator it=CName.begin(); it!=CName.end(); it++) {
		CName[it->first]->Write();
		CName[it->first]->Print("rho_figures/"+it->first+".png"); 
	}
	
	

  c->Print("rho_figures/rho.pdf");
	


}
