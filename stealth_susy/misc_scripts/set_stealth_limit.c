/*
 *  set_limit.c
 *
 *
 *  Created by Benjamin Carlson on 1/15/14.
 *  Copyright 2014 Carnegie Mellon University. All rights reserved.
 *
 */

#include "set_stealth_limit.h"

void set_stealth_limit(){
	gROOT->SetBatch();
	//GetHistograms();
	bookGraphs();
	create_theory();
	read_data();
	//draw_sensitivity(300);
	//draw_sensitivity(400);
	//draw_sensitivity(500);
	//draw_sensitivity(600);
    
	draw_plots();
	
	//model_Ind_Limit();
	//draw_model_independent();
	
	TFile *limit_plots = new TFile("limit_plots.root","RECREATE");
	for (std::map<TString,TCanvas*>::iterator it=CName.begin(); it!=CName.end(); it++) {
		CName[it->first]->Write();
		TString name="/uscms_data/d3/btcarlso/figtmp/";
		name=name+it->first+".pdf";
		CName[it->first]->Print(name);
        
	}
    
    
	
}

void CreateCanvas(TString Name,TString Title, int x, int y ){
	TCanvas *createC = new TCanvas(Name, Title,x,y);
	CName[Name]=createC;
}

void draw_sensitivity(int mass){
	CreateCanvas(Form("Limit_sensitivity_M%d",mass),"",600,600);
	CName[Form("Limit_sensitivity_M%d",mass)]->cd();
    
	hName[Form("h_sensitivity_M%d_exp",mass)]->SetLineStyle(kDashed);
	hName[Form("h_sensitivity_M%d_obs",mass)]->SetLineWidth(2);
    
	hName[Form("h_sensitivity_M%d_exp_1sigmaP",mass)]->SetLineStyle(kDashed);
	hName[Form("h_sensitivity_M%d_exp_1sigmaM",mass)]->SetLineStyle(kDashed);
	
	hName[Form("h_sensitivity_M%d_exp_1sigmaP",mass)]->SetLineColor(kRed);
	hName[Form("h_sensitivity_M%d_exp_1sigmaM",mass)]->SetLineColor(kRed);
    
	hName[Form("h_sensitivity_M%d_exp",mass)]->SetMinimum(0);
	double max=1.5*hName[Form("h_sensitivity_M%d_exp_1sigmaP",mass)]->GetMaximum();
	hName[Form("h_sensitivity_M%d_exp",mass)]->SetMaximum(max);
    
	
	hName[Form("h_sensitivity_M%d_exp",mass)]->Draw("histo");
	hName[Form("h_sensitivity_M%d_obs",mass)]->Draw("histo same");
	hName[Form("h_sensitivity_M%d_exp_1sigmaP",mass)]->Draw("histo same");
	hName[Form("h_sensitivity_M%d_exp_1sigmaM",mass)]->Draw("histo same");
    
	TLegend Leg(0.65,0.65,0.85,0.8);
	Leg.SetFillColor(10);
	Leg.SetBorderSize(0);
	Leg.SetLineColor(10);
	
	Leg.AddEntry(hName[Form("h_sensitivity_M%d_exp",mass)],"Expected","L" );
	Leg.AddEntry(hName[Form("h_sensitivity_M%d_exp_1sigmaP",mass)],"#pm 1 #sigma","L");
	Leg.AddEntry(hName[Form("h_sensitivity_M%d_obs",mass)],"Observed","L");
	Leg.DrawClone();
	
	TLatex txt;
	txt.SetNDC(kTRUE);
	txt.DrawLatex(0.62,0.85,Form("M_{#tilde{q}}=%d GeV",mass));
	
	/*
	 CreateGraph(Form("sensitivity_M%d_obs",mass) );
	 CreateGraph(Form("sensitivity_M%d_exp",mass) );
	 CreateGraph(Form("sensitivity_M%d_exp_1sigma",mass) );
	 CreateGraph(Form("sensitivity_M%d_exp_2sigma",mass) );
	 
     for(int i=0; i<grName[Form("sensitivity_M%d_exp",mass)]->GetN(); i++){
     cout << " x: " << grName[Form("sensitivity_M%d_exp",mass)]->GetX()[i] << " ";
     cout << " y: " << grName[Form("sensitivity_M%d_exp",mass)]->GetY()[i] << endl;
     
     }
     
     
     grName[Form("sensitivity_M%d_obs",mass)]->SetLineStyle(kDashed);
     grName[Form("sensitivity_M%d_obs",mass)]->SetLineColor(kBlack);
     grName[Form("sensitivity_M%d_exp",mass)]->SetFillStyle(1001);
     grName[Form("sensitivity_M%d_exp_1sigma",mass)]->SetFillColor(419);
     grName[Form("sensitivity_M%d_exp",mass)]->Draw("ap");
     
     int nJ=4;
     
     for(int i=0; i<15; i++){
     TString label=Form("%d",nJ);
     if(nJ==8)label=Form("#geq %d",nJ);
     grName[Form("sensitivity_M%d_exp",mass)]->GetXaxis()->SetBinLabel(i+1,label);
     
     nJ++;
     if(nJ>8)nJ=4;
     }
     
     //	grName[Form("sensitivity_M%d_exp",mass)]->GetXaxis()->SetRangeUser(4,24);
     //	grName["xs_limit_1sigma"]->Draw("a3");
     
     grName[Form("sensitivity_M%d_exp",mass)]->Draw("ap ");
     //grName["xs_limit_2sigma"]->Draw("a3 same");
     grName[Form("sensitivity_M%d_obs",mass)]->Draw("p same");
     */
    
}



void draw_plots(){
	CreateCanvas(Form("Limit_Plot"),"",600,600);
	CName["Limit_Plot"]->cd();
	gPad->SetLogy();
	
	grName["squark_xs"]->SetLineColor(kBlack);
    grName["squark_xs"]->SetLineStyle(kDotted);
    grName["squark_xs"]->SetFillColor(46);

	grName["xs_exp_combined"]->SetLineColor(kBlack);
	grName["xs_obs_combined"]->SetLineColor(kBlack);
	
	grName["xs_exp_combined"]->SetLineStyle(kDashed);
    
	grName["xs_exp_combined_1sigma"]->SetFillStyle(1001);
	grName["xs_exp_combined_2sigma"]->SetFillStyle(1001);
    
	grName["xs_exp_combined_1sigma"]->SetFillColor(3);
	grName["xs_exp_combined_2sigma"]->SetFillColor(5);
	
	grName["xs_exp_combined_1sigma"]->Draw("a3");
    grName["xs_exp_combined_2sigma"]->Draw("a3");
	grName["xs_exp_combined_1sigma"]->GetXaxis()->SetRangeUser(300,900);
    grName["xs_exp_combined_2sigma"]->GetXaxis()->SetRangeUser(300,900);
	grName["xs_exp_combined_1sigma"]->GetYaxis()->SetRangeUser(5,1.5*grName["squark_xs"]->Eval(300));
    grName["xs_exp_combined_2sigma"]->GetYaxis()->SetRangeUser(5,1.5*grName["squark_xs"]->Eval(300));

	grName["xs_exp_combined_1sigma"]->GetXaxis()->SetTitle("M_{#tilde{q}} (GeV)");
	grName["xs_exp_combined_1sigma"]->GetYaxis()->SetTitle("#sigma (fb)");
    
    grName["xs_exp_combined_2sigma"]->GetXaxis()->SetTitle("M_{#tilde{q}} (GeV)");
	grName["xs_exp_combined_2sigma"]->GetYaxis()->SetTitle("#sigma (fb)");
    
	//grName["xs_exp_combined_2sigma"]->Draw("a3 same");
    grName["xs_exp_combined_1sigma"]->Draw("a3");
	grName["xs_exp_combined"]->Draw("l same");
	grName["xs_obs_combined"]->Draw("l same");
    
    grName["squark_xs"]->Draw("3 same");
    grName["squark_xs"]->Draw("cx0 same");

    grName["xs_exp_combined"]->Draw("l same");
	grName["xs_obs_combined"]->Draw("l same");

	TLegend Leg(0.6,0.5,0.85,0.67);
	Leg.SetFillColor(10);
	Leg.SetBorderSize(0);
	Leg.SetLineColor(10);
	
	Leg.AddEntry(grName["xs_exp_combined"],"Expected","L" );
	Leg.AddEntry(grName["xs_obs_combined"],"Observed","L");
	Leg.AddEntry(grName["xs_exp_combined_1sigma"],"#pm 1 #sigma_{ex}","f");
   // Leg.AddEntry(grName["xs_exp_combined_2sigma"],"#pm 2 SD","f");

	Leg.AddEntry(grName["squark_xs"],"#sigma_{#tilde{q}#tilde{q}} #pm 1 #sigma_{th} ","LF");
	Leg.DrawClone();
	
	//TLatex txt(0.4,0.85,"M_{#tilde{#chi}_{1}^{#pm}}=1/2#timesM_{#tilde{q}}, Br(#tilde{#chi}_{1}^{#pm}#rightarrowW^{#pm},#tilde{S})=1");
    TLatex txt(0.6,0.85,"pp #rightarrow #tilde{q}_{L}#tilde{q}_{L}");
    //TLatex txt1(0.2,0.3,"M_{ #tilde{#chi}_{1}^{#pm} }=1/2#timesM_{ #tilde{q} }");
    TLatex txt1(0.2,0.3,"M_{ #tilde{#chi}_{1}^{#pm} }=1/2#timesM_{ #tilde{q} }");
	TLatex txt2(0.2,0.2,"M_{#tilde{S}} = 100 GeV, M_{S} = 90 GeV");

    //, Br(#tilde{#chi}_{1}^{#pm}#rightarrowW^{#pm},#tilde{S})=1");
    txt.SetTextSize(0.045);
	txt.SetNDC(kTRUE);
	txt.DrawClone();
    txt.DrawLatex(0.6,0.77,"#tilde{q}_{L}#rightarrow q #tilde{#chi}_{1}^{#pm}");
    txt.DrawLatex(0.6,0.69,"#tilde{#chi}_{1}^{#pm}#rightarrow W^{#pm} #tilde{S}");
    txt.DrawLatex(0.2,0.4,"#tilde{q}_{L}: #tilde{u}, #tilde{d}, #tilde{s}, #tilde{c}");
    
    txt1.SetTextSize(0.045);
	txt1.SetNDC(kTRUE);
	txt2.SetTextSize(0.045);
	txt2.SetNDC(kTRUE);
	txt1.DrawClone();
	txt2.DrawClone();
    //txt.DrawLatex(0.4,0.75,"M_{#tilde{S}}=100 GeV, M_{S}=90 GeV");
    
    CMS_lumi(CName["Limit_Plot"],2,10);
	//draw_header();
    gPad->RedrawAxis();

    
}


void
CMS_lumi( TPad* pad, int iPeriod, int iPosX )
{
    //https://ghm.web.cern.ch/ghm/plots/
    bool outOfFrame    = false;
    if( iPosX/10==0 )
    {
        outOfFrame = true;
    }
    int alignY_=3;
    int alignX_=2;
    if( iPosX/10==0 ) alignX_=1;
    if( iPosX==0    ) alignY_=1;
    if( iPosX/10==1 ) alignX_=1;
    if( iPosX/10==2 ) alignX_=2;
    if( iPosX/10==3 ) alignX_=3;
    int align_ = 10*alignX_ + alignY_;
    
    float H = pad->GetWh();
    float W = pad->GetWw();
    float l = pad->GetLeftMargin();
    float t = pad->GetTopMargin();
    float r = pad->GetRightMargin();
    float b = pad->GetBottomMargin();
    float e = 0.025;
    
    pad->cd();
    
    TString lumiText;
    if( iPeriod==1 )
    {
        lumiText += lumi_7TeV;
        lumiText += " (7 TeV)";
    }
    else if ( iPeriod==2 )
    {
        lumiText += lumi_8TeV;
        lumiText += " (8 TeV)";
    }
    else if( iPeriod==3 )
    {
        lumiText = lumi_8TeV;
        lumiText += " (8 TeV)";
        lumiText += " + ";
        lumiText += lumi_7TeV;
        lumiText += " (7 TeV)";
    }
    else if ( iPeriod==4 )
    {
        lumiText += lumi_13TeV;
        lumiText += " (13 TeV)";
    }
    else if ( iPeriod==7 )
    {
        if( outOfFrame ) lumiText += "#scale[0.85]{";
        lumiText += lumi_13TeV;
        lumiText += " (13 TeV)";
        lumiText += " + ";
        lumiText += lumi_8TeV;
        lumiText += " (8 TeV)";
        lumiText += " + ";
        lumiText += lumi_7TeV;
        lumiText += " (7 TeV)";
        if( outOfFrame) lumiText += "}";
    }
    else if ( iPeriod==12 )
    {
        lumiText += "8 TeV";
    }
    
    //cout << lumiText << endl;
    
    TLatex latex;
    latex.SetNDC();
    latex.SetTextAngle(0);
    latex.SetTextColor(kBlack);
    
    float extraTextSize = extraOverCmsTextSize*cmsTextSize;
    
    latex.SetTextFont(42);
    latex.SetTextAlign(31);
    latex.SetTextSize(lumiTextSize*t);
    latex.DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText);
    
    if( outOfFrame )
    {
        latex.SetTextFont(cmsTextFont);
        latex.SetTextAlign(11);
        latex.SetTextSize(cmsTextSize*t);
        latex.DrawLatex(l,1-t+lumiTextOffset*t,cmsText);
    }
    
    pad->cd();
    
    float posX_;
    if( iPosX%10<=1 )
    {
        posX_ =   l + relPosX*(1-l-r);
    }
    else if( iPosX%10==2 )
    {
        posX_ =  l + 0.5*(1-l-r);
    }
    else if( iPosX%10==3 )
    {
        posX_ =  1-r - relPosX*(1-l-r);
    }
    float posY_ = 1-t - relPosY*(1-t-b);
    if( !outOfFrame )
    {
        if( drawLogo )
        {
            posX_ =   l + 0.045*(1-l-r)*W/H;
            posY_ = 1-t - 0.045*(1-t-b);
            float xl_0 = posX_;
            float yl_0 = posY_ - 0.15;
            float xl_1 = posX_ + 0.15*H/W;
            float yl_1 = posY_;
            TASImage* CMS_logo = new TASImage("CMS-BW-label.png");
            TPad* pad_logo = new TPad("logo","logo", xl_0, yl_0, xl_1, yl_1 );
            pad_logo->Draw();
            pad_logo->cd();
            CMS_logo->Draw("X");
            pad_logo->Modified();
            pad->cd();
        }
        else
        {
            latex.SetTextFont(cmsTextFont);
            latex.SetTextSize(cmsTextSize*t);
            latex.SetTextAlign(align_);
            latex.DrawLatex(posX_+0.05, posY_, cmsText);
            if( writeExtraText )
            {
                latex.SetTextFont(extraTextFont);
                latex.SetTextAlign(align_);
                latex.SetTextSize(extraTextSize*t);
                latex.DrawLatex(posX_, posY_- relExtraDY*cmsTextSize*t, extraText);
            }
        }
    }
    else if( writeExtraText )
    {
        if( iPosX==0)
        {
            posX_ =   l +  relPosX*(1-l-r);
            posY_ =   1-t+lumiTextOffset*t;
        }
        latex.SetTextFont(extraTextFont);
        latex.SetTextSize(extraTextSize*t);
        latex.SetTextAlign(align_);
        latex.DrawLatex(posX_, posY_, extraText);
    }
    return;
}


void draw_header(){
	
	TString cms_pre = "CMS Preliminary";
	TString cms_sim = "CMS Simulation";
	TString lumi = "#sqrt{s} = 8 TeV, L = 19.7 fb^{-1}";
	TLatex L1;
	TLatex L2;
	
	L1.SetNDC(kTRUE);
	L2.SetNDC(kTRUE);
	L1.SetTextSize(0.03);
	L2.SetTextSize(0.03);
	L1.DrawLatex(0.15,0.92, cms_pre);
	L2.DrawLatex(0.43,0.92, lumi);
	
}

void create_theory(){
	double mass[]={300,400,500,600,700,800,900,1000};
	int N=sizeof(mass)/sizeof(double);
	double stopxs[]={1.99608,0.35683,0.0855847,0.0248009,0.0081141,0.00289588,0.00109501,0.000435488};
	double squarkxs[]={19.8283,3.54338,0.847051,0.244862,0.0799667,0.0284146, 0.0106744, 0.00424173,0.0106744};
    double squarkxsE[]={15.8736,15.6882,16.3916,17.9829,19.7578,22.2285,25.4549};
    
    double scale=0.4;
    
	CreateGraph("stop_xs_st_Acc0p5");
	
	for(int iP=0; iP<N; iP++){
		grName["stop_xs"]->SetPoint(iP,mass[iP],stopxs[iP]*1000);
		grName["stop_xs_st_Acc0p5"]->SetPoint(iP,mass[iP]*2,stopxs[iP]*1000*0.005);
	}
	for(int iP=0; iP<N; iP++){
		grName["squark_xs"]->SetPoint(iP,mass[iP],scale*squarkxs[iP]*1000);
        grName["squark_xs"]->SetPointError(iP,0,0,scale*squarkxs[iP]*(squarkxsE[iP]/100)*1000,squarkxs[iP]*(squarkxsE[iP]/100)*1000);

	}
	
}

void CreateGraph(TString name){
	TGraphAsymmErrors *gr = new TGraphAsymmErrors();
	gr->SetName(name);
	grName[name]=gr;
}

void CreateHistogram(const char* name,   const char* title,
					 const char* xTitle, const char* yTitle,
					 Int_t       nBinsX, Double_t    xLow, Double_t xUp)
{
	TH1F* h = new TH1F(name, title, nBinsX, xLow, xUp);
	
	h->GetXaxis()->SetTitle(xTitle);
	h->GetYaxis()->SetTitle(yTitle);
	
	h->Sumw2();
	h->SetStats(kFALSE);
	hName[name] = h;
}

void bookGraphs(){
	CreateHistogram("Mass_Bins","","","",8,250,1050);
	
	for(int mass=300; mass<=600; mass+=100){
		CreateHistogram(Form("h_sensitivity_M%d_exp",mass),"","n-jets","r",15,0.5,15.5 );
		CreateHistogram(Form("h_sensitivity_M%d_obs",mass),"","n-jets","r",15,0.5,15.5 );
        
		CreateHistogram(Form("h_sensitivity_M%d_exp_1sigmaP",mass),"","n-jets","r",15,0.5,15.5 );
		CreateHistogram(Form("h_sensitivity_M%d_exp_1sigmaM",mass),"","n-jets","r",15,0.5,15.5 );
		
		int nJ=4;
		
		for(int i=0; i<15; i++){
			TString label=Form("%d",nJ);
			if(nJ==8)label=Form("#geq %d",nJ);
			hName[Form("h_sensitivity_M%d_exp",mass)]->GetXaxis()->SetBinLabel(i+1,label);
			
			nJ++;
			if(nJ>8)nJ=4;
		}
		
		CreateGraph(Form("sensitivity_M%d_obs",mass) );
		CreateGraph(Form("sensitivity_M%d_exp",mass) );
		CreateGraph(Form("sensitivity_M%d_exp_1sigma",mass) );
		CreateGraph(Form("sensitivity_M%d_exp_2sigma",mass) );
        
		for(int i=0;  i<15; i++){
			grName[Form("sensitivity_M%d_obs",mass)]->SetPoint(i,i,0);
			grName[Form("sensitivity_M%d_exp",mass)]->SetPoint(i,i,0);
			grName[Form("sensitivity_M%d_exp_1sigma",mass)]->SetPoint(i,i,0);
			grName[Form("sensitivity_M%d_exp_2sigma",mass)]->SetPoint(i,i,0);
		}
		
	}
	
	CreateGraph("xs_exp_combined_1sigma");
	CreateGraph("xs_exp_combined_2sigma");
    
	CreateGraph("xs_exp_combined");
	CreateGraph("xs_obs_combined");
    
	CreateGraph("stop_xs");
	CreateGraph("squark_xs");
    
	
}

void read_data(){
	ifstream input;
	input.open("datacards/limit_summary.txt");
	int mass=300;
	int nJ;
	int st;
	
	double obs;
	double exp;
	double expp1;
	double expp2;
	double expm1;
	double expm2;
	
	
	bool print=true;
	int iP=0;
	/*
	cout << "first get binned limits:  " << endl;
	while(!input.eof()){
		input >> obs >> exp >> expm1 >> expm2 >> expp1 >> expp2;
		string line;
		getline(input,line);
		
		if(mass<300 || mass > 800)continue;
		if(nJ>8 || nJ<4) continue;
		if(st > 2 || st<0) continue;
		
		
		if(print){
			//cout << "mass: " << mass << " nJ : " << nJ << " st: " << st;
			//cout << " obs: " << obs << " exp: " << exp << " expm1: " << expm1 << " expm2: " << expm2;
			//cout << " expp1: " << expp1 << " expp2: " << expp2 << endl;
		}
		
		if(exp>0){
			iP=nJ-3;
			if(st==1)iP=nJ+5-3;
			if(st==2)iP=nJ+2*5-3;
			double ex=0.5;
			grName[Form("sensitivity_M%d_obs",mass)]->SetPoint(iP,iP,obs);
            
			grName[Form("sensitivity_M%d_exp",mass)]->SetPoint(iP,iP,exp);
			grName[Form("sensitivity_M%d_exp_1sigma",mass)]->SetPoint(iP,iP,exp);
			grName[Form("sensitivity_M%d_exp_2sigma",mass)]->SetPoint(iP,iP,exp);
            
			double eyl_1S=(exp-expm1);
			double eyh_1S=(expp1-exp);
			
			double eyl_2S=(exp-expm2);
			double eyh_2S=(expp2-exp);
			grName[Form("sensitivity_M%d_exp",mass)]->SetPointError(iP,ex,ex,0,0);
			grName[Form("sensitivity_M%d_exp_1sigma",mass)]->SetPointError(iP,ex,ex,eyl_1S,eyh_1S);
			grName[Form("sensitivity_M%d_exp_2sigma",mass)]->SetPointError(iP,ex,ex,eyl_2S,eyh_2S);
			
			hName[Form("h_sensitivity_M%d_obs",mass)]->SetBinContent(iP,obs);
			hName[Form("h_sensitivity_M%d_exp",mass)]->SetBinContent(iP,exp);
			hName[Form("h_sensitivity_M%d_exp_1sigmaP",mass)]->SetBinContent(iP,expp1);
			hName[Form("h_sensitivity_M%d_exp_1sigmaM",mass)]->SetBinContent(iP,expm1);
            
			//cout << "mass: " << mass << " iP: "<< obs << endl;
            
		}
		
	}
	input.close();
         */
    input.close();
	cout << "combined limit summary: " << endl;
	//input.open("datacards/limit_summary_combined.txt");
    input.open("datacards/limit_summary_toys.txt");
	iP=0;
    mass=300;
	
	while(!input.eof()){

		double MassE=50;
		input >> obs >> exp >> expm1 >> expm2 >> expp1 >> expp2;
		string line;
		getline(input,line);
        if(mass<300 || mass > 900)continue;
        //cout << "mass: " << mass << " obs: " << obs << " exp: " << exp << endl;
        //cout << "mass bin center: " << hName["Mass_Bins"]->GetBinCenter(hName["Mass_Bins"]->FindBin(mass)) << " " <<hName["Mass_Bins"]->FindBin(mass)-1 << endl;
        double squarkXS=grName["squark_xs"]->GetY()[hName["Mass_Bins"]->FindBin(mass)-1];
        //cout << "squark xs: " << squarkXS << endl;
        double expXS=squarkXS*exp;
        double obsXS=squarkXS*obs;
        
        //cout << "squarkXS: " << squarkXS << " " <<obsXS << endl;
        
        double eyl_1S=squarkXS*(exp-expm1);
        double eyh_1S=squarkXS*(expp1-exp);
        
        //cout << "exp +1sigma: " << EXPP1[min] << " - " << EXPM1[min] << endl;
        //cout << "expxs: " << expXS << " eyl:" << eyl_1S << " ey_h: "<< eyh_1S << endl;
        
        double eyl_2S=squarkXS*(exp-expm2);
        double eyh_2S=squarkXS*(expp2-exp);
        //cout << "fill graphs exp: " << iP << endl;
        grName["xs_exp_combined"]->SetPoint(iP,mass,expXS);
        grName["xs_exp_combined_1sigma"]->SetPoint(iP,mass,expXS);
        grName["xs_exp_combined_2sigma"]->SetPoint(iP,mass,expXS);
        
        grName["xs_exp_combined_1sigma"]->SetPointError(iP,MassE,MassE,eyl_1S,eyh_1S);
        grName["xs_exp_combined_2sigma"]->SetPointError(iP,MassE,MassE,eyl_2S,eyh_2S);
        
        grName["xs_obs_combined"]->SetPoint(iP,mass,obsXS);
        
        cout << "mass: " << mass << " obs " << obs << " exp " << exp << " + " << expp1 << " - " << expm1 << endl;
        cout << "mass: " << mass<< " theory: " << squarkXS <<" obs " << obsXS << " exp: " << expXS << " + " <<expXS+eyh_1S << " - " << expXS-eyl_1S << endl;
        
        iP++;
        mass+=100;
 
	}
    
	
	input.close();
    /*
    for(int i=0; i<7; i++){
        cout << "mass: " << grName["xs_obs_combined"]->GetX()[i] << " " << grName["xs_obs_combined"]->GetY()[i] << endl;
    }
    */
    
    for(double mass=500; mass<=800; mass+=5.0){
        double deltay=grName["xs_exp_combined"]->Eval(mass)-grName["squark_xs"]->Eval(mass);
        cout << "mass: " << mass << " " <<deltay << endl;
    }
     
    
    
	
}

void draw_model_independent(){
    
	
	
	
	for(int nJets=4; nJets<=8; nJets++){
		for (int nb=0; nb<=2; nb++) {
			TLegend Leg(0.65,0.65,0.85,0.8);
			Leg.SetFillColor(10);
			Leg.SetBorderSize(0);
			Leg.SetLineColor(10);
			
			TString nameObs=Form("ModelIndependent_obs_nJets%d_nb%d",nJets,nb);
			TString nameExp=Form("ModelIndependent_exp_nJets%d_nb%d",nJets,nb);
			TString nameExp_1S=Form("ModelIndependent_exp_1sigma_nJets%d_nb%d",nJets,nb);
			TString nameExp_2S=Form("ModelIndependent_exp_2sigma_nJets%d_nb%d",nJets,nb);
			
			CreateCanvas(Form("ModelIndependent_Plots_nJets%d_%dbtags",nJets,nb), "", 600,600);
			CName[Form("ModelIndependent_Plots_nJets%d_%dbtags",nJets,nb)]->cd();
			gPad->SetLogy();
			grName[nameObs]->SetLineColor(kBlack);
			grName[nameExp]->SetLineColor(kBlack);
			grName["stop_xs_st_Acc0p5"]->SetLineColor(kRed);
			grName["stop_xs_st_Acc0p5"]->SetLineWidth(2);
			
			grName[nameExp]->SetLineStyle(kDashed);
			
			grName[nameExp_1S]->SetFillStyle(1001);
			grName[nameExp_2S]->SetFillStyle(1001);
			
			grName[nameExp_1S]->SetFillColor(3);
			grName[nameExp_2S]->SetFillColor(5);
			
			grName[nameExp_1S]->Draw("a3");
			grName[nameExp_1S]->GetXaxis()->SetRangeUser(300,2000);
            grName[nameExp_1S]->GetYaxis()->SetRangeUser(0.05,20);

			grName[nameExp_1S]->GetXaxis()->SetNdivisions(6);
			//grName[nameExp_1S]->GetYaxis()->SetRangeUser(5,1.5*grName["squark_xs"]->Eval(300));
			grName[nameExp_1S]->GetXaxis()->SetTitle("S_{T}^{min} (GeV)");
			grName[nameExp_1S]->GetYaxis()->SetTitle("#sigma*A (fb)");
			
			Leg.AddEntry(grName[nameObs],"observed","L");
			Leg.AddEntry(grName[nameExp],"expected","L");
			Leg.AddEntry(grName[nameExp_1S],"#pm 1SD","F");
            
			TLatex txt;
			txt.SetNDC(kTRUE);
			
			grName[nameExp_1S]->Draw("a3");
			grName[nameExp]->Draw("l same");
			grName[nameObs]->Draw("l same");
			//grName["stop_xs_st_Acc0p5"]->Draw("l same");
			Leg.DrawClone();
			txt.DrawLatex(0.5,0.8,Form("%d-jets, %d b-tags",nJets,nb));
			draw_header();
			gPad->RedrawAxis();
			
		}
	}
}

void model_Ind_Limit(){
    
    ifstream input;
    input.open("limit_summary_modelInd.txt");
    int nb;
    int nJ;
    int st;
    
    double obs;
    double exp;
    double expp1;
    double expp2;
    double expm1;
    double expm2;
	int NJ=0;
	int NB=0;
	int iP=0;
	
	for(int nJets=4; nJets<=8; nJets++){
		for (int inb=0; inb<=2; inb++) {
			TString name=Form("ModelIndependent_obs_nJets%d_nb%d",nJets,inb);
			CreateGraph(name);
            
			name=Form("ModelIndependent_exp_nJets%d_nb%d",nJets,inb);
			CreateGraph(name);
            
			name=Form("ModelIndependent_exp_1sigma_nJets%d_nb%d",nJets,inb);
			CreateGraph(name);
            
			name=Form("ModelIndependent_exp_2sigma_nJets%d_nb%d",nJets,inb);
			CreateGraph(name);
            
		}
	}
	
	
	cout << "Get model ind. limits for each point. " << endl;
    while(!input.eof()){
        input >> nJ >> nb >> st >> obs >> exp >> expm1 >> expm2 >> expp1 >> expp2;
        if(NJ!=nJ || NB!=nb){
            NJ=nJ;
            NB=nb;
            iP=0;
        }
        
        TString nameObs=Form("ModelIndependent_obs_nJets%d_nb%d",nJ,nb);
        TString nameExp=Form("ModelIndependent_exp_nJets%d_nb%d",nJ,nb);
        TString nameExp_1S=Form("ModelIndependent_exp_1sigma_nJets%d_nb%d",nJ,nb);
        TString nameExp_2S=Form("ModelIndependent_exp_2sigma_nJets%d_nb%d",nJ,nb);
        
        string line;
        getline(input,line);
		
        if(exp>0){
            
            double lumi=19.6;
            obs=obs/lumi;
            exp=exp/lumi;
            expm1=expm1/lumi;
            expm2=expm2/lumi;
            expp1=expp1/lumi;
            expp2=expp2/lumi;
            
            cout << "nJ: " << nJ << " nb: " << nb << " st: " << st;
            cout << " iP: " << iP << endl;
            cout << "obs: "<< obs << "exp: " << exp	<< endl;
            
            double ex=50; // 50 GeV bins
            grName[nameObs]->SetPoint(iP,st,obs);
            
            grName[nameExp]->SetPoint(iP,st,exp);
            grName[nameExp_1S]->SetPoint(iP,st,exp);
            grName[nameExp_2S]->SetPoint(iP,st,exp);
            
            double eyl_1S=(exp-expm1);
            double eyh_1S=(expp1-exp);
            
            double eyl_2S=(exp-expm2);
            double eyh_2S=(expp2-exp);
            grName[nameObs]->SetPointError(iP,0,0,0,0);
            grName[nameExp_1S]->SetPointError(iP,ex,ex,eyl_1S,eyh_1S);
            grName[nameExp_2S]->SetPointError(iP,ex,ex,eyl_2S,eyh_2S);
            
            //cout << "mass: " << mass << " iP: "<< obs << endl;
            iP++; 
        }
        
    }
    input.close();
}

void GetHistograms(){
	TFile *input_file = new TFile("/eos/uscms/store/user/btcarlso/histogram_files/susy_histograms.root","READ"); 
	TH1F *TH1F_names=(TH1F*)input_file->FindObjectAny("TH1F_names"); 
	
	
	for (int i=1; i<=TH1F_names->GetNbinsX(); i++) {
		
		TString name=TH1F_names->GetXaxis()->GetBinLabel(i);
		if(name=="")continue; 
		if(name=="TH1F_names" || name=="TH2F_names" || name=="TProfile_names") continue; 
		//cout <<" filename: "<< file <<  " hist_name: "<< name << endl; 
		
		TH1F *h=(TH1F*)input_file->FindObjectAny(name); 
		h->SetStats(kFALSE); 
		h->SetName(h->GetName()); 
		hName[h->GetName()]=h; 
	}
	delete TH1F_names; 
	//cout << "end of GetHistograms() " << endl; 
}
