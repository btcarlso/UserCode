/*
 *  create_plots.c
 *  
 *
 *  Created by Benjamin Carlson on 8/9/13.
 *  Copyright 2013 Carnegie Mellon University. All rights reserved.
 *
 */

#include "compute_xs.h"

void compute_xs(){
	gROOT->SetBatch();
	create_plots();

}

void create_plots(){
	hist_file = new TFile("hist_file.root","READ"); 
	yields_file=new TFile("yields_file_test.root","READ");
	canvas_file=new TFile("canvas_plots.root","RECREATE"); 
	
	eff1S = new TFile("../plot_efficiencies/output_file_efficiencies1S.root","READ"); 
	eff2S = new TFile("../plot_efficiencies/output_file_efficiencies2S.root","READ"); 
	eff3S = new TFile("../plot_efficiencies/output_file_efficiencies3S.root","READ");

	get_atlas_ratios();
	readHisto();
	draw_muon_pt();
	draw_zeta();
	draw_zeta_m();
	draw_m();
	draw_pt();
	draw_dR();
	draw_zeta_eta();
	draw_zeta_profile();
	fit_nuisance(0);
	fit_nuisance(1); 
	plot_orderDist();
	
	rapidity_dep_atlas();
	rap_dep2010();
	for(int ups=1; ups<=3; ups++){
		plot_cms2010(ups);
		plot_atlas_pt(ups);
		for(int iy=0; iy<=2; iy++)read_NLO_xs(ups,iy);
	}
	
	for(int iy=0; iy<3; iy++)load_sg(iy);
	load_acceptance();
	rebin_acceptance(); //rebin acceptance to analysis bins
	
	draw_efficiency(); // makes some efficiency plots, and computes some ratios of efficiency 
    for(int iy=0; iy<3; iy++)efficiency_table(iy);
	for(int ups=1; ups<=3; ups++)combine_rapidity_bins(ups); // add yields in 2 bins
	
	
	for(int iy=0; iy<fNy+1; iy++){
		acceptance_ratio(iy);//compute acceptance ratio 
		
		for(int ups=1; ups<=3; ups++){
			uncertainty(iy, ups);
			dsigma_pt(iy,ups); 
			
			xs_fit(iy, ups); 
			xs_fit_po_study(iy,ups);
			bin_center(iy,ups);
			
		}
		for(int num=2; num<=3; num++){
			uncertaintyR(iy,num); 
			RN1(iy,num); //Ratio nS-1
		}
	}

	for(int iy=0; iy<=2; iy++){
		for(int ups=1; ups<=3; ups++)weighted_acceptance(iy, ups);
	}
	
	for(int iy=0; iy<=2; iy++){
		fit_table(iy);
		xs_table(iy); 
		xs_tablehepdata(iy);
		ratio_table( iy); 
		ratio_tablehepdata(iy); 
		for(int ups=1; ups<=3; ups++){
            acceptance_ratio( iy, ups, "");
            if(iy<2){
                acceptance_ratio( iy, ups, "Ep");
                acceptance_ratio( iy, ups, "Em");
                acceptance_ratio( iy, ups, "_longitudinal");
                acceptance_ratio( iy, ups, "_transverse");
            }

            acceptance_table( iy, ups);
            acceptance_table_summary(iy,ups);
			acceptance_table_summaryhepdata(iy,ups); 
        }
	}
    acceptance_ratio_plot();
	acceptance_table_smearing(0,1);
	
	draw_rapidity_ratio();
    table_rapidity_ratio();
	for(int iy=0; iy<=2; iy++){
		for(int ups=1; ups<=3; ups++){
			draw_xs_fits(iy,ups);
			draw_xs_fits_wNLO(iy,ups); 
		}
	}
	draw_ratio(3);
	draw_ratio(2);
    
	for(int iy=0; iy<fNy; iy++){
		draw_unc(iy); 
		draw_unc_ratio(iy); 
	}
	
    atlas_cmsratio();
    
	draw_xs(); 
	draw_xs_atlas();
	
    compare_ratios(3);
    compare_ratios(2);
    
	efficiency_closure();
	
	canvas_file->cd();
	yield_table();
	writeHisto(); 
	
	
}

void draw_headersim(){
	
	L1sim.SetNDC(kTRUE); 
	
	L1sim.DrawLatex(0.15,0.92, cms_sim); 
	
}

void CMS_lumi( TPad* pad, int iPeriod, int iPosX )
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
            latex.DrawLatex(posX_, posY_, cmsText);
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

	L1_header.SetNDC(kTRUE); 
	L2_header.SetNDC(kTRUE); 
	
	L1_header.DrawLatex(0.17,0.92, cmsText); 
	L2_header.DrawLatex(0.26,0.92, lumi_header); 
	
}

void read_NLO_xs(int ups, int iy){
	TGraphAsymmErrors *gr = new TGraphAsymmErrors();
	gr->SetName(Form("gr_NLOtheory_%dS_y%d",ups,iy));
	
	string name=Form("/eos/uscms/store/user/btcarlso/NLO_xs/upsilon%ds/",ups); 
	string rap_name;
	if(iy==0) rap_name=Form("cms-0.6-pt.dat"); 
	if(iy==1) rap_name=Form("cms-1.2-pt.dat");
	if(iy==2) rap_name=Form("cms-all-pt.dat");
	name=name+rap_name; 
	ifstream in; 
	in.open(name.c_str()); 
	double pt; double sigma; double sigmaU;
	int ip=0; 
	while(!in.eof()){
		in >> pt >> sigma >> sigmaU; 
		string tmp; 
		getline(in,tmp); 
		//    cout << pt << " " << sigma << " " << sigmaU << endl; 
		gr->SetPoint(ip,pt,sigma*1000);
		gr->SetPointError(ip,0,0,sigmaU*1000,sigmaU*1000); 
		ip++;
	}
	gr->GetYaxis()->SetTitle("d#sigma/dp_{T} #times Br(#mu#mu) [pb / GeV]"); 
	gr->GetXaxis()->SetTitle("p_{T} [GeV]"); 
	
	in.close();
	
	grName[gr->GetName()]=gr; 
	
	
}

void readHisto(){
	cout << "read Histo: " << endl; 
	TH1F *TH1F_names = (TH1F*)hist_file->FindObjectAny("TH1F_names"); 
	TH1F *TH2F_names = (TH1F*)hist_file->FindObjectAny("TH2F_names"); 
	TH1F *TProfile_names = (TH1F*)hist_file->FindObjectAny("TProfile_names"); 


	cout << "read Histogram: " << endl; 
	
	for (int i=1; i<=TH1F_names->GetNbinsX(); i++) {
		
		TString name=TH1F_names->GetXaxis()->GetBinLabel(i);
		if(name=="")continue; 
		if(name=="TH1F_names" || name=="TH2F_names" || name=="TProfile_names" || name.Contains("rho")) continue; 
		

		TH1F *h=(TH1F*)hist_file->FindObjectAny(name); 
		hName[h->GetName()]=h; 
	}
	
	for(int iy=0; iy<fNy; iy++){
		
		
		for(int ups=1; ups<=3; ups++){
			TString name=yield_histogram(iy,ups,"best","");
			TH1F *h=(TH1F*)yields_file->FindObjectAny(name); 
			hName[h->GetName()]=h; 
		}
		for(int ups=1; ups<=3; ups++){
			TString name=background_order_systematic(iy,Form("%dS",ups)); 
			
			TH1F *h=(TH1F*)yields_file->FindObjectAny(name); 
			hName[h->GetName()]=h; 
		}
		
		for(int num=2; num<=3; num++){
			TString name=Form("ratio_%d1_y%d_best",num,iy);
			TH1F *h=(TH1F*)yields_file->FindObjectAny(name); 
			hName[h->GetName()]=h; 
		}
		
		for(int num=2; num<=3; num++){
			TString name=background_order_systematic(iy,Form("R%d1",num));
		
			TH1F *h=(TH1F*)yields_file->FindObjectAny(name); 
			hName[h->GetName()]=h; 
		}
		TString name=Form("chi2_y%d",iy); 
		TH1F *h=(TH1F*)yields_file->FindObjectAny(name); 
		hName[h->GetName()]=h;
		
		name=Form("NDOF_y%d",iy); 
		TH1F *h2=(TH1F*)yields_file->FindObjectAny(name); 
		hName[h2->GetName()]=h2;
		
		TH1F *cw=(TH1F*)yields_file->FindObjectAny(Form("cw_y%d_best",iy)); 
		hName[cw->GetName()]=cw; 										   
		TH1F *deltaM=(TH1F*)yields_file->FindObjectAny(Form("deltaM_y%d_best",iy)); 
		hName[deltaM->GetName()]=deltaM; 
		
		TH1F *order=(TH1F*)yields_file->FindObjectAny(Form("Background_order_y%d_",iy)); 
		hName[order->GetName()]=order; 

	}
	
	cout << "read profiles: " << endl; 
	for (int i=1; i<=TProfile_names->GetNbinsX(); i++) {
		TString name=TProfile_names->GetXaxis()->GetBinLabel(i);
	//	cout << "getting name: " << name << endl; 
		if(name=="")continue; 
		if(name=="TH1F_names" || name=="TH2F_names" || name=="TProfile_names") continue; 
		TProfile *pr=(TProfile*)hist_file->FindObjectAny(name); 
		prName[pr->GetName()]=pr;
	}
	cout << "read TH2: "<< endl; 
	for (int i=1; i<=TH2F_names->GetNbinsX(); i++) {
		TString name=TH2F_names->GetXaxis()->GetBinLabel(i);
		if(name=="")continue; 
		if(name=="TH1F_names" || name=="TH2F_names"|| name=="TProfile_names") continue; 
		TH2F *h=(TH2F*)hist_file->FindObjectAny(name); 
		hName2D[h->GetName()]=h; 
	}
	cout << "end: " << endl; 
}

void writeHisto(){

	canvas_file->mkdir("eff"); 
	canvas_file->mkdir("mass"); 
	canvas_file->mkdir("zeta"); 
	
	for (std::map<TString,TCanvas*>::iterator it=CName.begin(); it!=CName.end(); it++) {
		canvas_file->cd();
		if(it->first.Contains("m_y")) canvas_file->cd("mass"); 
		if(it->first.Contains("h_eff")) canvas_file->cd("eff"); 
		if(it->first.Contains("zeta")) canvas_file->cd("zeta"); 

		CName[it->first]->Write();
		bool print = 1;
		if(print){
			TString plot_name=it->first+".pdf"; 
			if(it->first.Contains("zeta")) plot_name="zeta/"+it->first+".pdf";
			if(it->first.Contains("eff")){
				if(it->first.Contains("h_eff") && it->first.Contains("pt0")==0) continue; 
				plot_name="efficiency/"+it->first+".pdf";
			}
			//CName[it->first]->Print("/uscms/home/btcarlso/AN-UPS/tdr2/notes/AN-13-088/trunk/Fig/"+plot_name); 
			//CName[it->first]->Print("/uscms_data/d3/btcarlso/UpsAN/FigJPG/"+it->first+".jpg"); 
		}
		
	}
}


void fit_table(int iy){
	//F->GetChisquare()
	//F->GetNDF();
	ofstream output; 
	output.open(Form("/uscms/home/btcarlso/AN-UPS/tdr2/notes/AN-13-088/trunk/tables/fit_table_y%d.tex",iy)); 
	int N=4; 
	output <<"\\begin{tabular}{"; 
	for(int i=0; i<=N; i++) output << "c"; 
	output << "}\\\\ " << endl << " \\hline" << endl; 
	
	output << " & $\\Upsilon(1S)$ & $\\Upsilon(2S)$ & $\\Upsilon(3S)$ \\\\ \\hline " << endl; 
	
	for(int ipar=0; ipar<3; ipar++){
        TString name=f1Name[Form("xs_fit_power_y%d_%dS",iy,1)]->GetParName(ipar);
        if(name=="alpha") name="$\\alpha$";
		output << name << " & ";
		for(int ups=1; ups<=3; ups++){
			output << Form("%.2f $\\pm$ %.2f",f1Name[Form("xs_fit_power_y%d_%dS",iy,ups)]->GetParameter(ipar),f1Name[Form("xs_fit_power_y%d_%dS",iy,ups)]->GetParError(ipar)) << " & ";
			
		}
		output << "\\\\ " << endl; 

	}
	output << "$\\chi^{2}$ & "; 
	for(int ups=1; ups<=3; ups++){
		output <<setprecision(2)  << f1Name[Form("xs_fit_power_y%d_%dS",iy,ups)]->GetChisquare() 	<< " & ";
	}
	output << "\\\\" << endl; 
	output << " $n_{d}$ & "; 
	for(int ups=1; ups<=3; ups++){
		output << setprecision(2)  << f1Name[Form("xs_fit_power_y%d_%dS",iy,ups)]->GetNDF()	<< " & ";
	}
	output << "\\\\" << endl; 

	////// Exponential fit
	
	output << "\\hline \\hline " << endl; 
	for(int ipar=0; ipar<2; ipar++){
		output << f1Name[Form("xs_fit_expo_y%d_%dS",iy,1)]->GetParName(ipar) << " & ";
		for(int ups=1; ups<=3; ups++){
			output << Form("%.2f $\\pm$ %.2f",f1Name[Form("xs_fit_expo_y%d_%dS",iy,ups)]->GetParameter(ipar),f1Name[Form("xs_fit_expo_y%d_%dS",iy,ups)]->GetParError(ipar)) << " & ";
			
		}
		output << "\\\\ " << endl; 
	}
	output << "$\\chi^{2}$ & "; 
	for(int ups=1; ups<=3; ups++){
		output <<setprecision(2) << f1Name[Form("xs_fit_expo_y%d_%dS",iy,ups)]->GetChisquare() 	<< " & ";
	}
	output << "\\\\" << endl; 
	output << " $N_d$ & "; 
	for(int ups=1; ups<=3; ups++){
		output <<setprecision(2) << f1Name[Form("xs_fit_expo_y%d_%dS",iy,ups)]->GetNDF()	<< " & ";
	}
	
	output << "\\\\ \\hline" << endl;
	
	output << "\\end{tabular} " << endl; 
	output.close();
	//short form for paper. 
	
	output.open(Form("/uscms/home/btcarlso/AN-UPS/tdr2/notes/AN-13-088/trunk/tables/fit_table_y%d_paper.tex",iy)); 
	N=4; 
	output <<"\\begin{tabular}{"; 
	for(int i=0; i<=N; i++) output << "c"; 
	output << "}\\\\ " << endl << " \\hline" << endl; 
	
	output << " & $\\Upsilon(1S)$ & $\\Upsilon(2S)$ & $\\Upsilon(3S)$ \\\\ \\hline " << endl; 
	
	for(int ipar=0; ipar<3; ipar++){
        TString name=f1Name[Form("xs_fit_power_y%d_%dS",iy,1)]->GetParName(ipar);
        if(name=="alpha") name="$\\alpha$";
		output << name << " & ";
		for(int ups=1; ups<=3; ups++){
			output << Form("%.2f $\\pm$ %.2f",f1Name[Form("xs_fit_power_y%d_%dS",iy,ups)]->GetParameter(ipar),f1Name[Form("xs_fit_power_y%d_%dS",iy,ups)]->GetParError(ipar)) << " & ";
			
		}
		output << "\\\\ " << endl; 
		
	}
	output << "$\\chi^{2}$ & "; 
	for(int ups=1; ups<=3; ups++){
		output <<setprecision(2)  << f1Name[Form("xs_fit_power_y%d_%dS",iy,ups)]->GetChisquare() 	<< " & ";
	}
	output << "\\\\" << endl; 
	output << " $n_{d}$ & "; 
	for(int ups=1; ups<=3; ups++){
		output << setprecision(2)  << f1Name[Form("xs_fit_power_y%d_%dS",iy,ups)]->GetNDF()	<< " & ";
	}
	output << "\\\\ \\hline" << endl; 
	
	output << "\\end{tabular} " << endl; 
	output.close();
	
}

void xs_table(int iy){
	ofstream output; 
	output.open(Form("/uscms/home/btcarlso/AN-UPS/tdr2/notes/AN-13-088/trunk/tables/xs_table_y%d.tex",iy));
    output << "\\vspace*{.2cm}" << endl;
    output << "\\small" << endl;
    output << "\\setlength{\\extrarowheight}{5pt}" << endl;
    output <<"\\setlength{\\tabcolsep}{4pt}" << endl;
    output<<"\\begin{tabular}{cc|ccc|ccc|ccc}\\hline" << endl;
    
    output << "& &  \\multicolumn{3}{c}{$\\Upsilon(1S)$} & \\multicolumn{3}{c}{$\\Upsilon(2S)$} & \\multicolumn{3}{c}{$\\Upsilon(3S)$} \\\\ \\hline" << endl;
    output << "$p_{\\rm T}$ & $\\langle p_{\\rm T} \\rangle$ & $\\frac{\\text{d}\\sigma}{\\text{dp}_{\\rm T}} \\cdot B$  & $\\frac{\\sigma_{\\rm stat}}{d \\sigma/d p_{\\rm T}}$   & $\\frac{\\sigma_{\\rm syst}}{d\\sigma/d p_{\\rm T}}$ & $\\frac{\\text{d}\\sigma}{\\text{dp}_{\\rm T}} \\cdot B$  & $\\frac{\\sigma_{\\rm stat}}{d \\sigma/d p_{\\rm T}}$   & $\\frac{\\sigma_{\\rm syst}}{d\\sigma/d p_{\\rm T}}$ & $\\frac{\\text{d}\\sigma}{\\text{dp}_{\\rm T}} \\cdot B$  & $\\frac{\\sigma_{\\rm stat}}{d \\sigma/d p_{\\rm T}}$   & $\\frac{\\sigma_{\\rm syst}}{d\\sigma/d p_{\\rm T}}$    \\\\ \\hline " << endl;
   output <<" GeV & GeV & (fb/GeV) & (\\%) & (\\%) & (fb/GeV) & (\\%) & (\\%) & (fb/GeV) & (\\%) & (\\%) \\\\ \\hline" << endl;
    /*
	int N=5; 
	output <<"\\begin{tabular}{"; 
	for(int i=0; i<=N; i++) output << "c"; 
	output << "}\\hline" << endl; 
	output << "$p_{T}$ & $<p_{T}[1S]>$ & $d\\sigma[1S]/dp_{T}$ & $d\\sigma[2S]/dp_{T}$ &  $d\\sigma[3S]/dp_{T}$  \\\\ \\hline" << endl; 
	output << "GeV & GeV & fb/GeV  & fb/GeV  & fb/GeV  \\\\ \\hline" << endl; 
	output << "\\hline" << endl; 
	*/
    /*
	for(int ipt=0; ipt<fNpt; ipt++){
		output << setprecision(2) << fPTbin[ipt] << "--" << setprecision(3)<<  fPTbin[ipt+1] << " & ";
		output  << setprecision(3) << grName[xs_graph(iy,1)]->GetX()[ipt] << " & "; 
		for(int ups=1; ups<=3; ups++){
			int precision=set_precision(grName[xs_graph(iy,ups)]->GetY()[ipt]*1000);
			int precisionE=set_precisionE(grName[xs_graph(iy,ups)]->GetErrorYhigh(ipt)*1000 ); 
			
			output  << setprecision(precision) << grName[xs_graph(iy,ups)]->GetY()[ipt]*1000 << " & $\\pm$ ";
			output  << setprecision(2) << grName[xs_graph_stat(iy,ups)]->GetErrorYhigh(ipt)*100/grName[xs_graph(iy,ups)]->GetY()[ipt]; 
			
			
			double sysP=100*TMath::Sqrt(TMath::Power(grName[xs_graph(iy,ups)]->GetErrorYhigh(ipt),2)-TMath::Power(grName[xs_graph_stat(iy,ups)]->GetErrorYhigh(ipt),2));
			sysP=sysP/grName[xs_graph(iy,ups)]->GetY()[ipt];

			double sysM=100*TMath::Sqrt(TMath::Power(grName[xs_graph(iy,ups)]->GetErrorYlow(ipt),2)-TMath::Power(grName[xs_graph_stat(iy,ups)]->GetErrorYhigh(ipt),2)); 
			sysM=sysM/grName[xs_graph(iy,ups)]->GetY()[ipt];
			
			//cout << "xs: " << grName[xs_graph(iy,ups)]->GetY()[ipt]*1000  << endl; 
			
			//cout << "systematic: " << sysP << "  stat " << grName[xs_graph_stat(iy,ups)]->GetErrorYhigh(ipt)*1000 << " tot "  << 1000*grName[xs_graph(iy,ups)]->GetErrorYhigh(ipt) << endl; 
			
			//cout << "tot uncertainty: " << hName[tot_uncertainty(iy,ups,"Em")]->GetBinContent(ipt+1) << endl; 
			 
			output  << setprecision(2)<<" & $\\pm$ " << sysP << " (" << sysM <<  ") ";
            if(ups<3) output <<" & ";
		}
		output << "\\\\" << endl; 
	}*/
    for(int ipt=0; ipt<fNpt; ipt++){
        output << Form("%.0f--%.0f & %.1f & ",fPTbin[ipt],fPTbin[ipt+1],grName[xs_graph(iy,1)]->GetX()[ipt] );
		//output << setprecision(2) << fPTbin[ipt] << "--" << setprecision(3)<<  fPTbin[ipt+1] << " & ";
		//output  << setprecision(3) << grName[xs_graph(iy,1)]->GetX()[ipt] << " & ";
		for(int ups=1; ups<=3; ups++){
			int precision=set_precision(grName[xs_graph(iy,ups)]->GetY()[ipt]*1000);
			int precisionE=set_precisionE(grName[xs_graph(iy,ups)]->GetErrorYhigh(ipt)*1000 );
			
			//output  << setprecision(precision) << grName[xs_graph(iy,ups)]->GetY()[ipt]*1000 << " & $\\pm$ ";
			//output  << setprecision(2) << grName[xs_graph_stat(iy,ups)]->GetErrorYhigh(ipt)*100/grName[xs_graph(iy,ups)]->GetY()[ipt];
			
            if(grName[xs_graph(iy,ups)]->GetY()[ipt]*1000>=10) output << Form("%.0f & ",grName[xs_graph(iy,ups)]->GetY()[ipt]*1000);
            else output << Form("%.1f & ",grName[xs_graph(iy,ups)]->GetY()[ipt]*1000);

            double stat=grName[xs_graph_stat(iy,ups)]->GetErrorYhigh(ipt)*100/grName[xs_graph(iy,ups)]->GetY()[ipt];
            
            output << Form("%.1f & ",stat);
            //if(stat >=1 && stat<=10) output << Form("%.1f & $\\pm$ ",stat);
            //if(stat >=10) output << Form("%.0f & $\\pm$ ",stat);

			
			double sysP=100*TMath::Sqrt(TMath::Power(grName[xs_graph(iy,ups)]->GetErrorYhigh(ipt),2)-TMath::Power(grName[xs_graph_stat(iy,ups)]->GetErrorYhigh(ipt),2));
			sysP=sysP/grName[xs_graph(iy,ups)]->GetY()[ipt];
            
			double sysM=100*TMath::Sqrt(TMath::Power(grName[xs_graph(iy,ups)]->GetErrorYlow(ipt),2)-TMath::Power(grName[xs_graph_stat(iy,ups)]->GetErrorYhigh(ipt),2));
			sysM=sysM/grName[xs_graph(iy,ups)]->GetY()[ipt];
			
			//cout << "xs: " << grName[xs_graph(iy,ups)]->GetY()[ipt]*1000  << endl;
			
			//cout << "systematic: " << sysP << "  stat " << grName[xs_graph_stat(iy,ups)]->GetErrorYhigh(ipt)*1000 << " tot "  << 1000*grName[xs_graph(iy,ups)]->GetErrorYhigh(ipt) << endl;
			
			//cout << "tot uncertainty: " << hName[tot_uncertainty(iy,ups,"Em")]->GetBinContent(ipt+1) << endl;
            
			//output  << setprecision(2)<<" & $\\pm$ " << sysP << " (" << sysM <<  ") ";
            
            double sys=TMath::Max(sysP,sysM);
            
            output << Form("%.1f (%.1f)",sysP,sysM);
            //if(sys >=1 && sys<=10) output << Form("%.1f (%.1f)",sysP,sysM);
            //if(sys >=10) output << Form("%.0f (%.0f)",sysP,sysM);
            
            if(ups<3) output <<" & ";
		}
		output << "\\\\" << endl;
	}
    
	
	output << "\\hline" << endl; 
	output << "\\end{tabular}" << endl; 
	output.close();

}

void xs_tablehepdata(int iy){
	ofstream output; 
	output.open(Form("/uscms/home/btcarlso/hepdata/xs_table_y%d.tex",iy));
	output << "*reackey: P P --> UPSI(1S) X " << endl; 
	output << "*obskey: DSIGMA/DPT" << endl; 
	output << "*dserror: +2.2,-2.2 PCT : luminosity uncertainty" << endl; 
	if(iy==0) output << "*qual: ABS(YRAP) : < 0.6" << endl;   
	if(iy==1) output << "*qual: ABS(YRAP) : 0.6 TO 1.2" << endl;   
	if(iy==2) output << "*qual: ABS(YRAP) : < 1.2" << endl;   
	output << "*qual: RE : P P --> UPSI(nS)  < MU+ MU -> X" << endl; 
	output << "*qual: SQRT(S) IN GEV : 7000.0" << endl;
	output << "*qual: . : UPSI(1S) : UPSI(2S) : UPSI(3S)" << endl; 
	output << "*yheader: SIG IN FB/GEV" << endl; 
	output << "*xheader: PT IN GEV : $\\langle p_{T}\\rangle$ IN GeV" << endl;  
	
	
	//output << "& &  \\multicolumn{3}{c}{$\\Upsilon(1S)$} & \\multicolumn{3}{c}{$\\Upsilon(2S)$} & \\multicolumn{3}{c}{$\\Upsilon(3S)$} \\\\ \\hline" << endl;
    //output << "$p_{\\rm T}$ & $\\langle p_{\\rm T} \\rangle$ & $\\frac{\\text{d}\\sigma}{\\text{dp}_{\\rm T}} \\cdot B$  & $\\frac{\\sigma_{\\rm stat}}{d \\sigma/d p_{\\rm T}}$   & $\\frac{\\sigma_{\\rm syst}}{d\\sigma/d p_{\\rm T}}$ & $\\frac{\\text{d}\\sigma}{\\text{dp}_{\\rm T}} \\cdot B$  & $\\frac{\\sigma_{\\rm stat}}{d \\sigma/d p_{\\rm T}}$   & $\\frac{\\sigma_{\\rm syst}}{d\\sigma/d p_{\\rm T}}$ & $\\frac{\\text{d}\\sigma}{\\text{dp}_{\\rm T}} \\cdot B$  & $\\frac{\\sigma_{\\rm stat}}{d \\sigma/d p_{\\rm T}}$   & $\\frac{\\sigma_{\\rm syst}}{d\\sigma/d p_{\\rm T}}$    \\\\ \\hline " << endl;
	//output <<" GeV & GeV & (fb/GeV) & (\\%) & (\\%) & (fb/GeV) & (\\%) & (\\%) & (fb/GeV) & (\\%) & (\\%) \\\\ \\hline" << endl;
	output << "*data: x : x: y : y : y " << endl; 
       for(int ipt=0; ipt<fNpt; ipt++){
        output << Form("%.1f TO %.1f ; %.1f ; ",fPTbin[ipt],fPTbin[ipt+1],grName[xs_graph(iy,1)]->GetX()[ipt] );
		for(int ups=1; ups<=3; ups++){
			if(grName[xs_graph(iy,ups)]->GetY()[ipt]*1000>=10) output << Form("%.0f ",grName[xs_graph(iy,ups)]->GetY()[ipt]*1000);
            else output << Form("%.1f ",grName[xs_graph(iy,ups)]->GetY()[ipt]*1000);
			
            double stat=grName[xs_graph_stat(iy,ups)]->GetErrorYhigh(ipt)*100/grName[xs_graph(iy,ups)]->GetY()[ipt];
            
            output << Form("+- %.1f PCT",stat);
			double sysP=100*TMath::Sqrt(TMath::Power(grName[xs_graph(iy,ups)]->GetErrorYhigh(ipt),2)-TMath::Power(grName[xs_graph_stat(iy,ups)]->GetErrorYhigh(ipt),2));
			sysP=sysP/grName[xs_graph(iy,ups)]->GetY()[ipt];
            
			double sysM=100*TMath::Sqrt(TMath::Power(grName[xs_graph(iy,ups)]->GetErrorYlow(ipt),2)-TMath::Power(grName[xs_graph_stat(iy,ups)]->GetErrorYhigh(ipt),2));
			sysM=sysM/grName[xs_graph(iy,ups)]->GetY()[ipt];
            
            output << Form("(DSYS=+%.1f,-%.1f PCT) ; ",sysP,sysM);

		}
		output << endl;
	}
    
	output << "*dataend" << endl;
	output.close();
	
}


void ratio_table(int iy){
	ofstream output; 
	output.open(Form("/uscms/home/btcarlso/AN-UPS/tdr2/notes/AN-13-088/trunk/tables/ratio_table_y%d.tex",iy)); 
	/*
    int N=5;
	output <<"\\begin{tabular}{"; 
	for(int i=0; i<=N; i++) output << "c"; 
	output << "}\\hline" << endl; 
	output << "$p_{T}$ & $<p_{T}(2S)>$ & $R_{21}$ & $<p_{T}(3S)>$ & $R_{31}$  \\\\ \\hline" << endl; 
	output << "GeV & GeV & GeV & ratio & ratio  \\\\ \\hline" << endl; 
	output << "\\hline" << endl; 
	*/
    //output << "\\vspace*{.2cm}" << endl;
    //output << "\\setlength{\\extrarowheight}{5pt}" << endl;
    output << "\\begin{tabular}{c|ccc|ccc}\\hline" << endl;
    output << "$p_{\\rm T}$ & $R_{21}$ & $\\frac{\\sigma_{\\rm stat}}{R_{21}}$   & $\\frac{\\sigma_{\\rm syst}}{R_{21}}$ & $R_{31}$  & $\\frac{\\sigma_{\\rm stat}}{R_{31}}$   & $\\frac{\\sigma_{\\rm syst}}{R_{31}}$ \\\\ \\hline " << endl;
    output << "GeV &  & (\\%) & (\\%) &  & (\\%) & (\\%)  \\\\ \\hline" << endl;
    
	for(int ipt=0; ipt<fNpt; ipt++){
		//output << setprecision(2) << fPTbin[ipt] << "--" << setprecision(3)<<  fPTbin[ipt+1] << " & ";
        output << Form("%.0f--%.0f & ",fPTbin[ipt],fPTbin[ipt+1]);
		for(int ups=2; ups<=3; ups++){
            //if(ups==2)output << Form("%.1f & ",grName[xs_graph(iy,ups)]->GetX()[ipt]);
			//if(ups==2)output  << setprecision(3) << grName[xs_graph(iy,ups)]->GetX()[ipt] << " & ";
			
			int precision=set_precision(grName[xs_ratio(iy,ups)]->GetY()[ipt] );
			int precisionE=set_precisionE(grName[xs_ratio_stat(iy,ups)]->GetErrorYhigh(ipt)); 

            double stat=grName[xs_ratio_stat(iy,ups)]->GetErrorYhigh(ipt)/grName[xs_ratio(iy,ups)]->GetY()[ipt];
            double totP=grName[xs_ratio(iy,ups)]->GetErrorYhigh(ipt)/grName[xs_ratio(iy,ups)]->GetY()[ipt];
            double totM=grName[xs_ratio(iy,ups)]->GetErrorYhigh(ipt)/grName[xs_ratio(iy,ups)]->GetY()[ipt];

            double tot=TMath::Max(totP,totM);
            double sys=TMath::Sqrt(tot*tot-stat*stat);

            output << Form("%.2f & %.1f & %.1f ",grName[xs_ratio(iy,ups)]->GetY()[ipt],100*stat,100*sys);
            
			//output  << setprecision(3) << grName[xs_ratio(iy,ups)]->GetY()[ipt] << " & $\\pm$ ";
			//output  << setprecision(2) << 100*grName[xs_ratio_stat(iy,ups)]->GetErrorYhigh(ipt)/grName[xs_ratio(iy,ups)]->GetY()[ipt] ;
			//output  << setprecision(2) << " & $\\pm$ " << 100*grName[xs_ratio(iy,ups)]->GetErrorYhigh(ipt)/grName[xs_ratio(iy,ups)]->GetY()[ipt] << " (" << 100*grName[xs_ratio(iy,ups)]->GetErrorYlow(ipt)/grName[xs_ratio(iy,ups)]->GetY()[ipt] << ") ";
            if(ups==2) output << " & ";
		}
		output << "\\\\" << endl; 
	}
	
	output << "\\hline" << endl; 
	output << "\\end{tabular}" << endl; 
	output.close();
	
}


void ratio_tablehepdata(int iy){
	ofstream output; 
	output.open(Form("/uscms/home/btcarlso/hepdata/ratio_table_y%d.tex",iy)); 
	
    //output << "\\begin{tabular}{c|ccc|ccc}\\hline" << endl;
    //output << "$p_{\\rm T}$ & $R_{21}$ & $\\frac{\\asigma_{\\rm stat}}{R_{21}}$   & $\\frac{\\sigma_{\\rm syst}}{R_{21}}$ & $R_{31}$  & $\\frac{\\sigma_{\\rm stat}}{R_{31}}$   & $\\frac{\\sigma_{\\rm syst}}{R_{31}}$ \\\\ \\hline " << endl;
    //output << "GeV &  & (\\%) & (\\%) &  & (\\%) & (\\%)  \\\\ \\hline" << endl;
	output << "*reackey: P P --> UPSI(1S) X " << endl; 
	output << "*obskey: RN1/DPT" << endl; 
	output << "*dserror: +2.2,-2.2 PCT : luminosity uncertainty" << endl; 
	if(iy==0) output << "*qual: ABS(YRAP) : < 0.6" << endl;   
	if(iy==1) output << "*qual: ABS(YRAP) : 0.6 TO 1.2" << endl;   
	if(iy==2) output << "*qual: ABS(YRAP) : < 1.2" << endl;   
	output << "*qual: RE : P P --> UPSI(nS)  < MU+ MU -> X" << endl; 
	output << "*qual: SQRT(S) IN GEV : 7000.0" << endl;
	//output << "*qual: . : R21 : R31 " << endl; 
	output << "*yheader: R21 : R31" << endl; 
	output << "*xheader: PT IN GEV " << endl; 
    output << "*data: x : y : y " << endl; 
	for(int ipt=0; ipt<fNpt; ipt++){
        output << Form("%.1f TO %.1f ; ",fPTbin[ipt],fPTbin[ipt+1]);
		for(int ups=2; ups<=3; ups++){
					
			int precision=set_precision(grName[xs_ratio(iy,ups)]->GetY()[ipt] );
			int precisionE=set_precisionE(grName[xs_ratio_stat(iy,ups)]->GetErrorYhigh(ipt)); 
			
            double stat=grName[xs_ratio_stat(iy,ups)]->GetErrorYhigh(ipt)/grName[xs_ratio(iy,ups)]->GetY()[ipt];
            double totP=grName[xs_ratio(iy,ups)]->GetErrorYhigh(ipt)/grName[xs_ratio(iy,ups)]->GetY()[ipt];
            double totM=grName[xs_ratio(iy,ups)]->GetErrorYhigh(ipt)/grName[xs_ratio(iy,ups)]->GetY()[ipt];
			
            double tot=TMath::Max(totP,totM);
            double sys=TMath::Sqrt(tot*tot-stat*stat);
			
            output << Form("%.2f +- %.1f PCT (DSYS=+-%.1f PCT) ; ",grName[xs_ratio(iy,ups)]->GetY()[ipt],100*stat,100*sys);
                      
		}
		output << endl; 
	}
	output << "*dataend" << endl; 

	output.close();
	
}


int set_precision(double x){
	int precision=2;
    if(x>=100000) precision=6;
	if(x>=10000 && x<100000) precision=5;
	if(x>=1000 && x<10000) precision=4; 
	if(x>=100 && x<1000) precision=3;
	if(x>=10 && x<100) precision=2; 
	if(x>=1 && x<10) precision=2; 
	if(x>=0.1 && x<1) precision=1; 
	
	//cout << "x: " << x << " pre: "<< precision << endl; 
	
	return precision; 
}

int set_precisionE(double x){
	int precision=2; 
	if(x>=10000) precision=6; 
	if(x>=1000 && x<10000) precision=5; 
	if(x>=100 && x<1000) precision=3;
	if(x>=10 && x<100) precision=2; 
	if(x>=1 && x<10) precision=2; 
	if(x<=1) precision=1; 
	
	//cout << "xE: " << x << " pre: "<< precision << endl; 

	
	return precision; 
}

void yield_table(){
	
	
	for(int iy=0; iy<fNy; iy++){
		
		ofstream output; 
		output.open(Form("/uscms/home/btcarlso/AN-UPS/tdr2/notes/AN-13-088/trunk/tables/yields_table_y%d.tex",iy));
		int N=6; 
		output <<"\\begin{tabular}{"; 
		for(int i=0; i<=N; i++) output << "c"; 
		output << "}\\hline" << endl; 
		output << "$p_{T}$ & Yield(1S) & Yield(2S) & Yield(3S) & $\\chi^{2}$ & NDOF \\\\" << endl; 
		
		for(int ipt=1; ipt<=hName["yield_1S_y0_best"]->GetNbinsX(); ipt++){
			int pre=2; 
			if(hName["yield_1S_y0_best"]->GetBinCenter(ipt)>75) pre=3; 
			output << setprecision(pre) << hName["yield_1S_y0_best"]->GetBinCenter(ipt) << " & "; 
			for(int ups=1; ups<=3;ups++){
				
				TString name=yield_histogram(iy,ups,"best","");
				TString namesys=background_order_systematic(iy,Form("%dS",ups)); 
				int precision=set_precision(hName[name]->GetBinContent(ipt)); 
				int precisionE=set_precisionE(hName[name]->GetBinError(ipt)); 
				int precisionEsys=set_precisionE(hName[namesys]->GetBinContent(ipt)); 
				//output << setprecision(precision) << hName[name]->GetBinContent(ipt) << " $\\pm$ " << setprecision(precisionE) << hName[name]->GetBinError(ipt) << "[" << setprecision(3) << hName[namesys]->GetBinContent(ipt) << "] & ";
                if(hName[namesys]->GetBinContent(ipt) > 10) output << Form("%.0f $\\pm$ %.0f [%.0f] & ", hName[name]->GetBinContent(ipt), hName[name]->GetBinError(ipt), hName[namesys]->GetBinContent(ipt) );
                else output << Form("%.0f $\\pm$ %.0f [%.1f] & ", hName[name]->GetBinContent(ipt), hName[name]->GetBinError(ipt), hName[namesys]->GetBinContent(ipt) );
                
			}//ups
            if(iy<2) output << Form("%.0f & %.0f", hName[Form("chi2_y%d",iy)]->GetBinContent(ipt), hName[Form("NDOF_y%d",iy)]->GetBinContent(ipt)) << endl;
			//if(iy<2)output << setprecision(3) << hName[Form("chi2_y%d",iy)]->GetBinContent(ipt)<< " & " << hName[Form("NDOF_y%d",iy)]->GetBinContent(ipt);
			else output << "- & - "; 
			output << " \\\\ " << endl; 
		}//pt
		output << "\\hline" << endl << "\\end{tabular}" << endl; 
		output.close();
		}//iy 
					
}

void load_acceptance(){
	string method1S[]={"","Em","Ep","_longitudinal","_transverse","_unPol","_reco","_smeared"};

	string method[]={"","Em","Ep","_longitudinal","_transverse","_unPol"};
	int N=sizeof(method)/sizeof(string); 
	int N1=sizeof(method1S)/sizeof(string); 
	TH1F *hh=(TH1F*)eff1S->FindObjectAny("efftot_pt_TH1");
	hName[hh->GetName()]=hh;
	
	for(int i=0; i<N; i++){
		//if(i>2) continue;
		cout << Form("acceptance%s",method[i].c_str()) << endl; 
		TH1F *h1S=(TH1F*)eff1S->FindObjectAny(Form("acceptance%s",method[i].c_str())); 
		TString name=h1S->GetName();
		name=name+"_1S";
		h1S->SetName(name); 
		hName[h1S->GetName()]=h1S; 
		
		TH1F *h2S=(TH1F*)eff2S->FindObjectAny(Form("acceptance%s",method[i].c_str())); 
		name=h2S->GetName();
		name=name+"_2S";
		h2S->SetName(name); 
		hName[h2S->GetName()]=h2S; 
		
		TH1F *h3S=(TH1F*)eff3S->FindObjectAny(Form("acceptance%s",method[i].c_str())); 
		name=h3S->GetName();
		name=name+"_3S";
		h3S->SetName(name); 
		hName[h3S->GetName()]=h3S; 
	}
	
	for(int iy=0; iy<fNy; iy++){
		for(int i=0; i<N1;i++){
			cout << Form("acceptance%s_y%d",method[i].c_str(),iy) << endl; 
			TH1F *h1S=(TH1F*)eff1S->FindObjectAny(Form("acceptance%s_y%d",method1S[i].c_str(),iy)); 
			TString name=h1S->GetName();
			name=name+"_1S";
			h1S->SetName(name); 
			hName[h1S->GetName()]=h1S; 
		}
		
		for(int i=0; i<N; i++){
			TH1F *h2S=(TH1F*)eff2S->FindObjectAny(Form("acceptance%s_y%d",method[i].c_str(),iy)); 
			TString name=h2S->GetName();
			name=name+"_2S";
			h2S->SetName(name); 
			hName[h2S->GetName()]=h2S; 
			
			TH1F *h3S=(TH1F*)eff3S->FindObjectAny(Form("acceptance%s_y%d",method[i].c_str(),iy)); 
			name=h3S->GetName();
			name=name+"_3S";
			h3S->SetName(name); 
			hName[h3S->GetName()]=h3S; 
		}
		
	}
}

void load_sg(int iy){
	TString hist_name=Form("sg_y%d",iy);
	if(iy>1) hist_name="sg";
	hist_name=hist_name+"_TH1";

	
	TH1F *sg1S=(TH1F*)eff1S->FindObjectAny(hist_name); 
	TString name=sg1S->GetName(); 
	name=name+"_1S"; 
	sg1S->SetName(name); 
	hName[sg1S->GetName()]=sg1S; 
	Rebin(name); 
	
	TH1F *sg2S=(TH1F*)eff2S->FindObjectAny(hist_name); 
	name=sg2S->GetName(); 
	name=name+"_2S"; 
	sg2S->SetName(name); 
	hName[sg2S->GetName()]=sg2S; 
	Rebin(name);
	
	
	TH1F *sg3S=(TH1F*)eff3S->FindObjectAny(hist_name); 
	name=sg3S->GetName(); 
	name=name+"_3S"; 
	sg3S->SetName(name); 
	hName[sg3S->GetName()]=sg3S;
	
	Rebin(name); 
}

void Rebin(TString name){
	
	for(int i=1; i<=hName["bin_width_pt"]->GetNbinsX(); i++){
		hName["bin_width_pt"]->SetBinContent(i,hName["bin_width_pt"]->GetBinWidth(i)); 
	}
	//cout << "binwidth: " << hName["bin_width_pt"]->GetBinContent(13) << endl; 
	
	TH1F *h=(TH1F*)hName[name]->Rebin(fNpt,name+"_rebin",fPTbinD); 
	h->Divide(hName["bin_width_pt"]); 
	cout << "Rebin: " << h->GetName() << endl; 
	hName[h->GetName()]=h; 
}

void rebin_acceptance(){
	string method[]={"","Em","Ep","_longitudinal","_transverse","_unPol"};
	int N=sizeof(method)/sizeof(string); 
	for(int iy=0; iy<fNy; iy++){
		for(int ups=1; ups<=3; ups++){
			for(int i=0; i<N; i++){
				Rebin(Form("acceptance%s_y%d_%dS",method[i].c_str(),iy,ups)); 
			}
		}
	}
	
	Rebin("acceptance_reco_y0_1S"); 
	Rebin("acceptance_smeared_y0_1S"); 

	Rebin("acceptance_reco_y1_1S"); 
	Rebin("acceptance_smeared_y1_1S"); 
	
	cout << "|y| < 1.2 " << endl; 
	for(int ups=1; ups<=3; ups++){
		for(int i=0; i<N; i++){
			//if(i>2) continue;
			cout << Form("acceptance%s_%dS",method[i].c_str(),ups) << endl; 
			Rebin(Form("acceptance%s_%dS",method[i].c_str(),ups)); 
		}
	}
	
}

void weighted_acceptance(int iy, int ups){
	//cout << "iy:" << iy << " ups: " << ups << endl; 
	CreateHistogram(acc_name(iy,ups,"_reweighted"),"p_{T} re-weighted acceptance", "p_{T} [GeV]", "Events",fNpt,fPTbin); 
	//cout << "reweight acceptance: "<< endl; 
	for(int ipt=1; ipt<=hName[acc_name(iy,ups,"_reweighted")]->GetNbinsX(); ipt++){
		double ptlow=hName[acc_name(iy,ups,"_reweighted")]->GetBinCenter(ipt)-hName[acc_name(iy,ups,"_reweighted")]->GetBinWidth(ipt)/2; 
	//	cout << "pt start: " << ptlow << endl; 
		
		TString Acc1GeV=Form("acceptance_y%d_%dS",iy,ups);
		if(iy==2) Acc1GeV=Form("acceptance_%dS",ups);
		int start_bin=hName[Acc1GeV]->FindBin(ptlow); 
		int nSteps = static_cast<int>(hName[acc_name(iy,ups,"_reweighted")]->GetBinWidth(ipt));

	//	cout << "start bin: " << start_bin << " nSteps: "<< nSteps << endl; 
		
		double weightedAN=0;
		double weightedAD=0; 
		
		for(int i=start_bin; i<start_bin+nSteps; i++){
			double pt_center=hName[Acc1GeV]->GetBinCenter(i); 
			double w=f1Name[Form("xs_fit_expo_y%d_%dS",iy,ups)]->Eval(pt_center); 
			if(pt_center>20) w=f1Name[Form("xs_fit_power_y%d_%dS",iy,ups)]->Eval(pt_center);
			//cout << "pt: " << pt_center << " " << hName[Acc1GeV]->GetBinContent(i) << " w " << w << endl; 
			weightedAN+=hName[Acc1GeV]->GetBinContent(i)*w; 
			weightedAD+=w;
			
		}
		//cout << "acceptance: " << hName[acc_name(iy,ups,"")]->GetBinContent(ipt) << " corr: " << weightedAN/weightedAD << endl; 
		
		hName[acc_name(iy,ups,"_reweighted")]->SetBinContent(ipt,weightedAN/weightedAD); 
	}//end ipt loop 
}

void acceptance_table_smearing(int iy, int ups){
	cout << "acceptance table: " << iy << " ups: "<< ups << endl; 
	
	ofstream output; 
	output.open(Form("/uscms/home/btcarlso/AN-UPS/tdr2/notes/AN-13-088/trunk/tables/acceptance_smearing_y%d_%dS.tex",iy,ups)); 
	if(iy<2) {
		output << "\\begin{tabular}{cccc}" << endl; 
		output << "$p_{T}$ [GeV] &"<< Form(" Acceptance(%dS)",ups) << " & $\\delta_{Reco}$ [\\%] & $\\delta_{Smeared}$ [\\%] \\\\ \\hline" << endl;	
	}
	else  {
		output << "\\begin{tabular}{ccc}" << endl; 
		output << "$p_{T}$ [GeV] &"<< Form(" Acceptance(%dS)",ups) << " & $\\delta_{reweighted}$ [\\%] \\\\ \\hline" << endl;	
	}
	
	for(int ipt=1; ipt<=hName[acc_name(iy,ups,"")]->GetNbinsX(); ipt++){
		double Acc=hName[acc_name(iy,ups,"")]->GetBinContent(ipt); 
		double AccReco=hName[acc_name(iy,ups,"_reco")]->GetBinContent(ipt); 
		double AccSmeared=hName[acc_name(iy,ups,"_smeared")]->GetBinContent(ipt); 

		output << setprecision(3) << fPTbin[ipt-1] << "-" << fPTbin[ipt] << " & "; 
		output << setprecision(2) << Acc << " & ";
		output << setprecision(2) << 100*(AccReco-Acc)/Acc << " & ";
		output << setprecision(2) << 100*(AccSmeared-Acc)/Acc << " & "; 
		
		output << "\\\\" << endl; 
	}
	output << "\\end{tabular} " << endl; 
	output.close(); 
}


void acceptance_ratio(int iy, int ups, string mode){
    cout << "acceptance ratio: " << iy << " " << ups << "S " << mode << endl;
    string NN="RatioPol"+mode+"_unPol";
    TH1F *h=(TH1F*)hName[acc_name(iy,ups,mode)]->Clone(acc_name(iy,ups,NN));
    h->SetName(acc_name(iy,ups,NN));
    h->Divide(hName[acc_name(iy,ups,"_unPol")]);
    h->SetMinimum(0.);
    h->SetMaximum(2);
    hName[h->GetName()]=h;
    
    hName[acc_name(iy,ups,NN)]->GetYaxis()->SetTitle("A(pol)/A(un-pol)");
    
/*
    CreateCanvas(Form("Acceptance_Ratio_Pol_unPol_y%d_%dS",iy,ups),"",600,600);
    CName[Form("Acceptance_Ratio_Pol_unPol_y%d_%dS",iy,ups)]->cd();
    hName[acc_name(iy,ups,"RatioPol_unPol")]->Draw("histo");
    hName[acc_name(iy,ups,"RatioPol_unPol")]->GetYaxis()->SetTitleOffset(2);
    hName[acc_name(iy,ups,"RatioPol_unPol")]->Draw("histo");

    TString bt;
	if(iy==0)bt=	"|y| < 0.6";
	if(iy==1) bt="0.6 < |y| < 1.2";
	if(iy==2) bt="|y| < 1.2";
    TString upsT=Form("  #Upsilon(%dS)",ups);
    bt+=upsT;
	
	TLatex bin_t(0.2,0.8,bt);
	bin_t.SetNDC(kTRUE);
	bin_t.DrawClone();
 */

}
void acceptance_ratio_plot(){
    CreateCanvas(Form("Acceptance_Ratio_Pol_unPol_y%d_summary",0),"",600,600);

    CName[Form("Acceptance_Ratio_Pol_unPol_y%d_summary",0)]->cd();
    
    hName[acc_name(0,1,"RatioPolEp_unPol")]->SetLineColor(kBlack);
    hName[acc_name(0,1,"RatioPolEm_unPol")]->SetLineColor(kBlack);
    hName[acc_name(0,1,"RatioPolEp_unPol")]->SetLineStyle(kDashed);
    hName[acc_name(0,1,"RatioPolEm_unPol")]->SetLineStyle(kDashed);
    
    hName[acc_name(0,1,"RatioPol_longitudinal_unPol")]->SetLineColor(kRed);
    hName[acc_name(0,1,"RatioPol_transverse_unPol")]->SetLineColor(kBlue);

    hName[acc_name(0,2,"RatioPol_unPol")]->SetLineColor(kOrange);

    hName[acc_name(0,2,"RatioPolEp_unPol")]->SetLineColor(kRed);
    hName[acc_name(0,2,"RatioPolEm_unPol")]->SetLineColor(kRed);
    hName[acc_name(0,2,"RatioPolEp_unPol")]->SetLineStyle(kDashed);
    hName[acc_name(0,2,"RatioPolEm_unPol")]->SetLineStyle(kDashed);
    
    hName[acc_name(0,3,"RatioPol_unPol")]->SetLineColor(kGreen);
    hName[acc_name(0,3,"RatioPolEp_unPol")]->SetLineColor(kBlue);
    hName[acc_name(0,3,"RatioPolEm_unPol")]->SetLineColor(kBlue);
    hName[acc_name(0,3,"RatioPolEp_unPol")]->SetLineStyle(kDashed);
    hName[acc_name(0,3,"RatioPolEm_unPol")]->SetLineStyle(kDashed);

    
    hName[acc_name(0,1,"RatioPol_unPol")]->Draw("histo");
    hName[acc_name(0,1,"RatioPol_unPol")]->GetYaxis()->SetTitleOffset(2);
    hName[acc_name(0,1,"RatioPol_unPol")]->SetMinimum(0.5);
    hName[acc_name(0,1,"RatioPol_unPol")]->SetMaximum(1.5);

    hName[acc_name(0,1,"RatioPol_unPol")]->Draw("histo");
    hName[acc_name(0,1,"RatioPolEp_unPol")]->Draw("histo same");
    hName[acc_name(0,1,"RatioPolEm_unPol")]->Draw("histo same");


    hName[acc_name(0,2,"RatioPol_unPol")]->Draw("histo same");
    hName[acc_name(0,1,"RatioPol_longitudinal_unPol")]->Draw("histo same");
    hName[acc_name(0,1,"RatioPol_transverse_unPol")]->Draw("histo same");

    //hName[acc_name(0,2,"RatioPolEp_unPol")]->Draw("histo same");
    //hName[acc_name(0,2,"RatioPolEm_unPol")]->Draw("histo same");

    hName[acc_name(0,3,"RatioPol_unPol")]->Draw("histo same");
    //hName[acc_name(0,3,"RatioPolEp_unPol")]->Draw("histo same");
    //hName[acc_name(0,3,"RatioPolEm_unPol")]->Draw("histo same");

    
    TLegend L(0.5,0.7,0.8,0.85);
	L.SetFillColor(10);
	L.SetLineColor(10);
    L.SetBorderSize(0);
    L.AddEntry(hName[acc_name(0,1,"RatioPol_unPol")],"#Upsilon(1S)");
    L.AddEntry(hName[acc_name(0,1,"RatioPolEp_unPol")],"#pm 1 #sigma");
    L.AddEntry(hName[acc_name(0,1,"RatioPol_longitudinal_unPol")],"#Upsilon(1S)-L");
    L.AddEntry(hName[acc_name(0,1,"RatioPol_transverse_unPol")],"#Upsilon(1S)-T");

    L.AddEntry(hName[acc_name(0,2,"RatioPol_unPol")],"#Upsilon(2S)");
    L.AddEntry(hName[acc_name(0,3,"RatioPol_unPol")],"#Upsilon(3S)");

    L.DrawClone();

}


void acceptance_table(int iy, int ups){
	cout << "acceptance table: " << iy << " ups: "<< ups << endl; 

	ofstream output; 
	output.open(Form("/uscms/home/btcarlso/AN-UPS/tdr2/notes/AN-13-088/trunk/tables/acceptance_y%d_%dS.tex",iy,ups)); 
	if(iy<2) {
		output << "\\begin{tabular}{cccc}" << endl; 
		output << "$p_{T}$ [GeV] &"<< Form(" Acceptance(%dS)",ups) << " & $\\delta_{unpolarized}$ [\\%] & $\\delta_{reweighted}$ [\\%] \\\\ \\hline" << endl;	
	}
	else  {
		output << "\\begin{tabular}{ccc}" << endl; 
		output << "$p_{T}$ [GeV] &"<< Form(" Acceptance(%dS)",ups) << " & $\\delta_{reweighted}$ [\\%] \\\\ \\hline" << endl;	
	}

	for(int ipt=1; ipt<=hName[acc_name(iy,ups,"")]->GetNbinsX(); ipt++){
		double Acc=hName[acc_name(iy,ups,"")]->GetBinContent(ipt); 
		double deltaEp=100*TMath::Abs(hName[acc_name(iy,ups,"Ep")]->GetBinContent(ipt)-Acc)/Acc;
		double deltaEm=100*TMath::Abs(hName[acc_name(iy,ups,"Ep")]->GetBinContent(ipt)-Acc)/Acc; 

		double delta_reweighted=100*TMath::Abs(Acc-hName[acc_name(iy,ups,"_reweighted")]->GetBinContent(ipt))/Acc; 
	
		output << setprecision(3) << fPTbin[ipt-1] << "-" << fPTbin[ipt] << " & "; 
		output << setprecision(2) << Acc << " & ";
		if(iy<2){
			double delta_unPol=100*TMath::Abs(Acc-hName[acc_name(iy,ups,"_unPol")]->GetBinContent(ipt))/Acc;
			output << setprecision(2) << delta_unPol << " & ";
		}
		output << setprecision(2) << delta_reweighted;
	
		output << "\\\\" << endl; 
	}
	output << "\\end{tabular} " << endl; 
	output.close(); 
}

void acceptance_table_summary(int iy, int ups){
	cout << "acceptance table summary: " << iy << " ups: "<< ups << endl;
    cout << acc_name(iy,ups,"_unPol") << endl;

	ofstream output;
	output.open(Form("/uscms/home/btcarlso/AN-UPS/tdr2/notes/AN-13-088/trunk/tables/acceptance_summary_y%d_%dS.tex",iy,ups));
    output << "\\begin{tabular}{ccccccccc}" << endl;
    output << "\\hline \\hline" << endl;
    output << "$p_{\\rm T}$ [GeV] & $\\mathcal{A}$ & $\\mathcal{A}(\\sigma^{+})$ & $\\mathcal{A}(\\sigma^{-})$ & $\\mathcal{A}(\\text{unpol})$ & $\\mathcal{A}(T)$ & $\\mathcal{A}(L)$ \\\\ \\hline" << endl;
    
	   
	for(int ipt=1; ipt<=hName[acc_name(iy,ups,"")]->GetNbinsX(); ipt++){
		double Acc=hName[acc_name(iy,ups,"")]->GetBinContent(ipt);
		double AccEp=hName[acc_name(iy,ups,"Ep")]->GetBinContent(ipt);
        double AccEm=hName[acc_name(iy,ups,"Em")]->GetBinContent(ipt);
        
        double AccUnpol=hName[acc_name(iy,ups,"_unPol")]->GetBinContent(ipt);
        double AccT=hName[acc_name(iy,ups,"_transverse")]->GetBinContent(ipt);
        double AccL=hName[acc_name(iy,ups,"_longitudinal")]->GetBinContent(ipt);
        
        output << Form("%.0f--%.0f & ",fPTbin[ipt-1],fPTbin[ipt]);
        output << Form(" %.2f & %.2f & %.2f &", Acc, AccEp, AccEm);
        output << Form(" %.2f & %.2f & %.2f ", AccUnpol, AccT,AccL);
        
		output << "\\\\" << endl;
	}
    output << "\\hline\\hline" << endl;
	output << "\\end{tabular} " << endl;
	output.close(); 
}

void acceptance_table_summaryhepdata(int iy, int ups){
	//cout << "acceptance table summary: " << iy << " ups: "<< ups << endl;
    //cout << acc_name(iy,ups,"_unPol") << endl;

	ofstream output;
	output.open(Form("/uscms/home/btcarlso/hepdata/acceptance_summary_y%d_%dS.tex",iy,ups));
    //output << "\\begin{tabular}{ccccccccc}" << endl;
    //output << "\\hline \\hline" << endl;
	output << Form("*reackey: P P --> UPSI(%dS) X ",ups) << endl; 
	output << "*obskey: Acceptance" << endl; 
	if(iy==0) output << "*qual: ABS(YRAP) : < 0.6" << endl;   
	if(iy==1) output << "*qual: ABS(YRAP) : 0.6 TO 1.2" << endl;   
	if(iy==2) output << "*qual: ABS(YRAP) : < 1.2" << endl;   
	output << Form("*qual: RE : P P --> UPSI(%dS)  < MU+ MU -> X",ups) << endl; 
	output << "*qual: SQRT(S) IN GEV : 7000.0" << endl;
	//output << "*yheader: Acceptance" << endl; 
	output << "*yheader: $\\mathcal{A}$ : $\\mathcal{A}(\\sigma^{+})$ : $\\mathcal{A}(\\sigma^{-})$ : $\\mathcal{A}(\\text{unpol})$ : $\\mathcal{A}(T)$ : $\\mathcal{A}(L)$ " << endl;
	output << "*xheader: PT IN GEV " << endl; 
    
	output << "*data: x : y : y : y : y : y : y" << endl; 
	for(int ipt=1; ipt<=hName[acc_name(iy,ups,"")]->GetNbinsX(); ipt++){
		double Acc=hName[acc_name(iy,ups,"")]->GetBinContent(ipt);
		double AccEp=hName[acc_name(iy,ups,"Ep")]->GetBinContent(ipt);
        double AccEm=hName[acc_name(iy,ups,"Em")]->GetBinContent(ipt);
        
        double AccUnpol=hName[acc_name(iy,ups,"_unPol")]->GetBinContent(ipt);
        double AccT=hName[acc_name(iy,ups,"_transverse")]->GetBinContent(ipt);
        double AccL=hName[acc_name(iy,ups,"_longitudinal")]->GetBinContent(ipt);
        
        output << Form("%.1f TO %.1f ; ",fPTbin[ipt-1],fPTbin[ipt]);
        output << Form(" %.2f ; %.2f ; %.2f ;", Acc, AccEp, AccEm);
        output << Form(" %.2f ; %.2f ; %.2f ", AccUnpol, AccT,AccL);
        output << endl; 
		//output << "\\\\" << endl;
	}
    //output << "\\hline\\hline" << endl;
	//output << "\\end{tabular} " << endl;
	output << "*dataend" << endl; 
	output.close(); 
}

void scale_graph(TGraphAsymmErrors *gr, double SF){
	//Scale TGraphAsymmErrors 
	
	for (int i=0;i<gr->GetN();i++){
		gr->GetY()[i] *= SF;
		gr->GetEYhigh()[i]*=SF; 
		gr->GetEYlow()[i]*=SF;
	}
}

void plot_cms2010(int ups){
	//Table 38, BPH-11-001
	TGraphAsymmErrors *gr; 
	TGraphAsymmErrors *grRatio;
	
	double pTbins[]={0,2,5,8,10,13,16,18,22,38};
	
	int nR=9; 
	double ptR[9];
	double ptRE[9];
	
	//cout << "pT bins: " << endl; 
	
	for (int i=0; i<nR; i++) {
		ptR[i]=(pTbins[i+1]-pTbins[i])/2.+pTbins[i];
		ptRE[i]=pTbins[i+1]-ptR[i];
		//cout << pTbins[i+1] << "-" << pTbins[i] << "ptR: " << ptR[i] << " +/- " << ptRE[i] << endl;  
	}
	
	if(ups==1){
		
		//Table 58, page 97 from AN2011_107_v6
		//xs from http://arxiv.org/pdf/1303.5900v1.pdf
		// and exact numbers can be downloaded from the twiki: https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsBPH11001
		// though it's in a terrible format
		
		//xs measured for |y|<2.4, then scaled 
		double ptB[]={0,0.5,1,1.5,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,20,22,25,30,50};
		double pT[]={0.327,0.772,1.257,1.751,2.489,3.481,4.478,5.485,6.487,7.493,8.484,9.478,10.479,11.485,12.487,13.475,14.489,15.477,16.909,18.978,20.937,23.284,27.024,35.968};
		
		const int N=sizeof(ptB)/sizeof(double)-1;
		double ptP[N];
		double ptM[N];
		
		double xs[]={0.0859,	0.263,	0.374,	0.505,	1.16,	1.212,	1.084,	0.879,	0.68,	0.556,	0.419,	0.331,	0.238,	0.1793,	0.1451,	0.099,	0.075,	0.0595,	0.0732,0.05,0.0302,	0.0237,	0.0205,	0.0123};
		double xsEp[]={0.0094894678,0.0219317122,0.0297321375,0.0436806593,0.0921954446,0.0905207159,0.0680073525,0.0633245608,0.0481040539,0.0338378486,0.0256320112,0.0187882942,0.0134164079,0.0101355809,0.0090354856,0.0058830264,0.0054203321,0.0042720019,0.0048877398,0.0034669872,0.0024207437,0.0020615528,0.0019209373,0.001421267};
		
		double xsEm[]={0.0084599054,0.021023796,0.0278567766,0.0323109888,0.1118033989,0.0696419414,0.0596154342,0.071063352,0.0362353419,0.0273130006,0.0201246118,0.017,0.0116619038,0.0101355809,0.0078230429,0.0059682493,0.004580393,0.0038600518,0.0043011626,0.0032202484,0.0024207437,0.0020615528,0.0019209373,0.0013453624};
		
		for (int i=0; i<N; i++) {
			ptP[i]=(ptB[i+1]-pT[i]); 
			ptM[i]=(pT[i]-ptB[i]); 
			double BW=ptB[i+1]-ptB[i]; 
			xs[i]=xs[i]/BW; 
			//subtract lumi uncertainty
			double sysTot=xsEm[i];
			xsEm[i]=TMath::Sqrt(sysTot*sysTot-TMath::Power(0.04*sysTot,2)); 
			sysTot=xsEp[i];
			xsEp[i]=TMath::Sqrt(sysTot*sysTot-TMath::Power(0.04*sysTot,2)); 
			
			xsEm[i]=xsEm[i]/BW; 
			xsEp[i]=xsEp[i]/BW;
		}
		
		gr = new TGraphAsymmErrors(N,pT,xs,ptM,ptP,xsEp,xsEm);
		scale_graph(gr,1000*(0.5)); // scale from nb, then correction for rapidity scaling. 
		
		
	}//ups 1S
	
	
	if(ups==2){
		
		double ptB[]={0,1,2.5,4,5.5,7,8.5,10,11.5,13,14.5,16,18,19.5,22,26,42};
		double pT[]={0.656,1.785,3.214,4.708,6.217,7.711,9.206,10.691,12.213,13.698,15.222,16.875,18.757,20.647,23.688,31.3};
		
		
		const int N=sizeof(ptB)/sizeof(double)-1;
		double ptP[N];
		double ptM[N];
		
		
		double xs[]={0.0829,0.331,0.409,0.362,0.286,0.212,0.146,0.1123,0.0765,0.0519,0.0376,0.0373,0.0159,0.0204,0.0158,0.0126};
		double xsEp[]={0.0092962358,0.0386005181,0.0382753184,0.0304630924,0.0223606798,0.0187882942,0.0125299641,0.0097529483,0.0063513778,0.0043600459,0.0034205263,0.0033600595,0.0017691806,0.0021260292,0.0017029386,0.0015620499};
		//{13,12,12,12,13,13,13,15,14};
		double xsEm[]={0.0087920419,0.0245967478,0.0299666481,0.0259422435,0.0214709106,0.0144222051,0.01,0.0074027022,0.0056222771,0.0041231056,0.0039962482,0.0032015621,0.0017691806,0.0021260292,0.0017029386,0.0014866069};
		//{10,10,10,11,12,12,12,14,13};
		
		for (int i=0; i<N; i++) {
			ptP[i]=(ptB[i+1]-pT[i]); 
			ptM[i]=(pT[i]-ptB[i]); 
			double BW=ptB[i+1]-ptB[i]; 
			xs[i]=xs[i]/BW; 
			
			double sysTot=xsEm[i];
			xsEm[i]=TMath::Sqrt(sysTot*sysTot-TMath::Power(0.04*sysTot,2)); 
			sysTot=xsEp[i];
			xsEp[i]=TMath::Sqrt(sysTot*sysTot-TMath::Power(0.04*sysTot,2)); 
			
			xsEm[i]=xsEm[i]/BW; 
			xsEp[i]=xsEp[i]/BW;
		}
		
		gr = new TGraphAsymmErrors(N,pT,xs,ptM,ptP,xsEp,xsEm);
		scale_graph(gr,1000*(0.5)); // scale from nb, then correction for rapidity scaling. 
		
		double R[]={0.22,0.23,0.27,0.27,0.33,0.37,0.50,0.45,0.47};
		double Rp[]={0.03,0.03,0.04,0.04,0.04,0.04,0.06,0.05,0.09};
		double Rm[]={0.03,0.03,0.03,0.03,0.03,0.04,0.06,0.06,0.13};
		
		grRatio = new TGraphAsymmErrors(nR,ptR,R,ptRE,ptRE,Rp,Rm); 
		
		grRatio->GetYaxis()->SetLabelSize(0.025); 
		grRatio->GetYaxis()->SetTitle(Form("#sigma #timesBr(%dS)/#sigma #timesBr(%dS)",2,1)); 
		grRatio->GetXaxis()->SetRangeUser(0,100);
	}//ups 2S
	
	if(ups==3){
		double ptB[]={0,2.5,5,7.5,10,13,16,18,22,38};
		double pT[]={1.539,3.62,6.147,8.618,11.306,14.298,16.943,19.702,26.514};
		
		const int N=sizeof(ptB)/sizeof(double)-1;
		double ptP[N];
		double ptM[N];
		
		double xs[]={0.203,0.287,0.227,0.157,0.1127,0.0617,0.0227,0.0229,0.0185};
		double xsEp[]={0.0210950231,0.0281780056,0.0232594067,0.0196977156,0.0099849887,0.0051400389,0.0026248809,0.0024083189,0.0021260292};
		double xsEm[]={0.0194164878,0.0317804972,0.0205912603,0.0136014705,0.007981228,0.0052201533,0.0026248809,0.0026400758,0.0030413813};
		
		for (int i=0; i<N; i++) {
			ptP[i]=(ptB[i+1]-pT[i]); 
			ptM[i]=(pT[i]-ptB[i]); 
			double BW=ptB[i+1]-ptB[i]; 
			xs[i]=xs[i]/BW; 
			
			double sysTot=xsEm[i];
			xsEm[i]=TMath::Sqrt(sysTot*sysTot-TMath::Power(0.04*sysTot,2)); 
			sysTot=xsEp[i];
			xsEp[i]=TMath::Sqrt(sysTot*sysTot-TMath::Power(0.04*sysTot,2)); 
			
			xsEm[i]=xsEm[i]/BW; 
			xsEp[i]=xsEp[i]/BW;
		}
		
		gr = new TGraphAsymmErrors(N,pT,xs,ptM,ptP,xsEp,xsEm);
		scale_graph(gr,1000*(0.5)); // scale from nb, then correction for rapidity scaling. 
		
		//|y|<2.4 acceptance corrected 
		
		double R[]={0.11,0.10,0.12,0.15,0.20,0.26,0.31,0.28,0.32};
		double Rp[]={0.02,0.02,0.02,0.02,0.03,0.030,0.04,0.04,0.07};
		double Rm[]={0.02,0.01,0.02,0.02,0.02,0.03,0.04,0.04,0.09};
		
		grRatio=  new TGraphAsymmErrors(nR,ptR,R,ptRE,ptRE,Rp,Rm); 
		
		grRatio->GetYaxis()->SetLabelSize(0.025); 
		grRatio->GetYaxis()->SetTitle(Form("#sigma #timesBr(%dS)/#sigma #timesBr(%dS)",3,1)); 
		grRatio->GetXaxis()->SetRangeUser(0,100);
		
	}//ups 3S
	gr->GetXaxis()->SetTitle(x_label); 
	gr->GetYaxis()->SetTitle(xs_y); 
	gr->SetName(Form("xs_cms2010_%dS",ups)); 
	if(ups==2 || ups==3){
		grRatio->GetXaxis()->SetTitle("p_{T}(#mu#mu)");
		grRatio->SetMarkerStyle(8);
		grRatio->SetMarkerSize(0.5);
		grRatio->SetName(Form("ratio_cms2010_%dS-1S",ups)); 
		set_CMS2010(grRatio);
		grName[grRatio->GetName()]=grRatio; 

	}
	
	//store graphs of ratio and xs 
	set_CMS2010(gr); 
	grName[gr->GetName()]=gr; 
}

void draw_unc_ratio(int iy){
//hName[Rtot_uncertainty(iy,num,"Ep")]
	CreateCanvas(Form("uncertainty_ratio_y%d",iy),"",Cx,Cy); 
	CName[Form("uncertainty_ratio_y%d",iy)]->Divide(1,2); 
	CName[Form("uncertainty_ratio_y%d",iy)]->cd(1); 

	TLegend L(0.2,0.6,0.4,0.8);
	L.SetFillColor(10);
	L.SetLineColor(10);

	int num=2; 
	
	L.AddEntry(hName[Rstat_uncertainty(iy,num)],hName[Rstat_uncertainty(iy,num)]->GetTitle(),"L"); 
	L.AddEntry(hName[Reff_uncertainty(iy,num)],hName[Reff_uncertainty(iy,num)]->GetTitle(),"L"); 
	L.AddEntry(hName[Racc_uncertainty(iy,num)],hName[Racc_uncertainty(iy,num)]->GetTitle(),"L"); 
 
	
	hName[Rstat_uncertainty(iy,num)]->SetTitle(""); 
	
	hName[Rstat_uncertainty(iy,num)]->Draw(); 
	hName[Reff_uncertainty(iy,num)]->Draw("same"); 
	hName[Racc_uncertainty(iy,num)]->Draw("same"); 
	L.DrawClone("same");
	
	TLatex txt1(0.2,0.85,Form("ratio #Upsilon(%dS)/#Upsilon(1S)",num)); 
	txt1.SetNDC(kTRUE);
	txt1.DrawClone(); 
	draw_header();

	TString bt;
	if(iy==0)bt=	"|y| < 0.6"; 
	if(iy==1) bt="0.6 < |y| < 1.2"; 
	if(iy==2) bt="|y| < 1.2";
    
	TLatex bin_t(0.2,0.8,bt); 
	bin_t.SetNDC(kTRUE);
	bin_t.DrawClone(); 
	

	 num=3;
	
	CName[Form("uncertainty_ratio_y%d",iy)]->cd(2); 
	hName[Rstat_uncertainty(iy,num)]->SetTitle("");
	
	hName[Rstat_uncertainty(iy,num)]->Draw(); 
	hName[Reff_uncertainty(iy,num)]->Draw("same"); 
	hName[Racc_uncertainty(iy,num)]->Draw("same"); 
	
	TLatex txt2(0.2,0.85,Form("ratio #Upsilon(%dS)/#Upsilon(1S)",num)); 
	txt2.SetNDC(kTRUE);
	txt2.DrawClone(); 
	
}

void draw_unc(int iy){

	CreateCanvas(Form("uncertainty_y%d",iy),"",Cx,Cy);
	CName[Form("uncertainty_y%d",iy)]->Divide(2,2); 
	
	TLegend L(0.2,0.2,0.8,0.8);
	L.SetFillColor(10);
	L.SetLineColor(10);
	
	for(int ups=1; ups<=3; ups++){
			
		CName[Form("uncertainty_y%d",iy)]->cd(ups);
		
		if(ups==1){
			L.AddEntry(hName[stat_uncertainty(iy,ups)],hName[stat_uncertainty(iy,ups)]->GetTitle(),"L"); 
			L.AddEntry(hName[eff_uncertainty(iy,ups)],hName[eff_uncertainty(iy,ups)]->GetTitle(),"L"); 
			L.AddEntry(hName[acc_uncertainty(iy,ups)],hName[acc_uncertainty(iy,ups)]->GetTitle(),"L"); 
		}
		
		hName[stat_uncertainty(iy,ups)]->SetTitle(""); 
		hName[stat_uncertainty(iy,ups)]->Draw(); 
		hName[eff_uncertainty(iy,ups)]->Draw("same"); 
		hName[acc_uncertainty(iy,ups)]->Draw("same"); 

		TLatex txt(0.6,0.8,Form("#Upsilon(%dS)",ups)); 
		txt.DrawClone(); 
		draw_header();
	
	}
	CName[Form("uncertainty_y%d",iy)]->cd(4);
	L.DrawClone(); 
	
	TString bt;
	if(iy==0)bt=	"|y| < 0.6"; 
	if(iy==1) bt="0.6 < |y| < 1.2"; 
	if(iy==2) bt="|y| < 1.2"; 
	
	TLatex bin_t(0.2,0.8,bt); 
	bin_t.SetNDC(kTRUE);
	bin_t.DrawClone(); 
	
}

void compare_ratios(int num){
    
    int den=1;
	TLegend Leg(0.2,0.65,0.35,0.85);
	Leg.SetFillColor(10);
	Leg.SetLineColor(10);
	Leg.SetBorderSize(0);
    
    CreateCanvas(Form("RatioComparison_%dS1S",num),"",Cx,Cy);
    CName[Form("RatioComparison_%dS1S",num)]->cd();
    
    grName[xs_ratio(2,num)]->Draw("ap");
	grName[xs_ratio(2,num)]->GetYaxis()->SetRangeUser(0,1);
	grName[xs_ratio(2,num)]->GetXaxis()->SetRangeUser(10,100); // set to 0 for CMS 2010
	grName[xs_ratio(2,num)]->SetMarkerSize(1.5);
	grName[xs_ratio(2,num)]->Draw("ap");
    Leg.AddEntry(grName[xs_ratio(2,num)],"CMS 4.9 fb^{-1}","p");
    
    grName[Form("ratio_cms2010_%dS-1S",num)]->SetMarkerStyle(24);
    
	grName[Form("ratio_cms2010_%dS-1S",num)]->SetFillStyle(3944);//3013
    grName[Form("ratio_cms2010_%dS-1S",num)]->Draw("p same");
	grName[Form("ratio_cms2010_%dS-1S",num)]->Draw("e2 same");
    Leg.AddEntry(grName[Form("ratio_cms2010_%dS-1S",num)],"CMS 36 pb^{-1}","efp");
    
    grName[Form("AtlasRatio_%dS1S_1p2",num)]->SetMarkerStyle(25);
    grName[Form("AtlasRatio_%dS1S_1p2",num)]->SetMarkerColor(kRed);
    grName[Form("AtlasRatio_%dS1S_1p2",num)]->DrawClone("p same");
    
    Leg.AddEntry(grName[Form("AtlasRatio_%dS1S_1p2",num)],"ATLAS","p");
    
    Leg.DrawClone();

    
}

void draw_ratio(int num){
	//num 2S or 3S state 
	//denominator is 1S
	int den=1; 
	//cout << "draw ratio: " << num << "S/" << den << "S" << endl; 

	
	TLegend LegC(0.5,0.2,0.8,0.5);
	LegC.SetFillColor(10);
	LegC.SetLineColor(10);
	LegC.SetBorderSize(0); 
	
	
	CreateCanvas(Form("R_%d-1S_ycomparison_plot",num),"",Cx,Cy); 
	CName[Form("R_%d-1S_ycomparison_plot",num)]->cd(); 
	
	
	grName[xs_ratio(0,num)]->SetMarkerStyle(22); 
	grName[xs_ratio(0,num)]->SetMarkerSize(1);
	grName[xs_ratio(0,num)]->SetMarkerColor(kRed); 
	
	LegC.AddEntry(grName[xs_ratio(0,num)],"|y| < 0.6","LEP");
	LegC.AddEntry(grName[xs_ratio(1,num)],"0.6 < |y| < 1.2","LEP");
	
	
	grName[xs_ratio(0,num)]->Draw("ap"); 
	grName[xs_ratio(0,num)]->GetYaxis()->SetRangeUser(0,1.25); 
	grName[xs_ratio(0,num)]->GetXaxis()->SetRangeUser(10,100); // to overlay CMS 2010 set to 0
	grName[xs_ratio(0,num)]->SetMarkerSize(1.5);
	grName[xs_ratio(0,num)]->Draw("ap");
	//grName[xs_ratio(0,num)]->Draw("|| same");
	
	grName[xs_ratio(1,num)]->Draw("p same"); 
	//grName[xs_ratio_stat(1,num)]->Draw("|| same"); 
	
	LegC.DrawClone(); 
	
	CMS_lumi(CName[Form("R_%d-1S_ycomparison_plot",num)],1,0); 
	//draw_header(); 
	

	CreateCanvas(Form("R_%d-1S_ycombined_plot",num),"",Cx,Cy); 
	CName[Form("R_%d-1S_ycombined_plot",num)]->cd(); 
	
	TF1 *f = new TF1("f", "expo",10,20); 
	
	double p0=f1Name[Form("xs_fit_expo_y%d_%dS",2,num)]->GetParameter(0)-f1Name[Form("xs_fit_expo_y%d_%dS",2,den)]->GetParameter(0);
	double p1=f1Name[Form("xs_fit_expo_y%d_%dS",2,num)]->GetParameter(1)-f1Name[Form("xs_fit_expo_y%d_%dS",2,den)]->GetParameter(1);
	
	f->SetParameters(p0,p1); 
	
	TF1 *ratio=new TF1("ratio",ratio_function,20,100,6); 
	//Par[0] = Anum, Par[1]=Cnum, Par[2]=Alphanum
	//Par[3] = Aden, Par[4]=Cden, Par[5]=Alphaden
	
	ratio->SetParNames("Anum","alphanum","Cnum","Aden","alphaden","Cden"); 
	
	ratio->SetParameters(f1Name[Form("xs_fit_power_y%d_%dS",2,num)]->GetParameter(0),f1Name[Form("xs_fit_power_y%d_%dS",2,num)]->GetParameter(1),
						 f1Name[Form("xs_fit_power_y%d_%dS",2,num)]->GetParameter(2), f1Name[Form("xs_fit_power_y%d_%dS",2,den)]->GetParameter(0),
						 f1Name[Form("xs_fit_power_y%d_%dS",2,den)]->GetParameter(1),f1Name[Form("xs_fit_power_y%d_%dS",2,den)]->GetParameter(2)); 
	
	//ratio->Print();
	
	ratio->SetLineColor(kBlack); 
	ratio->SetLineWidth(2); 
	
	grName[xs_ratio(2,num)]->Draw("ap"); 
	grName[xs_ratio(2,num)]->GetYaxis()->SetRangeUser(0,1); 
	grName[xs_ratio(2,num)]->GetXaxis()->SetRangeUser(10,100); // set to 0 for CMS 2010
	grName[xs_ratio(2,num)]->SetMarkerSize(1.5);
	grName[xs_ratio(2,num)]->Draw("ap"); 
	//grName[xs_ratio(2,num)]->Draw("|| same"); 
	
	
	TLegend Leg(0.18,0.72,0.52,0.89);
	Leg.SetFillColor(10);
	Leg.SetLineColor(10);
	Leg.SetBorderSize(0); 
	Leg.SetTextSize(0.05);
	
	TLegend LegP(0.33,0.16,0.74,0.35);
	LegP.SetFillColor(10);
	LegP.SetLineColor(10);
	LegP.SetBorderSize(0);
	LegP.SetTextSize(0.05);
    
    f->SetLineWidth(2);
	f->SetLineStyle(kDashed);
    //: 10 < p_{T} < 20 GeV"
	LegP.AddEntry((TF1*)f->Clone("tmp"),"Ratio of exponential fits", "L"); 
	LegP.AddEntry(ratio,"Ratio of power-law fits", "L"); 
	Leg.AddEntry(grName[xs_ratio(2,num)],"CMS 2011 |y| < 1.2","LEP"); 
	//Leg.AddEntry(grName[Form("ratio_cms2010_%dS-1S",num)],"CMS 2010 |y| < 2.4","f");

	TString labelAB="(a)";
	if(num==3) labelAB="(b)";

	TLatex aLabel(0.39,0.82,labelAB); 
	aLabel.SetNDC(kTRUE);
	aLabel.SetTextSize(0.08);
	aLabel.DrawClone();
	

	f->DrawClone("same");
	f->SetRange(0,10); 
	f->DrawClone("same");
	ratio->DrawClone("same"); 
	
	//TLatex RatioLabel(0.67,0.835,Form("#Upsilon(%dS) / #Upsilon(1S)",num));
	TLatex RatioLabel(0.6,0.81,"|y| < 1.2");
	RatioLabel.SetNDC(kTRUE); 
	RatioLabel.SetTextSize(0.08); 
	
	grName[Form("ratio_cms2010_%dS-1S",num)]->SetFillStyle(3944);//3013 
	//grName[Form("ratio_cms2010_%dS-1S",num)]->Draw("e2 same");
	
	CMS_lumi(CName[Form("R_%d-1S_ycombined_plot",num)],1,10); 
	//draw_header(); 
	//Leg.DrawClone();
	LegP.DrawClone();
	RatioLabel.DrawClone();
	
	delete f; 
	delete ratio; 
}

TString divide_graphs(TString num, TString den){
	TString ratio_name=num+"_"+den; 
	TGraphAsymmErrors *gr=(TGraphAsymmErrors*)grName[num]->Clone(ratio_name); 
	grName[ratio_name]=gr; 
	for(int ipt=0; ipt<grName[num]->GetN(); ipt++){
		double ratio=grName[num]->GetY()[ipt]/grName[den]->GetY()[ipt]; 
		double numE=grName[num]->GetErrorYhigh(ipt)/grName[num]->GetY()[ipt];
		double denE=grName[den]->GetErrorYhigh(ipt)/grName[den]->GetY()[ipt];

		double ratioEhigh=TMath::Sqrt(numE*numE+denE*denE)*ratio;
		
		 numE=grName[num]->GetErrorYlow(ipt)/grName[num]->GetY()[ipt];
		 denE=grName[den]->GetErrorYlow(ipt)/grName[den]->GetY()[ipt];
		double ratioElow=TMath::Sqrt(numE*numE+denE*denE)*ratio;
		
		grName[ratio_name]->GetY()[ipt]=ratio;
		grName[ratio_name]->SetPointEYhigh(ipt,ratioEhigh);
		grName[ratio_name]->SetPointEYlow(ipt,ratioElow);

	}
	
	return ratio_name; 
	
}

void bin_center(int iy,int ups){
	cout << "iy, ups: " << iy << " " << ups << endl; 
	for (int ipt=0; ipt<grName[xs_graph(iy,ups)]->GetN(); ipt++) {
		double pt1=grName[xs_graph(iy,ups)]->GetX()[ipt]-grName[xs_graph(iy,ups)]->GetEXlow()[ipt];
		double pt2=grName[xs_graph(iy,ups)]->GetX()[ipt]+grName[xs_graph(iy,ups)]->GetEXhigh()[ipt];

		double deltapt=pt2-pt1; 
		
		double Int=f1Name[Form("xs_fit_expo_y%d_%dS",iy,ups)]->Integral(pt1,pt2)/deltapt;
		if(pt1>=20) Int=f1Name[Form("xs_fit_power_y%d_%dS",iy,ups)]->Integral(pt1,pt2)/deltapt;
		
		double pttmp=pt1; 
		double step_size=deltapt/100; 
		
		double delta[100];
		
		for(int i=0; i<100; i++){
			delta[i] = TMath::Abs(f1Name[Form("xs_fit_expo_y%d_%dS",iy,ups)]->Eval(pttmp)-Int); 
			if(pt1>=20) delta[i] = TMath::Abs(f1Name[Form("xs_fit_power_y%d_%dS",iy,ups)]->Eval(pttmp)-Int);  
			pttmp+=step_size; 
		}		
		int min=TMath::LocMin(100,delta); 
		double pt_mean=pt1+step_size*static_cast<double>(min); 
		
		double DPT=grName[xs_graph(iy,ups)]->GetX()[ipt]-pt_mean; //shift between the unweighted mean and the weighted mean
		
		grName[xs_graph(iy,ups)]->GetX()[ipt]=pt_mean; 
		grName[xs_graph(iy,ups)]->SetPointEXhigh(ipt,grName[xs_graph(iy,ups)]->GetEXhigh()[ipt]+DPT);//shift the bin edges
		grName[xs_graph(iy,ups)]->SetPointEXlow(ipt,grName[xs_graph(iy,ups)]->GetEXlow()[ipt]-DPT);
		
		grName[xs_graph_stat(iy,ups)]->GetX()[ipt]=pt_mean; 
		grName[xs_graph_stat(iy,ups)]->SetPointEXhigh(ipt,grName[xs_graph(iy,ups)]->GetEXhigh()[ipt]+DPT);//shift the bin edges
		grName[xs_graph_stat(iy,ups)]->SetPointEXlow(ipt,grName[xs_graph(iy,ups)]->GetEXlow()[ipt]-DPT);
		
		cout << "min: " << min << " <pT> " << DPT << endl;  
		
		//cout << "Mean pt: " << (pt1+pt2)/2 << " " << pt_mean << endl; 
		
	}
	
	
	
}

void plot_orderDist(){
	CreateCanvas("Order_Distribution","",Cx,Cy); 
	CName["Order_Distribution"]->Divide(2,2); 
	
	CreateHistogram("Order_Distribution_y0","","Order","N_{bins}",9,-0.5,9.5); 
	CreateHistogram("Order_Distribution_y1","","Order","N_{bins}",9,-0.5,9.5); 
	int iy=0; 
	for(int i=1; i<=hName[Form("Background_order_y%d_",iy)]->GetNbinsX();i++){
		hName[Form("Order_Distribution_y%d",iy)]->Fill(hName[Form("Background_order_y%d_",iy)]->GetBinContent(i)); 
	}
	iy=1; 
	for(int i=1; i<=hName[Form("Background_order_y%d_",iy)]->GetNbinsX();i++){
		hName[Form("Order_Distribution_y%d",iy)]->Fill(hName[Form("Background_order_y%d_",iy)]->GetBinContent(i)); 
	}
    
    TLatex txt;
    txt.SetNDC(kTRUE);
    
    for(int i=1; i<=hName["Background_order_y0_"]->GetNbinsX(); i++){
        hName["Background_order_y0_"]->SetBinError(i,0);
        hName["Background_order_y1_"]->SetBinError(i,0);
    }
    hName["Background_order_y0_"]->SetMarkerColor(kRed);
    hName["Background_order_y1_"]->SetMarkerColor(kBlue);
    
	
    hName["Background_order_y0_"]->SetMarkerStyle(22);
    hName["Background_order_y1_"]->SetMarkerStyle(22);
    
    CName["Order_Distribution"]->cd(1);
	hName["Background_order_y0_"]->SetMinimum(-0.5);
    hName["Background_order_y0_"]->DrawClone("histo ");
	hName["Background_order_y0_"]->DrawClone("p same");

    txt.DrawLatex(0.65,0.8, "|y| < 0.6");
    
	CName["Order_Distribution"]->cd(2); 
	hName["Background_order_y1_"]->SetMinimum(-0.5);
    hName["Background_order_y1_"]->DrawClone("histo ");
	hName["Background_order_y1_"]->DrawClone("p same");

    txt.DrawLatex(0.65,0.8, "0.6 < |y| < 1.2");
	
	CName["Order_Distribution"]->cd(3);
	hName["Order_Distribution_y0"]->DrawClone("histo");

	CName["Order_Distribution"]->cd(4); 
	hName["Order_Distribution_y1"]->DrawClone("histo"); 
	
}

void fit_nuisance(int iy){
	TF1 *fcw=new TF1("fcw","pol1",10,100); 
	TF1 *fdeltaM=new TF1("fdeltaM","pol1",10,100); 

	fcw->SetParameters(1.0,0);
	fdeltaM->SetParameters(0,0); 
	
	hName[Form("cw_y%d_best",iy)]->Fit("fcw","R0"); 	
	hName[Form("deltaM_y%d_best",iy)]->Fit("fdeltaM","R0"); 


	TGraphAsymmErrors *cw_pull = new TGraphAsymmErrors(hName[Form("cw_y%d_best",iy)]); 
	cw_pull->SetName(Form("cw_y%d_pull",iy)); 
	grName[cw_pull->GetName()]=cw_pull; 
	TGraphAsymmErrors *deltaM_pull = new TGraphAsymmErrors(hName[Form("deltaM_y%d_best",iy)]); 
	deltaM_pull->SetName(Form("deltaM_y%d_pull",iy)); 
	grName[deltaM_pull->GetName()]=deltaM_pull; 
	
	TString name=Form("cw_y%d_pull",iy); 
	for(int ipt=0; ipt<cw_pull->GetN(); ipt++){
		double pt=grName[name]->GetX()[ipt]; 
		double data=grName[name]->GetY()[ipt]; 
		double dataEh=grName[name]->GetErrorYhigh(ipt);
		double dataEl=grName[name]->GetErrorYlow(ipt);
		
		double dataE=(dataEh+dataEl)/2;
		
		double fit=fcw->Eval(pt); 
				
		double pull=(data-fit)/dataE;
		double pullEh=pull*(dataEh/data); 
		double pullEl=pull*(dataEl/data); 
		
		grName[name]->GetY()[ipt]=pull; 
		grName[name]->SetPointEYhigh(ipt,pullEh); 
		grName[name]->SetPointEYlow(ipt,pullEl); 
	}
	
	CreateCanvas(Form("cw_y%d_fit",iy),"",Cx,Cy); 
	CName[Form("cw_y%d_fit",iy)]->Divide(1,2); 
	CName[Form("cw_y%d_fit",iy)]->cd(1); 
	hName[Form("cw_y%d_best",iy)]->Draw("E1"); 
	fcw->DrawClone("same"); 
	CName[Form("cw_y%d_fit",iy)]->cd(2); 
	grName[name]->Draw("AP");
	grName[name]->GetXaxis()->SetRangeUser(10,100);
	grName[name]->GetYaxis()->SetTitle("Pull"); 
	grName[name]->DrawClone("ap"); 
	
	name=Form("deltaM_y%d_pull",iy); 
	for(int ipt=0; ipt<grName[name]->GetN(); ipt++){
		double pt=grName[name]->GetX()[ipt]; 
		double data=grName[name]->GetY()[ipt]; 
		double dataEh=grName[name]->GetErrorYhigh(ipt);
		double dataEl=grName[name]->GetErrorYlow(ipt);
		
		double dataE=(dataEh+dataEl)/2;
		
		double fit=fdeltaM->Eval(pt); 
		
		double pull=(data-fit)/dataE;
		double pullEh=pull*(dataEh/data); 
		double pullEl=pull*(dataEl/data); 
		
		grName[name]->GetY()[ipt]=pull; 
		grName[name]->SetPointEYhigh(ipt,pullEh); 
		grName[name]->SetPointEYlow(ipt,pullEl); 
	}
	

	
	CreateCanvas(Form("deltaM_y%d_fit",iy),"",Cx,Cy); 
	CName[Form("deltaM_y%d_fit",iy)]->Divide(1,2); 
	CName[Form("deltaM_y%d_fit",iy)]->cd(1);
	hName[Form("deltaM_y%d_best",iy)]->Draw("E1"); 
	fdeltaM->DrawClone("same"); 
	CName[Form("deltaM_y%d_fit",iy)]->cd(2);
	grName[name]->Draw("AP");
	grName[name]->GetXaxis()->SetRangeUser(10,100); 
	grName[name]->GetYaxis()->SetTitle("Pull"); 
	grName[name]->DrawClone("ap"); 
	
	delete fcw;
	delete fdeltaM;
}


void xs_pull(int iy, int ups){

	TGraphAsymmErrors *gr = (TGraphAsymmErrors*)grName[xs_graph(iy,ups)]->Clone(Form("xs_fit_pull_y%d_%dS",iy,ups)); 
	for(int ipt=0; ipt<gr->GetN(); ipt++){
		double pt=grName[xs_graph(iy,ups)]->GetX()[ipt]; 
		double data=grName[xs_graph(iy,ups)]->GetY()[ipt]; 
		double dataEh=grName[xs_graph(iy,ups)]->GetErrorYhigh(ipt);
		double dataEl=grName[xs_graph(iy,ups)]->GetErrorYlow(ipt);

		double dataE=(dataEh+dataEl)/2;
		
		double fit=f1Name[Form("xs_fit_expo_y%d_%dS",iy,ups)]->Eval(pt); 
		if(pt>20) fit=f1Name[Form("xs_fit_power_y%d_%dS",iy,ups)]->Eval(pt);
		
		
		double pull=(data-fit)/dataE;
		double pullEh=pull*(dataEh/data); 
		double pullEl=pull*(dataEl/data); 

		gr->GetY()[ipt]=pull; 
		gr->SetPointEYhigh(ipt,pullEh); 
		gr->SetPointEYlow(ipt,pullEl); 
	}
	gr->GetYaxis()->SetTitle("Pull"); 
	gr->GetYaxis()->SetRangeUser(-3,3); 
	gr->GetYaxis()->SetNdivisions(6);
	gr->GetYaxis()->SetLabelSize(0.1);
	grName[gr->GetName()]=gr; 

}

void draw_xs_fits(int iy,int ups){
	//int iy=2; 
	//cout << "draw xs fits " << endl; 
	xs_pull(iy, ups); 
	cout << "Create canvas: " << endl; 
	CreateCanvas(Form("xs_fit_y%d_%dS",iy,ups),"",Cx,Cy); 
	CName[Form("xs_fit_y%d_%dS",iy,ups)]->cd();

	cout << "Create pads: " << endl; 
	double eps=0.05;
	TPad *fit_pad = new TPad(Form("fit_pad_y%d_%dS",iy,ups),"xs-fit",0.0,0.3,1,1);
	TPad *pull_pad = new TPad(Form("pullpad_y%d_%dS",iy,ups),"pull",0,0.,1,0.3+eps);
	 
	fit_pad->Draw();
	pull_pad->Draw();
	fit_pad->SetTopMargin(0.1);
	fit_pad->SetBottomMargin(0.1);
	pull_pad->SetTopMargin(0.1);
	pull_pad->SetBottomMargin(0.3); 
	
	fit_pad->cd();
	
	gPad->SetLogy();
	//cout << "draw fits: "<< endl; 
	
	TGraphAsymmErrors *gr = (TGraphAsymmErrors*)grName[xs_graph(iy,ups)]->Clone("tmp");
	
	gr->Draw("ap");
	gr->GetXaxis()->SetRangeUser(10,100);
	double min=f1Name[Form("xs_fit_power_y%d_%dS",iy,ups)]->Eval(97);
	double max=2.5*f1Name[Form("xs_fit_expo_y%d_%dS",iy,ups)]->Eval(10);
	gr->GetYaxis()->SetRangeUser(min,max);
	gr->GetYaxis()->SetTitleOffset(0.85);
	gr->GetYaxis()->SetTitleSize(0.07);
	gr->GetYaxis()->SetLabelSize(0.06); 
	gr->GetXaxis()->SetLabelSize(0); 

	gr->DrawClone("ap");

	f1Name[Form("xs_fit_power_y%d_%dS",iy,ups)]->DrawClone("same"); 
	f1Name[Form("xs_fit_expo_y%d_%dS",iy,ups)]->SetLineStyle(kDashed);
	f1Name[Form("xs_fit_expo_y%d_%dS",iy,ups)]->DrawClone("same"); 

	
	TLegend Leg(0.38,0.49,0.65,0.69);
	Leg.SetFillColor(10);
	Leg.SetLineColor(10);
	Leg.SetBorderSize(0); 
	Leg.SetTextSize(0.055);
	
	Leg.AddEntry(f1Name[Form("xs_fit_power_y%d_%dS",iy,ups)],"power-law fit:   p_{T} > 20 GeV","L");
	Leg.AddEntry(f1Name[Form("xs_fit_expo_y%d_%dS",iy,ups)],"exponential fit:   10 < p_{T} < 20 GeV","L");
	Leg.DrawClone(); 
	
	TString title=Form("#Upsilon(%dS), |y| < 0.6",ups);
	if(iy==1)title = Form("#Upsilon(%dS), 0.6 < |y| < 1.2",ups); 
	if(iy==2) title = Form("#Upsilon(%dS), |y| < 1.2",ups);
	
	TLatex text(0.6,0.8,title); 
	text.SetNDC(kTRUE);
	text.SetTextSize(0.08);
	text.DrawClone(); 
	
	TString labelAB="(a)";
	if(iy==1 || iy==2) labelAB = "(b)";
	
	TLatex aLabel(0.39,0.79,labelAB); 
	aLabel.SetNDC(kTRUE);
	aLabel.SetTextSize(0.08);
	if(ups==1) aLabel.DrawClone(); 
	
	f1Name[Form("xs_fit_power_y%d_%dS",iy,ups)]->SetRange(10,20);
	f1Name[Form("xs_fit_power_y%d_%dS",iy,ups)]->SetLineStyle(kDashed);
	
	f1Name[Form("xs_fit_expo_y%d_%dS",iy,ups)]->SetRange(20,35);
	
	//f1Name[Form("xs_fit_power_y%d_%dS",iy,ups)]->DrawClone("same"); 
	f1Name[Form("xs_fit_expo_y%d_%dS",iy,ups)]->DrawClone("same"); 
	
	/*
	TLatex L1_;
	TLatex L2_; 
	L1_.SetNDC(kTRUE);
	L2_.SetNDC(kTRUE);
	L1_.SetTextSize(0.07);
	L2_.SetTextSize(0.07);
	L1_.DrawLatex(0.17,0.92, cms_pre); 
	L2_.DrawLatex(0.26,0.92, lumi); 
	*/
	//draw_header();
	//cout << "draw pulls: "<< endl; 
	pull_pad->cd();
	grName[Form("xs_fit_pull_y%d_%dS",iy,ups)]->Draw("ap"); 
	grName[Form("xs_fit_pull_y%d_%dS",iy,ups)]->GetXaxis()->SetRangeUser(10,100); 
	grName[Form("xs_fit_pull_y%d_%dS",iy,ups)]->GetXaxis()->SetTitleSize(0.15);
	grName[Form("xs_fit_pull_y%d_%dS",iy,ups)]->GetXaxis()->SetLabelSize(0.2);
	grName[Form("xs_fit_pull_y%d_%dS",iy,ups)]->GetXaxis()->SetTickLength(0.05);
	grName[Form("xs_fit_pull_y%d_%dS",iy,ups)]->GetXaxis()->SetLabelSize(0.12);

	
	grName[Form("xs_fit_pull_y%d_%dS",iy,ups)]->GetYaxis()->SetTitleOffset(0.45);
	grName[Form("xs_fit_pull_y%d_%dS",iy,ups)]->GetYaxis()->SetTitleSize(0.15);
	grName[Form("xs_fit_pull_y%d_%dS",iy,ups)]->Draw("ap"); 
	lumiTextSize*=1.5;
	cmsTextSize*=1.5; 
	CMS_lumi(fit_pad,1,10); 
	//CMS_lumi(CName[Form("xs_fit_y%d_%dS",iy,ups)],1,10); 
	cmsTextSize*=0.666666;
	lumiTextSize*=0.666666;
	delete gr; 
	
}


void draw_xs_fits_wNLO(int iy,int ups){
	//int iy=2; 
	//cout << "draw xs fits " << endl; 
	xs_pull(iy, ups); 
	cout << "Create canvas: " << endl; 
	CreateCanvas(Form("xs_fit_y%d_%dS_wNLO",iy,ups),"",Cx,Cy); 
	CName[Form("xs_fit_y%d_%dS_wNLO",iy,ups)]->cd();
	
	cout << "Create pads: " << endl; 
	double eps=0.05;
	TPad *fit_pad = new TPad(Form("fit_pad_y%d_%dS_wNLO",iy,ups),"xs-fit",0.0,0.3,1,1);
	TPad *pull_pad = new TPad(Form("pullpad_y%d_%dS_wNLO",iy,ups),"pull",0,0.,1,0.3+eps);
	
	fit_pad->Draw();
	pull_pad->Draw();
	fit_pad->SetTopMargin(0.1);
	fit_pad->SetBottomMargin(0.1);
	pull_pad->SetTopMargin(0.1);
	pull_pad->SetBottomMargin(0.3); 
	
	fit_pad->cd();
	
	gPad->SetLogy();
	//cout << "draw fits: "<< endl; 
	
	TGraphAsymmErrors *gr = (TGraphAsymmErrors*)grName[xs_graph(iy,ups)]->Clone("tmp");
	
	gr->Draw("ap");
	gr->GetXaxis()->SetRangeUser(10,100);
	double min=f1Name[Form("xs_fit_power_y%d_%dS",iy,ups)]->Eval(97);
	double max=2.5*f1Name[Form("xs_fit_expo_y%d_%dS",iy,ups)]->Eval(10);
	gr->GetYaxis()->SetRangeUser(min,max);
	gr->GetYaxis()->SetTitleOffset(0.85);
	gr->GetYaxis()->SetTitleSize(0.07);
	gr->GetYaxis()->SetLabelSize(0.06); 
	gr->GetXaxis()->SetLabelSize(0); 
	
	gr->DrawClone("ap");
	grName[Form("gr_NLOtheory_%dS_y%d",ups,iy)]->SetLineStyle(kDotted);
	grName[Form("gr_NLOtheory_%dS_y%d",ups,iy)]->DrawClone("l same"); 

	f1Name[Form("xs_fit_power_y%d_%dS",iy,ups)]->SetRange(20,100);
	f1Name[Form("xs_fit_power_y%d_%dS",iy,ups)]->SetLineStyle(1); 
 	f1Name[Form("xs_fit_power_y%d_%dS",iy,ups)]->DrawClone("same"); 

	f1Name[Form("xs_fit_expo_y%d_%dS",iy,ups)]->SetRange(10,20);
	f1Name[Form("xs_fit_expo_y%d_%dS",iy,ups)]->SetLineStyle(1);
	f1Name[Form("xs_fit_expo_y%d_%dS",iy,ups)]->DrawClone("same"); 
	
	TLegend Leg(0.38,0.49,0.65,0.69);
	Leg.SetFillColor(10);
	Leg.SetLineColor(10);
	Leg.SetBorderSize(0); 
	Leg.SetTextSize(0.055);
	
	Leg.AddEntry(f1Name[Form("xs_fit_power_y%d_%dS",iy,ups)],"power-law fit:   p_{T} > 20 GeV","L");
	Leg.AddEntry(f1Name[Form("xs_fit_expo_y%d_%dS",iy,ups)],"exponential fit:   10 < p_{T} < 20 GeV","L");
	Leg.AddEntry(grName[Form("gr_NLOtheory_%dS_y%d",ups,iy)],"NLO theory", "L"); 

	Leg.DrawClone(); 
	
	TString title=Form("#Upsilon(%dS), |y| < 0.6",ups);
	if(iy==1)title = Form("#Upsilon(%dS), 0.6 < |y| < 1.2",ups); 
	if(iy==2) title = Form("#Upsilon(%dS), |y| < 1.2",ups);
	
	TLatex text(0.6,0.8,title); 
	text.SetNDC(kTRUE);
	text.SetTextSize(0.08);
	text.DrawClone(); 
	
	TString labelAB="(a)";
	if(iy==1 || iy==2) labelAB = "(b)";
	
	TLatex aLabel(0.39,0.79,labelAB); 
	aLabel.SetNDC(kTRUE);
	aLabel.SetTextSize(0.08);
	if(ups==1) aLabel.DrawClone(); 
	
	//f1Name[Form("xs_fit_power_y%d_%dS",iy,ups)]->SetRange(10,20);
	//f1Name[Form("xs_fit_power_y%d_%dS",iy,ups)]->SetLineStyle(kDashed);
	
	f1Name[Form("xs_fit_expo_y%d_%dS",iy,ups)]->SetRange(20,35);
	f1Name[Form("xs_fit_expo_y%d_%dS",iy,ups)]->SetLineStyle(kDashed);
	//f1Name[Form("xs_fit_power_y%d_%dS",iy,ups)]->DrawClone("same"); 
	f1Name[Form("xs_fit_expo_y%d_%dS",iy,ups)]->DrawClone("same"); 
	/*
	 TLatex L1_;
	 TLatex L2_; 
	 L1_.SetNDC(kTRUE);
	 L2_.SetNDC(kTRUE);
	 L1_.SetTextSize(0.07);
	 L2_.SetTextSize(0.07);
	 L1_.DrawLatex(0.17,0.92, cms_pre); 
	 L2_.DrawLatex(0.26,0.92, lumi); 
	 */
	//draw_header();
	//cout << "draw pulls: "<< endl; 
	pull_pad->cd();
	grName[Form("xs_fit_pull_y%d_%dS",iy,ups)]->Draw("ap"); 
	grName[Form("xs_fit_pull_y%d_%dS",iy,ups)]->GetXaxis()->SetRangeUser(10,100); 
	grName[Form("xs_fit_pull_y%d_%dS",iy,ups)]->GetXaxis()->SetTitleSize(0.15);
	grName[Form("xs_fit_pull_y%d_%dS",iy,ups)]->GetXaxis()->SetLabelSize(0.2);
	grName[Form("xs_fit_pull_y%d_%dS",iy,ups)]->GetXaxis()->SetTickLength(0.05);
	grName[Form("xs_fit_pull_y%d_%dS",iy,ups)]->GetXaxis()->SetLabelSize(0.12);
	
	
	grName[Form("xs_fit_pull_y%d_%dS",iy,ups)]->GetYaxis()->SetTitleOffset(0.45);
	grName[Form("xs_fit_pull_y%d_%dS",iy,ups)]->GetYaxis()->SetTitleSize(0.15);
	grName[Form("xs_fit_pull_y%d_%dS",iy,ups)]->Draw("ap"); 
	lumiTextSize*=1.5;
	cmsTextSize*=1.5; 
	CMS_lumi(fit_pad,1,10); 
	//CMS_lumi(CName[Form("xs_fit_y%d_%dS",iy,ups)],1,10); 
	cmsTextSize*=0.666666;
	lumiTextSize*=0.666666;
	delete gr; 
	
}



void draw_rapidity_ratio(){
	CreateCanvas("rapidity_ratios","",Cx,Cy); 
	CName["rapidity_ratios"]->cd(); 
	
	TF1 *f=new TF1("f","expo",10,20); 
	int ups=1; 
	double p0=f1Name[Form("xs_fit_expo_y%d_%dS",0,ups)]->GetParameter(0)-f1Name[Form("xs_fit_expo_y%d_%dS",1,ups)]->GetParameter(0);
	double p1=f1Name[Form("xs_fit_expo_y%d_%dS",0,ups)]->GetParameter(1)-f1Name[Form("xs_fit_expo_y%d_%dS",1,ups)]->GetParameter(1);

	f->SetParameters(p0,p1); 
	f->SetLineStyle(kDotted); 
	
	TF1 *ratio=new TF1("ratio",ratio_function,20,100,6); 
	
	ratio->SetParameters(f1Name[Form("xs_fit_power_y%d_%dS",0,ups)]->GetParameter(0),f1Name[Form("xs_fit_power_y%d_%dS",0,ups)]->GetParameter(1),
						 f1Name[Form("xs_fit_power_y%d_%dS",0,ups)]->GetParameter(2), f1Name[Form("xs_fit_power_y%d_%dS",1,ups)]->GetParameter(0),
						 f1Name[Form("xs_fit_power_y%d_%dS",1,ups)]->GetParameter(1),f1Name[Form("xs_fit_power_y%d_%dS",1,ups)]->GetParameter(2)); 
	
	ratio->SetLineColor(kBlue); 
	ratio->SetLineWidth(1.5);
	
	TString ratio_name=divide_graphs(xs_graph_stat(0,1),xs_graph_stat(1,1));
	TGraphAsymmErrors *gr=(TGraphAsymmErrors*)grName[ratio_name]->Clone(ratio_name+"_stat");
    grName[ratio_name+"_stat"]=gr;
    
    //set systematic to 1.4%
    
    for(int i=0; i<grName[ratio_name]->GetN(); i++){
        int ipt=i+1;
        double y=grName[ratio_name]->GetY()[i];
        double eyh=grName[ratio_name]->GetEYhigh()[i];
        double eyl=grName[ratio_name]->GetEYlow()[i];
        
        double An=hName[acc_name(0,1,"")]->GetBinContent(ipt);
        double Ad=hName[acc_name(1,1,"")]->GetBinContent(ipt);

        double AnEp=TMath::Abs(hName[acc_name(0,1,"Ep")]->GetBinContent(ipt)-An);
        double AdEp=TMath::Abs(hName[acc_name(1,1,"Ep")]->GetBinContent(ipt)-Ad);
        
        double AnEm=TMath::Abs(hName[acc_name(0,1,"Em")]->GetBinContent(ipt)-An);
        double AdEm=TMath::Abs(hName[acc_name(1,1,"Em")]->GetBinContent(ipt)-Ad);
    
        double Ansigma = TMath::Max(AnEp,AnEm);
        double Adsigma = TMath::Max(AdEp,AdEm);

        double sigA2=TMath::Power(Ansigma/An,2)+TMath::Power(Adsigma/Ad,2); // uncertainty due to propagating unceranties on A, (dA/A)^2
        
        //cout << "A: " << An/Ad << " sigA: " << sigA2 << endl;
        
        //cout << " stat high: " <<eyh << " " << eyl << endl;
        
        eyh=y*TMath::Sqrt(TMath::Power(eyh/y,2)+0.014*0.014 + sigA2);
        eyl=y*TMath::Sqrt(TMath::Power(eyl/y,2)+0.014*0.014 + sigA2);
        
        //cout << "Eyh: "<< eyh << " " << eyl << endl;
        
        grName[ratio_name]->GetEYhigh()[i]=eyh;
        grName[ratio_name]->GetEYlow()[i]=eyl;
    }
    
    
    //compute weighted mean
    /*
    double *y=grName[ratio_name]->GetY();
    double *w=grName[ratio_name]->GetEYhigh();
    double *w2=grName[ratio_name]->GetEYlow();
    int N=sizeof(y)/sizeof(double);
    

    double sigma=0;
    
    for(int i=0; i<N; i++){
        double ww1=TMath::Power(w[i],2);
        double ww2=TMath::Power(w2[i],2);

        w[i]=TMath::Power(TMath::Max(ww1,ww2),-1); //set w=1/sigma^2
        sigma+=w[i];
    }
    cout << "check to see if this number has systematics or not!" << endl;
    cout << "mean: " << TMath::Mean(N,y,w) << " +/- " << TMath::Power(sigma,-0.5) <<  endl;
  
*/
    
	grName[ratio_name]->SetLineColor(kBlack);
	grName[ratio_name]->Draw("ap"); 
	grName[ratio_name]->GetYaxis()->SetRangeUser(0.75,1.75); 
	grName[ratio_name]->GetXaxis()->SetRangeUser(10,100);
	grName[ratio_name]->SetMarkerSize(1.5);
	
	grName[ratio_name]->GetYaxis()->SetTitle("#scale[0.8]{#frac{d#sigma}{dp_{T}}[|y| < 0.6]/#frac{d#sigma}{dp_{T}}[0.6 < |y| < 1.2]}");
	grName[ratio_name]->DrawClone("ap");	
	
	TLegend Leg(0.2,0.67,0.55,0.87);
	Leg.SetFillColor(10);
	Leg.SetLineColor(10);
	Leg.SetBorderSize(0); 
	Leg.SetTextSize(0.035); 
	Leg.AddEntry(grName[ratio_name],"CMS 2011","LEP"); 
	Leg.AddEntry((TF1*)f->Clone("tmp"),"exponential fit:   10 < p_{T} < 20 GeV","L"); 
	Leg.AddEntry(ratio,"power-law fit:   p_{T} > 20 GeV","L"); 
	
	//TLatex text(0.2,0.65,"|y| < 0.6 / 0.6 < |y| < 1.2"); 
	TLatex text(0.22,0.6,"#Upsilon(1S)"); 
	text.SetNDC(kTRUE);
	
	text.DrawClone();
	
	
	//grName["Ratio_CMS2010_y0p8_y0p4"]->Draw("p same"); 
	//grName["Ratio_CMS2010_y1p2_y0p4"]->Draw("p same"); 

	f->DrawClone("same");
	f->SetRange(0,10); 
	f->SetLineStyle(kDashed); 
	f->DrawClone("same");
	ratio->DrawClone("same"); 
	
	double atlasp0=f1Name["xs_atlas_y0_y1"]->GetParameter(0);
	double atlasp0E=f1Name["xs_atlas_y0_y1"]->GetParError(0);
	/*
	f1Name["xs_atlas_y0_y1"]->SetParameter(0,atlasp0+atlasp0E); 
	f1Name["xs_atlas_y0_y1"]->SetLineStyle(kDashed); 
	f1Name["xs_atlas_y0_y1"]->DrawClone("same"); 
	f1Name["xs_atlas_y0_y1"]->SetParameter(0,atlasp0-atlasp0E); 
	f1Name["xs_atlas_y0_y1"]->DrawClone("same"); 
	f1Name["xs_atlas_y0_y1"]->SetLineStyle(1);
	f1Name["xs_atlas_y0_y1"]->SetParameter(0,atlasp0); 
	f1Name["xs_atlas_y0_y1"]->DrawClone("same"); 
*/

	//Leg.AddEntry(grName["Ratio_CMS2010_y0p8_y0p4"],"CMS 2010, 36 pb^{-1}, |y|<0.4/0.4<|y|<0.8", "LEP");
	//Leg.AddEntry(grName["Ratio_CMS2010_y1p2_y0p4"],"CMS 2010, 36 pb^{-1}, |y|<0.4/0.8<|y|<1.2", "LEP");
//	Leg.AddEntry(f1Name["xs_atlas_y0_y1"],"ATLAS 2011, |y| < 0.6/0.6 < |y| < 1.2, 4 < p_{T} < 70 GeV", "L");
	draw_header();
	
	Leg.DrawClone();
	delete f; 
	delete ratio;
}

void table_rapidity_ratio(){
    ofstream output;
	output.open("/uscms/home/btcarlso/AN-UPS/tdr2/notes/AN-13-088/trunk/tables/rapidity_ratio.tex");
    int N=2;
	output <<"\\begin{tabular}{";
   	for(int i=0; i<=N; i++) output << "c";
    output << "}";
    output << " $p_{T}$ (GeV) & radio \\\\ \\hline" << endl;

    
    TString name=xs_graph_stat(0,1)+"_"+xs_graph_stat(1,1)+"_stat";
    for(int ipt=0; ipt<22; ipt++){
        output << Form("%.1f -- %.1f", fPTbin[ipt],fPTbin[ipt+1]) << " & ";
        output << Form("%.2f $\\pm$ %.2f ",grName[name]->GetY()[ipt],TMath::Max(grName[name]->GetEYhigh()[ipt],grName[name]->GetEYlow()[ipt]));
        output << " \\\\" << endl;
    }
    output << "\\end{tabular}";
}

void set_CMS2010(TGraphAsymmErrors *gr){
	
	gr->SetMarkerStyle(5); 
	gr->SetMarkerSize(1); //0.2
	gr->SetLineWidth(0.025); 
	gr->SetMarkerColor(kBlue); 
	gr->SetFillStyle(3005);
	gr->SetFillColor(9); 
	gr->SetLineColor(kBlue); 
}

void set_CMS2011(TGraphAsymmErrors *gr,string MODE, int ups){
	
	gr->SetMarkerStyle(20); 
	
	gr->SetMarkerSize(0.75); 
	gr->SetLineWidth(0.1); 
	gr->SetLineColor(kBlack); 
	if(MODE=="ratio") gr->SetMarkerSize(0.75); 
	
	if(MODE=="stat"){
		gr->SetLineColor(kRed); 
		//clear_points(gr,10,40);
	}
}

void atlas_cmsratio(){
    int iy=2;
    //2 GeV bins from 10-30
    //10-30GeV
    //>=40 GeV
    
    TH1F *h1S = new TH1F("h1S_atlas","",grName[Form("Atlas_xs_%dS",1)]->GetN()-1,grName[Form("Atlas_xs_%dS",1)]->GetX());
    TH1F *h2S = new TH1F("h2S_atlas","",grName[Form("Atlas_xs_%dS",2)]->GetN()-1,grName[Form("Atlas_xs_%dS",2)]->GetX());
    TH1F *h3S = new TH1F("h3S_atlas","",grName[Form("Atlas_xs_%dS",3)]->GetN()-1,grName[Form("Atlas_xs_%dS",3)]->GetX());

    TH1F *h1Scms = new TH1F("h1S_cms","",grName[xs_graph(iy,1)]->GetN()-1,grName[xs_graph(iy,1)]->GetX());
    TH1F *h2Scms = new TH1F("h2S_cms","",grName[xs_graph(iy,2)]->GetN()-1,grName[xs_graph(iy,2)]->GetX());
    TH1F *h3Scms = new TH1F("h3S_cms","",grName[xs_graph(iy,3)]->GetN()-1,grName[xs_graph(iy,3)]->GetX());
    
    TGraphAsymmErrors *gr1ScmsSys = (TGraphAsymmErrors*)grName[xs_graph(iy,1)]->Clone("gr1S_sys");
    TGraphAsymmErrors *gr2ScmsSys = (TGraphAsymmErrors*)grName[xs_graph(iy,2)]->Clone("gr2S_sys");
    TGraphAsymmErrors *gr3ScmsSys = (TGraphAsymmErrors*)grName[xs_graph(iy,3)]->Clone("gr3S_sys");

    /*
    for(int i=1; i<=h1S->GetNbinsX(); i++){
        cout << "ATLAS: " << grName[Form("Atlas_xs_%dS",1)]->GetY()[i-1] << endl;
        
    }
    */
    
    for(int i=1; i<=h1S->GetNbinsX(); i++){
        double BW=h1S->GetBinWidth(i);
        h1S->SetBinContent(i,grName[Form("Atlas_xs_%dS",1)]->GetY()[i-1]*BW);

    }
    
    for(int i=1; i<=h2S->GetNbinsX(); i++){
        double BW=h2S->GetBinWidth(i);
        h2S->SetBinContent(i,grName[Form("Atlas_xs_%dS",2)]->GetY()[i-1]*BW);
    }
    
    for(int i=1; i<=h3S->GetNbinsX(); i++){
        double BW=h3S->GetBinWidth(i);
        h3S->SetBinContent(i,grName[Form("Atlas_xs_%dS",3)]->GetY()[i-1]*BW);
    }
    //cms
    for(int i=1; i<=h1Scms->GetNbinsX(); i++){
        double BW=h1Scms->GetBinWidth(i);
        h1Scms->SetBinContent(i,grName[xs_graph(iy,1)]->GetY()[i-1]*BW);
        double unc=TMath::Max(grName[xs_graph(iy,1)]->GetEYhigh()[i-1]/grName[xs_graph(iy,1)]->GetY()[i-1],grName[xs_graph(iy,1)]->GetEYlow()[i-1]/grName[xs_graph(iy,1)]->GetY()[i-1]);
        gr1ScmsSys->GetY()[i-1]=1;
        gr1ScmsSys->GetEYhigh()[i-1]=unc;
        gr1ScmsSys->GetEYlow()[i-1]=unc;

        h2Scms->SetBinContent(i,grName[xs_graph(iy,2)]->GetY()[i-1]*BW);
        unc=TMath::Max(grName[xs_graph(iy,2)]->GetEYhigh()[i-1]/grName[xs_graph(iy,2)]->GetY()[i-1],grName[xs_graph(iy,2)]->GetEYlow()[i-1]/grName[xs_graph(iy,2)]->GetY()[i-1]);
        gr2ScmsSys->GetY()[i-1]=1;
        gr2ScmsSys->GetEYhigh()[i-1]=unc;
        gr2ScmsSys->GetEYlow()[i-1]=unc;
        
        h3Scms->SetBinContent(i,grName[xs_graph(iy,3)]->GetY()[i-1]*BW);
        unc=TMath::Max(grName[xs_graph(iy,3)]->GetEYhigh()[i-1]/grName[xs_graph(iy,3)]->GetY()[i-1],grName[xs_graph(iy,3)]->GetEYlow()[i-1]/grName[xs_graph(iy,3)]->GetY()[i-1]);
        gr3ScmsSys->GetY()[i-1]=1;
        gr3ScmsSys->GetEYhigh()[i-1]=unc;
        gr3ScmsSys->GetEYlow()[i-1]=unc;
    }
    
    double ptbins[]={10,12,14,16,18,20,22,24,26,28,30,40,70};
    int Nbins=sizeof(ptbins)/sizeof(double)-1;
    
    h1S=(TH1F*)h1S->Rebin(Nbins, "h1S_atlas",ptbins);
    h2S=(TH1F*)h2S->Rebin(Nbins, "h2S_atlas",ptbins);
    h3S=(TH1F*)h3S->Rebin(Nbins, "h3S_atlas",ptbins);
    
    h1Scms=(TH1F*)h1Scms->Rebin(Nbins, "h1S_cms",ptbins);
    h2Scms=(TH1F*)h2Scms->Rebin(Nbins, "h2S_cms",ptbins);
    h3Scms=(TH1F*)h3Scms->Rebin(Nbins, "h3S_cms",ptbins);
    /*
    for(int i=1; i<=h1S->GetNbinsX(); i++){
        double BW=h1S->GetBinWidth(i);
        h1S->SetBinContent(i,h1S->GetBinContent(i)/BW);
        h2S->SetBinContent(i,h2S->GetBinContent(i)/BW);
        h3S->SetBinContent(i,h3S->GetBinContent(i)/BW);

    }
    
    for(int i=1; i<=h1S->GetNbinsX(); i++){
        cout << "CMS: " << h1Scms->GetBinContent(i) << " ATLAS: " << h1S->GetBinContent(i) << endl;
    
    }
*/
    h1Scms->Divide(h1S);
    h2Scms->Divide(h2S);
    h3Scms->Divide(h3S);

    CreateCanvas("AtlasCMSRatio_1S","",Cx,Cy);
    CreateCanvas("AtlasCMSRatio_2S","",Cx,Cy);
    CreateCanvas("AtlasCMSRatio_3S","",Cx,Cy);

    gr1ScmsSys->SetFillColor(kBlue);
    gr1ScmsSys->SetFillStyle(3004);
    
    gr2ScmsSys->SetFillColor(kBlue);
    gr2ScmsSys->SetFillStyle(3004);
    
    gr3ScmsSys->SetFillColor(kBlue);
    gr3ScmsSys->SetFillStyle(3004);
    
    CName["AtlasCMSRatio_1S"]->cd();
    h1Scms->Draw("histo p");
    gr1ScmsSys->Draw("E2 same");
    
    CName["AtlasCMSRatio_2S"]->cd();
    h2Scms->Draw("histo E0");
//    h2Scms->Draw("EX0 same");
    gr2ScmsSys->Draw("E2 same");


    CName["AtlasCMSRatio_3S"]->cd();
    h3Scms->Draw("histo p");
  //  h3Scms->Draw("EX0 same");
    gr3ScmsSys->Draw("E2 same");

    
}

void draw_xs_atlas(){
	CreateCanvas("xs_atlas_comp","",Cx,Cy); 
	CName["xs_atlas_comp"]->cd();
	gPad->SetLogy();
	
	double min=0.01*0.3*grName[xs_graph(2,3)]->GetY()[grName[xs_graph(2,1)]->GetN()-1];
	double max=6*grName[Form("Atlas_xs_%dS",1)]->Eval(4);
	
	//xs_graph(a,b) where a is the rapidity bin, 0, 1, or 2 and b is the Upsilon state 1-3
	
	TGraphAsymmErrors *gr2S=(TGraphAsymmErrors*)grName[xs_graph(2,2)]->Clone("gr2S");
	TGraphAsymmErrors *gr3S=(TGraphAsymmErrors*)grName[xs_graph(2,3)]->Clone("gr3S");
	/*
    for(int i=0; i<gr3S->GetN(); i++){
        cout << "ATLAS: " << " pt: " << grName[Form("Atlas_xs_%dS",3)]->GetX()[i+10] << " xs " << grName[Form("Atlas_xs_%dS",3)]->GetY()[i+10];
        cout <<  " CMS pt " << grName[xs_graph(2,3)]->GetX()[i] << " xs " <<grName[xs_graph(2,3)]->GetY()[i];
        cout << " Ratio(1S): " << grName[xs_graph(2,1)]->GetY()[i]/grName[Form("Atlas_xs_%dS",1)]->GetY()[i+10] << " ";
        cout << " Ratio(2S): " << grName[xs_graph(2,2)]->GetY()[i]/grName[Form("Atlas_xs_%dS",2)]->GetY()[i+10] << " ";
        cout << " Ratio(3S): " << grName[xs_graph(2,3)]->GetY()[i]/grName[Form("Atlas_xs_%dS",3)]->GetY()[i+10] << endl;
    }
    */
	scale_graph(gr2S,0.1);
	scale_graph(gr3S,0.01);

	
	grName[Form("Atlas_xs_%dS",1)]->SetMarkerSize(0.25);
	grName[Form("Atlas_xs_%dS",1)]->SetFillStyle(3002);
	grName[Form("Atlas_xs_%dS",1)]->SetFillColor(kRed); 
	
	grName[Form("Atlas_xs_%dS",2)]->SetMarkerSize(0.25);
	grName[Form("Atlas_xs_%dS",2)]->SetFillStyle(3002);
	grName[Form("Atlas_xs_%dS",2)]->SetFillColor(kRed); 
	
	grName[Form("Atlas_xs_%dS",3)]->SetMarkerSize(0.25);
	grName[Form("Atlas_xs_%dS",3)]->SetFillStyle(3002);
	grName[Form("Atlas_xs_%dS",3)]->SetFillColor(kRed); 
	
	
	grName[xs_graph(2,1)]->DrawClone("ap");
	grName[xs_graph(2,1)]->GetYaxis()->SetTitleSize(0.08);
	grName[xs_graph(2,1)]->GetYaxis()->SetTitleOffset(0.9);
	grName[xs_graph(2,1)]->GetYaxis()->SetRangeUser(min,max); 
	grName[xs_graph(2,1)]->GetXaxis()->SetRangeUser(0,100); 
	grName[xs_graph(2,1)]->DrawClone("ap");

	grName[Form("Atlas_xs_%dS",1)]->Draw("e2 same"); 
//	grName[Form("Atlas_xs_%dS",1)]->Draw("p0 same");


	scale_graph(grName[Form("Atlas_xs_%dS",2)],0.1); 
	grName[Form("Atlas_xs_%dS",2)]->DrawClone("e2 same");
 //   grName[Form("Atlas_xs_%dS",2)]->DrawClone("p0 same");

	scale_graph(grName[Form("Atlas_xs_%dS",3)],0.01);
	grName[Form("Atlas_xs_%dS",3)]->DrawClone("e2 same");
//    grName[Form("Atlas_xs_%dS",3)]->DrawClone("p0 same");

	gr2S->DrawClone("p same");
	gr3S->DrawClone("p same");
	
	TLatex label1(0.583,0.505,"#Upsilon(1S)"); 
	label1.SetNDC(kTRUE); 
	
	TLatex label2(0.707,0.31,"#Upsilon(2S) #times 0.1"); 
	label2.SetNDC(kTRUE);
	
	TLatex label3(0.311,0.27,"#Upsilon(3S) #times 0.01"); 
	label3.SetNDC(kTRUE); 
	
	
	label1.SetTextSize(0.04);
	label2.SetTextSize(0.04);
	label3.SetTextSize(0.04);
	
	label1.DrawClone();
	label2.DrawClone();
	label3.DrawClone(); 
	
	TLegend Leg(0.51,0.66,0.71,0.86);
	Leg.SetFillColor(10);
	Leg.SetLineColor(10);
	Leg.SetBorderSize(0); 
	Leg.SetEntrySeparation(0);
	Leg.SetTextSize(0.035); 
	Leg.AddEntry(grName[xs_graph(2,1)],"CMS 2011","P"); 
	Leg.AddEntry(grName[Form("Atlas_xs_%dS",1)],"ATLAS", "f"); 
	Leg.DrawClone();
	
	TLatex rap_text(0.3,0.8,"|y| < 1.2"); 
	rap_text.SetNDC(kTRUE);
	rap_text.DrawClone();
	draw_header();
	
}

void draw_xs(){
	CreateCanvas("xs_can","",Cx,Cy); 
	CName["xs_can"]->cd();
	gPad->SetLogy();

	double min=0.01*0.3*grName[xs_graph(2,3)]->GetY()[grName[xs_graph(2,1)]->GetN()-1];
	double max=10000; 
	//2*grName[Form("xs_cms2010_%dS",1)]->Eval(4);

	//xs_graph(a,b) where a is the rapidity bin, 0, 1, or 2 and b is the Upsilon state 1-3
	
	grName[xs_graph(2,1)]->Draw("ap");
	grName[xs_graph(2,1)]->GetYaxis()->SetRangeUser(min,max); 
	grName[xs_graph(2,1)]->GetXaxis()->SetRangeUser(0,100); 
	grName[xs_graph(2,1)]->DrawClone("ap");
	//grName[xs_graph_stat(2,1)]->DrawClone("|| same");

	grName[Form("xs_cms2010_%dS",1)]->SetFillStyle(3002); 
	grName[Form("xs_cms2010_%dS",2)]->SetFillStyle(3002); 
	grName[Form("xs_cms2010_%dS",3)]->SetFillStyle(3002); 

	
	grName[Form("xs_cms2010_%dS",1)]->Draw("e2 same"); 
	grName[Form("gr_NLOtheory_%dS_y%d",1,2)]->Draw("l same"); 
	
	TLatex label1(0.583,0.505,"#Upsilon(1S)"); 
	label1.SetNDC(kTRUE); 
	
	TLatex label2(0.70,0.30,"#Upsilon(2S) #times 0.1"); 
	label2.SetNDC(kTRUE);

	
	TLatex label3(0.311,0.27,"#Upsilon(3S) #times 0.01"); 
	label3.SetNDC(kTRUE); 

	
	label1.SetTextSize(0.04);
	label2.SetTextSize(0.04);
	label3.SetTextSize(0.04);
	
	TLatex aLabel(0.39,0.82,"(a)"); 
	aLabel.SetNDC(kTRUE);
	aLabel.SetTextSize(0.06);

	
	TGraphAsymmErrors *gr2S=(TGraphAsymmErrors*)grName[xs_graph(2,2)]->Clone("gr2S");
	TGraphAsymmErrors *gr3S=(TGraphAsymmErrors*)grName[xs_graph(2,3)]->Clone("gr3S");
	
	TGraphAsymmErrors *gr2Sstat=(TGraphAsymmErrors*)grName[xs_graph_stat(2,2)]->Clone("gr2S");
	TGraphAsymmErrors *gr3Sstat=(TGraphAsymmErrors*)grName[xs_graph_stat(2,3)]->Clone("gr3S");
	

	TGraphAsymmErrors *grcms20102S=(TGraphAsymmErrors*)grName[Form("xs_cms2010_%dS",2)]->Clone("grcms20102S");
	TGraphAsymmErrors *grcms20103S=(TGraphAsymmErrors*)grName[Form("xs_cms2010_%dS",3)]->Clone("grcms20103S");
	
	TGraphAsymmErrors *grNLO_2S=(TGraphAsymmErrors*)grName[Form("gr_NLOtheory_%dS_y%d",2,2)]->Clone("grNLO_2S");
	TGraphAsymmErrors *grNLO_3S=(TGraphAsymmErrors*)grName[Form("gr_NLOtheory_%dS_y%d",3,2)]->Clone("grNLO_3S");
	
	
	scale_graph(gr2S,0.1);
	scale_graph(gr3S,0.01); 
	
	scale_graph(gr2Sstat,0.1);
	scale_graph(gr3Sstat,0.01); 
	
	scale_graph(grcms20102S,0.1);
	scale_graph(grcms20103S,0.01); 
	
	scale_graph(grNLO_2S,0.1);
	scale_graph(grNLO_3S,0.01); 
	
	
	gr2S->DrawClone("p same");
	//gr2Sstat->Draw("|| same");
	grcms20102S->Draw("e2 same");
	
	grNLO_2S->Draw("l same"); 
	
	gr3S->Draw("p same");
	//gr3Sstat->Draw("|| same"); 
	grcms20103S->Draw("e2 same"); 
	grNLO_3S->Draw("l same"); 
	//draw_header();
	CMS_lumi(CName["xs_can"],1,10); 
	
	label1.DrawClone();
	label2.DrawClone();
	label3.DrawClone(); 
	aLabel.DrawClone();
	
	TLegend Leg(0.45,0.6,0.87,0.84);
	Leg.SetFillColor(10);
	Leg.SetLineColor(10);
	Leg.SetBorderSize(0); 
	Leg.SetEntrySeparation(0);
	Leg.SetTextSize(0.035); 
	
	gr2S->SetMarkerSize(1.5); 
	
	Leg.AddEntry(gr2S,"CMS 2011 |y| < 1.2", "LEP"); 
	//Leg.AddEntry((TObject*)0, "", "");
	Leg.AddEntry(grName[Form("xs_cms2010_%dS",1)],"CMS 2010 |y| < 2.4 #times 0.5", "f"); 
	//Leg.AddEntry((TObject*)0, "", "");
	Leg.AddEntry(grName[Form("gr_NLOtheory_%dS_y%d",1,2)],"NLO theory","l");

	Leg.DrawClone();
}

void uncertaintyR(int iy, int num){
	//uncertainty in ratios 3/1 and 2/1
	CreateHistogram(Rstat_uncertainty(iy,num),"statistical","p_{T} [GeV]", "fractional uncertainty",fNpt, fPTbin); 
	CreateHistogram(Rbkg_sys(iy,num),"bkg order systematic","p_{T} [GeV]", "fractional uncertainty",fNpt, fPTbin); 

	CreateHistogram(Reff_uncertainty(iy,num),"#epsilon #times 5","p_{T} [GeV]", "fractional #sigma_{sys}", fNpt, fPTbin); 
	CreateHistogram(Racc_uncertainty(iy,num),"Acceptance","p_{T} [GeV]", "fractional #sigma_{sys}", fNpt, fPTbin); 
	CreateHistogram(Rtot_sys_uncertainty(iy,num),"total uncertainty","p_{T} [GeV]", "fractional #sigma_{sys}", fNpt, fPTbin); 
	CreateHistogram(Rtot_uncertainty(iy,num,"Ep"),"total uncertainty","p_{T} [GeV]", "fractional #sigma_{sys}", fNpt, fPTbin); 
	CreateHistogram(Rtot_uncertainty(iy,num,"Em"),"total uncertainty","p_{T} [GeV]", "fractional #sigma_{sys}", fNpt, fPTbin); 
	
	hName[Rstat_uncertainty(iy,num)]->SetLineColor(kRed);
	hName[Reff_uncertainty(iy,num)]->SetLineColor(kGreen);
	hName[Racc_uncertainty(iy,num)]->SetLineColor(kBlue); 
	
	ofstream output; 
	output.open(Form("/uscms/home/btcarlso/AN-UPS/tdr2/notes/AN-13-088/trunk/tables/ratio_uncertainty_table_y%d_%d-1S.tex",iy,num)); 
	int N=5; 
	output <<"\\begin{tabular}{"; 
	for(int i=0; i<=N; i++) output << "c"; 
	output << "}\\hline" << endl; 
	output << "$p_{T}$ & fit sys & $\\epsilon_{\\mu\\mu}$ & A & stat & total   \\\\ \\hline" << endl;
	
	for(int ipt=1; ipt<=hName[ratio_histogram(iy,num,"best","")]->GetNbinsX(); ipt++){
		double stat=hName[ratio_histogram(iy,num,"best","")]->GetBinError(ipt)/hName[ratio_histogram(iy,num,"best","")]->GetBinContent(ipt);
		hName[Rstat_uncertainty(iy,num)]->SetBinContent(ipt,stat); 
		hName[Rstat_uncertainty(iy,num)]->SetBinError(ipt,0); 
		
		double bkg_order=hName[background_order_systematic(iy,Form("R%d1",num))]->GetBinContent(ipt)/hName[ratio_histogram(iy,num,"best","")]->GetBinContent(ipt); 
		
		hName[Rbkg_sys(iy,num)]->SetBinContent(ipt,bkg_order); 
		hName[Rbkg_sys(iy,num)]->SetBinError(ipt,0);

		
		double Eff=hName[Form("eff_r1%d_y%d",num,iy)]->GetBinContent(ipt);
		
		double effP=TMath::Abs(hName[Form("effEp_r1%d_y%d",num,iy)]->GetBinContent(ipt)-Eff)/Eff;
		double effM=TMath::Abs(hName[Form("effEm_r1%d_y%d",num,iy)]->GetBinContent(ipt)-Eff)/Eff;
		double avgE_U=(effP+effM)/2; 
		
		hName[Reff_uncertainty(iy,num)]->SetBinContent(ipt,avgE_U*5); 
		hName[Reff_uncertainty(iy,num)]->SetBinError(ipt,0); 
		
		double Acc=hName[Form("Acc_r1%d_y%d",num,iy)]->GetBinContent(ipt); 
		
		double accP=TMath::Abs(hName[Form("AccEp_r1%d_y%d",num,iy)]->GetBinContent(ipt)-Acc)/Acc; 
		double accM=TMath::Abs(hName[Form("AccEm_r1%d_y%d",num,iy)]->GetBinContent(ipt)-Acc)/Acc;
		double avgA_U=(accP+accM)/2;
		
		hName[Racc_uncertainty(iy,num)]->SetBinContent(ipt,avgA_U); 
		hName[Racc_uncertainty(iy,num)]->SetBinError(ipt,0); 
		
		double totUP=TMath::Sqrt(effP*effP+accP*accP+stat*stat+bkg_order*bkg_order);
		double totUM=TMath::Sqrt(effM*effM+accM*accM+stat*stat+bkg_order*bkg_order);
		
		double totSys=TMath::Sqrt(avgE_U*avgE_U+avgA_U*avgA_U+bkg_order*bkg_order);		
		
		hName[Rtot_uncertainty(iy,num,"Ep")]->SetBinContent(ipt,totUP); 
		hName[Rtot_uncertainty(iy,num,"Ep")]->SetBinError(ipt,0);
		
		hName[Rtot_uncertainty(iy,num,"Em")]->SetBinContent(ipt,totUM); 
		hName[Rtot_uncertainty(iy,num,"Em")]->SetBinError(ipt,0); 

		double pt1=hName[yield_histogram(iy,ups,"best","")]->GetBinCenter(ipt)-hName[yield_histogram(iy,ups,"best","")]->GetBinWidth(ipt)/2; 
		double pt2=hName[yield_histogram(iy,ups,"best","")]->GetBinCenter(ipt)+hName[yield_histogram(iy,ups,"best","")]->GetBinWidth(ipt)/2; 
		
		output << Form("%.0f-%.0f",pt1,pt2) << " & "; 
		output << Form("%.1f &",bkg_order*100); 
		output << Form("%.1f(%.1f)",100*effP,100*effM) << " & "; 
		output << Form("%.1f(%.1f)",100*accP,100*accM) << " & "; 
		output << Form("%.1f",100*stat) << " & "; 
		output << Form("%.1f(%.1f)",100*totUP,100*totUM) << "\\\\ " <<  endl; 
	}		
	
	output << "\\hline" << endl; 
	output << "\\end{tabular}" << endl; 
	output.close();
	
}

void uncertainty(int iy, int ups){
    cout << "uncertainty " << iy << " " << ups << endl;
	CreateHistogram(stat_uncertainty(iy,ups),"statistical","p_{T} [GeV]", "fractional uncertainty",fNpt, fPTbin); 
	CreateHistogram(fit_sys(iy,ups),"fit systematic","p_{T} [GeV]", "fractional uncertainty",fNpt, fPTbin); 
	CreateHistogram(eff_uncertainty(iy,ups),"#epsilon","p_{T} [GeV]", "fractional #sigma_{sys}", fNpt, fPTbin); 
	CreateHistogram(acc_uncertainty(iy,ups),"Acceptance","p_{T} [GeV]", "fractional #sigma_{sys}", fNpt, fPTbin); 
	CreateHistogram(tot_sys_uncertainty(iy,ups),"total uncertainty","p_{T} [GeV]", "fractional #sigma_{sys}", fNpt, fPTbin); 
	CreateHistogram(tot_uncertainty(iy,ups,"Ep"),"total uncertainty","p_{T} [GeV]", "fractional #sigma_{sys}", fNpt, fPTbin); 
	CreateHistogram(tot_uncertainty(iy,ups,"Em"),"total uncertainty","p_{T} [GeV]", "fractional #sigma_{sys}", fNpt, fPTbin); 

	hName[stat_uncertainty(iy,ups)]->SetLineColor(kRed);
	hName[eff_uncertainty(iy,ups)]->SetLineColor(kGreen);
	hName[acc_uncertainty(iy,ups)]->SetLineColor(kBlue); 
	
	ofstream output; 
	output.open(Form("/uscms/home/btcarlso/AN-UPS/tdr2/notes/AN-13-088/trunk/tables/xs_uncertainty_table_y%d_ups%d.tex",iy,ups)); 
	int N=5; 
	output <<"\\begin{tabular}{"; 
	for(int i=0; i<=N; i++) output << "c"; 
	output << "}\\hline" << endl; 
	output << "$p_{T}$ & fit sys & $\\epsilon_{\\mu\\mu}$ & A & stat & total   \\\\" << endl; 
	
    double Effmass_dep=0.01;
	
	for(int ipt=1; ipt<=hName[yield_histogram(iy,ups,"best","")]->GetNbinsX(); ipt++){
		double stat=hName[yield_histogram(iy,ups,"best","")]->GetBinError(ipt)/hName[yield_histogram(iy,ups,"best","")]->GetBinContent(ipt);
		
		double bkg_order=hName[background_order_systematic(iy,Form("%dS",ups))]->GetBinContent(ipt)/hName[yield_histogram(iy,ups,"best","")]->GetBinContent(ipt);

        bkg_order=TMath::Sqrt(bkg_order*bkg_order+0.01*0.01); // 1% lineshape uncertainty
        
		double Eff=prName[eff_name(iy, ups, "")]->GetBinContent(ipt);
		
		double effP=TMath::Abs(prName[eff_name(iy, ups, "Ep")]->GetBinContent(ipt)-Eff)/Eff;
		double effM=TMath::Abs(prName[eff_name(iy, ups, "Em")]->GetBinContent(ipt)-Eff)/Eff;
		double avgE_U=(effP+effM)/2; 
		
		double Acc=hName[acc_name(iy,ups,"")]->GetBinContent(ipt); 
		
		double accP=TMath::Abs(hName[acc_name(iy,ups,"Ep")]->GetBinContent(ipt)-Acc)/Acc; 
		double accM=TMath::Abs(hName[acc_name(iy,ups,"Em")]->GetBinContent(ipt)-Acc)/Acc;
		double avgA_U=(accP+accM)/2;
		
		
		double totUP=TMath::Sqrt(effP*effP+accP*accP+stat*stat+bkg_order*bkg_order+Effmass_dep*Effmass_dep);
		double totUM=TMath::Sqrt(effM*effM+accM*accM+stat*stat+bkg_order*bkg_order+Effmass_dep*Effmass_dep);

		double totSys=TMath::Sqrt(avgE_U*avgE_U+avgA_U*avgA_U+bkg_order*bkg_order+Effmass_dep*Effmass_dep);
		
		
		hName[eff_uncertainty(iy,ups)]->SetBinContent(ipt,avgE_U); 
		hName[eff_uncertainty(iy,ups)]->SetBinError(ipt,0);
		
		hName[acc_uncertainty(iy,ups)]->SetBinContent(ipt,avgA_U); 
		hName[acc_uncertainty(iy,ups)]->SetBinError(ipt,0);
		
		hName[stat_uncertainty(iy,ups)]->SetBinContent(ipt,stat);
		hName[stat_uncertainty(iy,ups)]->SetBinError(ipt,0); 
		
		hName[fit_sys(iy,ups)]->SetBinContent(ipt,bkg_order);
		hName[fit_sys(iy,ups)]->SetBinError(ipt,0); 
		
		hName[tot_sys_uncertainty(iy,ups)]->SetBinContent(ipt,totSys);
		hName[tot_sys_uncertainty(iy,ups)]->SetBinError(ipt,0);

		hName[tot_uncertainty(iy,ups,"Ep")]->SetBinContent(ipt,totUP);
		hName[tot_uncertainty(iy,ups,"Ep")]->SetBinError(ipt,0);
		
		hName[tot_uncertainty(iy,ups,"Em")]->SetBinContent(ipt,totUM);
		hName[tot_uncertainty(iy,ups,"Em")]->SetBinError(ipt,0);
		
		double pt1=hName[yield_histogram(iy,ups,"best","")]->GetBinCenter(ipt)-hName[yield_histogram(iy,ups,"best","")]->GetBinWidth(ipt)/2; 
		double pt2=hName[yield_histogram(iy,ups,"best","")]->GetBinCenter(ipt)+hName[yield_histogram(iy,ups,"best","")]->GetBinWidth(ipt)/2; 
		
		output << Form("%.0f-%.0f",pt1,pt2) << " & "; 
		output << Form("%.2f &",bkg_order*100);
		output << Form("%.1f(%.1f)",100*effP,100*effM) << " & "; 
		output << Form("%.1f(%.1f)",100*accP,100*accM) << " & "; 
		output << Form("%.1f",100*stat) << " & "; 
		output << Form("%.1f(%.1f)",100*totUP,100*totUM) << "\\\\ " <<  endl; 
		
		
	}
	output << "\\hline" << endl; 
	output << "\\end{tabular}" << endl; 
	output.close();
}

void acceptance_ratio(int iy){
	
	int num=1; 
	
	for(int den=2; den<=3; den++){
		CreateHistogram(Form("Acc_r1%d_y%d",den,iy),"","p_{T} [GeV]",Form("A(1S)/A(%dS)",den),fNpt,fPTbin); 
		CreateHistogram(Form("AccEp_r1%d_y%d",den,iy),"","p_{T} [GeV]",Form("A(1S)/A(%dS)",den),fNpt,fPTbin); 
		CreateHistogram(Form("AccEm_r1%d_y%d",den,iy),"","p_{T} [GeV]",Form("A(1S)/A(%dS)",den),fNpt,fPTbin); 

		hName[Form("Acc_r1%d_y%d",den,iy)]->SetLineColor(kBlack);
		hName[Form("AccEp_r1%d_y%d",den,iy)]->SetLineColor(kRed);
		hName[Form("AccEm_r1%d_y%d",den,iy)]->SetLineColor(kGreen); 
		
		for(int ipt=1; ipt<=hName[Form("Acc_r1%d_y%d",den,iy)]->GetNbinsX(); ipt++){
			float acc_num=hName[acc_name(iy,num,"")]->GetBinContent(ipt); 
			float acc_den=hName[acc_name(iy,den,"")]->GetBinContent(ipt); 
			
            float acc_numP=hName[acc_name(iy,num,"Ep")]->GetBinContent(ipt);
			float acc_denP=hName[acc_name(iy,den,"Ep")]->GetBinContent(ipt);
            
            float acc_numM=hName[acc_name(iy,num,"Ep")]->GetBinContent(ipt);
			float acc_denM=hName[acc_name(iy,den,"Ep")]->GetBinContent(ipt);
			
            float sigmaNP=TMath::Abs(acc_numP-acc_num);
            float sigmaNM=TMath::Abs(acc_numM-acc_num);
            float sigmaN=TMath::Max(sigmaNP,sigmaNM);
            
            float sigmaDP=TMath::Abs(acc_denP-acc_den);
            float sigmaDM=TMath::Abs(acc_denM-acc_den);
            float sigmaD=TMath::Max(sigmaDP,sigmaDM); // maximum of 2 choices
            
            float sigmaR=TMath::Power(sigmaN/acc_num,2)+TMath::Power(sigmaD/acc_den,2);
            sigmaR=TMath::Sqrt(sigmaR);
            
			hName[Form("Acc_r1%d_y%d",den,iy)]->SetBinContent(ipt,acc_num/acc_den);
			hName[Form("Acc_r1%d_y%d",den,iy)]->SetBinError(ipt,0);
            
            hName[Form("AccEp_r1%d_y%d",den,iy)]->SetBinContent(ipt,(acc_num/acc_den)*(1+sigmaR));
			hName[Form("AccEp_r1%d_y%d",den,iy)]->SetBinError(ipt,0);
            
            hName[Form("AccEm_r1%d_y%d",den,iy)]->SetBinContent(ipt,(acc_num/acc_den)*(1-sigmaR));
			hName[Form("AccEm_r1%d_y%d",den,iy)]->SetBinError(ipt,0);
			/*
             old way, vary acceptance coherently. too conservative!
			acc_num=hName[acc_name(iy,num,"Ep")]->GetBinContent(ipt); 
			acc_den=hName[acc_name(iy,den,"Ep")]->GetBinContent(ipt); 
			
			
			hName[Form("AccEp_r1%d_y%d",den,iy)]->SetBinContent(ipt,acc_num/acc_den); 
			hName[Form("AccEp_r1%d_y%d",den,iy)]->SetBinError(ipt,0); 
			
			acc_num=hName[acc_name(iy,num,"Em")]->GetBinContent(ipt); 
			acc_den=hName[acc_name(iy,den,"Em")]->GetBinContent(ipt); 
			
			hName[Form("AccEm_r1%d_y%d",den,iy)]->SetBinContent(ipt,acc_num/acc_den); 
			hName[Form("AccEm_r1%d_y%d",den,iy)]->SetBinError(ipt,0); 
             */
		}
		
		CreateCanvas(Form("Acceptance_Ratio_1-%dS_y%d",den,iy), "",Cx,Cy); 
		CName[Form("Acceptance_Ratio_1-%dS_y%d",den,iy)]->cd();
		
		hName[Form("Acc_r1%d_y%d",den,iy)]->SetLineColor(kBlack);
		hName[Form("AccEp_r1%d_y%d",den,iy)]->SetLineColor(kRed);
		hName[Form("AccEm_r1%d_y%d",den,iy)]->SetLineColor(kGreen); 
		
		hName[Form("Acc_r1%d_y%d",den,iy)]->SetMinimum(0.8);
		hName[Form("Acc_r1%d_y%d",den,iy)]->SetMaximum(1.2);
		
		hName[Form("Acc_r1%d_y%d",den,iy)]->Draw();
		hName[Form("AccEp_r1%d_y%d",den,iy)]->Draw("same");
		hName[Form("AccEm_r1%d_y%d",den,iy)]->Draw("same");

		draw_headersim(); 

		TLegend leg(0.3,0.25,0.6,0.4); 
		leg.SetFillColor(10);
		leg.SetLineColor(10);
		
		leg.AddEntry(hName[Form("Acc_r1%d_y%d",den,iy)],Form("Acceptance %d/1",den),"L"); 
		leg.AddEntry(hName[Form("AccEp_r1%d_y%d",den,iy)],Form("Acceptance %d/1 +#sigma HX",den),"L"); 
		leg.AddEntry(hName[Form("AccEm_r1%d_y%d",den,iy)],Form("Acceptance %d/1 -#sigma HX",den),"L"); 
		leg.DrawClone("same"); 
		
		TString bt;
		if(iy==0)bt=	"|y| < 0.6"; 
		if(iy==1) bt="0.6 < |y| < 1.2"; 
		if(iy==2) bt="|y| < 1.2"; 
		
		TLatex text(0.3,0.85,bt);
		text.SetNDC(kTRUE); 
		text.DrawClone();
	}
	
	if(iy!=0) return; 
	CreateHistogram("Acc_Ry1_y0","","p_{T} [GeV]","0.6 < |y| < 1.2/|y| < 0.6",fNpt,fPTbin);
	CreateHistogram("AccEp_Ry1_y0","","p_{T} [GeV]","0.6 < |y| < 1.2/|y| < 0.6",fNpt,fPTbin);
	CreateHistogram("AccEm_Ry1_y0","","p_{T} [GeV]","0.6 < |y| < 1.2/|y| < 0.6",fNpt,fPTbin);
	hName["Acc_Ry1_y0"]->SetLineColor(kBlack);
	hName["AccEp_Ry1_y0"]->SetLineColor(kRed);
	hName["AccEm_Ry1_y0"]->SetLineColor(kGreen);
	
	for(int ipt=1; ipt<=hName["Acc_Ry1_y0"]->GetNbinsX(); ipt++){
		float acc_num=hName[acc_name(1,1,"")]->GetBinContent(ipt); 
		float acc_den=hName[acc_name(0,1,"")]->GetBinContent(ipt); 
		
		hName["Acc_Ry1_y0"]->SetBinContent(ipt,acc_num/acc_den); 
		hName["Acc_Ry1_y0"]->SetBinError(ipt,0); 
		
		acc_num=hName[acc_name(1,1,"Ep")]->GetBinContent(ipt); 
		acc_den=hName[acc_name(0,1,"Ep")]->GetBinContent(ipt); 
		
		hName["AccEp_Ry1_y0"]->SetBinContent(ipt,acc_num/acc_den); 
		hName["AccEp_Ry1_y0"]->SetBinError(ipt,0); 
		
		acc_num=hName[acc_name(1,1,"Em")]->GetBinContent(ipt); 
		acc_den=hName[acc_name(0,1,"Em")]->GetBinContent(ipt); 
		
		hName["AccEm_Ry1_y0"]->SetBinContent(ipt,acc_num/acc_den); 
		hName["AccEm_Ry1_y0"]->SetBinError(ipt,0); 
	}
	
	CreateCanvas("Acceptance_ratio_y1-0","",Cx,Cy); 
	CName["Acceptance_ratio_y1-0"]->cd(); 
	
	hName["Acc_Ry1_y0"]->SetMinimum(0.8); 
	hName["Acc_Ry1_y0"]->SetMaximum(1.1); 
	
	hName["Acc_Ry1_y0"]->Draw();
	hName["AccEp_Ry1_y0"]->Draw("same");
	hName["AccEm_Ry1_y0"]->Draw("same");
	
	
	TLegend leg2(0.3,0.25,0.6,0.4); 
	leg2.SetFillColor(10);
	leg2.SetLineColor(10);
	
	leg2.AddEntry(hName["Acc_Ry1_y0"],"0.6 < |y| < 1.2/|y| < 0.6","L"); 
	leg2.AddEntry(hName["AccEp_Ry1_y0"],"0.6 < |y| < 1.2/|y| < 0.6 +#sigma HX","L"); 
	leg2.AddEntry(hName["AccEm_Ry1_y0"],"0.6 < |y| < 1.2/|y| < 0.6 -#sigma HX","L"); 
	leg2.DrawClone("same"); 
	
	TLatex text2(0.3,0.85,"#Upsilon(1S)");
	text2.SetNDC(kTRUE);  
	text2.DrawClone();
	

	
}

void draw_acc_ratio(int iy){
	for(int den=2; den<=3; den++){
		CreateCanvas(Form("AccRplot_1%d_y%d",den,iy),"",Cx,Cy); 
		CName[Form("AccRplot_1%d_y%d",den,iy)]->cd();
		
		TLegend L(0.2,0.2,0.5,0.5);
		L.SetFillColor(10);
		L.SetLineColor(10);
		
		hName[Form("Acc_r1%d_y%d",den,iy)]->Draw(); 
		hName[Form("AccEp_r1%d_y%d",den,iy)]->Draw("same");
		hName[Form("AccEm_r1%d_y%d",den,iy)]->Draw("same"); 
		
		L.AddEntry(hName[Form("Acc_r1%d_y%d",den,iy)],"Acc","L"); 
		L.AddEntry(hName[Form("AccEp_r1%d_y%d",den,iy)],"Acc +#sigma","L"); 
		L.AddEntry(hName[Form("AccEm_r1%d_y%d",den,iy)],"Acc -#sigma","L"); 

		L.DrawClone(); 
		draw_header();
	}
}

void combine_rapidity_bins(int ups){

	//Add yields, and add uncertainties on yield in quadrature
	//Do weighted average of ratios
	
	CreateHistogram(yield_histogram(2,ups,"best",""),"","p_{T} [GeV]","Events", fNpt, fPTbin);
	CreateHistogram(background_order_systematic(2,Form("%dS",ups)),"","p_{T} [GeV]","Events", fNpt, fPTbin);
	
	if(ups>=2) {
		CreateHistogram(ratio_histogram(2,ups,"best",""),"","p_{T} [GeV]", "Ratio",fNpt, fPTbin);
		CreateHistogram(background_order_systematic(2,Form("R%d1",ups)),"","p_{T} [GeV]","Events", fNpt, fPTbin);
	}
	//CreateHistogram(ratio_histogram(2,2,"best",""),"","p_{T} [GeV]", "Ratio",fNpt, fPTbin);

	for(int ipt=1; ipt<=hName[yield_histogram(2,ups,"best","")]->GetNbinsX(); ipt++){
		double yield=hName[yield_histogram(0,ups,"best","")]->GetBinContent(ipt)+hName[yield_histogram(1,ups,"best","")]->GetBinContent(ipt);
		double yieldE0=hName[yield_histogram(0,ups,"best","")]->GetBinError(ipt);
		double yieldE1=hName[yield_histogram(1,ups,"best","")]->GetBinError(ipt);
		
		double sys0=hName[background_order_systematic(0,Form("%dS",ups))]->GetBinContent(ipt); 
		double sys1=hName[background_order_systematic(1,Form("%dS",ups))]->GetBinContent(ipt); 

		
		double yieldE=TMath::Sqrt(yieldE0*yieldE0+yieldE1*yieldE1); 
		double sysE=TMath::Sqrt(sys0*sys0+sys1*sys1); 
		

		hName[yield_histogram(2,ups,"best","")]->SetBinContent(ipt,yield);
		hName[yield_histogram(2,ups,"best","")]->SetBinError(ipt,yieldE);
		
		hName[background_order_systematic(2,Form("%dS",ups))]->SetBinContent(ipt,sysE); 
		hName[background_order_systematic(2,Form("%dS",ups))]->SetBinError(ipt,0); 
		
		if(ups>=2){
			int num=ups; 
			double sigma0=hName[ratio_histogram(0,num,"best","")]->GetBinError(ipt);
			double sigma1=hName[ratio_histogram(1,num,"best","")]->GetBinError(ipt);

			double sys0=hName[background_order_systematic(0,Form("R%d1",ups))]->GetBinContent(ipt);
			double sys1=hName[background_order_systematic(1,Form("R%d1",ups))]->GetBinContent(ipt);
			
			double w0=1./TMath::Power(sigma0,2);
			double w1=1./TMath::Power(sigma1,2);
			
			double R0=hName[ratio_histogram(0,num,"best","")]->GetBinContent(ipt); 
			double R1=hName[ratio_histogram(1,num,"best","")]->GetBinContent(ipt); 

			double Rt=(R0*w0+R1*w1)/(w0+w1);
			//double sigmat=TMath::Sqrt(sigma0*sigma0+sigma1*sigma1);
			//double syst=TMath::Sqrt(sys0*sys0+sys1*sys1);
			
            
            double sigmat=TMath::Sqrt(1/(w0+w1));
            double syst=TMath::Sqrt(1/(TMath::Power(sys0,-2)+TMath::Power(sys1,-2)));

            //cout << "sumquad: " <<TMath::Sqrt(sys0*sys0+sys1*sys1) << " sigmat: " <<syst << endl;

            
			hName[ratio_histogram(2,num,"best","")]->SetBinContent(ipt,Rt);
			hName[ratio_histogram(2,num,"best","")]->SetBinError(ipt,sigmat);
			hName[background_order_systematic(2,Form("R%d1",ups))]->SetBinContent(ipt,syst); 
			hName[background_order_systematic(2,Form("R%d1",ups))]->SetBinError(ipt,0); 
		}
		
	}
	
}

void RN1(int iy, int num){
	//Ratio of Upsilon(nS)/Upsilon(1S)
	
	TH1F *R = (TH1F*)hName[ratio_histogram(iy,num,"best","")]->Clone(Form("R_%d1_y%d",num,iy));
	TGraphAsymmErrors *gr_stat = new TGraphAsymmErrors(R); 
	gr_stat->SetName(xs_ratio_stat(iy,num)); 
	
	TGraphAsymmErrors *gr = new TGraphAsymmErrors(R); 
	gr->SetName(xs_ratio(iy,num));
	
	for(int ipt=1; ipt<=R->GetNbinsX(); ipt++){
		double pt_=R->GetBinCenter(ipt); 
		pt_=grName[xs_graph(iy,num)]->GetX()[ipt-1];
		
		double eff_ratio=hName[Form("eff_r1%d_y%d",num,iy)]->GetBinContent(ipt); 
		//double effsg_ratio=hName[Form("eff_sg_r1%d_y%d",num,iy)]->GetBinContent(ipt); 
		double effsg_ratio=1; 
		double acc_ratio=hName[Form("Acc_r1%d_y%d",num,iy)]->GetBinContent(ipt); 
		double ratio_=R->GetBinContent(ipt)*eff_ratio*effsg_ratio*acc_ratio;
		
		double ratio_stat =ratio_*hName[Rstat_uncertainty(iy,num)]->GetBinContent(ipt); 
		double ratio_UP =ratio_*hName[Rtot_uncertainty(iy,num,"Ep")]->GetBinContent(ipt); 
		double ratio_UM =ratio_*hName[Rtot_uncertainty(iy,num,"Em")]->GetBinContent(ipt); 

		double pt_L=pt_-R->GetBinLowEdge(ipt); 
		double pt_H=R->GetBinLowEdge(ipt)+R->GetBinWidth(ipt)-pt_; 
		
		gr_stat->SetPoint(ipt-1,pt_, ratio_); 
		gr_stat->SetPointEXhigh(ipt-1,0);
		gr_stat->SetPointEXlow(ipt-1,0);
		
		gr_stat->SetPointEYhigh(ipt-1,ratio_stat);
		gr_stat->SetPointEYlow(ipt-1,ratio_stat);
		
		gr->SetPoint(ipt-1,pt_, ratio_); 
		gr->SetPointEXhigh(ipt-1,pt_H);
		gr->SetPointEXlow(ipt-1,pt_L);
		
		gr->SetPointEYhigh(ipt-1,ratio_UP);
		gr->SetPointEYlow(ipt-1,ratio_UM);
		
	}
	
	gr_stat->GetXaxis()->SetTitle(x_label);
	gr->GetXaxis()->SetTitle(x_label);

	TString laby;
	if (num==3) laby=Ratio_y3;
	if(num==2) laby=Ratio_y2; 
	
	gr_stat->GetYaxis()->SetTitle(laby); 
	gr->GetYaxis()->SetTitle(laby); 

	set_CMS2011(gr,"ratio",num); 
	set_CMS2011(gr_stat,"stat",num); 

	grName[gr_stat->GetName()]=gr_stat;
	grName[gr->GetName()]=gr; 

	
}
double power_law_variablepo(double *x, double*par){
	return par[0]/(par[2]+TMath::Power(x[0]/par[3],par[1])); 
}
double power_law(double *x, double*par){
	return par[0]/(par[2]+TMath::Power(x[0]/20,par[1])); 
}

double ratio_function(double *x, double *par){
	
	//Par[0] = Anum, Par[1]=Cnum, Par[2]=Alphanum
	//Par[3] = Aden, Par[4]=Cden, Par[5]=Alphaden
	
	double SN=par[0]/(par[2]+TMath::Power((x[0]/20),par[1])); 
	double SD=par[3]/(par[5]+TMath::Power((x[0]/20),par[4])); 
	
	double RR=SN/SD;
	
	return RR;
}

void xs_fit_po_study(int iy, int ups){
    cout << "xs fit po study: " << endl;
	CreateHistogram(Form("xsfit_postudy_y%d_%dS",iy,ups),"Fit(Po)","p_{o}","#alpha",10,20,30); 
	
	TH1F *h= new TH1F("h","",fNpt,fPTbin); 
	
	for (int ipt=1; ipt<=h->GetNbinsX(); ipt++) {
		
		double xs=grName[xs_graph(iy,ups)]->GetY()[ipt-1];
		
		double P=grName[xs_graph(iy,ups)]->GetErrorYhigh(ipt-1);; 
		double M=grName[xs_graph(iy,ups)]->GetErrorYlow(ipt-1);
		
		double err=(P+M)/2;
		
		h->SetBinContent(ipt,xs); 
		h->SetBinError(ipt,err); 
		
	}
	
	TF1 *xs_fit = new TF1(Form("xs_fit_power_y%d_%dS",iy,ups), power_law_variablepo,20,100,4); 
	xs_fit->SetLineColor(kBlack); 
	xs_fit->SetParNames("A","alpha","C","po");
	xs_fit->SetParameters(1,5,1,20);
	
	
	for(int po=20; po<=30; po++){
		double PO=static_cast<double>(po);
		xs_fit->FixParameter(3,PO);
		h->Fit(Form("xs_fit_power_y%d_%dS",iy,ups),"QR0");
		h->Fit(Form("xs_fit_power_y%d_%dS",iy,ups),"RI0");
		cout << "PO: " << PO << " " << xs_fit->GetParameter(1) << endl; 
		
		hName[Form("xsfit_postudy_y%d_%dS",iy,ups)]->SetBinContent(po-20,xs_fit->GetParameter(1)); 
		hName[Form("xsfit_postudy_y%d_%dS",iy,ups)]->SetBinError(po-20,xs_fit->GetParError(1)); 

	}
	delete xs_fit;
	delete h;

	CreateCanvas(Form("po_study_y%d_%dS",iy,ups),"",Cx,Cy);
	CName[Form("po_study_y%d_%dS",iy,ups)]->cd();
	hName[Form("xsfit_postudy_y%d_%dS",iy,ups)]->Draw("E1");

}

void xs_fit(int iy, int ups){
    bool print=true;
    if(print)cout << "xs fit " << endl;
	//create a histogram for the rapidity combined bin
	//fill the histogram with the average total uncertainty 
	//fit the histogram using the bin integral method
	
	TH1F *h= new TH1F("h","",fNpt,fPTbin); 
	
	for (int ipt=1; ipt<=h->GetNbinsX(); ipt++) {
		
		double xs=grName[xs_graph(iy,ups)]->GetY()[ipt-1];
		
		double P=grName[xs_graph(iy,ups)]->GetErrorYhigh(ipt-1);; 
		double M=grName[xs_graph(iy,ups)]->GetErrorYlow(ipt-1);
		
		double err=(P+M)/2;
		
		h->SetBinContent(ipt,xs); 
		h->SetBinError(ipt,err); 
		
	}
	if(print)cout << "define function: " << endl;
	TF1 *xs_fit = new TF1(Form("xs_fit_power_y%d_%dS",iy,ups), power_law,20,100,3); 
	xs_fit->SetLineColor(kBlack);
	f1Name[xs_fit->GetName()]=xs_fit; 
	f1Name[xs_fit->GetName()]->SetParNames("A","alpha","C");
	f1Name[xs_fit->GetName()]->SetParameters(1,5,1);
	
	TF1 *expo = new TF1(Form("xs_fit_expo_y%d_%dS",iy,ups), "expo",10,20);
	expo->SetLineColor(kRed); 
	f1Name[expo->GetName()]=expo; 
	f1Name[expo->GetName()]->SetParameters(15,5); 
	if(print)cout << "do fit: " << endl;
	h->Fit(f1Name[xs_fit->GetName()],"QR0"); //First do fit without option I to get the approximately correct parameters 
	TFitResultPtr fitRes=h->Fit(f1Name[xs_fit->GetName()],"QRIS0"); //Then do fit with option I
	
    TMatrixDSym cor=fitRes->GetCorrelationMatrix();
    CreateCanvas(Form("correlationmatrix_power_law_y%d_%dS",iy,ups),"",Cx,Cy);
    CName[Form("correlationmatrix_power_law_y%d_%dS",iy,ups)]->cd();
    cor.DrawClone("COLZ text");
    /*
    TH2F *hh=(TH2F*)cor.DrawClone("COLZ text");
    hh->GetXaxis()->SetBinLabel(1,"A");
    hh->GetXaxis()->SetBinLabel(2,"#alpha");
    hh->GetXaxis()->SetBinLabel(3,"C");

    hh->GetYaxis()->SetBinLabel(1,"A");
    hh->GetYaxis()->SetBinLabel(2,"#alpha");
    hh->GetYaxis()->SetBinLabel(3,"C");
    
    hh->DrawClone("COLZ text");
    */
	h->Fit(f1Name[Form("xs_fit_expo_y%d_%dS",iy,ups)],"QRI0"); 
	
	delete h; 
	
	if(print)cout << "xs fit done: " << endl;
}

void dsigma_pt(int iy, int ups){
	cout << "compute xs y:" << iy << " ups: " << ups << endl; 
	TH1F *xs = (TH1F*)hName[yield_histogram(iy,ups,"best","")]->Clone(Form("yieldtmp_y%d_%dS",iy,ups));
	TGraphAsymmErrors *gr = new TGraphAsymmErrors(xs); 
	gr->SetName(xs_graph(iy,ups)); 
	
	TGraphAsymmErrors *gr_stat = new TGraphAsymmErrors(xs); 
	gr_stat->SetName(xs_graph_stat(iy,ups)); 
	
	for(int ipt=1; ipt<=xs->GetNbinsX(); ipt++){
		double pt_=xs->GetBinCenter(ipt); 

		double y=xs->GetBinContent(ipt); 
		double A=hName[acc_name(iy,ups,"")]->GetBinContent(ipt);
		double eff=prName[eff_name(iy, ups, "")]->GetBinContent(ipt);
		TString SGName=Form("sg_y%d_TH1_%dS_rebin",iy,ups); 
		if(iy>1) SGName=Form("sg_TH1_%dS_rebin",ups);
		//double eff_sg=hName[SGName]->GetBinContent(ipt); 
		double eff_sg=0.5;
		double delta_pt=hName["bin_width_pt"]->GetBinContent(ipt); 
        double eff_track=0.99;
        
		double xs_=y/(A*eff*eff_sg*4900*delta_pt*eff_track*eff_track);
		
		double xs_stat=hName[stat_uncertainty(iy,ups)]->GetBinContent(ipt)*xs_; 
		double xs_totP=hName[tot_uncertainty(iy,ups,"Ep")]->GetBinContent(ipt)*xs_; 
		double xs_totM=hName[tot_uncertainty(iy,ups,"Em")]->GetBinContent(ipt)*xs_; 

		double pt_L=pt_-xs->GetBinLowEdge(ipt); 
		double pt_H=xs->GetBinLowEdge(ipt)+xs->GetBinWidth(ipt)-pt_; 
		
		gr->SetPoint(ipt-1,pt_, xs_); 
		gr->SetPointEXhigh(ipt-1,pt_H);
		gr->SetPointEXlow(ipt-1,pt_L);
		
		gr->SetPointEYhigh(ipt-1,xs_totP);
		gr->SetPointEYlow(ipt-1,xs_totM);
		
		gr_stat->SetPoint(ipt-1,pt_, xs_); 
		gr_stat->SetPointEXhigh(ipt-1,pt_H);
		gr_stat->SetPointEXlow(ipt-1,pt_L);
		
		gr_stat->SetPointEYhigh(ipt-1,xs_stat);
		gr_stat->SetPointEYlow(ipt-1,xs_stat);
		
		//cout << "pt: " << pt_ << " "  << xs_ << endl; 
		

	}
	set_CMS2011(gr,"",ups); 
	
	gr->GetXaxis()->SetTitle(x_label); 
	gr->GetYaxis()->SetTitle(xs_y); 
	grName[gr->GetName()]=gr; 
	
	set_CMS2011(gr_stat,"stat",ups); 

	gr_stat->GetXaxis()->SetTitle(x_label); 
	gr_stat->GetYaxis()->SetTitle(xs_y); 
	grName[gr_stat->GetName()]=gr_stat; 
	
	cout << "xs computed: " << endl;
}

void draw_dR(){
	CreateCanvas("dR_ptE_pt","",Cx,Cy); 
	CName["dR_ptE_pt"]->cd();
	hName2D["h_DeltaRPtE_pt"]->SetTitle("");
	hName2D["h_DeltaRPtE_pt"]->Draw("COLZ"); 
	draw_header(); 
	
}

void draw_pt(){
	//plot pT distributions 
	
	for(int ups=1; ups<=3; ups++){
		CreateCanvas(Form("pt_distribution_%dS",ups),"",Cx,Cy); 
		CName[Form("pt_distribution_%dS",ups)]->cd();
		gPad->SetLogy();
		
		TLegend Leg(0.35,0.65,0.8,0.85);
		Leg.SetFillColor(10);
		Leg.SetLineColor(10);
		Leg.SetBorderSize(0); 
		Leg.SetTextSize(0.035); 
		double max=0; 
		for(int iy=0; iy<fNy; iy++){
			cout << "iy: " << iy << " ups: " << ups << endl; 
			TString BT;
			if(iy==0 && ups==1) BT=Form("|y| < 0.6, %.2f < M_{#mu#mu} < %.2f GeV",y1m[0],y1m[1]); 
			if(iy==0 && ups==2) BT=Form("|y| < 0.6, %.2f < M_{#mu#mu} < %.2f GeV",y2m[0],y2m[1]); 
			if(iy==0 && ups==3) BT=Form("|y| < 0.6, %.2f < M_{#mu#mu} < %.2f GeV",y3m[0],y3m[1]); 
			
			if(iy==1 && ups==1) BT=Form("0.6 < |y| < 1.2, %.2f < M_{#mu#mu} < %.2f GeV",y1m[0],y1m[1]); 
			if(iy==1 && ups==2) BT=Form("0.6 < |y| < 1.2, %.2f < M_{#mu#mu} < %.2f GeV",y2m[0],y2m[1]); 
			if(iy==1 && ups==3) BT=Form("0.6 < |y| < 1.2, %.2f < M_{#mu#mu} < %.2f GeV",y3m[0],y3m[1]); 
			
	
			if(iy==1){
				hName[Form("h_pt_y%d_%dS",iy,ups)]->SetMarkerStyle(23); 
				hName[Form("h_pt_y%d_%dS",iy,ups)]->SetMarkerColor(kRed); 
			}
			
			Leg.AddEntry(hName[Form("h_pt_y%d_%dS",iy,ups)],BT,"LEP"); 

			if(iy==0)hName[Form("h_pt_y%d_%dS",iy,ups)]->Draw(); 
			else hName[Form("h_pt_y%d_%dS",iy,ups)]->Draw("same");

		}
		cout << "Draw legend: " << endl; 
		Leg.DrawClone(); 
		//draw_header(); 
		CMS_lumi(CName[Form("pt_distribution_%dS",ups)],1,10); 
		
		for(int i=fNpt-1; i>14; i--){
			max=5*hName[Form("h_pt_y%d_%dS",0,ups)]->GetBinContent(fPTbin[i]);
			TLine Lpt(fPTbin[i],0,fPTbin[i],max); 
			Lpt.SetLineColor(kBlue);
			Lpt.DrawClone(); 
		}
		TLatex KC; 
		KC.SetTextSize(0.03); 
		KC.SetNDC(kTRUE);
		float x1=0.55;
		float y1=0.6;
		KC.DrawLatex(x1,y1, "1.4 <|#eta(#mu)|< 1.6: p_{T}(#mu)> 3 GeV"); 
		y1-=0.035; 
		KC.DrawLatex(x1,y1, "1.2 <|#eta(#mu)|< 1.4: p_{T}(#mu)> 3.5 GeV"); 
		y1-=0.035; 
		KC.DrawLatex(x1,y1, "|#eta(#mu)|< 1.2: p_{T}(#mu)> 4.5 GeV"); 
		x1=0.2; y1=0.45;
		KC.SetTextSize(0.05); 
		KC.DrawLatex(x1,y1,"2 GeV bins"); 
		
	}
}

void draw_m(){
	/*
	for(int iy=0; iy<fNy;iy++){
		for(int ipt=0; ipt<fNpt; ipt++){
			TString can_name=h_m(iy,ipt,"") + "_can"; 
			CreateCanvas(can_name,"",Cx,Cy); 
			CName[can_name]->cd(); 
			hName[h_m(iy,ipt,"")]->Draw("E1"); 
			draw_header();
		}
	}
	*/
    bool print_muon_cuts=false;
    
	CreateCanvas("m_summary","",Cx,Cy); 
	CName["m_summary"]->cd();
	
	TLegend Leg(0.45,0.7,0.85,0.87); 
	Leg.SetFillColor(10);
	Leg.SetLineColor(10); 
	Leg.SetBorderSize(0); 
	
	hName["h_m_y0"]->SetMarkerSize(0.75); 
	hName["h_m_y1"]->SetMarkerColor(kRed); 
	hName["h_m_y1"]->SetMarkerStyle(23);
	hName["h_m_y1"]->SetMarkerSize(0.75); 
	
	hName["h_m_y0"]->GetYaxis()->SetTitleOffset(1.4);
	
	hName["h_m_y0"]->GetYaxis()->SetTitle("Events / 10 MeV"); 
	
	hName["h_m_y0"]->DrawClone("p");
	hName["h_m_y1"]->DrawClone("p same"); 
	
	hName["h_m_y0"]->SetMarkerSize(1.5);
	hName["h_m_y1"]->SetMarkerSize(1.5);
	
	hName["h_m_y0"]->SetMaximum(hName["h_m_y0"]->GetMaximum()*5);
	
	Leg.AddEntry(hName["h_m_y0"],"|y| < 0.6","P");// , 10 < p_{T} < 100 GeV","P"); 
	Leg.AddEntry(hName["h_m_y1"],"0.6 < |y| < 1.2", "P");// , 10 < p_{T} < 100 GeV","P"); 
	Leg.SetTextSize(0.045);
	Leg.DrawClone(); 

	double x1=0.4;
	double y1=0.65; 
	
	TLatex KC; 
	KC.SetTextSize(0.05); 
	KC.SetNDC(kTRUE);
	x1=0.48;
	y1=0.62;
    
	KC.DrawLatex(x1,y1,"10 < p_{T} < 100 GeV");
	if(print_muon_cuts)KC.DrawLatex(x1,y1, "1.4 < |#eta(#mu)| < 1.6: p_{T}(#mu) > 3 GeV");
	y1-=0.067; 
	if(print_muon_cuts)KC.DrawLatex(x1,y1, "1.2 < |#eta(#mu)| < 1.4: p_{T}(#mu) > 3.5 GeV");
	y1-=0.067; 
	if(print_muon_cuts)KC.DrawLatex(x1,y1, "|#eta(#mu)| < 1.2: p_{T}(#mu) > 4.5 GeV");
	
	
	KC.SetTextSize(0.045); 
	KC.DrawLatex(0.36,0.24,"1S"); 
	KC.DrawLatex(0.52,0.19,"2S"); 
	KC.DrawLatex(0.61,0.16,"3S"); 
	
	//draw_header();
	CMS_lumi(CName["m_summary"],1,10);
	//line Upsilon 1S low
	TLine y1L(y1m[0],0,y1m[0],0.7*hName["h_m_y0"]->GetBinContent(hName["h_m_y0"]->FindBin(PDG_mass[0]))); 
	y1L.DrawClone(); 
	
	TLine y1H(y1m[1],0,y1m[1],0.7*hName["h_m_y0"]->GetBinContent(hName["h_m_y0"]->FindBin(PDG_mass[0]))); 
	y1H.DrawClone(); 
	//2S
	TLine y2L(y2m[0],0,y2m[0],hName["h_m_y0"]->GetBinContent(hName["h_m_y0"]->FindBin(PDG_mass[1]))); 
	y2L.DrawClone(); 
	
	TLine y2H(y2m[1],0,y2m[1],hName["h_m_y0"]->GetBinContent(hName["h_m_y0"]->FindBin(PDG_mass[1]))); 
	y2H.DrawClone(); 
	//3S
	
	TLine y3L(y3m[0],0,y3m[0],hName["h_m_y0"]->GetBinContent(hName["h_m_y0"]->FindBin(PDG_mass[2]))); 
	y3L.DrawClone(); 
	
	TLine y3H(y3m[1],0,y3m[1],hName["h_m_y0"]->GetBinContent(hName["h_m_y0"]->FindBin(PDG_mass[2]))); 
	y3H.DrawClone(); 
	
}

void draw_zeta_eta(){
	CreateCanvas("zeta_eta_sections","",Cx,Cy); 
	CName["zeta_eta_sections"]->cd();
	TLegend L(0.55,0.5,0.75,0.85);
	L.SetFillColor(10);
	L.SetLineColor(10);
	L.SetTextSize(0.05);
	for(int ieta=3; ieta>=0; ieta--){
		TString name = Form("zeta_eta%d",ieta); 
		if(ieta==3) {
			hName[name]->SetLineColor(kBlack); 
			hName[name]->SetLineWidth(2);
			hName[name]->SetLineStyle(9); 
		}
		if(ieta==2){
			hName[name]->SetLineColor(kBlue);
			hName[name]->SetLineStyle(7);
			hName[name]->SetLineWidth(2);

		}
		if(ieta==1){
			hName[name]->SetLineColor(kRed);
			hName[name]->SetLineStyle(5);
			hName[name]->SetLineWidth(2);


		}
		if(ieta==0){
			hName[name]->SetLineColor(28); 
			hName[name]->SetLineStyle(2);
			hName[name]->SetLineWidth(2);
		}
		hName[name]->GetYaxis()->SetTitle(Form("Events / %.2f MeV",hName[name]->GetBinWidth(1))); 
		if(ieta==3) L.AddEntry(hName[name]->Clone(name+Form("_%d",ieta)),"|#eta(#mu)| < 1.6","L");
		if(ieta==2) L.AddEntry(hName[name]->Clone(name+Form("_%d",ieta)),"1.0 < |#eta(#mu)| < 1.6","L");
		if(ieta==1) L.AddEntry(hName[name]->Clone(name+Form("_%d",ieta)),"0.5 < |#eta(#mu)| < 1.0","L");
		if(ieta==0) L.AddEntry(hName[name]->Clone(name+Form("_%d",ieta)),"|#eta(#mu)| < 0.5","L");

		hName[name]->SetTitle(""); 
		hName[name]->SetMaximum(hName[name]->GetMaximum()*1.3);
		if(ieta==3)hName[name]->Draw("histo"); 
		else hName[name]->Draw("histo same"); 
	}
	L.DrawClone("same");
	//draw_header();
	CMS_lumi(CName["zeta_eta_sections"],1,10); 
		
	
}

void draw_zeta_m(){
	for(int iy=0; iy<fNy;iy++){
		for(int ipt=0; ipt<fNpt; ipt++){
			TString can_name=h_zeta_m(iy,ipt) + "_can"; 
			CreateCanvas(can_name,"",Cx,Cy); 
			CName[can_name]->cd(); 
			hName2D[h_zeta_m(iy,ipt)]->SetAxisRange(0,150,"Y"); 
			TLatex L(0.2,0.8,BT(iy,ipt));
			L.SetNDC(kTRUE);
			hName2D[h_zeta_m(iy,ipt)]->Draw("COLZ"); 
			draw_header();
			L.DrawClone(); 

		}
	}
}

void draw_zeta_profile(){
	for(int iy=0; iy<fNy;iy++){
		for(int ipt=0; ipt<fNpt; ipt++){
			TProfile *prof_dm=(TProfile*)hName2D[h_zeta_m(iy,ipt)]->ProfileX("prof_dm",0,hName2D[h_zeta_m(iy,ipt)]->GetNbinsY());
			prof_dm->GetYaxis()->SetTitle("#zeta Mean [MeV]"); 
			prof_dm->GetXaxis()->SetTitle("M_{#mu#mu} [GeV]"); 
			prof_dm->SetStats(kFALSE);
			prof_dm->SetMinimum(0);
			prof_dm->SetMarkerSize(0.5);
			TString can_name="profile_"+h_zeta_m(iy,ipt);
			CreateCanvas(can_name,"",Cx,Cy); 
			TLatex L(0.2,0.81,BT(iy,ipt)); 
			L.SetNDC(kTRUE);
			prof_dm->DrawCopy(); 
			L.DrawClone(); 
			draw_header();

//			delete prof_dm; 
		}
	}
}

void draw_zeta(){

	for(int iy=0; iy<fNy;iy++){
		for(int ipt=0; ipt<fNpt; ipt++){
			for(int ups=1; ups<=3; ups++){
				TString can_name=h_zeta(iy,ipt,ups)+"_can";
				CreateCanvas(can_name,"",Cx,Cy); 
				CName[can_name]->cd(); 
				hName[h_zeta(iy,ipt,ups)]->DrawCopy("histo"); 
				draw_header();
				
				TLatex bin_label(0.35,0.83,BT(iy,ipt)); 
				bin_label.SetNDC(kTRUE);
				bin_label.DrawClone(); 
				TPad *nPad=new TPad(Form("zetaPad_y%d_pt%d_ups%d",iy,ipt,ups),"",0.3,0.3,0.8,0.8); 
				nPad->Draw(); 
				nPad->cd();
				
				TH1F *z=(TH1F*)hName[h_zeta(iy,ipt,ups)]->Clone(h_zeta(iy,ipt,ups)+"_zoom"); 
				z->SetAxisRange(40,150); 
				z->DrawCopy("histo"); 
				delete z; 
			}//ups
		}//ipt
	}//iy
	
}

void plot_atlas_pt(int ups){
	double SF=1000/2.4; // 1000 for fb->pb, 0.6 for rapidity region
	
	// Note that Atlas plots d^2sigma/dPtdy... so we have to account for this
	
	if(ups==1){
		// Plot: p8322_d6x1y2
		double p8322_d6x1y1_xval[] = { 0.25, 0.75, 1.25, 1.75, 2.25, 2.75, 3.25, 3.75, 4.25, 
			4.75, 5.25, 5.75, 6.25, 6.75, 7.25, 7.75, 8.25, 8.75, 9.25, 
			9.75, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 
			19.5, 20.5, 21.5, 22.5, 23.5, 24.5, 25.5, 26.5, 27.5, 28.5, 
			29.5, 30.5, 31.5, 33.0, 35.0, 37.0, 39.0, 42.5, 47.5, 55.0, 
			65.0 };
		/*
		 double p8322_d6x1y1_xerrminus[] = { 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 
		 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 
		 0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
		 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
		 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 2.5, 2.5, 5.0, 
		 5.0 };
		 
		 double p8322_d6x1y1_xerrplus[] = { 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 
		 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 
		 0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
		 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
		 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 2.5, 2.5, 5.0, 
		 5.0 };
		 
		 */
		double p8322_d6x1y1_yval[] = { 35392.84, 98878.16, 148372.15, 182173.1, 206775.92, 231574.2, 242509.16, 240950.58, 229797.15, 
			208978.78, 184938.25, 179761.51, 157285.63, 136729.51, 122675.42, 108822.11, 94001.57, 82014.84, 71983.66, 
			62278.51, 49124.29, 35832.87, 27727.3, 20263.43, 16302.79, 12331.68, 9516.72, 7338.76, 5701.8, 
			4476.17, 3595.25, 2801.82, 2288.65, 1867.8, 1472.76, 1240.35, 988.08, 780.36, 662.63, 
			599.29, 488.3, 374.2, 309.67, 236.99, 164.84, 126.49, 74.05, 41.24, 18.07, 
			7.28 };
		double p8322_d6x1y1_yerrminus[] = { 10040.187744644021, 29165.59695075347, 44279.714260927656, 53244.39693433479, 58505.326055196034, 66398.84614706342, 67324.27643658638, 67963.08338160136, 65933.42573497679, 
			61006.88931563139, 54157.40316778953, 52706.095821221665, 46209.00046205717, 38514.3298287196, 35650.19575819044, 32012.60791104811, 26752.9963248792, 23204.191100292635, 20099.002776060806, 
			17090.635482351146, 13269.239589588395, 9396.921362808142, 7052.913343363862, 5006.5960821500275, 3853.3447551445483, 2811.3084500104214, 2090.0416370493676, 1552.4720581382455, 1173.3702550772282, 
			886.0189593343925, 684.7933779615571, 519.8145781911085, 413.08033613329985, 327.06105760239933, 251.20286105058597, 204.92675276790973, 158.94243014374732, 122.73967655163509, 103.61023260276951, 
			89.65957282967615, 69.64123275761278, 54.57913612361413, 42.612845481145705, 31.766877403987948, 21.940068368170596, 16.76881927864929, 8.930761445699913, 5.295431993709295, 2.098952119511067, 
			1.0630616162763098 };
		double p8322_d6x1y1_yerrplus[] = { 8699.38653125035, 32576.66833363719, 76805.90114189795, 236700.99427909675, 453550.8479124082, 624837.9455057994, 712761.5664403944, 718525.2938768687, 666355.8395618561, 
			575792.9900605836, 468661.24430696777, 419350.19024839095, 327654.08364405914, 253070.89968113083, 197740.74979823153, 147822.2229745988, 104222.93167987072, 75840.32422017855, 56664.81719483528, 
			43287.58068620376, 29173.177920019272, 18089.10823410596, 11864.652619925288, 7629.401854418208, 5552.797011849073, 3822.807598022689, 2722.015197992105, 1958.100781088655, 1400.19382069055, 
			1038.064565718337, 789.4447244107722, 586.6845467540456, 455.53080170280475, 355.2579786577636, 269.27636565432175, 218.1831874824456, 168.5934043787004, 128.73072088666325, 103.87868212487103, 
			93.38044441958927, 75.79567072069486, 56.37493592014096, 43.5152789259129, 32.44108660325668, 22.11096108268476, 16.809961332495682, 9.070622911355096, 5.281098370604358, 2.1842389979120873, 
			1.0938921336219583 };
		double p8322_d6x1y1_ystatminus[] = { 936.86, 2155.91, 2587.05, 2985.8, 3191.88, 2809.87, 3533.43, 3492.8, 3058.9, 
			2937.09, 2939.08, 2576.18, 2420.64, 1824.43, 1731.36, 1984.32, 1559.85, 1322.67, 1157.24, 
			746.27, 663.51, 426.62, 319.35, 332.39, 182.0, 142.19, 126.55, 107.58, 88.57, 
			76.01, 63.94, 56.97, 51.14, 46.03, 40.61, 37.7, 32.44, 30.2, 26.21, 
			24.61, 22.76, 19.91, 12.43, 11.26, 9.21, 7.96, 4.0, 2.96, 1.38, 
			0.89 };
		double p8322_d6x1y1_ystatplus[] = { 936.86, 2155.91, 2587.05, 2985.8, 3191.88, 2809.87, 3533.43, 3492.8, 3058.9, 
			2937.09, 2939.08, 2576.18, 2420.64, 1824.43, 1731.36, 1984.32, 1559.85, 1322.67, 1157.24, 
			746.27, 663.51, 426.62, 319.35, 332.39, 182.0, 142.19, 126.55, 107.58, 88.57, 
			76.01, 63.94, 56.97, 51.14, 46.03, 40.61, 37.7, 32.44, 30.2, 26.21, 
			24.61, 22.76, 19.91, 12.43, 11.26, 9.21, 7.96, 4.0, 2.96, 1.38, 
			0.89 };
		int p8322_d6x1y1_numpoints = 50;
		
		double p8322_d6x1y1_xerrminus[]={ 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 
			0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 
			0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
			0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
			0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 2.5, 2.5, 5.0, 
			5.0 };
		double p8322_d6x1y1_xerrplus[]={ 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 
			0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 
			0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
			0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
			0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 2.5, 2.5, 5.0, 
			5.0 };
		
		TGraphAsymmErrors *p8322_d6x1y1 =new TGraphAsymmErrors(p8322_d6x1y1_numpoints, p8322_d6x1y1_xval, p8322_d6x1y1_yval, p8322_d6x1y1_xerrminus, p8322_d6x1y1_xerrplus, p8322_d6x1y1_yerrminus, p8322_d6x1y1_yerrplus);
		scale_graph(p8322_d6x1y1,(1./SF));
		p8322_d6x1y1->SetName(Form("Atlas_xs_%dS",ups));
		grName[p8322_d6x1y1->GetName()]=p8322_d6x1y1; 
		
	} //ups 1S
	
	if(ups==2){
	double p8322_d7x1y1_xval[] = { 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 
		9.5, 11.0, 13.0, 15.0, 17.0, 19.0, 21.0, 23.0, 25.0, 27.0, 
		29.0, 31.0, 34.0, 38.0, 45.0, 60.0 };
	/*
	 double p8322_d7x1y1_xerrminus[] = { 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
	 0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
	 1.0, 1.0, 2.0, 2.0, 5.0, 10.0 };
	 double p8322_d7x1y1_xerrplus[] = { 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
	 0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
	 1.0, 1.0, 2.0, 2.0, 5.0, 10.0 };
	 */
	double p8322_d7x1y1_yval[] = { 14099.19, 37169.61, 43895.41, 53632.18, 51016.7, 41440.45, 35410.02, 29636.38, 22632.57, 
		19053.93, 12611.54, 7939.33, 5343.06, 3309.25, 2072.13, 1337.29, 897.91, 608.07, 418.21, 
		296.1, 209.3, 128.9, 74.37, 31.83, 7.51 };
	double p8322_d7x1y1_yerrminus[] = { 3915.9007974794254, 10660.241072564917, 12853.494024929563, 14293.457798353062, 13238.90561274987, 11292.63679040905, 9961.513892185263, 8653.921212323346, 6430.445462843145, 
		5246.714875557847, 3297.1416707960852, 1967.388458693402, 1217.0259453273788, 705.8831526534685, 410.33545277004765, 246.36233255106188, 157.04912002300426, 101.24744243683392, 67.16017048816954, 
		44.87941398904402, 31.301439264033853, 17.62770830255595, 10.05802167426577, 3.9284348028190563, 1.1434159348198711 };
	double p8322_d7x1y1_yerrplus[] = { 3584.398011884283, 14539.50576821647, 48442.20717914285, 101196.61661934454, 106239.46591437947, 80524.4152455173, 58953.63710351805, 39473.663513499734, 22059.308130478617, 
		13904.259904140888, 7055.801357018209, 3287.4933701225923, 1749.2841637652814, 920.9084363279554, 497.0995806676968, 286.2217689135472, 173.65078692594514, 108.14296371008149, 70.38218595639098, 
		47.09316192399911, 32.442011343318406, 17.764475787368454, 10.074691062260918, 3.9676315353117153, 1.162497311824849 };
	double p8322_d7x1y1_ystatminus[] = { 462.5, 901.0, 1013.67, 1219.81, 1122.06, 989.2, 818.33, 1227.45, 664.19, 
		465.43, 283.93, 153.6, 82.14, 59.1, 43.89, 33.56, 27.06, 21.34, 18.39, 
		14.47, 12.56, 6.82, 5.13, 2.21, 0.67 };
	double p8322_d7x1y1_ystatplus[] = { 462.5, 901.0, 1013.67, 1219.81, 1122.06, 989.2, 818.33, 1227.45, 664.19, 
		465.43, 283.93, 153.6, 82.14, 59.1, 43.89, 33.56, 27.06, 21.34, 18.39, 
		14.47, 12.56, 6.82, 5.13, 2.21, 0.67 };
	int p8322_d7x1y1_numpoints = 25;
	
	double p8322_d7x1y1_xerrminus[]={ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
		0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
		1.0, 1.0, 2.0, 2.0, 5.0, 10.0 };
	double p8322_d7x1y1_xerrplus[]= { 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
		0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
		1.0, 1.0, 2.0, 2.0, 5.0, 10.0 };
		
		TGraphAsymmErrors *p8322_d7x1y1 =new TGraphAsymmErrors(p8322_d7x1y1_numpoints, p8322_d7x1y1_xval, p8322_d7x1y1_yval, p8322_d7x1y1_xerrminus, p8322_d7x1y1_xerrplus, p8322_d7x1y1_yerrminus, p8322_d7x1y1_yerrplus);
		scale_graph(p8322_d7x1y1,1./SF); 
		p8322_d7x1y1->SetName(Form("Atlas_xs_%dS",ups));
		grName[p8322_d7x1y1->GetName()]=p8322_d7x1y1;
		
}//Ups 2S	
	
	if(ups==3){
		// Plot: p8322_d8x1y2
		double p8322_d8x1y1_xval[] = { 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 
			9.5, 11.0, 13.0, 15.0, 17.0, 19.0, 21.0, 23.0, 25.0, 27.0, 
			29.0, 31.0, 34.0, 38.0, 45.0, 60.0 };
		
		/*
		 double p8322_d8x1y1_xerrminus[] = { 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
		 0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
		 1.0, 1.0, 2.0, 2.0, 5.0, 10.0 };
		 double p8322_d8x1y1_xerrplus[] = { 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
		 0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
		 1.0, 1.0, 2.0, 2.0, 5.0, 10.0 };
		 */
		
		double p8322_d8x1y1_yval[] = { 6022.32, 14021.13, 20711.77, 20317.61, 19951.64, 18330.26, 17320.64, 14677.17, 12151.44, 
			9357.4, 6510.87, 4765.1, 3157.02, 1992.37, 1334.81, 946.0, 610.91, 437.57, 331.78, 
			208.81, 161.25, 85.68, 53.61, 19.27, 5.49 };
		double p8322_d8x1y1_yerrminus[] = { 1572.0742719413736, 3795.0331310016254, 5680.281874528763, 5146.731171403846, 4966.7803433008785, 4885.1679069403535, 5283.628766387737, 4121.8282818550315, 3441.0830164644385, 
			2577.152806102114, 1724.769824179447, 1227.0599322771484, 718.3371423642244, 425.4963911715351, 266.156395564713, 178.58454944367386, 107.63820371968309, 73.39877110687888, 55.09720319580659, 
			34.224321760993305, 24.980154122823183, 12.657551105960426, 8.39753535270915, 2.8718460961548757, 1.1015443704181869 };
		double p8322_d8x1y1_yerrplus[] = { 1487.2148013316703, 4576.14641170931, 14303.538135066443, 28691.751653766278, 33960.546032234226, 31241.88712183373, 25458.392024750898, 18944.708991259275, 11234.75530593791, 
			6797.286830199238, 3500.615076811502, 1952.7710345557668, 1037.0850676776713, 552.0350199942029, 322.970171068475, 204.0989519816307, 119.8237192712695, 79.175352225298, 57.12247368593205, 
			33.601421100899884, 25.62148707628033, 12.793592146070626, 7.988372800514507, 2.885186302476843, 1.1015443704181869 };
		double p8322_d8x1y1_ystatminus[] = { 357.98, 537.42, 639.56, 820.97, 750.21, 762.74, 698.45, 583.39, 587.97, 
			401.56, 233.58, 135.86, 73.7, 53.92, 40.56, 31.02, 24.73, 20.24, 17.8, 
			13.49, 12.29, 6.38, 4.99, 2.05, 0.71 };
		double p8322_d8x1y1_ystatplus[] = { 357.98, 537.42, 639.56, 820.97, 750.21, 762.74, 698.45, 583.39, 587.97, 
			401.56, 233.58, 135.86, 73.7, 53.92, 40.56, 31.02, 24.73, 20.24, 17.8, 
			13.49, 12.29, 6.38, 4.99, 2.05, 0.71 };
		const int p8322_d8x1y1_numpoints = 25;
		
		double p8322_d8x1y1_xerrminus[]={ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
			0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
			1.0, 1.0, 2.0, 2.0, 5.0, 10.0 };
		double p8322_d8x1y1_xerrplus[]={ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
			0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
			1.0, 1.0, 2.0, 2.0, 5.0, 10.0 };
		
		TGraphAsymmErrors *p8322_d8x1y1 =new TGraphAsymmErrors(p8322_d8x1y1_numpoints, p8322_d8x1y1_xval, p8322_d8x1y1_yval, p8322_d8x1y1_xerrminus, p8322_d8x1y1_xerrplus, p8322_d8x1y1_yerrminus, p8322_d8x1y1_yerrplus);
		
		scale_graph(p8322_d8x1y1,1./SF); 
		p8322_d8x1y1->SetName(Form("Atlas_xs_%dS",ups));
		grName[p8322_d8x1y1->GetName()]=p8322_d8x1y1;
		
	}

}

void get_atlas_ratios(){
    
    // Plot: p8322_d11x1y2
    double p8322_d11x1y1_xval[] = { 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5,
        9.5, 11.0, 13.0, 15.0, 17.0, 19.0, 21.0, 23.0, 25.0, 27.0,
        29.0, 31.0, 34.0, 38.0, 45.0, 60.0 };
    double p8322_d11x1y1_xerrminus[] = { 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
        0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 2.0, 2.0, 5.0, 10.0 };
    double p8322_d11x1y1_xerrplus[] = { 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
        0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 2.0, 2.0, 5.0, 10.0 };
    double p8322_d11x1y1_yval[] = { 0.09, 0.08, 0.1, 0.08, 0.09, 0.1, 0.12, 0.13, 0.14,
        0.14, 0.15, 0.2, 0.22, 0.24, 0.26, 0.3, 0.29, 0.32, 0.37,
        0.33, 0.38, 0.31, 0.36, 0.34, 0.45 };
    double p8322_d11x1y1_yerrminus[] = { 0.01414213562373095, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01, 0.01414213562373095, 0.01414213562373095,
        0.01414213562373095, 0.01414213562373095, 0.01414213562373095, 0.01, 0.01, 0.01, 0.01, 0.01, 0.02, 0.02,
        0.02, 0.03, 0.03, 0.04, 0.04123105625617661, 0.08 };
    double p8322_d11x1y1_yerrplus[] = { 0.01414213562373095, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01, 0.01414213562373095, 0.01414213562373095,
        0.01414213562373095, 0.01414213562373095, 0.01414213562373095, 0.01, 0.01, 0.01, 0.01, 0.01, 0.02, 0.02,
        0.02, 0.03, 0.03, 0.04, 0.04123105625617661, 0.08 };
    double p8322_d11x1y1_ystatminus[] = { 0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01, 0.01,
        0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.02, 0.02,
        0.02, 0.03, 0.03, 0.04, 0.04, 0.08 };
    double p8322_d11x1y1_ystatplus[] = { 0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01, 0.01,
        0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.02, 0.02,
        0.02, 0.03, 0.03, 0.04, 0.04, 0.08 };
    int p8322_d11x1y1_numpoints = 25;
    TGraphAsymmErrors *gr3S1S = new TGraphAsymmErrors(p8322_d11x1y1_numpoints, p8322_d11x1y1_xval, p8322_d11x1y1_yval, p8322_d11x1y1_xerrminus, p8322_d11x1y1_xerrplus, p8322_d11x1y1_yerrminus, p8322_d11x1y1_yerrplus);
    gr3S1S->SetName("AtlasRatio_3S1S_1p2");
    grName[gr3S1S->GetName()]=gr3S1S;
 
    // Plot: p8322_d10x1y2
    double p8322_d10x1y1_xval[] = { 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5,
        9.5, 11.0, 13.0, 15.0, 17.0, 19.0, 21.0, 23.0, 25.0, 27.0,
        29.0, 31.0, 34.0, 38.0, 45.0, 60.0 };
    double p8322_d10x1y1_xerrminus[] = { 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
        0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 2.0, 2.0, 5.0, 10.0 };
    double p8322_d10x1y1_xerrplus[] = { 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
        0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 2.0, 2.0, 5.0, 10.0 };
    double p8322_d10x1y1_yval[] = { 0.21, 0.22, 0.2, 0.22, 0.23, 0.23, 0.24, 0.26, 0.26,
        0.28, 0.3, 0.33, 0.37, 0.39, 0.41, 0.42, 0.43, 0.45, 0.47,
        0.47, 0.49, 0.47, 0.51, 0.54, 0.62 };
    double p8322_d10x1y1_yerrminus[] = { 0.01414213562373095, 0.01414213562373095, 0.01, 0.01, 0.01, 0.01414213562373095, 0.01414213562373095, 0.01414213562373095, 0.01414213562373095,
        0.01414213562373095, 0.01414213562373095, 0.01414213562373095, 0.01414213562373095, 0.01414213562373095, 0.01, 0.01414213562373095, 0.01414213562373095, 0.022360679774997897, 0.022360679774997897,
        0.03162277660168379, 0.03162277660168379, 0.03, 0.04123105625617661, 0.050990195135927854, 0.08246211251235322 };
    double p8322_d10x1y1_yerrplus[] = { 0.01414213562373095, 0.01414213562373095, 0.01, 0.01, 0.01, 0.01414213562373095, 0.01414213562373095, 0.01414213562373095, 0.01414213562373095,
        0.01414213562373095, 0.01414213562373095, 0.01414213562373095, 0.01414213562373095, 0.01414213562373095, 0.01, 0.01414213562373095, 0.01414213562373095, 0.022360679774997897, 0.022360679774997897,
        0.03162277660168379, 0.03162277660168379, 0.03, 0.04123105625617661, 0.050990195135927854, 0.08246211251235322 };
    double p8322_d10x1y1_ystatminus[] = { 0.01, 0.01, 0.0, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
        0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.02, 0.02,
        0.03, 0.03, 0.03, 0.04, 0.05, 0.08 };
    double p8322_d10x1y1_ystatplus[] = { 0.01, 0.01, 0.0, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
        0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.02, 0.02,
        0.03, 0.03, 0.03, 0.04, 0.05, 0.08 };
    int p8322_d10x1y1_numpoints = 25;
    TGraphAsymmErrors *gr2S1S = new TGraphAsymmErrors(p8322_d10x1y1_numpoints, p8322_d10x1y1_xval, p8322_d10x1y1_yval, p8322_d10x1y1_xerrminus, p8322_d10x1y1_xerrplus, p8322_d10x1y1_yerrminus, p8322_d10x1y1_yerrplus);
    gr2S1S->SetName("AtlasRatio_2S1S_1p2");
    grName[gr2S1S->GetName()]=gr2S1S;
}

void rapidity_dep_atlas(){
	// Plot: p8322_d5x1y3
	double xval[] = { 0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 
		0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 
		0.975, 1.025, 1.075, 1.125, 1.175, 1.225, 1.275, 1.325, 1.375, 1.425, 
		1.475, 1.525, 1.575, 1.625, 1.675, 1.725, 1.775, 1.825, 1.875, 1.925, 
		1.975, 2.025, 2.075, 2.125, 2.175, 2.225 };
	double xerrminus[] = { 0.025, 0.024999999999999994, 0.024999999999999994, 0.024999999999999994, 0.024999999999999994, 0.025000000000000022, 0.025000000000000022, 0.025000000000000022, 0.024999999999999967, 
		0.024999999999999967, 0.025000000000000022, 0.02499999999999991, 0.025000000000000022, 0.025000000000000022, 0.025000000000000022, 0.025000000000000022, 0.02499999999999991, 0.025000000000000022, 0.025000000000000022, 
		0.025000000000000022, 0.02499999999999991, 0.02499999999999991, 0.02499999999999991, 0.025000000000000133, 0.025000000000000133, 0.02499999999999991, 0.02499999999999991, 0.02499999999999991, 0.025000000000000133, 
		0.025000000000000133, 0.02499999999999991, 0.02499999999999991, 0.02499999999999991, 0.025000000000000133, 0.025000000000000133, 0.02499999999999991, 0.02499999999999991, 0.02499999999999991, 0.025000000000000133, 
		0.025000000000000133, 0.02499999999999991, 0.025000000000000355, 0.02499999999999991, 0.02499999999999991, 0.02499999999999991 };
	double xerrplus[] = { 0.025, 0.02500000000000001, 0.024999999999999994, 0.025000000000000022, 0.024999999999999994, 0.024999999999999967, 0.024999999999999967, 0.025000000000000022, 0.025000000000000022, 
		0.025000000000000022, 0.025000000000000022, 0.025000000000000022, 0.025000000000000022, 0.02499999999999991, 0.025000000000000022, 0.025000000000000022, 0.025000000000000022, 0.025000000000000022, 0.02499999999999991, 
		0.025000000000000022, 0.025000000000000133, 0.025000000000000133, 0.02499999999999991, 0.02499999999999991, 0.02499999999999991, 0.025000000000000133, 0.025000000000000133, 0.02499999999999991, 0.02499999999999991, 
		0.02499999999999991, 0.025000000000000133, 0.025000000000000133, 0.02499999999999991, 0.02499999999999991, 0.02499999999999991, 0.025000000000000133, 0.025000000000000133, 0.02499999999999991, 0.02499999999999991, 
		0.02499999999999991, 0.02499999999999991, 0.02499999999999991, 0.02499999999999991, 0.025000000000000355, 0.02499999999999991 };
	double yval[] = { 470.342, 467.979, 476.144, 467.122, 471.467, 471.223, 474.822, 469.293, 468.917, 
		489.283, 490.82, 463.538, 471.25, 469.363, 475.7, 462.838, 477.397, 463.457, 457.915, 
		452.526, 474.651, 467.965, 480.585, 500.732, 489.191, 493.106, 503.849, 511.277, 496.55, 
		482.684, 460.354, 469.869, 467.296, 472.721, 448.888, 435.157, 403.599, 394.87, 359.604, 
		317.949, 264.182, 233.939, 162.114, 120.417, 62.35 };
	double yerrminus[] = { 28.94770464130101, 22.733582779667618, 21.850906846170023, 22.210643304506064, 22.16478865227458, 23.105940578128386, 24.5578927638346, 22.100128234922078, 31.61897773489839, 
		23.54050691892594, 25.18560733831924, 28.758150462086395, 43.04854211236427, 25.072283761157458, 22.682841665011903, 35.543688398926754, 27.885868284849945, 22.98816880484394, 22.958584625363994, 
		26.06250527098268, 39.73932315729597, 57.76217143598395, 48.352155577595504, 55.46015186600196, 51.93031725302667, 35.08147274559607, 38.24690233208436, 38.919451691923925, 40.15135270448556, 
		35.66673902952161, 35.02917870861376, 36.389579524913444, 34.99792408129373, 41.04011653248563, 56.281983778115, 49.961942015898465, 49.14916313631393, 48.01041972322258, 49.09009735170628, 
		28.234585883274438, 24.997221125557136, 21.51739505609357, 14.021974789593655, 17.415829724707347, 10.640816697979531 };
	double yerrplus[] = { 28.94770464130101, 22.733582779667618, 21.850906846170023, 22.210643304506064, 22.16478865227458, 23.105940578128386, 24.5578927638346, 22.100128234922078, 31.61897773489839, 
		23.54050691892594, 25.18560733831924, 28.758150462086395, 43.04854211236427, 25.072283761157458, 22.682841665011903, 35.543688398926754, 27.885868284849945, 22.98816880484394, 22.958584625363994, 
		26.06250527098268, 39.73932315729597, 57.76217143598395, 48.352155577595504, 55.46015186600196, 51.93031725302667, 35.08147274559607, 38.24690233208436, 38.919451691923925, 40.15135270448556, 
		35.66673902952161, 35.02917870861376, 36.389579524913444, 34.99792408129373, 41.04011653248563, 56.281983778115, 49.961942015898465, 49.14916313631393, 48.01041972322258, 49.09009735170628, 
		28.234585883274438, 24.997221125557136, 21.51739505609357, 14.021974789593655, 17.415829724707347, 10.640816697979531 };
	double ystatminus[] = { 4.77, 5.035, 4.783, 4.826, 4.884, 4.651, 4.776, 5.588, 5.112, 
		5.365, 5.561, 6.547, 6.513, 8.233, 6.575, 5.747, 6.053, 6.456, 5.652, 
		5.985, 6.638, 6.232, 6.465, 6.902, 7.451, 6.463, 10.217, 7.954, 8.07, 
		7.712, 7.825, 6.027, 8.623, 9.067, 10.897, 8.375, 7.186, 9.609, 6.727, 
		9.228, 7.858, 7.683, 4.489, 4.713, 4.202 };
	double ystatplus[] = { 4.77, 5.035, 4.783, 4.826, 4.884, 4.651, 4.776, 5.588, 5.112, 
		5.365, 5.561, 6.547, 6.513, 8.233, 6.575, 5.747, 6.053, 6.456, 5.652, 
		5.985, 6.638, 6.232, 6.465, 6.902, 7.451, 6.463, 10.217, 7.954, 8.07, 
		7.712, 7.825, 6.027, 8.623, 9.067, 10.897, 8.375, 7.186, 9.609, 6.727, 
		9.228, 7.858, 7.683, 4.489, 4.713, 4.202 };
	int numpoints = 45;
	CreateCanvas("tmp","tmp",Cx,Cy); 
	CName["tmp"]->cd(); 
	TGraphAsymmErrors *p8322_d5x1y1 = new TGraphAsymmErrors(numpoints, xval, yval, xerrminus, xerrplus, yerrminus, yerrplus);
	p8322_d5x1y1->SetName("d5x1y1");
	p8322_d5x1y1->SetTitle("d5x1y1");
	p8322_d5x1y1->Draw("AP");
	
	double y0=0; 
	double y0E=0; 
	
	double y1=0; 
	double y1E=0; 
	
	for(int i=0; i<numpoints; i++){
		if(p8322_d5x1y1->GetX()[i]<0.6) {
			y0+=p8322_d5x1y1->GetY()[i]; 
			y0E+=yerrplus[i]*yerrplus[i];
		}
		
		if(p8322_d5x1y1->GetX()[i]>0.6 && p8322_d5x1y1->GetX()[i]<1.2){
			y1+=p8322_d5x1y1->GetY()[i]; 
			y1E+=yerrplus[i]*yerrplus[i];
		}
		
	}
	
	y0E=TMath::Sqrt(y0E);
	y1E=TMath::Sqrt(y1E);
	
	//cout << "y0: "<< y0 << " + " << y0E << endl;
	//cout << "y1: " << y1 << " + " << y1E << endl; 
	
	//cout << "Ratio: " << y0/y1 << endl; 
	
	double ratioUnc=(y0/y1)*TMath::Sqrt(TMath::Power(y0E/y0,2)+TMath::Power(y1E/y1,2)); 
	
	TF1 *F=new TF1("xs_atlas_y0_y1","pol0",4,70); 
	F->SetParameter(0,y0/y1);
	F->SetParError(0,ratioUnc); 
	F->SetLineColor(kBlack); 
	f1Name[F->GetName()]=F; 
	
}

void rap_dep2010(){
	
	double x[]={1,3,5,7,9.5,13,32.5};
	double xE[]={1,1,1,1,1.5,2,17.5}; 
	double y0[]={0.216,0.387,0.355,0.224,0.190,0.0914,0.059};
	double y0E[]={15,12,11,10,8,7,7};
	
	double y1[]={0.22,0.409,0.367,0.231,0.180,0.0915,0.0492};
	double y1E[]={14,11,9,9,8,8,7};
	
	double y2[]={0.198,0.426,0.331,0.217,0.174,0.0879,0.0482};
	double y2E[]={12,10,8,9,7,8,8};
	
	double R01[7]; 
	double R02[7]; 
	
	
	double R01E[7]; 
	double R02E[7]; 
	
	for(int i=0; i<7; i++) {
		R01E[i]=TMath::Sqrt(y0E[i]*y0E[i]+y1E[i]*y1E[i])/100;
		R02E[i]=TMath::Sqrt(y0E[i]*y0E[i]+y2E[i]*y2E[i])/100;
		
		
		y0E[i]=y0E[i]*y0[i]/100;
		y1E[i]=y1E[i]*y1[i]/100;
		y2E[i]=y2E[i]*y2[i]/100;
		
		R01[i]=y0[i]/y1[i]; 
		R02[i]=y0[i]/y2[i];
		
		R01E[i]=R01[i]*R01E[i]; 
		R02E[i]=R02[i]*R02E[i]; 
		
	}
	
	TGraphAsymmErrors *ratio_01=new TGraphAsymmErrors(7, x,R01,xE,xE,R01E,R01E);
	ratio_01->SetName("Ratio_CMS2010_y0p8_y0p4"); 
	ratio_01->SetMarkerColor(kRed);
	ratio_01->SetMarkerStyle(24);
	ratio_01->SetMarkerSize(1); 
	ratio_01->SetTitle(""); 
	
	
	TGraphAsymmErrors *ratio_02=new TGraphAsymmErrors(7,x,R02,xE,xE,R02E,R02E); 
	ratio_02->SetMarkerStyle(25);
	ratio_02->SetMarkerColor(kGreen); 
	ratio_02->SetLineColor(kGreen); 
	ratio_02->SetTitle(""); 
	ratio_02->SetName("Ratio_CMS2010_y1p2_y0p4"); 

	grName[ratio_01->GetName()]=ratio_01; 
	grName[ratio_02->GetName()]=ratio_02; 

	
		
}

void efficiency_closure(){
	int iy=2;
	int ups=1; 
	CreateCanvas("eff_closure","",Cx,Cy); 
	CName["eff_closure"]->cd();
	
	TProfile pr_scaled;
	TProfile pr_scaledEp;
	TProfile pr_scaledEm;
	
	prName[eff_name(iy, ups, "")]->Copy(pr_scaled);
	prName[eff_name(iy, ups, "Ep")]->Copy(pr_scaledEp);
	prName[eff_name(iy, ups, "Em")]->Copy(pr_scaledEm);

	pr_scaled.Scale(0.5);
	pr_scaledEm.Scale(0.5);
	pr_scaledEp.Scale(0.5);
	
	pr_scaled.SetMinimum(0);
	pr_scaled.SetMaximum(1); 
	pr_scaled.DrawCopy();
	pr_scaledEm.SetLineColor(kGreen);
	pr_scaledEp.SetLineColor(kRed); 
	pr_scaledEm.DrawCopy("histo same");
	pr_scaledEp.DrawCopy("histo same");
	hName["efftot_pt_TH1"]->SetLineColor(kBlack);
	hName["efftot_pt_TH1"]->SetLineStyle(kDashed);
	hName["efftot_pt_TH1"]->SetLineWidth(2);
	
	hName["efftot_pt_TH1"]->Draw("histo same");	
	
	TLegend Leg(0.22,0.57,0.62,0.87);
	Leg.SetFillColor(10);
	Leg.SetLineColor(10); 
	Leg.SetBorderSize(0);
	
	Leg.AddEntry(prName[eff_name(iy, ups, "")],"#epsilon_{#mu#mu}=#frac{1}{2} #epsilon_{1}(#mu) #epsilon_{1}(#mu) #rho: Data","P"); 
	Leg.AddEntry(prName[eff_name(iy, ups, "Ep")],"+1 #sigma","L"); 
	Leg.AddEntry(prName[eff_name(iy, ups, "Em")],"-1 #sigma","L"); 
	Leg.AddEntry(hName["efftot_pt_TH1"],"#epsilon_{#mu#mu}: MC","L"); 
	Leg.DrawClone();
	
}

void efficiency_table(int iy){
	ofstream output;
	output.open(Form("/uscms/home/btcarlso/AN-UPS/tdr2/notes/AN-13-088/trunk/tables/efficiency_table_y%d.tex",iy));
	int N=4;
	output <<"\\begin{tabular}{";
	for(int i=0; i<=N; i++) output << "c";
	output << "}\\hline\\hline" << endl;
	
	output << "$p_{\\rm T}$ & $\\Upsilon(1S)$ & $\\Upsilon(2S)$ & $\\Upsilon(3S)$ \\\\ \\hline " << endl;
    
    for(int i=1; i<=prName[eff_name(iy,ups,"")]->GetNbinsX(); i++){
        for(int ups=1; ups<=3; ups++){
            double BW=prName[eff_name(iy,ups,"")]->GetBinWidth(i);
            double pt=prName[eff_name(iy,ups,"")]->GetBinCenter(i);
            double eff=prName[eff_name(iy,ups,"")]->GetBinContent(i);
            double effUp=TMath::Abs(prName[eff_name(iy,ups,"Ep")]->GetBinContent(i)-eff);
            double effUm=TMath::Abs(prName[eff_name(iy,ups,"Em")]->GetBinContent(i)-eff);

            double effU=TMath::Max(effUm,effUp);
            
            if(ups==1)output << Form("%.0f--%.0f & ",pt-BW/2,pt+BW/2);
            if(ups<3) output << Form("%.2f $\\pm$ %.2f &",eff,effU);
        }
        output << "\\\\" << endl;
    }
    output << "\\hline\\hline" << endl; 
    output.close();
}

void draw_efficiency(){
	for(int iy=0; iy<fNy;iy++){
		for(int ipt=0; ipt<fNpt; ipt++){
			for(int ups=1; ups<=3; ups++){
				TString can_name=eff_dist(iy,ipt,ups)+"_can";
				CreateCanvas(can_name,"",Cx,Cy); 
				double mean=hName[eff_dist(iy,ups,ups)]->GetMean(); 
				TLine line_mean(mean,0,mean,hName[eff_dist(iy,ups,ups)]->GetMaximum()); 
				hName[eff_dist(iy,ups,ups)]->DrawCopy(); 
				draw_header();
				TString massT;
				if(ups==1) massT=Form("%.2f < M_{#mu#mu} < %.2f GeV",y1m[0],y1m[1]); 
				if(ups==2) massT=Form("%.2f < M_{#mu#mu} < %.2f GeV",y2m[0],y2m[1]); 
				if(ups==3) massT=Form("%.2f < M_{#mu#mu} < %.2f GeV",y3m[0],y3m[1]); 
				
				TLatex L(0.2,0.8,massT); 
				L.SetNDC(kTRUE);
				L.DrawClone();
				TLatex L2(0.2,0.7,BT(iy,ipt)); 
				L2.SetNDC(kTRUE);
				L2.DrawClone(); 
				
				TString eff_mean=Form("<#epsilon_{#mu#mu}>=%.2f",mean); 
				TLatex L3(mean-0.3,0.6,eff_mean); 
				L3.SetNDC(kTRUE);
				L3.DrawClone(); 
				line_mean.DrawClone("same"); 
			}
		}
	}

	
	for(int iy=0; iy<fNy+1;iy++){
		for(int ups=1; ups<=3; ups++){
			
			CreateCanvas(Form("eff_summary_y%d_%dS",iy,ups),"",Cx,Cy); 
			CName[Form("eff_summary_y%d_%dS",iy,ups)]->cd();
			prName[eff_name(iy, ups, "")]->SetMinimum(0);
			prName[eff_name(iy, ups, "")]->SetMaximum(1); 
			prName[eff_name(iy, ups, "")]->Draw();
			prName[eff_name(iy, ups, "Em")]->SetLineColor(kGreen);
			prName[eff_name(iy, ups, "Ep")]->SetLineColor(kRed); 
			prName[eff_name(iy, ups, "Em")]->Draw("histo same");
			prName[eff_name(iy, ups, "Ep")]->Draw("histo same");
	
			TString BT;
			if(iy==0 && ups==1) BT=Form("|y| < 0.6, %.2f < M_{#mu#mu} < %.2f GeV",y1m[0],y1m[1]); 
			if(iy==0 && ups==2) BT=Form("|y| < 0.6, %.2f < M_{#mu#mu} < %.2f GeV",y2m[0],y2m[1]); 
			if(iy==0 && ups==3) BT=Form("|y| < 0.6, %.2f < M_{#mu#mu} < %.2f GeV",y3m[0],y3m[1]); 

			if(iy==1 && ups==1) BT=Form("0.6 < |y| < 1.2, %.2f < M_{#mu#mu} < %.2f GeV",y1m[0],y1m[1]); 
			if(iy==1 && ups==2) BT=Form("0.6 < |y| < 1.2, %.2f < M_{#mu#mu} < %.2f GeV",y2m[0],y2m[1]); 
			if(iy==1 && ups==3) BT=Form("0.6 < |y| < 1.2, %.2f < M_{#mu#mu} < %.2f GeV",y3m[0],y3m[1]); 
			
			TLatex L(0.2,0.65,BT); 
			L.SetNDC(kTRUE);
			L.DrawClone(); 
			draw_header();
			
			TLegend Leg(0.2,0.2,0.6,0.5);
			Leg.SetFillColor(10);
			Leg.SetLineColor(10); 
			Leg.SetBorderSize(0);
			Leg.AddEntry(eff_name(iy, ups, ""),"#epsilon_{#mu#mu}","LEP"); 
			Leg.AddEntry(eff_name(iy, ups, "Ep"),"#epsilon_{#mu#mu} +#sigma","L"); 
			Leg.AddEntry(eff_name(iy, ups, "Em"),"#epsilon_{#mu#mu} -#sigma","L"); 
			Leg.DrawClone("same");
		}
	}
	
	CreateHistogram("eff_y1_y0","","p_{T} [GeV]","#epsilon_{#mu#mu}(0.6 < |y| < 1.2)/#epsilon_{#mu#mu}(|y| < 0.6)",fNpt,fPTbin); 
	CreateHistogram("effEp_y1_y0","","p_{T} [GeV]","#epsilon_{#mu#mu}(0.6 < |y| < 1.2)/#epsilon_{#mu#mu}(|y| < 0.6)",fNpt,fPTbin); 
	CreateHistogram("effEm_y1_y0","","p_{T} [GeV]","#epsilon_{#mu#mu}(0.6 < |y| < 1.2)/#epsilon_{#mu#mu}(|y| < 0.6)",fNpt,fPTbin); 
	
	CreateHistogram("effsg_y1_y0","","p_{T} [GeV]","#epsilon_{sg}(0.6 < |y| < 1.2)/#epsilon_{sg}(|y| < 0.6)",fNpt,fPTbin); 
	
	
	hName["eff_y1_y0"]->SetLineColor(kBlack); 
	hName["effEp_y1_y0"]->SetLineColor(kRed); 
	hName["effEm_y1_y0"]->SetLineColor(kGreen); 
	
	for(int ipt=1; ipt<=fNpt; ipt++){
		float eff_num=prName[eff_name(1, 1, "")]->GetBinContent(ipt); 
		float eff_den=prName[eff_name(0, 1, "")]->GetBinContent(ipt); 
		
		hName["eff_y1_y0"]->SetBinContent(ipt,eff_num/eff_den); 
		hName["eff_y1_y0"]->SetBinError(ipt,0); 
		
		eff_num=prName[eff_name(1, 1, "Ep")]->GetBinContent(ipt); 
		eff_den=prName[eff_name(0, 1, "Ep")]->GetBinContent(ipt); 
		
		hName["effEp_y1_y0"]->SetBinContent(ipt,eff_num/eff_den); 
		hName["effEp_y1_y0"]->SetBinError(ipt,0); 
		
		eff_num=prName[eff_name(1, 1, "Em")]->GetBinContent(ipt); 
		eff_den=prName[eff_name(0, 1, "Em")]->GetBinContent(ipt); 
		
		hName["effEm_y1_y0"]->SetBinContent(ipt,eff_num/eff_den); 
		hName["effEm_y1_y0"]->SetBinError(ipt,0); 
		
		TString SG_name_num=Form("sg_y%d_TH1_%dS_rebin",1,1);
		TString SG_name_den=Form("sg_y%d_TH1_%dS_rebin",0,1); 
		
		eff_num=hName[SG_name_num]->GetBinContent(ipt); 
		eff_den=hName[SG_name_den]->GetBinContent(ipt); 
		
		hName["effsg_y1_y0"]->SetBinContent(ipt,eff_num/eff_den); 
		hName["effsg_y1_y0"]->SetBinError(ipt,0); 
		
	}
	
	CreateCanvas("eff_y1_y0_plot","",Cx,Cy);
	CName["eff_y1_y0_plot"]->cd(); 
	hName["eff_y1_y0"]->Draw();
	hName["effEp_y1_y0"]->Draw("same"); 
	hName["effEm_y1_y0"]->Draw("same"); 
	
	CreateCanvas("eff_sg_y1_y0_plot","",Cx,Cy);
	CName["eff_sg_y1_y0_plot"]->cd(); 
	hName["effsg_y1_y0"]->Draw();
	
	
	for(int iy=0; iy<fNy+1;iy++){
		int num=1; 
		for(int den=2; den<=3; den++){
			

			CreateHistogram(Form("eff_r1%d_y%d",den,iy),"","p_{T} [GeV]",Form("#epsilon_{#mu#mu}^{1S}/#epsilon_{#mu#mu}^{%dS}",den),fNpt,fPTbin); 
			CreateHistogram(Form("effEp_r1%d_y%d",den,iy),"","p_{T} [GeV]",Form("#epsilon_{#mu#mu}^{1S}/#epsilon_{#mu#mu}^{%dS}",den),fNpt,fPTbin); 
			CreateHistogram(Form("effEm_r1%d_y%d",den,iy),"","p_{T} [GeV]",Form("#epsilon_{#mu#mu}^{1S}/#epsilon_{#mu#mu}^{%dS}",den),fNpt,fPTbin); 
			
			CreateHistogram(Form("eff_sg_r1%d_y%d",den,iy),"","p_{T} [GeV]",Form("#epsilon_{sg}^{1S}/#epsilon_{sg}^{%dS}",den),fNpt,fPTbin); 
			
			for(int ipt=1; ipt<=fNpt; ipt++){
				float eff_num=prName[eff_name(iy, num, "")]->GetBinContent(ipt); 
				float eff_den=prName[eff_name(iy, den, "")]->GetBinContent(ipt); 
			
				hName[Form("eff_r1%d_y%d",den,iy)]->SetBinContent(ipt,eff_num/eff_den); 
				hName[Form("eff_r1%d_y%d",den,iy)]->SetBinError(ipt,0); 

				eff_num=prName[eff_name(iy, num, "Ep")]->GetBinContent(ipt); 
				eff_den=prName[eff_name(iy, den, "Ep")]->GetBinContent(ipt); 
				
				hName[Form("effEp_r1%d_y%d",den,iy)]->SetBinContent(ipt,eff_num/eff_den); 

				eff_num=prName[eff_name(iy, num, "Em")]->GetBinContent(ipt); 
				eff_den=prName[eff_name(iy, den, "Em")]->GetBinContent(ipt); 
				
				hName[Form("effEm_r1%d_y%d",den,iy)]->SetBinContent(ipt,eff_num/eff_den); 
				
				TString SG_name_num=Form("sg_y%d_TH1_%dS_rebin",iy,num);
				TString SG_name_den=Form("sg_y%d_TH1_%dS_rebin",iy,den); 
				if(iy>1){
					SG_name_num=Form("sg_TH1_%dS_rebin",num); 
					SG_name_den=Form("sg_TH1_%dS_rebin",den); 
				}
				float eff_sg_num=hName[SG_name_num]->GetBinContent(ipt); 
				float eff_sg_den=hName[SG_name_den]->GetBinContent(ipt); 
	
				hName[Form("eff_sg_r1%d_y%d",den,iy)]->SetBinContent(ipt,eff_sg_num/eff_sg_den); 
				
				
			}
			
			CreateCanvas(Form("eff_R1-%d_y%d",den,iy),"",Cx,Cy); 
			CName[Form("eff_R1-%d_y%d",den,iy)]->cd();
			
			hName[Form("effEp_r1%d_y%d",den,iy)]->SetLineColor(kRed); 
			hName[Form("effEm_r1%d_y%d",den,iy)]->SetLineColor(kGreen); 
		
			TString BT;
			if(iy==0) BT=Form("|y| < 0.6, #epsilon_{#mu#mu}^{1S}/#epsilon_{#mu#mu}^{%dS}",den); 
			if(iy==1) BT=Form("0.6 < |y| < 1.2, #epsilon_{#mu#mu}^{1S}/#epsilon_{#mu#mu}^{%dS}",den);
			if(iy==2) BT=Form("|y| < 1.2, #epsilon_{#mu#mu}^{1S}/#epsilon_{#mu#mu}^{%dS}",den);
			TLatex L(0.2,0.45,BT); 
			L.SetNDC(kTRUE);
			
			hName[Form("eff_r1%d_y%d",den,iy)]->SetMinimum(0.9);
			hName[Form("eff_r1%d_y%d",den,iy)]->SetMaximum(1.03);

			hName[Form("eff_r1%d_y%d",den,iy)]->Draw(); 
			hName[Form("effEp_r1%d_y%d",den,iy)]->Draw("histo same"); 
			hName[Form("effEm_r1%d_y%d",den,iy)]->Draw("histo same"); 
			L.DrawClone(); 
			draw_header();
			TLegend Leg(0.2,0.2,0.6,0.4);
			Leg.SetFillColor(10);
			Leg.SetLineColor(10); 
			Leg.SetBorderSize(0);
			
			Leg.AddEntry(eff_name(iy, den, ""),Form("#epsilon_{#mu#mu}^{1S}/#epsilon_{#mu#mu}^{%dS}",den),"L"); 
			Leg.AddEntry(hName[Form("effEp_r1%d_y%d",den,iy)],Form("#epsilon_{#mu#mu}^{1S}/#epsilon_{#mu#mu}^{%dS} +1#sigma",den),"L"); 
			Leg.AddEntry(hName[Form("effEm_r1%d_y%d",den,iy)],Form("#epsilon_{#mu#mu}^{1S}/#epsilon_{#mu#mu}^{%dS} +1#sigma",den),"L"); 

			Leg.DrawClone("same");

		}
		
	}
	CreateCanvas("eff_mass_dependence","",Cx,Cy); 
	CName["eff_mass_dependence"]->cd(); 
	int iy=0; int ipt=0; 
	
	TLegend Leg(0.26,0.2,0.65,0.4);
	Leg.SetFillColor(10);
	Leg.SetLineColor(10); 
	Leg.SetBorderSize(0); 
	
	prName[eff_name_mass(iy,ipt,"")]->SetMinimum(0.65);
	prName[eff_name_mass(iy,ipt,"")]->SetMaximum(0.9); 
	prName[eff_name_mass(iy,ipt,"")]->Draw(); 
	Leg.AddEntry(prName[eff_name_mass(iy,ipt,"")],BT(iy,ipt),"LEP");
	
	iy=0; ipt=5; 
	prName[eff_name_mass(iy,ipt,"")]->SetMarkerColor(kRed); 
	prName[eff_name_mass(iy,ipt,"")]->Draw("same"); 
	Leg.AddEntry(prName[eff_name_mass(iy,ipt,"")],BT(iy,ipt),"LEP");
	
	iy=0; ipt=17; 
	prName[eff_name_mass(iy,ipt,"")]->SetMarkerColor(kBlue); 
	prName[eff_name_mass(iy,ipt,"")]->Draw("same"); 
	Leg.AddEntry(prName[eff_name_mass(iy,ipt,"")],BT(iy,ipt),"LEP");
	

	iy=1; ipt=0; 
	prName[eff_name_mass(iy,ipt,"")]->SetMarkerStyle(23); 
	prName[eff_name_mass(iy,ipt,"")]->Draw("same"); 
	Leg.AddEntry(prName[eff_name_mass(iy,ipt,"")],BT(iy,ipt),"LEP");
	
	iy=1; ipt=5; 
	prName[eff_name_mass(iy,ipt,"")]->SetMarkerStyle(23);
	prName[eff_name_mass(iy,ipt,"")]->SetMarkerColor(kRed); 
	prName[eff_name_mass(iy,ipt,"")]->Draw("same"); 
	Leg.AddEntry(prName[eff_name_mass(iy,ipt,"")],BT(iy,ipt),"LEP");
	
	iy=1; ipt=17; 
	prName[eff_name_mass(iy,ipt,"")]->SetMarkerStyle(23);
	prName[eff_name_mass(iy,ipt,"")]->SetMarkerColor(kBlue); 
	prName[eff_name_mass(iy,ipt,"")]->Draw("same"); 
	Leg.AddEntry(prName[eff_name_mass(iy,ipt,"")],BT(iy,ipt),"LEP");
	
	Leg.DrawClone();
	draw_header();
	
	CreateCanvas("sg_eff_plot","",Cx,Cy);
	CName["sg_eff_plot"]->cd();
	
	hName[Form("sg_y%d_TH1_%dS_rebin",0,1)]->Draw(); 
	
	
}

void draw_muon_pt(){
	
	
	float use_bins[]={10,20,60};
	
	for (int iy=0; iy<fNy; iy++) {
		TLegend *L = new TLegend(0.4,0.5,0.8,0.8); 
		L->SetFillColor(10);
		CreateCanvas(Form("muonpt_overlay_y%d",iy),"",Cx,Cy); 
		CName[Form("muonpt_overlay_y%d",iy)]->cd(); 
				
		for (int iB=0; iB<3; iB++) {
			int ipt=hName["bin_width_pt"]->FindBin(use_bins[iB])-1;
			TH1F *y=(TH1F*)hName[h_muon(iy,ipt)]->Clone(h_muon(iy,ipt)+"_overlay"); 
			if(iB==0) {
				y->GetYaxis()->SetTitle("Arb. Units"); 
				y->SetLineColor(kRed);
				y->SetLineStyle(kDashed);
			}
			if(iB==1) {
				y->SetLineColor(kBlack);
				y->SetLineStyle(kDotted); 
				y->SetLineWidth(1.5);
			}
			
			if(iB==2){
				y->SetLineColor(kBlue);
			}
			y->Scale(1./y->Integral()); 
			if(iB==0)y->DrawCopy("histo");
			else y->DrawCopy("histo same");
			L->AddEntry(y,BT(iy,ipt),"L"); 
		}//loop over bins...
		L->DrawClone("same"); 
		draw_header();
		delete L; 

	}
}

void CreateHistogram(const char* name,   const char* title,
					 const char* xTitle, const char* yTitle,
					 Int_t       nBinsX, Double_t    xLow, Double_t xUp)
{
	TH1F* h = new TH1F(name, title, nBinsX, xLow, xUp);
	
	h->GetXaxis()->SetTitle(xTitle);
	h->GetYaxis()->SetTitle(yTitle);
	
	h->Sumw2();
	
	hName[name] = h;
}


void CreateHistogram(const char* name,   const char* title,
					 const char* xTitle, const char* yTitle,
					 Int_t       nBinsX, const Float_t* xBins)
{
	TH1F* h = new TH1F(name, title, nBinsX, xBins);
	
	h->GetXaxis()->SetTitle(xTitle);
	h->GetYaxis()->SetTitle(yTitle);
	
	h->Sumw2();
	
	hName[name] = h;
}

void CreateCanvas(TString Name,TString Title, int x, int y ){
	TCanvas *createC = new TCanvas( Name,  Title,x,y);
	CName[Name]=createC;
}