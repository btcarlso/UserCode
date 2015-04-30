//
//  plots.c
//  
//
//  Created by Benjamin Carlson on 11/24/14.
//
//

#include <stdio.h>
#include "plots.h"


void plots(){
    gROOT->SetBatch();
    open_file("ttFullLept");
    open_file("dy");
    open_file("singleTop");
    open_file("diboson");
    open_file("singleMu");
    //open_file("allMC");
    open_file(Form("RPV_LQD221_M%d",M));
    rescale();
    combine_all();
    ratio();

    plot_distribution();
    
    TFile *output = new TFile("summary_plots.root","RECREATE");
    
    for (std::map<TString,TCanvas*>::iterator it=CName.begin(); it!=CName.end(); it++) {
        CName[it->first]->Write();
        CName[it->first]->Print("plots/"+it->first+".jpg");
    }
    
    
}

void CreateGraph(TH1F *h,TString graphName){
    TGraphAsymmErrors *gr = new TGraphAsymmErrors(h);
    gr->SetName(graphName);
    gr->SetFillStyle(3244);
    gr->SetFillColor(38);
    
    grName[gr->GetName()]=gr;
}

double poisson(double n){
    
    const double alpha = 1 - 0.6827;
    double l = 0.5*TMath::ChisquareQuantile(alpha/2,2*n);
    double h = (n==0) ? ( 0.5*TMath::ChisquareQuantile(1-alpha,2*(n+1)) ) : ( 0.5*TMath::ChisquareQuantile(1-alpha/2,2*(n+1)) );
    
    return TMath::Max(h-n,n-l);
    
    //cout << "n: " << n << " +/- " << h << " " << l << endl;
    //if(n>0)cout << "n: " << n << " +/- [%] " << 100*(h-n)/n << " " << 100*(n-l)/n << endl;
}
void combine_all(){
    std::vector<TString> variable_names;
    for(int nJ=2; nJ<=7; nJ++){
        variable_names.push_back(Form("h_nJets%d_MuE_st",nJ));
        variable_names.push_back(Form("h_nJets%d_st",nJ));

    }
    for(int nB=0; nB<=2; nB++){
        variable_names.push_back(Form("h_nJets_%dbtags",nB));
        variable_names.push_back(Form("h_nJets_%dbtags_lowM",nB));
        variable_names.push_back(Form("h_nJets_%dbtags_MuE",nB));
    }
    
    std::vector<TString> names;
    names.push_back("_diboson");
    names.push_back("_dy");
    names.push_back("_singleTop");
    names.push_back("_ttFullLept");
    for(int i=0; i<variable_names.size(); i++){
        combine_histograms(variable_names.at(i),names,"_allMC");
    }

}

void ratio(){
 
    
    for(int nJ=2; nJ<=7; nJ++){
        cout << "nJ: " << nJ << endl;
       // combine_histograms(Form("h_nJets%d_st",nJ),names,"_allMC");
        CreateGraph(hName[Form("h_nJets%d_st_allMC",nJ)],Form("h_nJets%d_st_sys",nJ));
        for(int ist=1; ist<=hName[Form("h_nJets%d_MuE_st_ttFullLept",nJ)]->GetNbinsX(); ist++){
            double MCtt=hName[Form("h_nJets%d_MuE_st_ttFullLept",nJ)]->GetBinContent(ist);
            double MCdy=hName[Form("h_nJets%d_MuE_st_dy",nJ)]->GetBinContent(ist);
            double MCt=hName[Form("h_nJets%d_MuE_st_singleTop",nJ)]->GetBinContent(ist);
            double MCdb=hName[Form("h_nJets%d_MuE_st_diboson",nJ)]->GetBinContent(ist);
            double Ndata=hName[Form("h_nJets%d_MuE_st_singleMu",nJ)]->GetBinContent(ist);
            double R=(Ndata-MCdy-MCt-MCdb)/MCtt;
            double sigmaR=R*poisson(Ndata)/Ndata;
            
            double y=grName[Form("h_nJets%d_st_sys",nJ)]->GetY()[ist-1];
            /*
            grName[Form("h_nJets%d_st_sys",nJ)]->GetEXhigh()[ist-1]=hName[Form("h_nJets%d_st_ttFullLept",nJ)]->GetBinWidth(ist)/2;
            grName[Form("h_nJets%d_st_sys",nJ)]->GetEXlow()[ist-1]=hName[Form("h_nJets%d_st_ttFullLept",nJ)]->GetBinWidth(ist)/2;
            grName[Form("h_nJets%d_st_sys",nJ)]->GetEYhigh()[ist-1]=y*poisson(Ndata)/Ndata;
            grName[Form("h_nJets%d_st_sys",nJ)]->GetEYlow()[ist-1]=y*poisson(Ndata)/Ndata;
             */
            //cout << "Nj: " << nJ << " st: " << hName[Form("h_nJets%d_MuE_st_ttFullLept",nJ)]->GetBinCenter(ist);
            //cout << " sys: " << y*poisson(Ndata)/Ndata << endl;
            if(Ndata>0 && MCtt>0){
                double BW=hName[Form("h_nJets%d_MuE_st_ttFullLept",nJ)]->GetBinWidth(ist)/2;
                cout << hName[Form("h_nJets%d_MuE_st_ttFullLept",nJ)]->GetBinCenter(ist)-BW;
                cout << "--" << hName[Form("h_nJets%d_MuE_st_ttFullLept",nJ)]->GetBinCenter(ist)+BW;
                cout << Form(" %.2f $\\pm$ %.2f ",R,sigmaR) << endl;
                //cout << " R: " << R << " +/- " << sigmaR << endl;
            }
            int ICntrlMax=hName[Form("h_nJets%d_MuE_st_ttFullLept",nJ)]->GetNbinsX();
            int ISigMax=hName[Form("h_nJets%d_st_ttFullLept",nJ)]->GetNbinsX();
           // for(int ist=ICntrlMax; ist<=ISigMax; ist++){
            for(int ist=1; ist<=ISigMax; ist++){
               // double Ndata=hName[Form("h_nJets%d_st_singleMu",nJ)]->GetBinContent(ICntrlMax);
                double y=grName[Form("h_nJets%d_st_sys",nJ)]->GetY()[ist-1];
                double unc=1;
                if(nJ==3) unc=0.03;
                if(nJ==4) unc=0.04;
                if(nJ==5) unc=0.05;
                if(nJ==6) unc=0.20;
                if(nJ==7) unc=0.46;

                grName[Form("h_nJets%d_st_sys",nJ)]->GetEXhigh()[ist-1]=hName[Form("h_nJets%d_st_ttFullLept",nJ)]->GetBinWidth(ist)/2;
                grName[Form("h_nJets%d_st_sys",nJ)]->GetEXlow()[ist-1]=hName[Form("h_nJets%d_st_ttFullLept",nJ)]->GetBinWidth(ist)/2;
                grName[Form("h_nJets%d_st_sys",nJ)]->GetEYhigh()[ist-1]=y*unc;
                grName[Form("h_nJets%d_st_sys",nJ)]->GetEYlow()[ist-1]=y*unc;
            }
        }
    }
    float dy[]={1,1,1.24,1.24,1.35,1.16,1.74,1.74};
    
    CreateGraph(hName["h_nJets_0btags_lowM_allMC"],"h_nJets_0btags_lowM_sys");
    CreateGraph(hName["h_nJets_1btags_lowM_allMC"],"h_nJets_1btags_lowM_sys");

    cout << "fill graph: " << endl;
    for(int nB=0; nB<=1; nB++){
    for(int nJ=2; nJ<=7; nJ++){
        cout << "nJ : " << nJ << endl;
        double BW=0.5;
        double y=grName[Form("h_nJets_%dbtags_lowM_sys",nB)]->GetY()[nJ-1];
        grName[Form("h_nJets_%dbtags_lowM_sys",nB)]->GetEYhigh()[nJ-1]=y*(dy[nJ]-1);
        grName[Form("h_nJets_%dbtags_lowM_sys",nB)]->GetEYlow()[nJ-1]=y*(dy[nJ]-1);
        grName[Form("h_nJets_%dbtags_lowM_sys",nB)]->GetEXhigh()[nJ-1]=0.5;
        grName[Form("h_nJets_%dbtags_lowM_sys",nB)]->GetEXlow()[nJ-1]=0.5;
    }
    }
    
    
}

void rescale(){
    float tt[]={1,1,1,1.01,1.03,1.05,1.01,0.81};
    for(int nJ=3; nJ<=7; nJ++){
        //cout << "tt: " << tt[nJ] << endl;
        hName[Form("h_nJets%d_st_ttFullLept",nJ)]->Scale(tt[nJ]);
    }
    float dy[]={1,1,1.24,1.24,1.35,1.16,1.74,1.74};
    for(int nJ=3; nJ<=7; nJ++){
        //cout << "dy: " << dy[nJ] << endl;
        hName[Form("h_nJets%d_st_dy",nJ)]->Scale(dy[nJ]);
    }
    for(int nJ=2; nJ<=7; nJ++){
        //cout << "dy: " << dy[nJ] << endl;
        double Ndy=hName["h_nJets_0btags_lowM_dy"]->GetBinContent(nJ);
        Ndy*=dy[nJ];
        hName["h_nJets_0btags_lowM_dy"]->SetBinContent(nJ,Ndy);
    }
    
    for(int nJ=2; nJ<=7; nJ++){
        //cout << "dy: " << dy[nJ] << endl;
        double Ndy=hName["h_nJets_1btags_lowM_dy"]->GetBinContent(nJ);
        Ndy*=dy[nJ];
        hName["h_nJets_1btags_lowM_dy"]->SetBinContent(nJ,Ndy);
    }
  
}

void plot_distribution(){
    vector<TString> names;
    names.push_back("_diboson");
    names.push_back("_dy");
    names.push_back("_singleTop");
    names.push_back("_ttFullLept");

    std::vector<TString> signal_names;
    signal_names.push_back(Form("_RPV_LQD221_M%d",M));
    
	vector<TString> legend_names;
    legend_names.push_back("Dibsoon");
	legend_names.push_back("Drell-Yan");
	legend_names.push_back("Single t");
	legend_names.push_back("t#bar{t}");
	legend_names.push_back(Form("M_{#tilde{t}} = %d GeV",M));
    
    right_=false;
    for(int nB=0; nB<=2; nB++)fill_stack(Form("h_nJets_%dbtags",nB),"_singleMu",signal_names,names,legend_names,"_allMC");
    for(int nB=0; nB<=2; nB++){
       if(nB==0 || nB==1)plotSystematic=true;
      fill_stack(Form("h_nJets_%dbtags_lowM",nB),"_singleMu",signal_names,names,legend_names,"_allMC");
      plotSystematic=false;
    }
    right_=true;
    
    for(int nB=0; nB<=2; nB++)fill_stack(Form("h_nJets_%dbtags_MuE",nB),"_singleMu",signal_names,names,legend_names,"_allMC");

    for(int nJ=2; nJ<=7; nJ++)fill_stack(Form("h_nJets%d_MuE_st",nJ),"_singleMu",signal_names,names,legend_names,"_allMC");
    plotSystematic=true;
    for(int nJ=2; nJ<=7; nJ++){
        if(nJ==5||nJ==6)right_=false;
        fill_stack(Form("h_nJets%d_st",nJ),"_singleMu",signal_names,names,legend_names,"_allMC");
        right_=true;
    }
}

void fill_stack(TString variable_name, TString data_name, std::vector<TString> signal_names, std::vector<TString> names, std::vector<TString> legend_names, TString MC_name){
    //gStyle->SetErrorX(0);
	bool print=false;
	if(print) cout << variable_name << "  " << data_name << endl;
    TString canvas_name=variable_name+"_stack_canvas";
    int Cx=800;
    int Cy=800;
    CreateCanvas(canvas_name,"",Cx,Cy);
    TString stack_name=variable_name+"_stack";
    CreateStack(stack_name,"");
	
    int Nbkg=names.size();
  
    //8,30,46,
    vector<int> colors;
    for(int i=0; i<Nbkg; i++){
        if(names.at(i)=="_dy")colors.push_back(593);
        if(names.at(i)=="_ttFullLept")colors.push_back(625);
        if(names.at(i).Contains("singleTop"))colors.push_back(411);
        if(names.at(i).Contains("diboson"))colors.push_back(38);
    }
	
	
	if(Nbkg+signal_names.size()>legend_names.size() ){
		cout << "Legend names need more entries. " << endl;
		return;
	}
	
	
    
    //colors.push_back(7);
    if(colors.size()<Nbkg) {
        cout << "please specify colors" << endl;
        return;
    }
    TString labelTxt="#mu#mu";
    for(int nJ=2; nJ<=7; nJ++){
        if(variable_name.Contains(Form("nJets%d",nJ)) && nJ<7) labelTxt+=Form(" N_{jets} = %d",nJ);
        if(variable_name.Contains(Form("nJets%d",nJ)) && nJ==7) labelTxt+=Form(" N_{jets} #geq %d",nJ);

    }
    
    TLatex binlabel(0.6,0.85,labelTxt);
    
    
    TLegend *L;
    if(right_)L= new TLegend(0.18,0.51,0.41,0.81);
    else L = new TLegend(0.6,0.55,0.89,0.85);
	L->SetFillColor(10);
	L->SetLineColor(10);
	L->SetLineWidth(0);
    L->SetTextSize(0.04);
    
    L->AddEntry(hName[variable_name+data_name], "Data","p");
    
    for(int ibkg=0; ibkg<Nbkg; ibkg++){
        TString histName=variable_name+names.at(ibkg);
        hName[histName]->SetFillColor(colors.at(ibkg));
        stackName[stack_name]->Add(hName[histName]);
        
    }
    
    for(int ibkg=Nbkg-1; ibkg>=0; ibkg--){
        TString histName=variable_name+names.at(ibkg);
        L->AddEntry(hName[histName], legend_names.at(ibkg),"F");
    }
    
    /*
	TString eventspad_name="eventspad_"+variable_name;
	TString ratiopad_name="ratiopad_"+variable_name;
	
	TPad *events_pad = new TPad(eventspad_name,"Events",0.0,0.3,1,1);
	TPad *ratio_pad = new TPad(ratiopad_name,"Ratio",0,0.,1,0.3);
	
	events_pad->SetTopMargin(0.1);
	events_pad->SetBottomMargin(0.0);
	ratio_pad->SetTopMargin(0.0);
	ratio_pad->SetBottomMargin(0.4);
	
	CName[canvas_name]->cd();
	events_pad->Draw();
	ratio_pad->Draw();
	
	events_pad->cd();
     */
    CName[canvas_name]->cd();
    gPad->SetLogy();
    
    hName[variable_name+data_name]->GetYaxis()->SetTitleOffset(0.65);
    hName[variable_name+data_name]->GetYaxis()->SetTitleSize(0.08);
    hName[variable_name+data_name]->GetYaxis()->SetLabelSize(0.05);
    //hName[variable_name+data_name]->GetXaxis()->SetLabelSize(0);
    //hName[variable_name+data_name]->GetXaxis()->SetTitleSize(0);
    
    hName[variable_name+data_name]->SetMarkerSize(0.75);
    hName[variable_name+data_name]->SetMarkerStyle(20);
    hName[variable_name+data_name]->SetMinimum(0.1);
    double max=hName[variable_name+MC_name]->GetMaximum();
    for(int isig=0; isig<signal_names.size();isig++){
        TString sigName=variable_name+signal_names.at(isig);
        max = TMath::Max(max,hName[sigName]->GetMaximum());
    }
    
    hName[variable_name+data_name]->SetMaximum(100*max);
    
    hName[variable_name+data_name]->DrawCopy("E1");
   // stackName[stack_name]->Draw("histo ");
   // stackName[stack_name]->SetMaximum(75);
   // stackName[stack_name]->SetMinimum(0.1);
    stackName[stack_name]->DrawClone("histo same");
    if(plotSystematic)grName[variable_name+"_sys"]->DrawClone("E2 same");
    hName[variable_name+data_name]->DrawCopy("E1 same");// draw data on top of MC and MC on top of data, so that data will always be visible
	
    for(int isig=0; isig<signal_names.size();isig++){
        TString sigName=variable_name+signal_names.at(isig);
		if(print) cout << "signal drawing settings. " << endl;
        hName[sigName]->SetLineWidth(3);
        hName[sigName]->SetLineStyle(kDashed);
        hName[sigName]->SetLineColor(1+isig);
        //L->AddEntry((TObject*)0,"","");
        L->AddEntry(hName[sigName], legend_names.at(Nbkg+isig),"L");
		if(print) cout << "draw signal: " << isig << endl;
        hName[sigName]->DrawCopy("hist same");
    }
    L->DrawClone("same");
    binlabel.SetNDC(kTRUE);
    binlabel.DrawClone();
    gPad->RedrawAxis();
	if(print) cout << "set ratio stuff. " << endl;
    //draw_header();
    CMS_lumi(CName[canvas_name],2,10);
    
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



void CreateStack(TString Name,TString Title){
	THStack *hstack = new THStack(Name,Title);
	stackName[Name]=hstack;
}

void CreateCanvas(TString Name,TString Title, int x, int y ){
	TCanvas *createC = new TCanvas(Name, Title,x,y);
	CName[Name]=createC;
}

void read_histograms(TString file){
    
    TIter nextkey(fName[file]->GetListOfKeys());
    TKey *key;
    
    while((key=(TKey*)nextkey())){
        TString className=key->ReadObj()->ClassName();
        TString keyName=key->GetName();
        //cout << keyName << endl;
        //if(keyName.Contains("st"))cout << "histogram: " << keyName << endl;
        if(className=="TH1F" ){
            TH1F *h=(TH1F*)key->ReadObj();
            h->SetStats(kFALSE);
            h->SetBinErrorOption(TH1::kPoisson);
            hName[h->GetName()+file]=h;
            TString NN=h->GetName()+file;
            //if(NN.Contains("st")==0) continue;
        
        }
    }
    
}



void open_file(TString name){
	cout << "open file:" << name << endl;
    
	TString file_name="output_file_"+name+".root";
    TFile *f = new TFile(file_name,"READ");
	if(f->IsOpen()!=1) {
		cout << "File: " << file_name << " Failed to open!" << endl;
		return;
	}
	name="_"+name;
	fName[name]=f;
	read_histograms(name);
	
}

void combine_histograms(TString variable_name,std::vector<TString> names, TString outname){
    int N=names.size();
    outname=variable_name+outname;
	//cout << "combine histograms: " << variable_name << endl;
	//for(int i=0; i<N; i++) cout << names.at(i) << endl;
	//cout << "clone histogram: " << endl;
    clone_histogram(variable_name+names.at(0),outname);
	//  cout << "adding histograms: " << endl;
    for(int i=1; i<N;i++){
		//        cout << "adding histogram: " << names.at(i) << endl;
        hName[outname]->Add(hName[variable_name+names.at(i)]);
    }
    
}

void clone_histogram(TString name1, TString clone_name){
    //cout << "clone histogram: " << name1 << " " << clone_name << endl;
    TH1F *h=(TH1F*)hName[name1]->Clone(clone_name);
    //h->Sumw2();
    h->SetName(clone_name);
    hName[clone_name]=h;
}

