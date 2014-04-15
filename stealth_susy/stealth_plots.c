//
//  stealth_plots.c
//  
//
//  Created by Benjamin Carlson on 4/10/14.
//
//

#include <stealth_plots.h>

void stealth_plots(){
    gROOT->SetBatch();
 	open_file("dy");
	open_file("ttbar");
	open_file("singleTop");
	open_file("diboson");
	open_file("singleMu");
    open_file("stealth_300_200");
    open_file("stealth_500_400");


    TString scales[]={"_SF","_TopCor"};
    
    int Nscales=sizeof(scales)/sizeof(TString);
    
    std::vector<TString> signal_names;
    signal_names.push_back("_stealth_300_200");
    signal_names.push_back("_stealth_500_400");

    for(int ib=0; ib<=2; ib++){
        for(int i=0; i<Nscales; i++){
            TString variable_name=Form("EventCategories_1Mu_1El_%dbtag_tight",ib)+scales[i];
            vector<TString> names;
            names.push_back("_dy");
            names.push_back("_diboson");
            names.push_back("_singleTop");
            names.push_back("_ttbar");
            
            combine_histograms(variable_name,names,"_allMC");
            blind_categoryplots(ib, scales[i]);
            compute_ratio(variable_name, "_allMC");
            fill_stack_eventcategories(variable_name,"_singleMu",signal_names,names);
        }
        
    }
    
    
    
    TFile *output_file = new TFile("stealth_plots.root","RECREATE");
    output_file->cd();
    
	for (std::map<TString,TCanvas*>::iterator it=CName.begin(); it!=CName.end(); it++) {
		CName[it->first]->Write();
        TString name="/uscms/home/btcarlso/AN-13-083/tdr2/notes/AN-13-083/trunk/";
        name=name+it->first+".pdf";
        CName[it->first]->Print(name);
    }
    output_file->Close();
    
    
}

void blind_categoryplots(int nb, TString scale){
    int nJetmax=8;
    TString name=Form("EventCategories_1Mu_1El_%dbtag_tight",nb)+scale;
    name=name+"_singleMu";
    int nJ=1;
    for(int ibin=1; ibin<=hName[name]->GetNbinsX(); ibin++){
        if(nb==0 && nJ>3)hName[name]->SetBinContent(ibin,0);
        if(nb==1 && nJ>6)hName[name]->SetBinContent(ibin,0);
        nJ++;
        if(nJ>nJetmax)nJ=1;
    }

}

void compute_ratio(TString variable_name, TString MC_name){
    
    clone_histogram(variable_name+"_singleMu",variable_name+"_RatioDataMC");
    for (int i=1; i<=hName[variable_name+"_RatioDataMC"]->GetNbinsX(); i++){
        double R=hName[variable_name+"_RatioDataMC"]->GetBinContent(i);
        if(R<0.01) hName[variable_name+"_RatioDataMC"]->SetBinContent(i,0);
        if(R<0.01) hName[variable_name+"_RatioDataMC"]->SetBinError(i,0);

    }
    hName[variable_name+"_RatioDataMC"]->GetYaxis()->SetTitle("Data/MC");
    hName[variable_name+"_RatioDataMC"]->SetAxisRange(0.5,2.0,"Y");
    hName[variable_name+"_RatioDataMC"]->SetMarkerStyle(20);
    hName[variable_name+"_RatioDataMC"]->SetMarkerSize(0.5);
    hName[variable_name+"_RatioDataMC"]->Divide(hName[variable_name+MC_name]);
}

void fill_stack_eventcategories(TString variable_name, TString data_name, std::vector<TString> signal_names, std::vector<TString> names){
    TString canvas_name=variable_name+"_stack_canvas";
    CreateCanvas(canvas_name,"",900,600);
    TString stack_name=variable_name+"_stack";
    CreateStack(stack_name,"");

    int Nbkg=names.size();
    
    vector<int> colors;
    colors.push_back(4);
    colors.push_back(9);
    colors.push_back(5);
    colors.push_back(419);
    
    colors.push_back(7);
    if(colors.size()<Nbkg) {
        cout << "please specify colors" << endl;
        return;
    }
    
    TLegend *L = new TLegend(0.6,0.57,0.89,0.89);
	L->SetFillColor(10);
	L->SetLineColor(10);
	L->SetLineWidth(0);
    
    for(int ibkg=0; ibkg<Nbkg; ibkg++){
        TString histName=variable_name+names.at(ibkg);
        hName[histName]->SetFillColor(colors.at(ibkg));
        stackName[stack_name]->Add(hName[histName]);
        L->AddEntry(hName[histName], names.at(ibkg));
    }
    
    
	TString eventspad_name="eventspad_"+variable_name;
	TString ratiopad_name="ratiopad_"+variable_name;
	
	TPad *events_pad = new TPad(eventspad_name,"Events",0.0,0.3,1,1);
	TPad *ratio_pad = new TPad(ratiopad_name,"Ratio",0,0.,1,0.3);
	
	events_pad->SetTopMargin(0.1);
	events_pad->SetBottomMargin(0.05);
	ratio_pad->SetTopMargin(0.05);
	ratio_pad->SetBottomMargin(0.3);
	
	CName[canvas_name]->cd();
	events_pad->Draw();
	ratio_pad->Draw();
	
	events_pad->cd();
    gPad->SetLogy();
    
    hName[variable_name+data_name]->GetYaxis()->SetTitleOffset(0.65);
    hName[variable_name+data_name]->GetYaxis()->SetTitleSize(0.08);
    hName[variable_name+data_name]->GetYaxis()->SetLabelSize(0.05);
    hName[variable_name+data_name]->GetXaxis()->SetLabelSize(0);
    hName[variable_name+data_name]->GetXaxis()->SetTitleSize(0);
    
    hName[variable_name+data_name]->SetMarkerSize(0.5);
    hName[variable_name+data_name]->SetMarkerStyle(20);
    hName[variable_name+data_name]->DrawCopy("E1");
    stackName[stack_name]->DrawClone("histo same");
    hName[variable_name+data_name]->DrawCopy("E1 same");// draw data on top of MC and MC on top of data, so that data will always be visible

    
    TLatex txt;
    txt.SetNDC(kTRUE);
    txt.DrawLatex(0.2,0.85,"S_{T} > 300 GeV");
    txt.DrawLatex(0.4,0.75,"S_{T} > 500 GeV");
    txt.DrawLatex(0.65,0.42,"S_{T} > 1000 GeV");

    
    for(int isig=0; isig<signal_names.size();isig++){
        TString sigName=variable_name+signal_names.at(isig);
        hName[sigName]->SetLineWidth(3);
        hName[sigName]->SetLineStyle(kDashed);
        hName[sigName]->SetLineColor(1+isig);
        L->AddEntry(hName[sigName], signal_names.at(isig));
        hName[sigName]->DrawCopy("hist same");
    }
    L->DrawClone("same");
    ratio_pad->cd();
    TString dataMC=variable_name+"_RatioDataMC";

    hName[dataMC]->GetYaxis()->SetTitleOffset(0.45);
	hName[dataMC]->GetYaxis()->SetTitleSize(0.15);
	hName[dataMC]->GetYaxis()->SetTitleOffset(0.45);
	hName[dataMC]->GetYaxis()->SetLabelSize(0.08);
	
	
	hName[dataMC]->GetXaxis()->SetTitleSize(0.15);
	hName[dataMC]->GetXaxis()->SetLabelSize(0.12);
	hName[dataMC]->GetXaxis()->SetTitleSize(0.15);
    
    hName[dataMC]->DrawCopy("E1");
    
}

void open_file(TString name){
	cout << "open file:" << name << endl;
    
	TString file_name="/uscms_data/d3/btcarlso/TEST_DIR/CMSSW_5_3_8_patch3/src/SusyAnalysis/SusyNtuplizer/macro/output_file_"+name+".root";
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
    clone_histogram(variable_name+names.at(0),outname);
    
    for(int i=1; i<N;i++){
        hName[outname]->Add(hName[variable_name+names.at(i)]);
    }
    
}

void clone_histogram(TString name1, TString clone_name){
    TH1F *h=(TH1F*)hName[name1]->Clone(clone_name);
    //h->Sumw2();
    hName[h->GetName()]=h;
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
void CreateHistogram(const char* name,   const char* title,
					 const char* xTitle, const char* yTitle,
					 Int_t nBinsX, Double_t xLow, Double_t xUp,
					 Int_t nBinsY,Double_t yLow, Double_t yUp)
{
	TH2F* h = new TH2F(name, title, nBinsX, xLow,xUp,nBinsY, yLow,yUp);
	
	h->GetXaxis()->SetTitle(xTitle);
	h->GetYaxis()->SetTitle(yTitle);
	
	h->Sumw2();
	
	hName2D[name] = h;
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
        if(className=="TH1F" && (keyName.Contains("EventCategories") )){
            TH1F *h=(TH1F*)key->ReadObj();
            h->SetStats(kFALSE);
            hName[h->GetName()+file]=h;
        }
    }
    
}

