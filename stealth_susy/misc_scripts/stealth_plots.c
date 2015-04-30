//
//  stealth_plots.c
//  
//
//  Created by Benjamin Carlson on 4/10/14.
//
//

#include <stealth_plots.h>

void stealth_plots(){
    bool printPDF=true;
    gROOT->SetBatch();
    
 	open_file("dy");
	open_file("ttbar");
	
	open_file("ttJets_scaledown");
	open_file("ttJets_scaleup");

	open_file("ttJets_matchingdown");
	open_file("ttJets_matchingup");
	
	open_file("singleTop");
	open_file("diboson");
	open_file("singleMu");
    open_file("stealth_300_200");
	open_file("stealth_400_200");
    open_file("stealth_500_300");
	open_file("stealth_600_300");
	open_file("stealth_700_400");
    open_file("stealth_800_400");
    open_file("stealth_900_500");


//    open_file("stealth_600_300_compressed_ptmu20_pte20");
//	open_file("stealth_800_300");

	//open_file("UDD300"); 
	
	
	//  plot_jet_distributions();

    plot_distributions();
    gStyle->SetErrorX(0.);
    ErrorX=false;
    sum_histograms();
    compute_QCD(0);
    compute_QCD(1);

    compute_DY(0);
    compute_DY(1);

    plot_categories();
    plot_categories_expected();
	draw_correction(0,0);
    draw_correction(0,1);
	draw_correction(0,2);
    draw_correction(0);
    shape_difference(0);
    table(0,0);//table 0 btags, st0
    table(0,1);
    table(0,2);
    //2 btag table
    table_simulation(2,0);
    table_simulation(2,1);
    table_simulation(2,2);

    plot_MC_exp(0);
    acceptance_table();
	/*
    for(int nJ=4; nJ<=nJetmax; nJ++)predict_systematic(0, nJ,""); // predict systematic for exclusive bins
    predict_systematic(0, 4,"_Inclusive");

    plot_zbi(0, 4, "_Inclusive", "_stealth_300_200");
    plot_zbi(0, 4, "_Inclusive", "_stealth_400_200");
    plot_zbi(0, 4, "_Inclusive", "_stealth_500_300");
    plot_zbi(0, 4, "_Inclusive", "_stealth_600_300");
    plot_zbi(0, 4, "_Inclusive", "_stealth_700_400");

    for(int nJ=4; nJ<=nJetmax; nJ++){
        plot_zbi(0, nJ, "", "_stealth_300_200");
        plot_zbi(0, nJ, "", "_stealth_400_200");
        plot_zbi(0, nJ, "", "_stealth_500_300");
        plot_zbi(0, nJ, "", "_stealth_600_300");
        plot_zbi(0, nJ, "", "_stealth_700_400");
    }
    
    summarize_zbi(0,"");
    
	for(int nJets=4; nJets<=nJetmax; nJets++){
		for(int nbtag=0; nbtag<=2; nbtag++){
			TString variable="stInclusive_Inclusive"; 
			model_independent_datacard(variable,nJets,nbtag); 
		}
	}
	*/
   // plot_mass();
    //print_fractions();
    
    TFile *output_file = new TFile("stealth_plots.root","RECREATE");
    output_file->cd();
    
	for (std::map<TString,TCanvas*>::iterator it=CName.begin(); it!=CName.end(); it++) {
		//CName[it->first]->Write();
        //TString name="/uscms/home/btcarlso/AN-13-083/tdr2/notes/AN-13-083/trunk/figures/";
        
        TString name="/uscms_data/d3/btcarlso/figtmp/";
        name=name+it->first+".pdf";
       if(printPDF)CName[it->first]->Print(name);
        CName[it->first]->Write();
    }
    /*
    for (std::map<TString,TGraphAsymmErrors*>::iterator it=grName.begin(); it!=grName.end(); it++) {
//		grName[it->first]->Write();
    }
     */
    output_file->Close();

    cout << "close program: " << endl;
    gSystem->Exit(1,1);
    
}

void compute_DY(int nbtag){
    TString name= Form("EventCategories_1Mu_1El_%dbtag_SF_dy",nbtag);
    TString nameDD= Form("EventCategories_1Mu_1El_%dbtag_SF_dy_datadriven",nbtag);

    clone_histogram(name,nameDD);
    int nJ=1;
    int ist=0;
    
    double R[]={0,0,1.03,1.096,1.11,1.23,1.65,1.65,1.65}; //0-8 jets
    //use overflow for >6 jets
    
    for(int ibin=1; ibin<=hName[name]->GetNbinsX();ibin++){
        double Ndy=hName[name]->GetBinContent(ibin);
        double NdyE=hName[name]->GetBinError(ibin);
        hName[nameDD]->SetBinContent(ibin, Ndy*R[nJ]);
        hName[nameDD]->SetBinError(ibin, (NdyE/Ndy)*Ndy*R[nJ]);

        nJ++;
        if(nJ>nJetmax){
            nJ=1;
            ist++;
        }
    }
    clone_histogram(nameDD,Form("EventCategories_1Mu_1El_%dbtag_TopCor_dy_datadriven",nbtag));
    
    name= Form("EventCategories_1Mu_1El_%dbtag_SF_diboson",nbtag);
    nameDD= Form("EventCategories_1Mu_1El_%dbtag_SF_diboson_datadriven",nbtag);
    
    clone_histogram(name,nameDD);
    nJ=1;
    ist=0;
    
    //use overflow for >6 jets
    
    for(int ibin=1; ibin<=hName[name]->GetNbinsX();ibin++){
        double Ndy=hName[name]->GetBinContent(ibin);//actually diboson here, but I'm too lazy to make Ndy->Ndiboson.
        double NdyE=hName[name]->GetBinError(ibin);
        hName[nameDD]->SetBinContent(ibin, Ndy*R[nJ]);
        hName[nameDD]->SetBinError(ibin, (NdyE/Ndy)*Ndy*R[nJ]);
        
        nJ++;
        if(nJ>nJetmax){
            nJ=1;
            ist++;
        }
    }
    
    
}

void compute_QCD(int nbtag){
    bool print=false;
    TString name= Form("EventCategories_SS_1Mu_1El_%dbtag_SF_QCD",nbtag);
    
    for(int ist=0; ist<3; ist++){
        CreateHistogram(Form("h_njets_%dbtag_QCD_st%d",nbtag,ist),"","N_{jets}","Events",nJetmax-1,1.5,nJetmax+0.5);
        for(int nJ=1; nJ<nJetmax-1; nJ++)hName[Form("h_njets_%dbtag_QCD_st%d",nbtag,ist)]->GetXaxis()->SetBinLabel(nJ,Form("%d", nJ+1));
        hName[Form("h_njets_%dbtag_QCD_st%d",nbtag,ist)]->GetXaxis()->SetBinLabel(nJetmax-1,Form("#geq%d", nJetmax));
    }
    
    for(int ist=0; ist<3; ist++){
        for(int nJ=2; nJ<=nJetmax; nJ++){
            int ibin=category_bin(nJ,ist);
            if(print)cout << "nJ: " << nJ << " ist: " << ist << " ibin: " << ibin << " bin content: " << hName[name]->GetBinContent(ibin) << endl;
            
            hName[Form("h_njets_%dbtag_QCD_st%d",nbtag,ist)]->SetBinContent(nJ-1,hName[name]->GetBinContent(ibin));
            hName[Form("h_njets_%dbtag_QCD_st%d",nbtag,ist)]->SetBinError(nJ-1,hName[name]->GetBinError(ibin));
        }
            
    }
    
    
    RooRealVar njets("njets","njets",1.5,nJetmax+0.5);
    RooDataHist dh_qcd("dh_qcd","dh_qcd",njets,Import(*hName[Form("h_njets_%dbtag_QCD_st0",nbtag)]));
    
    RooDataHist dh_qcd1("dh_qcd","dh_qcd",njets,Import(*hName["h_njets_QCD_st1"]));
    RooDataHist dh_qcd2("dh_qcd","dh_qcd",njets,Import(*hName["h_njets_QCD_st2"]));

    RooRealVar p1("p1","p1",0.01,-2,1);
    RooExponential expo("expo","expo",njets,p1);
    if(print)cout << "st>300 fit: " << endl;
    expo.fitTo(dh_qcd,SumW2Error(kTRUE));
    
    if(print)cout << "p1: " << p1.getVal() << " +/- " << p1.getError() << endl;
    
    TF1 *exp=new TF1("exp","expo",1.5,7.5);
    exp->SetParameters(0,p1.getVal());
    double Norm=exp->Integral(1.5,7.5);
   
    if(nbtag==0 && print){
        cout<< "st>300: " << p1.getVal() << " +/- " << p1.getError() << endl << endl;
        expo.fitTo(dh_qcd1,SumW2Error(kTRUE));
        cout<< "st>700: " << p1.getVal() << " +/- " << p1.getError() << endl << endl;
        expo.fitTo(dh_qcd2,SumW2Error(kTRUE));
        cout<< "st>1200: " << p1.getVal() << " +/- " << p1.getError() << endl << endl;
    }
    
    
	TString variable[]={"EventCategories","EventCategories_SS"};
	int Nvariable=sizeof(variable)/sizeof(TString);
    
    int nJ=1;
    int ist=0;
    
    double Ntot[]={0,0,0};
    
    for(int ibin=1; ibin<=hName[Form("EventCategories_SS_1Mu_1El_%dbtag_SF_QCD",nbtag)]->GetNbinsX();ibin++){
        Ntot[ist]+=hName[Form("EventCategories_SS_1Mu_1El_%dbtag_SF_QCD",nbtag)]->GetBinContent(ibin);
        
        nJ++;
        if(nJ>nJetmax){
            nJ=1;
            ist++;
        }
    }
    
    if(print)cout << "Ntot: " << Ntot[0] << " " << Ntot[1] << " " << Ntot[2] << endl;
    if(print) cout << "Norm: " << Norm << endl;
    if(print)cout << "fill QCD: "<< endl;
    
    TString scales[]={"_noSF","_SF","_SFPbc","_SFMbc","_SFPl","_SFMl","_TopCor"};
    int Nscales=sizeof(scales)/sizeof(TString);
    
    for(int iscale=0; iscale<Nscales;iscale++){
        for (int ivar=0; ivar< Nvariable; ivar++) {
            nJ=1;
            ist=0;
            TString variable_name=variable[ivar]+Form("_1Mu_1El_%dbtag",nbtag)+scales[iscale]+"_QCD";
            //cout << variable_name << endl;
            for(int ibin=1; ibin<=hName[variable_name]->GetNbinsX();ibin++){
                double Ni=Ntot[ist]*exp->Integral(nJ-0.5,nJ+0.5)/Norm;
                //cout << "nJ "<< nJ << " st: "<< ist << " " << Ni << endl;
                if(nJ>=2)hName[variable_name]->SetBinContent(ibin,Ni);
                
                nJ++;
                if(nJ>nJetmax){
                    nJ=1;
                    ist++;
                }
            }//loop over bins
        }
    }
    
    if(nbtag==0){
        for(int ist=0; ist<3; ist++){
            double A=TMath::Log(Ntot[ist]/Norm);
            if(print)cout << "Norm: " << Norm << endl;
            if(print)cout << "A: " << A << endl;
            exp->SetParameter(0,A);
            CreateCanvas(Form("SameSign_Valdiation_st%d",ist),"",600,600);
            CName[Form("SameSign_Valdiation_st%d",ist)]->cd();
            gPad->SetLogy();
            hName[Form("h_njets_%dbtag_QCD_st%d",nbtag,ist)]->SetMaximum(150);
            hName[Form("h_njets_%dbtag_QCD_st%d",nbtag,ist)]->Draw("");
            exp->DrawClone("same");
            CMS_lumi(CName[Form("SameSign_Valdiation_st%d",ist)],2,10);
            //draw_header();
            int st[]={300,700,1200};
            TLatex txt1(0.4,0.7,"SS e,#mu");
            txt1.SetNDC(kTRUE);
            TLatex txt(0.4,0.65,Form("%d b-tag, S_{T} > %d GeV",nbtag,st[ist]));
            txt.SetNDC(kTRUE);
            txt1.DrawClone();
            txt.DrawClone();
            
        }
    }
    delete exp;
}

void summarize_zbi(int nbsig,TString stType){
    TString signal="_stealth_600_300";
    
    ofstream output;
	output.open("ZBI_nJets.tex");
	output << "\\begin{tabular}{";
	int N=3;
	for(int i=0; i<=N; i++) output << "c";
	output << "}" << endl;
    
    output << "n-jets & $Z_{BI}$ & $S_{T}^{min}$" << "\\\\" << endl;

    for(int nJ=4; nJ<=nJetmax; nJ++){
        TString variable_name="stInclusive"+stType+Form("_nJets%d_1Mu_1El_%dbtag_SF",nJ,nbsig)+signal+"_ZBI";
        
        double stmin=hName[variable_name]->GetBinCenter(hName[variable_name]->GetMaximumBin());
        double zbi=hName[variable_name]->GetBinContent(hName[variable_name]->GetMaximumBin());
        
        output << nJ << " & ";
        output << zbi << " & ";
        output << stmin << "\\\\" << endl;
        
    }
    output << "\\end{tabular}" << endl;
    output.close();

    int nJ=5;

	output.open("ZBI_nJets5.tex");
	output << "\\begin{tabular}{";
    N=5;
	for(int i=0; i<=N; i++) output << "c";
	output << "}" << endl;
    
    output << "$M_{\\tilde{q}}$ & optimum: $Z_{BI}$ & $S_{T}^{min}$ & used: $Z_{BI}$ & $S_{T}^{min}$" << "\\\\" << endl;
    
    TString variable_name="stInclusive"+stType+Form("_nJets%d_1Mu_1El_%dbtag_SF",nJ,nbsig);
    output << print_zbi(300,variable_name,"_stealth_300_200",300) << endl;
    output << print_zbi(400,variable_name,"_stealth_400_200",700) << endl;
    output << print_zbi(500,variable_name,"_stealth_500_300",700) << endl;
    output << print_zbi(600,variable_name,"_stealth_600_300",1200) << endl;
    output << print_zbi(700,variable_name,"_stealth_700_400",1200) << endl;
    output << "\\end{tabular}";
    output.close();
}

TString print_zbi(double mass, TString variable_name, TString signal, double st_used){
    variable_name=variable_name+signal+"_ZBI";
    double stmin=hName[variable_name]->GetBinCenter(hName[variable_name]->GetMaximumBin())-50;
    double zbi=hName[variable_name]->GetBinContent(hName[variable_name]->GetMaximumBin());
    double zbi_used=hName[variable_name]->GetBinContent(hName[variable_name]->FindBin(st_used));
    
    
    return Form("%.0f & %.1f & %.0f & %.1f & %.0f \\\\",mass,zbi,stmin,zbi_used,st_used);
}

void plot_categories_expected(){
    expected_bkg=true;
    int nbsig=0;
    normalization(nbsig, "_SF","");
    normalization(nbsig, "_TopCor","");

    predict_background(nbsig,"_TopCor","");
    
    predict_background(nbsig,"_SF","");
    predict_background(nbsig,"_SF","_scaleup");
    predict_background(nbsig,"_SF","_scaledown");

    predict_background(nbsig,"_SFPbc","");
    predict_background(nbsig,"_SFPl","");
	
    predict_background(nbsig,"_SFMbc","");
    predict_background(nbsig,"_SFMl","");
    //perform for 1 btag bin
    nbsig=1;
    predict_background(nbsig,"_TopCor","");

    predict_background(nbsig,"_SF","");
    predict_background(nbsig,"_SFPbc","");
    predict_background(nbsig,"_SFPl","");
    
    predict_background(nbsig,"_SFMbc","");
    predict_background(nbsig,"_SFMl","");
	
	
    overlay_stlabels=true;
    plot_systematics=true;
	drawQCD=true;
    
    TString scales[]={"_SF","_TopCor"};
    
    int Nscales=sizeof(scales)/sizeof(TString);
    
    std::vector<TString> signal_names;
    //signal_names.push_back("_stealth_300_200");
    signal_names.push_back("_stealth_600_300");
    
    vector<TString> names;
   names.push_back("_QCD");
    names.push_back("_diboson");
    names.push_back("_dy_datadriven");
    //names.push_back("_singleTop");
    names.push_back("_ttbarExp");
	vector<TString> legend_names;
    legend_names.push_back("Non-prompt");
    legend_names.push_back("Diboson");
    legend_names.push_back("Drell-Yan");
	legend_names.push_back("Single t + t#bar{t}");
	
	legend_names.push_back("M_{#tilde{q}} = 600 GeV");
	
    
    for(int ib=0; ib<=1; ib++){
        for(int i=0; i<Nscales; i++){
            TString variable_name=Form("EventCategories_1Mu_1El_%dbtag",ib)+scales[i];
            combine_histograms(variable_name,names,"_allExp");
        }
    }
    cout << "prediction systematic: " << endl;
    prediction_systematic(0,"_SF");
    prediction_systematic(1,"_SF");
    
    cout << "prediction systematic TopCor:" << endl;

    //prediction_systematic(0,"_TopCor");
    //prediction_systematic(1,"_TopCor");
    
    scaleUpDown=false;
    drawQCD=true;
    cout << "plots: " << endl;
    for(int ib=0; ib<=1; ib++){
        for(int i=0; i<Nscales; i++){
            TString variable_name=Form("EventCategories_1Mu_1El_%dbtag",ib)+scales[i];
            TString variable_name_cntrl=Form("EventCategories_1Mu_1El_%dbtag",2)+scales[i];

            // combine_histograms(variable_name,names,"_allMC");
            //blind_categoryplots(ib, scales[i]);
            plot_systematics=true;
            if(scales[i]=="_TopCor") plot_systematics=false;
            compute_ratio(variable_name, "_allExp");
        
            fill_stack_eventcategories(variable_name,"_singleMu",signal_names,names,legend_names,"_allExp");
			if(ib==0 && scales[i]=="_SF"){
				datacards_new(variable_name,variable_name_cntrl,"_stealth_300_200");
                datacards_new(variable_name,variable_name_cntrl,"_stealth_400_200");
				datacards_new(variable_name,variable_name_cntrl,"_stealth_500_300");
				datacards_new(variable_name,variable_name_cntrl,"_stealth_600_300");
				datacards_new(variable_name,variable_name_cntrl,"_stealth_700_400");
                datacards_new(variable_name,variable_name_cntrl,"_stealth_800_400");
                datacards_new(variable_name,variable_name_cntrl,"_stealth_900_500");


				//datacards(variable_name,"_singleMu", "_stealth_800_400","_allExp");

			}

        }
        
    }
    
    print_loose_events();
    
    plot_systematics=false;
	scaleUpDown=false; 
    expected_bkg=false;
    drawQCD=false;
}


void plot_MC_exp(int nbtag){
    
    TString variable_name=Form("EventCategories_1Mu_1El_%dbtag_SF",nbtag);
    TString variable_nameCtrl=Form("EventCategories_1Mu_1El_%dbtag_SF_singleMu",2);
    
    CreateCanvas(Form("MC_Predicted_Comparison_%dbtag",nbtag),"",600,600);
    CName[Form("MC_Predicted_Comparison_%dbtag",nbtag)]->cd();
    
    gPad->SetLogy();
    
    hName[variable_name+"_singleTopttbar"]->SetAxisRange(2,7);
    hName[variable_name+"_singleTopttbar"]->SetLineWidth(2);
    hName[variable_name+"_singleTopttbar"]->DrawCopy("histo");
    
    hName[variable_name+"_ttbarExp"]->SetLineWidth(2);
    hName[variable_name+"_ttbarExp"]->SetLineColor(kRed);
    hName[variable_name+"_ttbarExp"]->SetLineStyle(kDashed);
    
    hName[variable_name+"_ttbarExp"]->SetFillStyle(0);
    
    hName[variable_name+"_ttbarExp"]->DrawCopy("same");
    grName[variable_name+"_TopOnlysystematic"]->SetFillStyle(3244);
    grName[variable_name+"_TopOnlysystematic"]->SetFillColor(38);

    grName[variable_name+"_TopOnlysystematic"]->DrawClone("E2 same");
    
    TLegend *Leg = new TLegend(0.25,0.25,0.5,0.4);
    Leg->SetTextSize(0.025);
    Leg->SetFillColor(10);
	Leg->SetLineColor(10);
	Leg->SetLineWidth(0);
    Leg->SetBorderSize(0);
    Leg->AddEntry(hName[variable_name+"_singleTopttbar"],"Top MC", "L");
    Leg->AddEntry(hName[variable_name+"_ttbarExp"],"Data driven top background", "L");
    Leg->AddEntry(grName[variable_name+"_TopOnlysystematic"],"Data driven top uncertainty", "F");
    Leg->DrawClone();
    //draw_header();
    CMS_lumi(CName[Form("MC_Predicted_Comparison_%dbtag",nbtag)],2,10);
                        
}

void shape_difference(int nbtag){
    TString variable_name=Form("EventCategories_1Mu_1El_%dbtag_SF_ttbarExp",nbtag);
    TString sysName=Form("EventCategories_1Mu_1El_%dbtag_SF_allMC_Rsystematic",nbtag);
    
    CreateHistogram("ShapeSystematicU","","N_{jets}","(Top^{Scaleup}-Top)/Top",nJetmax,0.5,nJetmax+0.5);
    CreateHistogram("ShapeSystematicD","","N_{jets}","(Top^{Scaledown}-Top)/Top",nJetmax,0.5,nJetmax+0.5);
    CreateHistogram("NominalTop","","N_{jets}","Events",nJetmax,0.5,nJetmax+0.5);
    CreateHistogram("TopU","","N_{jets}","Events",nJetmax,0.5,nJetmax+0.5);
    CreateHistogram("TopD","","N_{jets}","Events",nJetmax,0.5,nJetmax+0.5);

    hName["ShapeSystematicU"]->GetXaxis()->SetBinLabel(2,"2");
    hName["ShapeSystematicU"]->GetXaxis()->SetBinLabel(3,"3");
    hName["ShapeSystematicU"]->GetXaxis()->SetBinLabel(4,"4");
    hName["ShapeSystematicU"]->GetXaxis()->SetBinLabel(5,"#geq 5");
    
    for(int nJ=2; nJ<=7; nJ++){
        int ibin=category_bin(nJ,0);
        double top=hName[variable_name]->GetBinContent(ibin);
        double topU=hName[variable_name+"_scaleup"]->GetBinContent(ibin);
        double topD=hName[variable_name+"_scaledown"]->GetBinContent(ibin);
        hName["NominalTop"]->SetBinContent(nJ,top);
        hName["TopU"]->SetBinContent(nJ,topU);
        hName["TopD"]->SetBinContent(nJ,topD);

        double topsys=grName[sysName]->GetEYhigh()[ibin-1];
    
        double Nsu=hName["Acceptance_nEvents_ttJets_scaleup"]->GetBinContent(nJ);
        double Nsd=hName["Acceptance_nEvents_ttJets_scaledown"]->GetBinContent(nJ);

        double uncU=TMath::Sqrt(topsys*topsys+2*TMath::Power(poisson(Nsu)/Nsu,2));
        double uncD=TMath::Sqrt(topsys*topsys+2*TMath::Power(poisson(Nsd)/Nsd,2));

        hName["NominalTop"]->SetBinError(nJ,topsys*top);
        hName["TopU"]->SetBinError(nJ,uncU*topU);
        hName["TopD"]->SetBinError(nJ,uncD*topD);

        double deltaR_U=(topU/top)*TMath::Sqrt(2*topsys*topsys+TMath::Power(poisson(Nsu)/Nsu,2));
        double deltaR_D=(topD/top)*TMath::Sqrt(2*topsys*topsys+TMath::Power(poisson(Nsd)/Nsd,2));

        hName["ShapeSystematicU"]->SetBinContent(nJ,(topU-top)/top);//(top-topU)/top);
        hName["ShapeSystematicU"]->SetBinError(nJ,deltaR_U);
        
        hName["ShapeSystematicD"]->SetBinContent(nJ,(topD-top)/top);//(top-topD)/top);
        hName["ShapeSystematicD"]->SetBinError(nJ,deltaR_D);

    }
    double topU=hName[variable_name+"_scaleup"]->GetBinContent(5)+hName[variable_name+"_scaleup"]->GetBinContent(6)+hName[variable_name+"_scaleup"]->GetBinContent(7);
    double topD=hName[variable_name+"_scaledown"]->GetBinContent(5)+hName[variable_name+"_scaledown"]->GetBinContent(6)+hName[variable_name+"_scaledown"]->GetBinContent(7);
    double top=hName[variable_name]->GetBinContent(5)+hName[variable_name]->GetBinContent(6)+hName[variable_name]->GetBinContent(7);
    hName["ShapeSystematicU"]->SetBinContent(5,(topU-top)/top);
    hName["ShapeSystematicD"]->SetBinContent(5,(topD-top)/top);
    
    
    hName["ShapeSystematicU"]->SetLineColor(kRed);
    hName["TopU"]->SetLineColor(kRed);

    hName["NominalTop"]->SetLineColor(kBlack);
    
    hName["ShapeSystematicD"]->SetLineColor(kBlue);
    hName["TopD"]->SetLineColor(kBlue);
    
    CreateCanvas(Form("ShapeComparison_%dbtag",nbtag),"",800,600);
    CName[Form("ShapeComparison_%dbtag",nbtag)]->Divide(2,1);
    CName[Form("ShapeComparison_%dbtag",nbtag)]->cd(1);
    
    double max=1.5*TMath::Max( hName["TopU"]->GetBinContent(2), hName["TopD"]->GetBinContent(2));
    
    hName["TopU"]->SetMaximum(max);
    
    gPad->SetLogy();
    hName["TopU"]->Draw("histo E1");
    hName["NominalTop"]->Draw("histo E1 same");
    hName["TopD"]->Draw("histo E1 same");

    CName[Form("ShapeComparison_%dbtag",nbtag)]->cd(2);
    gPad->SetLogy();
    hName[Form("EventCategories_1Mu_1El_%dbtag_SF_ttbar",nbtag)]->SetLineColor(kBlack);
    hName[Form("EventCategories_1Mu_1El_%dbtag_SF_ttbar",nbtag)]->SetFillStyle(0);
    hName[Form("EventCategories_1Mu_1El_%dbtag_SF_ttJets_scaleup",nbtag)]->SetLineColor(kRed);
    hName[Form("EventCategories_1Mu_1El_%dbtag_SF_ttJets_scaledown",nbtag)]->SetLineColor(kBlue);
    

    hName[Form("EventCategories_1Mu_1El_%dbtag_SF_ttbar",nbtag)]->SetAxisRange(1,7);
    hName[Form("EventCategories_1Mu_1El_%dbtag_SF_ttbar",nbtag)]->SetMaximum(max);
    hName[Form("EventCategories_1Mu_1El_%dbtag_SF_ttbar",nbtag)]->Draw("histo");
    hName[Form("EventCategories_1Mu_1El_%dbtag_SF_ttJets_scaleup",nbtag)]->Draw("histo E1 same");
    hName[Form("EventCategories_1Mu_1El_%dbtag_SF_ttJets_scaledown",nbtag)]->Draw("histo E1 same");

    CreateCanvas(Form("ShapeSystematic_%dbtag",nbtag),"",600,600);
    
    CName[Form("ShapeSystematic_%dbtag",nbtag)]->cd();
    
    hName["ShapeSystematicU"]->SetAxisRange(2,5);
    hName["ShapeSystematicU"]->SetMinimum(-0.5);
    hName["ShapeSystematicU"]->SetMaximum(0.5);
    
    hName["ShapeSystematicU"]->Draw("histo");
    hName["ShapeSystematicD"]->Draw("histo same");

    
    
}

void predict_systematic(int nbsig, int nJ,TString stType){
    bool print=false;
    int nbcontrol=2;
//stInclusive_nJets4_1Mu_1El_0btag_SF_stack_canvas
    TString variable_name="stInclusive"+stType+Form("_nJets%d_1Mu_1El_%dbtag_SF",nJ,nbsig);
    TString variable_nameCtrl="stInclusive"+stType+Form("_nJets%d_1Mu_1El_%dbtag_SF",nJ,nbcontrol);

    TString allMC="_allMC";
    TString data="_singleMu";
    
    if(print)cout << variable_name << endl;
    
    
    CreateGraph(hName[variable_name+allMC],variable_name+allMC+"_systematic");
    if(print) cout << variable_name+"Pbc"+allMC << endl;
    
    for(int ibin=1; ibin<=hName[variable_name+allMC]->GetNbinsX();ibin++){
        double y=hName[variable_name+allMC]->GetBinContent(ibin);
        double x=hName[variable_name+allMC]->GetBinCenter(ibin);
        
        double ybcP=hName[variable_name+"Pbc"+allMC]->GetBinContent(ibin)-y;
        double ylP=hName[variable_name+"Pl"+allMC]->GetBinContent(ibin)-y;
  
        double ybcM=hName[variable_name+"Pbc"+allMC]->GetBinContent(ibin)-y;
        double ylM=hName[variable_name+"Pl"+allMC]->GetBinContent(ibin)-y;
        
        double yE=y*hName[variable_nameCtrl+data]->GetBinError(ibin)/hName[variable_nameCtrl+data]->GetBinContent(ibin);//number of events
     
        double sysP=TMath::Sqrt(ybcP*ybcP+ylP*ylP+yE*yE);
        double sysM=TMath::Sqrt(ybcM*ybcM+ylM*ylM+yE*yE);
        
        double sys=TMath::Max(sysP,sysM);
        
        if(print)cout << "N events: " << y << " sys: " << sys << endl;
        
        
        grName[variable_name+allMC+"_systematic"]->SetPoint(ibin-1,x,sys);
    }
    
    
}

void plot_zbi(int nbsig, int nJ, TString stType, TString signal){
    TString variable_name="stInclusive"+stType+Form("_nJets%d_1Mu_1El_%dbtag_SF",nJ,nbsig);
    TString allMC="_allMC";
    
    clone_histogram(variable_name+allMC,variable_name+signal+"_ZBI");
    hName[variable_name+signal+"_ZBI"]->GetYaxis()->SetTitle("Z_{BI}");
    
    for(int ibin=1; ibin<=hName[variable_name+allMC]->GetNbinsX();ibin++){
        double x=hName[variable_name+allMC]->GetBinCenter(ibin);
        double bkg=hName[variable_name+allMC]->GetBinContent(ibin);
        double sig=hName[variable_name+signal]->GetBinContent(ibin);
        double bkgE=grName[variable_name+allMC+"_systematic"]->GetY()[ibin-1];
        
        double zbi=getZBI(sig,bkg,bkgE);
        //cout << "stmin: " << x << " zbi: "<< zbi << endl;
        if(zbi>0)hName[variable_name+signal+"_ZBI"]->SetBinContent(ibin,zbi);
    }
    
    TString canvas_name=variable_name+signal+"_ZBI";
    CreateCanvas(canvas_name,"",600,600);
    CName[canvas_name]->cd();
    hName[variable_name+signal+"_ZBI"]->Draw("histo");
}

void prediction_systematic(int nbsig, TString tag){
    cout << "prediction systematic: " << nbsig << "  " << tag << endl;
    //tag nominal tag
    bool print=false;
    TString exp="_ttbarExp";
    TString allExp="_allExp";
    TString data="_singleMu";
	//   int nbsig=0;
    int nbcontrol=2;
    
    TString variable_name=Form("EventCategories_1Mu_1El_%dbtag",nbsig)+tag;
    TString variable_nameDY=Form("EventCategories_1Mu_1El_%dbtag_SF",nbsig); // for dy, do not reweight top pt as that does not matter.

    TString variable_nameMuMu=Form("EventCategories_2Mu_onZ_%dbtag",nbsig)+tag+"_singleMu";
    
    cout << variable_name << endl;
    cout << hName[variable_name+allExp] << " " << hName[variable_name+"_ttbarExp"] << " ";
    cout << hName[variable_name+"_singleTopttbar"] << " " << hName[variable_nameDY+"_dy_datadriven"] << endl;

    CreateGraph(hName[variable_name+allExp],variable_name+allExp+"_systematic");
    CreateGraph(hName[variable_name+"_ttbarExp"],variable_name+"_TopOnlysystematic");

    CreateGraph(hName[variable_name+allExp],variable_name+allExp+"_Rsystematic");

    CreateGraph(hName[variable_name+allExp],variable_name+allExp+"_Rsystematic_stat");
    grName[variable_name+allExp+"_Rsystematic_stat"]->SetFillColor(kRed);
    
    TString variable_nameDataCtrl=Form("EventCategories_1Mu_1El_%dbtag",nbcontrol)+tag+data;//data control
	
	if(print)cout <<  "nbtag: " << nbsig << " tag: " << tag << endl;
	
    TString PlusM[]={"P","M"};
    
    for(int nJ=1; nJ<=nJetmax; nJ++){
        for(int ist=0; ist<3; ist++){
            int ibin=category_bin(nJ,ist);
            int ibin_ctr=category_bin(nJ,0);
            
            //for(int ibin=1; ibin<=hName[variable_name+exp]->GetNbinsX(); ibin++){
            double nom=hName[variable_name+exp]->GetBinContent(ibin);
            double nomTop=hName[variable_name+"_singleTopttbar"]->GetBinContent(ibin);
            double Nsys[]={0,0};
            double NsysTop[]={0,0};
            double y=grName[variable_name+allExp+"_systematic"]->GetY()[ibin-1]; //hist index starts at 1, graph at 0
            
            for(int pm=0; pm<2;  pm++){
                
                if(print)cout << " bin: " << ibin << " pm: "<< PlusM[pm] << endl;
                TString variable_nameExp_bc=variable_name+PlusM[pm]+"bc"+exp;
                TString variable_nameExp_l=variable_name+PlusM[pm]+"l"+exp;
                
                double sysbc=hName[variable_nameExp_bc]->GetBinContent(ibin);
                double sysl=hName[variable_nameExp_l]->GetBinContent(ibin);
                
                // cout << "ttbc: " << sysbc << " nom : " << nom << endl;
                // cout << "ttl: " << sysl << " nom : " << nom << endl;
                
                Nsys[pm]+=TMath::Power(sysbc-nom,2)+TMath::Power(sysl-nom,2);
                
                double sysBtag=TMath::Sqrt(TMath::Power(sysbc-nom,2)+TMath::Power(sysl-nom,2))/y; //sysBtag
                
                double Nstat=0;
                Nstat+=TMath::Power(y*poisson(hName[variable_nameDataCtrl]->GetBinContent(ibin_ctr))/hName[variable_nameDataCtrl]->GetBinContent(ibin_ctr),2);
                //cout <<"systematic bin: "<< ibin <<" " <<  TMath::Sqrt(Nstat) << " " << TMath::Sqrt(Nsys[pm]) << endl;
                // Nsys[pm]+=TMath::Power(hName[variable_nameDataCtrl]->GetBinError(ibin), 2);
                NsysTop[pm]+=TMath::Power(nomTop*poisson(hName[variable_nameDataCtrl]->GetBinContent(ibin_ctr))/hName[variable_nameDataCtrl]->GetBinContent(ibin_ctr),2);;
                //data-driven QCD + DY
                
               //Nsys[pm]+=TMath::Power(hName[variable_name+"_dy"]->GetBinContent(ibin)*poisson(hName[variable_nameMuMu]->GetBinContent(ibin))/hName[variable_nameMuMu]->GetBinContent(ibin),2);
                
                double nqcd=hName[variable_name+"_QCD"]->GetBinContent(ibin)/2;
                double NdyCor=(hName[variable_nameDY+"_dy_datadriven"]->GetBinContent(ibin)-hName[variable_name+"_dy"]->GetBinContent(ibin))/2;
                if(ist==2)nqcd=hName[variable_name+"_QCD"]->GetBinContent(ibin)*1.2;
                double NdibosonCor=(hName[variable_name+"_diboson_datadriven"]->GetBinContent(ibin)-hName[variable_name+"_diboson"]->GetBinContent(ibin))/2;
                
                
                Nsys[pm]+=TMath::Power(nqcd,2);
                Nsys[pm]+=TMath::Power(NdyCor,2);
                Nsys[pm]+=TMath::Power(NdibosonCor,2);
                Nsys[pm]+=TMath::Power(nom*0.15,2); //R0 uncertainty
                Nsys[pm]+=TMath::Power(nom*0.1,2); //Top shape uncertainty
                NsysTop[pm]+=TMath::Power(nomTop*0.15,2);
                NsysTop[pm]+=TMath::Power(nomTop*0.1,2);


                /*
                 if(ibin<=7) Nsys[pm]+=TMath::Power(nom*0.02,2);// R0 uncertainty on the three bins
                 if(ibin>7 && ibin<=14) Nsys[pm]+=TMath::Power(nom*0.12,2);
                 if(ibin>14 && ibin<=21) Nsys[pm]+=TMath::Power(nom*0.6,2);
                 */
                Nsys[pm]+=Nstat;
                
                NsysTop[pm]=TMath::Sqrt(NsysTop[pm]);
                Nsys[pm]=TMath::Sqrt(Nsys[pm]);
                if(print)cout << " stat: " << TMath::Sqrt(Nstat)/y << " b-tag: " << TMath::Sqrt(TMath::Power(sysbc-nom,2)+TMath::Power(sysl-nom,2))/y << " dy " << TMath::Sqrt(NdyCor)/y;
                if(print)cout << " diboson: " << TMath::Sqrt(NdibosonCor)/y << " Non-prompt: " << TMath::Sqrt(nqcd)/y << endl;
                
            }//+/0
            grName[variable_name+allExp+"_systematic"]->GetEXhigh()[ibin-1]=0.5;
            grName[variable_name+allExp+"_systematic"]->GetEXlow()[ibin-1]=0.5;
            
            grName[variable_name+"_TopOnlysystematic"]->GetEXhigh()[ibin-1]=0.5;
            grName[variable_name+"_TopOnlysystematic"]->GetEXlow()[ibin-1]=0.5;
            
            grName[variable_name+allExp+"_Rsystematic"]->GetEXhigh()[ibin-1]=0.5;
            grName[variable_name+allExp+"_Rsystematic"]->GetEXlow()[ibin-1]=0.5;
            
            
            grName[variable_name+allExp+"_Rsystematic_stat"]->GetEXhigh()[ibin-1]=0.5;
            grName[variable_name+allExp+"_Rsystematic_stat"]->GetEXlow()[ibin-1]=0.5;
            
            grName[variable_name+"_TopOnlysystematic"]->GetEYhigh()[ibin-1]=NsysTop[0];
            grName[variable_name+"_TopOnlysystematic"]->GetEYlow()[ibin-1]=NsysTop[1];
            
            grName[variable_name+allExp+"_systematic"]->GetEYhigh()[ibin-1]=Nsys[0];
            grName[variable_name+allExp+"_systematic"]->GetEYlow()[ibin-1]=Nsys[1];
            
            grName[variable_name+allExp+"_Rsystematic"]->GetY()[ibin-1]=1;
            grName[variable_name+allExp+"_Rsystematic"]->GetEYhigh()[ibin-1]=Nsys[0]/y;
            grName[variable_name+allExp+"_Rsystematic"]->GetEYlow()[ibin-1]=Nsys[1]/y;
            
            //cout << "ibin: " << ibin << grName[variable_name+allExp+"_Rsystematic"]->GetY()[ibin-1] << " " << grName[variable_name+allExp+"_Rsystematic"]->GetEYhigh()[ibin-1] << endl;
            
            grName[variable_name+allExp+"_Rsystematic_stat"]->GetY()[ibin-1]=1;
            grName[variable_name+allExp+"_Rsystematic_stat"]->GetEYhigh()[ibin-1]=hName[variable_nameDataCtrl]->GetBinError(ibin)/y;
            grName[variable_name+allExp+"_Rsystematic_stat"]->GetEYlow()[ibin-1]=hName[variable_nameDataCtrl]->GetBinError(ibin)/y;
			
        }//ist
    }//nj bin
    
    
}

double poisson(double n){
    
    const double alpha = 1 - 0.6827;
    double l = 0.5*TMath::ChisquareQuantile(alpha/2,2*n);
    double h = (n==0) ? ( 0.5*TMath::ChisquareQuantile(1-alpha,2*(n+1)) ) : ( 0.5*TMath::ChisquareQuantile(1-alpha/2,2*(n+1)) );
    
    return TMath::Max(h-n,n-l);
    
    //cout << "n: " << n << " +/- " << h << " " << l << endl;
    //if(n>0)cout << "n: " << n << " +/- [%] " << 100*(h-n)/n << " " << 100*(n-l)/n << endl;
}


double compute_R0(double Ndata, double *N){
    
    double Ntt=N[0];
    double Ntop=N[1];
    double Ndy=N[2];
    double Ndb=N[3];
    double Nnp=N[4];
    
    return (Ndata-Ntop-Ndy-Ndb-Nnp)/Ntt;
    
}

void normalization(int nbsig, TString tag, TString scale){
    TString variable_name=Form("EventCategories_1Mu_1El_%dbtag",nbsig)+tag;
    double Ntt=0;
    double Ntop=0;
    double Ndy=0;
    double Ndb=0;
    double Nnp=0;
    double Ndata=0;
    
    double NdyE=0;
    double NdbE=0;
    double NnpE=0;
    
    int ist=0;
    int istP=1;
    //subtract next bin
    for(int nJ=2; nJ<=3; nJ++){
        int ibin=category_bin(nJ,ist);
        int ibinP=category_bin(nJ,istP);
        Ndata+=hName[variable_name+"_singleMu"]->GetBinContent(ibin)-hName[variable_name+"_singleMu"]->GetBinContent(ibinP);
        Ntt+=hName[variable_name+"_ttbar"]->GetBinContent(ibin)-hName[variable_name+"_ttbar"]->GetBinContent(ibinP);
        Ntop+=hName[variable_name+"_singleTop"]->GetBinContent(ibin)-hName[variable_name+"_singleTop"]->GetBinContent(ibinP);
        Ndy+=hName[variable_name+"_dy_datadriven"]->GetBinContent(ibin)-hName[variable_name+"_dy_datadriven"]->GetBinContent(ibinP);
        double Ndytmp=(hName[variable_name+"_dy_datadriven"]->GetBinContent(ibin)-hName[variable_name+"_dy"]->GetBinContent(ibin))/2;
        Ndytmp-=(hName[variable_name+"_dy_datadriven"]->GetBinContent(ibinP)-hName[variable_name+"_dy"]->GetBinContent(ibinP))/2;
        
        NdyE+=TMath::Power(Ndytmp,2);
        
        Ndb+=hName[variable_name+"_diboson"]->GetBinContent(ibin)-hName[variable_name+"_diboson"]->GetBinContent(ibinP);
        NdbE+=TMath::Power((hName[variable_name+"_diboson"]->GetBinContent(ibin)-hName[variable_name+"_diboson"]->GetBinContent(ibinP))*0.3,2);
        
        Nnp+=hName[variable_name+"_QCD"]->GetBinContent(ibin)-hName[variable_name+"_QCD"]->GetBinContent(ibinP);
    }
    NdbE=TMath::Sqrt(NdbE);
    NdyE=TMath::Sqrt(NdyE);
    NnpE=Nnp*0.5;
    
    cout << "NdyE: " << NdyE << "  " << Ndy*0.02 << endl;
    cout << "NdbE: " << NdbE << " " << Ndb*0.3 << endl;
    
    
    double Nst0[]={Ntt,Ntop,Ndy,Ndb,Nnp};
    double Nst0EdyP[]={Ntt,Ntop,Ndy+NdyE,Ndb,Nnp};
    double Nst0EdbP[]={Ntt,Ntop,Ndy,Ndb+NdbE,Nnp};
    double Nst0EnpP[]={Ntt,Ntop,Ndy,Ndb,Nnp+NnpE};
    
    double Nst0EdyM[]={Ntt,Ntop,Ndy-NdyE,Ndb,Nnp};
    double Nst0EdbM[]={Ntt,Ntop,Ndy,Ndb-NdbE,Nnp};
    double Nst0EnpM[]={Ntt,Ntop,Ndy,Ndb,Nnp-NnpE};
    
    cout << "R st0: " << compute_R0(Ndata, Nst0) << endl;
    cout << "R st0 sigma: " << compute_R0(Ndata+poisson(Ndata), Nst0)-compute_R0(Ndata, Nst0) << endl;
    cout << "R st0 dy sys sigma+: " << compute_R0(Ndata, Nst0EdyP)-compute_R0(Ndata, Nst0) << endl;
    cout << "R st0 db sys sigma+: " << compute_R0(Ndata, Nst0EdbP)-compute_R0(Ndata, Nst0) << endl;
    cout << "R st0 np sys sigma+: " << compute_R0(Ndata, Nst0EnpP)-compute_R0(Ndata, Nst0) << endl;

    cout << "R st0 dy sys sigma-: " << compute_R0(Ndata, Nst0EdyM)-compute_R0(Ndata, Nst0) << endl;
    cout << "R st0 db sys sigma-: " << compute_R0(Ndata, Nst0EdbM)-compute_R0(Ndata, Nst0) << endl;
    cout << "R st0 np sys sigma-: " << compute_R0(Ndata, Nst0EnpM)-compute_R0(Ndata, Nst0) << endl;
    
    Ntt=0;
    Ntop=0;
    Ndy=0;
    Ndb=0;
    Nnp=0;
    Ndata=0;
    
    NdyE=0;
    NdbE=0;
    NnpE=0;
    
    for(int nJ=2; nJ<=3; nJ++){
        int ibin=category_bin(nJ,ist);
        int ibinP=category_bin(nJ,istP);
        Ndata+=hName[variable_name+"_singleMu"]->GetBinContent(ibinP);
        Ntt+=hName[variable_name+"_ttbar"]->GetBinContent(ibinP);
        Ntop+=hName[variable_name+"_singleTop"]->GetBinContent(ibinP);
        Ndy+=hName[variable_name+"_dy_datadriven"]->GetBinContent(ibinP);
        Ndb+=hName[variable_name+"_diboson"]->GetBinContent(ibinP);
        Nnp+=hName[variable_name+"_QCD"]->GetBinContent(ibinP);
        
        NdyE+=TMath::Power(hName[variable_name+"_dy_datadriven"]->GetBinContent(ibinP)-hName[variable_name+"_dy"]->GetBinContent(ibinP),2);

    }
    NdyE=TMath::Sqrt(NdyE);
    NdbE=Ndb*0.3;
    NnpE=Nnp*0.5;
    
    double Nst1[]={Ntt,Ntop,Ndy,Ndb,Nnp};
    double Nst1Edy[]={Ntt,Ntop,Ndy+NdyE,Ndb,Nnp};
    double Nst1Edb[]={Ntt,Ntop,Ndy,Ndb+NdbE,Nnp};
    double Nst1Enp[]={Ntt,Ntop,Ndy,Ndb,Nnp+NnpE};
    cout << "R st1: " << compute_R0(Ndata,Nst1) << endl;
    cout << "R st1 sigma: " << compute_R0(Ndata+poisson(Ndata), Nst1)-compute_R0(Ndata, Nst1) << endl;
    cout << "R st0 dy sys sigma: " << compute_R0(Ndata, Nst1Edy)-compute_R0(Ndata, Nst1) << endl;
    cout << "R st0 db sys sigma: " << compute_R0(Ndata, Nst1Edb)-compute_R0(Ndata, Nst1) << endl;
    cout << "R st0 np sys sigma: " << compute_R0(Ndata, Nst1Enp)-compute_R0(Ndata, Nst1) << endl;

}

void predict_background(int nbsig,TString tag, TString scale){
	//    int nbsig=0;
    int nbcontrol=2;
    
    TString data="_singleMu";
    TString ttbarSample="_singleTopttbar"+scale;
    TString exp="_ttbarExp"+scale;
    
    bool print=false; 
	//if(tag=="_SF")print=true; 
    
    TString TopName="_ttbar";
    if(scale=="_scaleup")TopName="_ttJets_scaleup";
    if(scale=="_scaledown")TopName="_ttJets_scaledown";

    
    TString variable_nameCtrl=Form("EventCategories_1Mu_1El_%dbtag",nbcontrol)+tag+"_allMC";
    TString variable_nameSig=Form("EventCategories_1Mu_1El_%dbtag",nbsig)+tag+"_allMC";
    
    TString VNSig=Form("EventCategories_1Mu_1El_%dbtag",nbsig)+tag;
    
    TString variable_nameCtrldy=Form("EventCategories_2Mu_onZ_%dbtag",nbsig)+tag;
    
    TString variable_nameCtrlTT=Form("EventCategories_1Mu_1El_%dbtag",nbcontrol)+tag+ttbarSample;
    TString variable_nameSigTT=Form("EventCategories_1Mu_1El_%dbtag",nbsig)+tag+ttbarSample;
    
    TString variable_nameExp=Form("EventCategories_1Mu_1El_%dbtag",nbsig)+tag+exp;
    
    TString variable_nameDataSig=Form("EventCategories_1Mu_1El_%dbtag",nbsig)+tag+data; //data signal region
    TString variable_nameDataCtrl=Form("EventCategories_1Mu_1El_%dbtag",nbcontrol)+tag+data;//data control
	
    clone_histogram(variable_nameCtrl,variable_nameExp);
    
	// cout << "predict background: " << variable_nameExp<< endl;
    
    int ist=0;
    int nJ=1;
    
    double Ndatasig[]={0,0,0};//ist bins
    double NdataCtrl[]={0,0,0};
    
    double NMCsig[]={0,0,0};//ist bins
    double NMCCtrl[]={0,0,0};
    
    double NMCothersig[]={0,0,0};//ist bins, 0 btags, 2-3 jets

    double NMCttsig[]={0,0,0};//ist bins, 0 btags, 2-3 jets
    double Nsysother[]={0,0,0};
    
    double Ntop_uwMC[]={0,0,0};
    double Ndy_uwMC[]={0,0,0};
    double Ndb_uwMC[]={0,0,0};
    
    double Ndy[]={0,0,0};
    double Ndb[]={0,0,0};
    
    double deltaDY[]={0,0,0};
    double DYCntrl[]={0,0,0};

    double DY[]={0,0,0};
    
    double deltaDB[]={0,0,0};
    double deltaQCD[]={0,0,0};

    double R0[]={0,0,0}; //normalize 0 btag
    double R2[]={0,0,0}; //normalize 2 btag
    
    double R[]={0,0,1.03,1.096,1.11,1.23,1.56,1.56,1.56}; //dy correction 0-8 jets

	
    for(int ibin=1; ibin<=hName[variable_nameCtrl]->GetNbinsX();ibin++){
        
     
        if(nJ>=2 && nJ<=3)NMCsig[ist]+=hName[variable_nameSig]->GetBinContent(ibin);
        if(nJ>=2 && nJ<=nJetmax)NMCCtrl[ist]+=hName[variable_nameCtrl]->GetBinContent(ibin);
        
        if(nJ>=2 && nJ<=nJetmax)NdataCtrl[ist]+=hName[variable_nameDataCtrl]->GetBinContent(ibin);
        if(nJ>=2 && nJ<=3){
            
            double dydata=hName[variable_nameCtrldy+"_singleMu"]->GetBinContent(ibin);
            DYCntrl[ist]+=dydata;
            DY[ist]+=hName[VNSig+"_dy"]->GetBinContent(ibin)*(R[nJ]);
            deltaDY[ist]+=hName[VNSig+"_dy"]->GetBinContent(ibin)*(R[nJ]-1)/2;
            deltaDB[ist]+=hName[VNSig+"_diboson"]->GetBinContent(ibin)*(R[nJ]-1)/2;
            if(ist<2)deltaQCD[ist]+=hName[VNSig+"_QCD"]->GetBinContent(ibin)/2;
            else deltaQCD[ist]+=hName[VNSig+"_QCD"]->GetBinContent(ibin)*1.2;

            Ntop_uwMC[ist]+=hName["Acceptance_nEvents_singleTop"]->GetBinContent(ibin);
            Ntop_uwMC[ist]+=hName["Acceptance_nEvents_ttbar"]->GetBinContent(ibin);
            Ndy_uwMC[ist]+=hName["Acceptance_nEvents_dy"]->GetBinContent(ibin);
            Ndb_uwMC[ist]+=hName["Acceptance_nEvents_diboson"]->GetBinContent(ibin);

            Ndy[ist]+=R[nJ]*hName[VNSig+"_dy"]->GetBinContent(ibin);
            Ndb[ist]+=R[nJ]*hName[VNSig+"_diboson"]->GetBinContent(ibin);

            Ndatasig[ist]+=hName[variable_nameDataSig]->GetBinContent(ibin);
            NMCothersig[ist]+=R[nJ]*hName[VNSig+"_dy"]->GetBinContent(ibin);
            NMCothersig[ist]+=hName[VNSig+"_diboson"]->GetBinContent(ibin);
            NMCothersig[ist]+=hName[VNSig+"_QCD"]->GetBinContent(ibin);
            NMCttsig[ist]+=hName[VNSig+TopName]->GetBinContent(ibin);
            NMCttsig[ist]+=hName[VNSig+"_singleTop"]->GetBinContent(ibin);
            
            if(nbsig==0 && ist==0){
                /*
                cout << tag << " nj " << nJ << " ist: " << ist;
                cout << " dy " << hName[VNSig+"_dy"]->GetBinContent(ibin) <<  " +/- " << hName[VNSig+"_dy"]->GetBinContent(ibin)*poisson(hName["Acceptance_nEvents_dy"]->GetBinContent(ibin))/hName["Acceptance_nEvents_dy"]->GetBinContent(ibin);
                cout << " db " << hName[VNSig+"_diboson"]->GetBinContent(ibin) << " +/- " << hName[VNSig+"_diboson"]->GetBinContent(ibin)*poisson(hName["Acceptance_nEvents_diboson"]->GetBinContent(ibin))/hName["Acceptance_nEvents_diboson"]->GetBinContent(ibin);
                 
                cout << " np " << hName[VNSig+"_QCD"]->GetBinContent(ibin);
                cout << " top " << hName[VNSig+"_singleTop"]->GetBinContent(ibin);
                cout << " tt " << hName[VNSig+"_ttbar"]->GetBinContent(ibin);
                cout << " data " << hName[VNSig+"_singleMu"]->GetBinContent(ibin) << endl;
                */
            }
            
        }
        
        nJ++;
        if(nJ>nJetmax){
            nJ=1;
            ist++;
        }
    }
  
    
    for(int ist=0; ist<3; ist++) {
        //cout << "Ndata: " << Ndatasig[ist] << " +/- " <<poisson(Ndatasig[ist]) << endl;
     
        R0[ist]=(Ndatasig[ist]-NMCothersig[ist])/NMCttsig[ist];
        Nsysother[ist]+=TMath::Power(poisson(Ndatasig[ist]),2);
        /*
        Nsysother[ist]+=TMath::Power(R0[ist]*NMCttsig[ist]*poisson(Ntop_uwMC[ist])/Ntop_uwMC[ist],2);
        Nsysother[ist]+=TMath::Power(Ndy[ist]*poisson(Ndy_uwMC[ist])/Ndy_uwMC[ist],2);
        Nsysother[ist]+=TMath::Power(Ndb[ist]*poisson(Ndb_uwMC[ist])/Ndb_uwMC[ist],2);
        */
        
        double finite_dyunc=DY[ist]*poisson(DYCntrl[ist])/DYCntrl[ist];
        
     //   Nsysother[ist]+=TMath::Power(finite_dyunc,2);
        Nsysother[ist]+=TMath::Power(deltaDY[ist],2);
        Nsysother[ist]+=TMath::Power(deltaDB[ist],2);
        Nsysother[ist]+=TMath::Power(deltaQCD[ist],2);
        Nsysother[ist]=TMath::Sqrt(Nsysother[ist]);
        double R0E=Nsysother[ist]/NMCttsig[ist];
        R2[ist]=NMCCtrl[ist]/NdataCtrl[ist];
        
        if(print)cout << "R0: " << R0[ist] << " R2: " << R2[ist] << endl;
        if(nbsig==0 && ist==0){
            /*
            cout << "relative importance of each term in sys: " << endl;
            cout << "data: " <<poisson(Ndatasig[ist])/NMCttsig[ist] << endl;
            cout << "dy: " << finite_dyunc/NMCttsig[ist] << endl;
            cout << "db: " << deltaDB[ist]/NMCttsig[ist] << endl;
            cout << "qcd: " <<deltaQCD[ist]/NMCttsig[ist] << endl;
            
            cout << "dy fit uncertainty: " << deltaDY[ist] <<" data: " << DY[ist]*poisson(DYCntrl[ist])/DYCntrl[ist] << endl;
            cout << "relative uncertainty on data: " << poisson(Ndatasig[ist])/Ndatasig[ist] << endl;
            
            cout <<"tag: " << tag << " data: " << Ndatasig[ist] << " mc other: " << NMCothersig[ist] << " top " << NMCttsig[ist] << endl;
            cout << "R0: " << ist << " "<< R0[ist] << " +/- : " << R0E << endl;
            */
        }
    }
    
    double R2exc[]={0,0,0};
    R2exc[0]=(NMCCtrl[0]-NMCCtrl[1])/(NdataCtrl[0]-NdataCtrl[1]);
    R2exc[1]=(NMCCtrl[1]-NMCCtrl[2])/(NdataCtrl[1]-NdataCtrl[2]);
    R2exc[2]=NMCCtrl[2]/NdataCtrl[2];
    
    double R0exc[]={0,0,0};
    R0exc[0]=((Ndatasig[0]-Ndatasig[1])-(NMCothersig[0]-NMCothersig[1]))/(NMCttsig[0]-NMCttsig[1]);
    R0exc[1]=(Ndatasig[1]-NMCothersig[1])/NMCttsig[1];
    R0exc[2]=(Ndatasig[2]-NMCothersig[2])/NMCttsig[2];

    
   // cout << "tag: " << tag << endl;
    cout << "scale: " << scale <<  " " << tag << endl;
    if(nbsig==0 && (tag=="_SF" || tag=="_TopCor")){
        cout << "R0exc0: " << R0exc[0] << " " << R0exc[0]*poisson(Ndatasig[0]-Ndatasig[1])/(Ndatasig[0]-Ndatasig[1]) << endl;
        cout << "R0exc1: " << R0exc[1] << " " << R0exc[1]*poisson(Ndatasig[1]-NMCothersig[1])/(Ndatasig[1]-NMCothersig[1]) << endl;
        cout << "R0exc2: " << R0exc[2] << endl;
    }
    
    R0[1]=R0[0];
    R0[2]=R0[0];

    nJ=1;
    ist=0;
	
    for(int ist=0; ist<3; ist++){
        TString correctionName=Form("Correction_%dbcntl_%dbsig_%dst",nbcontrol,nbsig,ist)+tag;
		if(scale=="")CreateHistogram(correctionName,"","N_{jets}","Shape Correction",nJetmax,0.5,nJetmax+0.5);
        
        correctionName=Form("Correction_%dbcntl_%dbsig_excl%dst",nbcontrol,nbsig,ist)+tag;
		if(scale=="")CreateHistogram(correctionName,"","N_{jets}","Shape Correction",nJetmax,0.5,nJetmax+0.5);
	}
	
    
    for(int ist=0; ist<3; ist++){
        for(int nJ=2; nJ<=nJetmax; nJ++){
            int ibin=category_bin(nJ,ist);
            int ibin_ctr=category_bin(nJ,0);
            int ibin_st1=category_bin(nJ,1);
            int ibin_st2=category_bin(nJ,2);
            
            TString correctionNameExcl_=Form("Correction_%dbcntl_%dbsig_excl%dst",nbcontrol,nbsig,ist)+tag;
            TString correctionName_=Form("Correction_%dbcntl_%dbsig_%dst",nbcontrol,nbsig,ist)+tag;

            double MCsig=hName[variable_nameSigTT]->GetBinContent(ibin);
            double MCCtrl=hName[variable_nameCtrlTT]->GetBinContent(ibin); // always use the lowest st bins...
            double dataCtrl=hName[variable_nameDataCtrl]->GetBinContent(ibin); // always use the lowest est bins
            
            double MCCtrl_exc=0;
            if(ist==0)MCCtrl_exc=hName[variable_nameCtrlTT]->GetBinContent(ibin_ctr)-hName[variable_nameCtrlTT]->GetBinContent(ibin_st1); //exclusive samples
            if(ist==1)MCCtrl_exc=hName[variable_nameCtrlTT]->GetBinContent(ibin_st1);
            if(ist==2)MCCtrl_exc=hName[variable_nameCtrlTT]->GetBinContent(ibin);
            
            double dataCtrl_exc=0;
            if(ist==0) dataCtrl_exc=hName[variable_nameDataCtrl]->GetBinContent(ibin_ctr)-hName[variable_nameDataCtrl]->GetBinContent(ibin_st1); //exclusive samples
            if(ist==1) dataCtrl_exc=hName[variable_nameDataCtrl]->GetBinContent(ibin_st1);
            if(ist==2) dataCtrl_exc=hName[variable_nameDataCtrl]->GetBinContent(ibin);
            
            double dataCtrlE_exc=poisson(dataCtrl_exc);
            
            double dataCtrlE=poisson(hName[variable_nameDataCtrl]->GetBinContent(ibin));//hName[variable_nameDataCtrl]->GetBinError(ibin);
            
            double MCCtrl_=hName[variable_nameCtrlTT]->GetBinContent(ibin_ctr);
            double dataCtrl_=hName[variable_nameDataCtrl]->GetBinContent(ibin_ctr); // always use the lowest st bins
            
            double exp=MCsig*R0[ist]*R2[0]*(dataCtrl_/MCCtrl_);
            
            double Corr=R2[ist]*(dataCtrl_/MCCtrl_);
            double CorrE=Corr*(poisson(dataCtrl_)/dataCtrl_);
            
            double CorrExcl=R2exc[ist]*(dataCtrl_exc/MCCtrl_exc);
            double CorrExclE=(dataCtrlE_exc/dataCtrl_exc)*R2exc[ist]*(dataCtrl_exc/MCCtrl_exc);

            if(nJ>1 && scale=="")hName[correctionName_]->SetBinContent(nJ,Corr);
            if(nJ>1 && scale=="")hName[correctionName_]->SetBinError(nJ,CorrE);
            
            if(nJ>1 && scale=="")hName[correctionNameExcl_]->SetBinContent(nJ,CorrExcl);
            if(nJ>1 && scale=="")hName[correctionNameExcl_]->SetBinError(nJ,CorrExclE);
        }
    }
    
    TF1 *f=new TF1("f_pol","pol1",1.5,7.5);
    
    double p0=1;
    double p1=0;
    
    
    TString correctionName_=Form("Correction_%dbcntl_%dbsig_%dst",nbcontrol,nbsig,0)+tag;
    hName[correctionName_]->Fit("f_pol","Q+");
    
    p0=f->GetParameter(0);
    p1=f->GetParameter(1);

    for(int ist=0; ist<3; ist++){
        for(int nJ=2; nJ<=nJetmax; nJ++){
            int ibin=category_bin(nJ,ist);
            int ibin_ctr=category_bin(nJ,0);
            int ibin_st1=category_bin(nJ,1);
            int ibin_st2=category_bin(nJ,2);

            double MCsig=hName[variable_nameSigTT]->GetBinContent(ibin);
            double MCCtrl=hName[variable_nameCtrlTT]->GetBinContent(ibin); // always use the lowest st bins...
            double dataCtrl=hName[variable_nameDataCtrl]->GetBinContent(ibin); // always use the lowest est bins
            
            double MCCtrl_exc=0;
            if(ist==0)MCCtrl_exc=hName[variable_nameCtrlTT]->GetBinContent(ibin_ctr)-hName[variable_nameCtrlTT]->GetBinContent(ibin_st1); //exclusive samples
            if(ist==1)MCCtrl_exc=hName[variable_nameCtrlTT]->GetBinContent(ibin_st1);
            if(ist==2)MCCtrl_exc=hName[variable_nameCtrlTT]->GetBinContent(ibin);

            double dataCtrl_exc=0;
            if(ist==0) dataCtrl_exc=hName[variable_nameDataCtrl]->GetBinContent(ibin_ctr)-hName[variable_nameDataCtrl]->GetBinContent(ibin_st1); //exclusive samples
            if(ist==1) dataCtrl_exc=hName[variable_nameDataCtrl]->GetBinContent(ibin_st1);
            if(ist==2) dataCtrl_exc=hName[variable_nameDataCtrl]->GetBinContent(ibin);

            double dataCtrlE_exc=poisson(dataCtrl_exc);
            
            double dataCtrlE=poisson(hName[variable_nameDataCtrl]->GetBinContent(ibin));//hName[variable_nameDataCtrl]->GetBinError(ibin);
            
            double MCCtrl_=hName[variable_nameCtrlTT]->GetBinContent(ibin_ctr);
            double dataCtrl_=hName[variable_nameDataCtrl]->GetBinContent(ibin_ctr); // always use the lowest st bins

            double exp=MCsig*R0[ist]*R2[0]*(dataCtrl_/MCCtrl_);
            
            //double exp=MCsig*R0[ist]*(p0+p1*static_cast<float>(nJ));
            double Corr=R2exc[ist]*(dataCtrl_exc/MCCtrl_exc);
            double CorrE=(dataCtrlE_exc/dataCtrl_exc)*R2exc[ist]*(dataCtrl_exc/MCCtrl_exc);
            
            if(nbsig==0 && ist==0) {
               // cout << "nJ: " << nJ << " " << Corr << " " << CorrE << endl;
                //cout << dataCtrl_/MCCtrl_ << " " << p0[0]+p1[0]*static_cast<float>(nJ) << endl;
            }
            
            if(print && nJ>1) cout << "ist: "<< ist << " nJ: " << nJ << " " << R2[ist]*(dataCtrl/MCCtrl) << " +/- " << (dataCtrlE/dataCtrl)*R2[ist]*(dataCtrl/MCCtrl) << endl;
            if(print && nJ>1) cout << "nJ: " <<(dataCtrl/MCCtrl) << endl;
            
            if(nJ==1) exp=0.000000000001;
            hName[variable_nameExp]->SetBinContent(ibin,exp);

        }
    }
  
}

void predict_background(int nbsig, int nJ,TString stType){
    //predict background for inclusive st plot
    //could be inclusive or exclusive jet multiplicity
    
    TString baseName="stInclusive"+stType+Form("_nJets%d_1Mu_1El",nJ);
    int nbcontrol=2;

    TString bTagControl=Form("_%dbtag",nbcontrol);
    TString bTagSig=Form("_%dbtag",nbsig);
    
    TString data="_singleMu";
    TString ttbarSample="_singleTopttbar";
    TString exp="_ttbarExp";
    bool print=false;
	//if(tag=="_SF")print=true;
    
    
    TString variable_nameCtrl=baseName+bTagControl+"_allMC";
    TString variable_nameSig=baseName+bTagSig+"_allMC";
    
    TString variable_nameCtrlTT=baseName+bTagControl+ttbarSample;
    TString variable_nameSigTT=baseName+bTagSig+ttbarSample;
    
    TString variable_nameExp=baseName+bTagSig+exp;
    
    TString variable_nameDataSig=baseName+bTagSig+data; //data signal region
    TString variable_nameDataCtrl=baseName+bTagControl+data;//data control
    
    //compute R0
    //compute R2

    
    for(int ibin=1; ibin<=hName[variable_nameCtrl]->GetNbinsX();ibin++){
        TString name2J="stInclusive"+stType+Form("_nJets%d_1Mu_1El_%dbtag",2,nbsig)+data;
        TString name3J="stInclusive"+stType+Form("_nJets%d_1Mu_1El_%dbtag",3,nbsig)+data;

        TString name2J_MC="stInclusive"+stType+Form("_nJets%d_1Mu_1El_%dbtag_allMC",2,nbsig);
        TString name3J_MC="stInclusive"+stType+Form("_nJets%d_1Mu_1El_%dbtag_allMC",3,nbsig);
        
        
        //Compute R0. Sum 2-3 jet contents, data/MC for 0 btags (or nbsig, could also be 1 I suppose).
        double N0=hName[name2J]->GetBinContent(ibin);
        N0+=hName[name3J]->GetBinContent(ibin);
        
        double D0=hName[name2J_MC]->GetBinContent(ibin);
        D0+=hName[name3J_MC]->GetBinContent(ibin);
        
        double R0=N0/D0;
        
        //Compute R2. Use all jet bins in control region - 2 btags
        
        double N2=0;
        double D2=0;
        double R2=0;
        
        for(int nJ=2; nJ<=nJetmax; nJ++){
            TString name="stInclusive"+stType+Form("_nJets%d_1Mu_1El_%dbtag",nJ,nbcontrol)+data;//data
            TString name_MC="stInclusive"+stType+Form("_nJets%d_1Mu_1El_%dbtag_allMC",nJ,nbcontrol);//MC

            N2+=hName[name]->GetBinContent(ibin);
            D2+=hName[name_MC]->GetBinContent(ibin);
            
        }
        R2=N2/D2;

        double MCsig=hName[variable_nameSigTT]->GetBinContent(ibin);
        double MCCtrl=hName[variable_nameCtrlTT]->GetBinContent(ibin);
        double dataCtrl=hName[variable_nameDataCtrl]->GetBinContent(ibin);
		double dataCtrlE=hName[variable_nameDataCtrl]->GetBinError(ibin);
        
        double exp=MCsig*R0*R2*(dataCtrl/MCCtrl);
		double Corr=R2*(dataCtrl/MCCtrl) ;
		double CorrE=(dataCtrlE/dataCtrl)*R2*(dataCtrl/MCCtrl);
    }//end for loop over st bins
    
}

void draw_correction(int nbsig){
	TString tag="_SF";
	int nbcontrol=2;
    
	CreateCanvas(Form("Correction_%dbsig",nbsig),"",600,600);
	CName[Form("Correction_%dbsig",nbsig)]->cd();
    
    TLegend Leg(0.25,0.55,0.55,0.85);
    Leg.SetFillColor(10);
	Leg.SetLineColor(10);
	Leg.SetLineWidth(0);
    Leg.SetBorderSize(0);
    
    int st[]={300,700,1200};
    
    double chi2=0;
    double NDF=0;
    
    TF1 *fit_fcn = new TF1("fit_fcn","pol1",2,7);
    
    for(int ist=0; ist<2; ist++) {
        TString correctionName=Form("Correction_%dbcntl_%dbsig_excl%dst",nbcontrol,nbsig,ist)+tag;

        if(ist==0){
            hName[correctionName]->SetMinimum(-0.5);
            hName[correctionName]->SetMaximum(2.5);
        }
        
        hName[correctionName]->Fit("fit_fcn","Q0");
        
        cout << "ist: " << ist << " p0 " << fit_fcn->GetParameter(0) << " +/- " << fit_fcn->GetParError(0) << endl;
        cout << "ist: " << ist << " p1 " << fit_fcn->GetParameter(1) << " +/- " << fit_fcn->GetParError(1) << endl;

        for(int nJ=2; nJ<=nJetmax; nJ++){
            double y=hName[correctionName]->GetBinContent(nJ);
            double yE=hName[correctionName]->GetBinError(nJ);
            chi2+=TMath::Power((y-1)/yE,2);
            NDF+=1.0;
        }
        
        
        hName[correctionName]->SetMarkerStyle(22+ist);
        if(ist==0)hName[correctionName]->SetLineColor(kRed);
        if(ist==1)hName[correctionName]->SetLineColor(kBlue);
        
        if(ist==0)hName[correctionName]->SetMarkerColor(kRed);
        if(ist==1)hName[correctionName]->SetMarkerColor(kBlue);
        
        if(ist==0)Leg.AddEntry(hName[correctionName],Form("%d < S_{T} < %d GeV",st[ist],st[ist+1]),"p");
        if(ist==1)Leg.AddEntry(hName[correctionName],Form("S_{T} > %d GeV",st[ist]),"p");

        if(ist==0)hName[correctionName]->Draw("E1");
        else hName[correctionName]->Draw("E1 same");
    }
    
    cout << "combined chi^2/NDF: " << chi2 << " / " << NDF << endl;
    
    TF1 *ff=new TF1("ff","pol0", 0.5,8);
    ff->SetParameter(0,1);
    ff->SetLineColor(kBlack);
    ff->DrawClone("same");
    Leg.DrawClone();
    CMS_lumi(CName[Form("Correction_%dbsig",nbsig)],2,10);
    //draw_header();
    delete fit_fcn;
	
}

void draw_correction(int nbsig,int ist){
	TString tag="_SF"; 
	int nbcontrol=2;
    int st[]={300,700,1200};
	TString correctionName=Form("Correction_%dbcntl_%dbsig_%dst",nbcontrol,nbsig,ist)+tag; 

	CreateCanvas(Form("Correction_%dbsig_st%d",nbsig,ist),"",600,600);
	CName[Form("Correction_%dbsig_st%d",nbsig,ist)]->cd();
	hName[correctionName]->Draw("E1");
    
    TF1 *ff=new TF1("ff","pol0", 0.5,8);
    ff->SetParameter(0,1);
    ff->SetLineColor(kBlack);
    ff->Draw("same");
    
    TLatex txt(0.2,0.7,Form("S_{T} > %d GeV, #geq%d b-tag",st[ist],2));
    txt.SetNDC(kTRUE);
    txt.DrawClone();
    //draw_header();
    CMS_lumi(CName[Form("Correction_%dbsig_st%d",nbsig,ist)],2,10);
	
}

void sum_histograms(){
    vector<TString> names;
    names.push_back("_dy");
    names.push_back("_diboson");
    names.push_back("_singleTop");
    names.push_back("_ttbar");
    
    
    TString scales[]={"_noSF","_SF","_SFPbc","_SFMbc","_SFPl","_SFMl","_TopCor"};
    
    int Nscales=sizeof(scales)/sizeof(TString);
    
	TString variable[]={"EventCategories","EventCategories_SS","st","met","bjet1pt","bjet2pt","ptl1","ptl2"};
	int Nvariable=sizeof(variable)/sizeof(TString);

	
    for(int ib=0; ib<=2; ib++){
        for(int i=0; i<Nscales; i++){
			for(int ivar=0; ivar<Nvariable; ivar++){
				
				TString variable_name=variable[ivar]+Form("_1Mu_1El_%dbtag",ib)+scales[i];
				combine_histograms(variable_name,names,"_allMC");
			}
        }
    }
    
    //loose btag study
    for(int ib=0; ib<=2; ib++){
        TString variable_name=variable[0]+Form("_1Mu_1El_%dbtag_loose_noSF",ib);
        combine_histograms(variable_name,names,"_allMC");
        
        variable_name=variable[0]+Form("_SS_1Mu_1El_%dbtag_loose_noSF",ib);
        combine_histograms(variable_name,names,"_allMC");
    }

    
    //estimate QCD
    
    for(int ib=0; ib<=2; ib++){
        for(int i=0; i<Nscales; i++){
			for(int ivar=0; ivar<Nvariable; ivar++){
				if(variable[ivar].Contains("EventCategories")==0) continue;
				TString variable_name=variable[ivar]+Form("_1Mu_1El_%dbtag",ib)+scales[i];
				clone_histogram(variable_name+"_singleMu",variable_name+"_QCD");
                hName[variable_name+"_QCD"]->Add(hName[variable_name+"_allMC"],-1);
			}
        }
    }
    
    vector<TString> names2;
    names2.push_back("_singleTop");
    names2.push_back("_ttbar");
	
    
    for(int ib=0; ib<=2; ib++){
        for(int i=0; i<Nscales; i++){
            TString variable_name=Form("EventCategories_1Mu_1El_%dbtag",ib)+scales[i];
            combine_histograms(variable_name,names2,"_singleTopttbar");
        }
    }
    
    //scaleup
    names2.erase(names2.begin()+1);
    names2.push_back("_ttJets_scaleup");
    
    for(int ib=0; ib<=2; ib++){
        for(int i=0; i<Nscales; i++){
            TString variable_name=Form("EventCategories_1Mu_1El_%dbtag",ib)+scales[i];
            combine_histograms(variable_name,names2,"_singleTopttbar_scaleup");
        }
    }
    
    //scaledown
    names2.erase(names2.begin()+1);
    names2.push_back("_ttJets_scaledown");
    
    for(int ib=0; ib<=2; ib++){
        for(int i=0; i<Nscales; i++){
            TString variable_name=Form("EventCategories_1Mu_1El_%dbtag",ib)+scales[i];
            combine_histograms(variable_name,names2,"_singleTopttbar_scaledown");
        }
    }
    
	//scaleup
	
	vector<TString> namesScaleup;
	namesScaleup.push_back("_dy");
    namesScaleup.push_back("_diboson");
    namesScaleup.push_back("_singleTop");
    namesScaleup.push_back("_ttJets_scaleup");
	
	
    for(int ib=0; ib<=2; ib++){
        for(int i=0; i<Nscales; i++){
            TString variable_name=Form("EventCategories_1Mu_1El_%dbtag",ib)+scales[i];
            combine_histograms(variable_name,namesScaleup,"_allMC_scaleup");
            //normalize_systematic_categories(variable_name, "_allMC_scaleup");

        }
    }
	//scale down ttbar
	vector<TString> namesScaledown;
	namesScaledown.push_back("_dy");
    namesScaledown.push_back("_diboson");
    namesScaledown.push_back("_singleTop");
    namesScaledown.push_back("_ttJets_scaledown");
	
    
    for(int ib=0; ib<=2; ib++){
        for(int i=0; i<Nscales; i++){
            TString variable_name=Form("EventCategories_1Mu_1El_%dbtag",ib)+scales[i];
            combine_histograms(variable_name,namesScaledown,"_allMC_scaledown");
			//normalize_systematic_categories(variable_name, "_allMC_scaledown"); 

        }
    }
	
	//matching up
	vector<TString> namesMatchingup;
	namesMatchingup.push_back("_dy");
    namesMatchingup.push_back("_diboson");
    namesMatchingup.push_back("_singleTop");
    namesMatchingup.push_back("_ttJets_matchingup");
	
    
    for(int ib=0; ib<=2; ib++){
        for(int i=0; i<Nscales; i++){
            TString variable_name=Form("EventCategories_1Mu_1El_%dbtag",ib)+scales[i];
            combine_histograms(variable_name,namesScaledown,"_allMC_matchingup");
			//normalize_systematic_categories(variable_name, "_allMC_matchingup"); 

        }
    }
	//matching down
	vector<TString> namesMatchingdown;
	namesMatchingdown.push_back("_dy");
    namesMatchingdown.push_back("_diboson");
    namesMatchingdown.push_back("_singleTop");
    namesMatchingdown.push_back("_ttJets_matchingup");
	
    
    for(int ib=0; ib<=2; ib++){
        for(int i=0; i<Nscales; i++){
            TString variable_name=Form("EventCategories_1Mu_1El_%dbtag",ib)+scales[i];
            combine_histograms(variable_name,namesMatchingdown,"_allMC_matchingdown");
			//normalize_systematic_categories(variable_name, "_allMC_matchingdown"); 

        }
    }
	
    
}

void normalize_systematic_categories(TString variable_name, TString sys){
	//cout << variable_name << " sys: " << sys << endl;
	TString data=variable_name+"_singleMu"; 
	TString MC=variable_name+sys; 

	int nJ=1; 
	int ist=0; 
	double scale[3]; 
	
	double Ndata=0; 
	double NMC=0; 
	for(int i=1; i<=hName[data]->GetNbinsX(); i++){
		Ndata+=hName[data]->GetBinContent(i); 
		NMC+=hName[MC]->GetBinContent(i); 

		if(nJ==3) scale[ist]=Ndata/NMC; 
		
		nJ++; 
		
		if(nJ>nJetmax){
			nJ=1; 
			ist++; 
			Ndata=0; 
			NMC=0; 
		}
	}
	
	cout << "scale: " << scale[0] << " " << scale[1] << " " << scale[2] << endl; 
	
	ist=0;
	nJ=1; 
	for(int i=1; i<=hName[data]->GetNbinsX(); i++){
		
		hName[MC]->SetBinContent(i,hName[MC]->GetBinContent(i)*scale[ist]); 
		nJ++; 
		
		if(nJ>nJetmax){
			nJ=1; 
			ist++; 
		}
	}
	
	
}

void plot_categories(){
    overlay_stlabels=true;
    TString scales[]={"_noSF","_SF"};
    
    int Nscales=sizeof(scales)/sizeof(TString);
    
    std::vector<TString> signal_names;
    //signal_names.push_back("_stealth_300_200");
    signal_names.push_back("_stealth_600_300");
	//signal_names.push_back("_UDD300"); 
    
    vector<TString> names;
    names.push_back("_diboson");
    names.push_back("_dy");
    names.push_back("_singleTop");
    names.push_back("_ttbar");
    
    plot_systematics=true;
	
	vector<TString> legend_names;
    legend_names.push_back("Diboson");
	legend_names.push_back("Drell-Yan"); 
	legend_names.push_back("Single t");
	legend_names.push_back("t#bar{t}");
	
	legend_names.push_back("M_{#tilde{q}} = 600 GeV");
	//legend_names.push_back("UDD 212, M_{#tilde{t}}=300 GeV"); 

	
    compute_systematic(names);
    compute_systematic_MC("EventCategories_1Mu_1El",names);
 
    
    for(int nJ=4; nJ<=nJetmax; nJ++){
        compute_systematic_MC(Form("stInclusive_Inclusive_nJets%d_1Mu_1El",nJ),names);
    }
    
    plot_systematics=false;
    scaleUpDown=true; 
    for(int ib=0; ib<=2; ib++){
        for(int i=0; i<Nscales; i++){
			plot_systematics=true;
            scaleUpDown=true;
           // notblind=false;
           // blind_categoryplots(ib, scales[i]);

			if(scales[i].Contains("_noSF") || scales[i].Contains("_TopCor")) {
                plot_systematics=false;
                scaleUpDown=false;
            }
            TString variable_name=Form("EventCategories_1Mu_1El_%dbtag",ib)+scales[i];

            compute_ratio(variable_name, "_allMC");
            compute_ratio(variable_name, "_allMC_scaleup");
			compute_ratio(variable_name, "_allMC_scaledown");

            fill_stack_eventcategories(variable_name,"_singleMu",signal_names,names,legend_names,"_allMC");
        }
        
    }
	scaleUpDown=false; 
    plot_systematics=false;
    
    //make SS plots
    
    
    for(int ib=0; ib<=2; ib++){
            TString variable_name=Form("EventCategories_SS_1Mu_1El_%dbtag_SF",ib);
            compute_ratio(variable_name, "_allMC");

            fill_stack_eventcategories(variable_name,"_singleMu",signal_names,names,legend_names,"_allMC");
    }
    
    cout << "loose btag plots: " << endl;
    //loose b-tags
    
    
    for(int ib=0; ib<=2; ib++){
        TString variable_name=Form("EventCategories_1Mu_1El_%dbtag_loose_noSF",ib);
        //cout << "bin content before blinding: " << hName[variable_name+"_allMC"]->GetBinContent(category_bin(4,2)) << endl;

        blind_categoryplots(ib,"_loose_noSF");
        //cout << "bin content: " << hName[variable_name+"_allMC"]->GetBinContent(category_bin(5,2)) << endl;
        
        compute_ratio(variable_name, "_allMC");
        //cout << "bin content after ratio: " << hName[variable_name+"_allMC"]->GetBinContent(category_bin(5,2)) << endl;

        
        fill_stack_eventcategories(variable_name,"_singleMu",signal_names,names,legend_names,"_allMC");
    }
    //SS background
    for(int ib=0; ib<=2; ib++){
        TString variable_name=Form("EventCategories_SS_1Mu_1El_%dbtag_loose_noSF",ib);
        compute_ratio(variable_name, "_allMC");
        
        fill_stack_eventcategories(variable_name,"_singleMu",signal_names,names,legend_names,"_allMC");
    }
    
    cout << "Top Cor plots: " << endl;

    for(int ib=0; ib<=2; ib++){
        TString variable_name=Form("EventCategories_1Mu_1El_%dbtag_TopCor",ib);
        compute_ratio(variable_name, "_allMC");
        
        fill_stack_eventcategories(variable_name,"_singleMu",signal_names,names,legend_names,"_allMC");
    }
    
    
  
    
}

void compute_systematic_MC(TString variable, std::vector<TString> names){
    cout << "Compute systematic using MC only: " << endl;
    //cout << variable << endl;
    TString tag[]={""};
    int Ntag=sizeof(tag)/sizeof(TString);
    
    int Nnames=names.size();
    if(Nnames==0){
        cout << "No MC samples given: " << endl;
        return;
    }
    
    
    TString signPM[]={"P","M"};
    int Nsign=sizeof(signPM)/sizeof(TString);
    
    
    for(int itag=0; itag<Ntag;itag++){
        for(int ib=0; ib<=2; ib++){
            TString variable_name=variable+Form("_%dbtag",ib)+tag[itag] +"_SF";
            double BW=hName[variable_name+"_allMC"]->GetBinWidth(1)/2;
            // cout << variable_name << endl;
            CreateGraph(hName[variable_name+"_allMC"],variable_name+"_allMC_MCsystematic");
            CreateGraph(hName[variable_name+"_allMC"],variable_name+"_allMC_RMCsystematic");
			
            //make graphs
            for(int ibin=1; ibin<=hName[variable_name+names.at(0)]->GetNbinsX(); ibin++){
                double Nsys[2];
                for(int pm=0; pm<Nsign; pm++){
                    double Nisys=0;
                    for(int isample=0; isample<Nnames; isample++){
                        TString bc=signPM[pm]+"bc";
                        TString l=signPM[pm]+"l";
                        
                        double N=hName[variable_name+names.at(isample)]->GetBinContent(ibin);
                        double Nbc=hName[variable_name+bc+names.at(isample)]->GetBinContent(ibin);
                        double Nl=hName[variable_name+l+names.at(isample)]->GetBinContent(ibin);
                        
                        Nisys+=TMath::Power(Nbc-N,2)+TMath::Power(Nl-N,2);
                        
                    }//sample loop
                    //cout << "btag sys done: " << endl;
                    double N=grName[variable_name+"_allMC_MCsystematic"]->GetY()[ibin-1]; //hist index starts at 1, graph at 0

                    double sysScaleU=hName[variable_name+"_allMC_scaleup"]->GetBinContent(ibin)-N;
                    double sysScaleD=hName[variable_name+"_allMC_scaledown"]->GetBinContent(ibin)-N;
                    //cout << "scaleup/down accessed. " << variable_name+"_allMC_scaleup" << endl;
                    double sys=TMath::Max(sysScaleU*sysScaleU,sysScaleD*sysScaleD); //max of sysU*sysU or down
                    
                    Nisys+=sys;
                    
                    Nisys=TMath::Sqrt(Nisys);
                    Nsys[pm]=Nisys;
                    //cout << " scaleup/down done " << endl;
					
                }//+/- systematic
                //set y errors
                double y=grName[variable_name+"_allMC_MCsystematic"]->GetY()[ibin-1]; //hist index starts at 1, graph at 0
                grName[variable_name+"_allMC_MCsystematic"]->GetEXhigh()[ibin-1]=BW;
                grName[variable_name+"_allMC_MCsystematic"]->GetEXlow()[ibin-1]=BW;
				
                grName[variable_name+"_allMC_RMCsystematic"]->GetEXhigh()[ibin-1]=BW;
                grName[variable_name+"_allMC_RMCsystematic"]->GetEXlow()[ibin-1]=BW;
                
                grName[variable_name+"_allMC_MCsystematic"]->GetEYhigh()[ibin-1]=Nsys[0];
                grName[variable_name+"_allMC_MCsystematic"]->GetEYlow()[ibin-1]=Nsys[1];
                
                grName[variable_name+"_allMC_RMCsystematic"]->GetY()[ibin-1]=1;
                if(y>0){
                    grName[variable_name+"_allMC_RMCsystematic"]->GetEYhigh()[ibin-1]=Nsys[0]/y;
                    grName[variable_name+"_allMC_RMCsystematic"]->GetEYlow()[ibin-1]=Nsys[1]/y;
                }
                else {
                    grName[variable_name+"_allMC_RMCsystematic"]->GetEYhigh()[ibin-1]=0;
                    grName[variable_name+"_allMC_RMCsystematic"]->GetEYlow()[ibin-1]=0;
                }
				//cout << "bin: " << ibin << " " << y << " " << Nsys[0] << " " << Nsys[0]/y << endl;
            }//loop over bins
        }//btag loop
    }//tag working point
    
    
    cout << "finished systematic calcs" << endl;
}

void compute_systematic(std::vector<TString> names){
    cout << "Compute systematic: " << endl;
    TString tag[]={""};
    int Ntag=sizeof(tag)/sizeof(TString);
    
    int Nnames=names.size();
    if(Nnames==0){
        cout << "No MC samples given: " << endl;
        return;
    }
    
    TString signPM[]={"P","M"};
    int Nsign=sizeof(signPM)/sizeof(TString);
    
    
    for(int itag=0; itag<Ntag;itag++){
        for(int ib=0; ib<=2; ib++){
            TString variable_name=Form("EventCategories_1Mu_1El_%dbtag",ib)+tag[itag] +"_SF";
            CreateGraph(hName[variable_name+"_allMC"],variable_name+"_allMC_systematic");
            CreateGraph(hName[variable_name+"_allMC"],variable_name+"_allMC_Rsystematic");
			
            //make graphs
            for(int ibin=1; ibin<=hName[variable_name+names.at(0)]->GetNbinsX(); ibin++){
                double Nsys[2];
                for(int pm=0; pm<Nsign; pm++){
                    double Nisys=0;
                    for(int isample=0; isample<Nnames; isample++){
                        TString bc=signPM[pm]+"bc";
                        TString l=signPM[pm]+"l";
                        
                        double N=hName[variable_name+names.at(isample)]->GetBinContent(ibin);
                        double Nbc=hName[variable_name+bc+names.at(isample)]->GetBinContent(ibin);
                        double Nl=hName[variable_name+l+names.at(isample)]->GetBinContent(ibin);
                        
                        Nisys+=TMath::Power(Nbc-N,2)+TMath::Power(Nl-N,2);
                        
                    }//sample loop
                    Nisys=TMath::Sqrt(Nisys);
                    Nsys[pm]=Nisys;
					
                }//+/- systematic
                //set y errors
                double y=grName[variable_name+"_allMC_systematic"]->GetY()[ibin-1]; //hist index starts at 1, graph at 0
                grName[variable_name+"_allMC_systematic"]->GetEXhigh()[ibin-1]=0.5;
                grName[variable_name+"_allMC_systematic"]->GetEXlow()[ibin-1]=0.5;
				
                grName[variable_name+"_allMC_Rsystematic"]->GetEXhigh()[ibin-1]=0.5;
                grName[variable_name+"_allMC_Rsystematic"]->GetEXlow()[ibin-1]=0.5;
                
                grName[variable_name+"_allMC_systematic"]->GetEYhigh()[ibin-1]=Nsys[0];
                grName[variable_name+"_allMC_systematic"]->GetEYlow()[ibin-1]=Nsys[1];
                
                grName[variable_name+"_allMC_Rsystematic"]->GetY()[ibin-1]=1;
                grName[variable_name+"_allMC_Rsystematic"]->GetEYhigh()[ibin-1]=Nsys[0]/y;
                grName[variable_name+"_allMC_Rsystematic"]->GetEYlow()[ibin-1]=Nsys[1]/y;
                
				//cout << "bin: " << ibin << " " << y << " " << y+Nsys[0] << " " << Nsys[0]/y << endl;
            }//loop over bins
        }//btag loop
    }//tag working point
}

void CreateGraph(TH1F *h,TString graphName){
    TGraphAsymmErrors *gr = new TGraphAsymmErrors(h);
    gr->SetName(graphName);
    gr->SetFillStyle(3244);
    gr->SetFillColor(38);
    
    grName[gr->GetName()]=gr;
}

void blind_categoryplots(int nb, TString scale){
    TString name=Form("EventCategories_1Mu_1El_%dbtag",nb)+scale;
    name=name+"_singleMu";
    int nJ=1;
    //cout << "blinding histogram: " << name << endl;
    for(int nJ=2; nJ<=nJetmax; nJ++){
        for(int ist=0; ist<3; ist++){
            int ibin=category_bin(nJ,ist);
            if(nb==0 && nJ>=4)hName[name]->SetBinContent(ibin,0);
        }
    }
	
}


void compute_ratio(TString variable_name, TString MC_name, int N, double *x){
	//rebin 
    bool print=false;
    if(print)cout << "Compute Rebinned Ratio: " << variable_name << endl;
    clone_histogram(variable_name+"_singleMu",variable_name+MC_name+"_RatioDataMC");
    for (int i=1; i<=hName[variable_name+MC_name+"_RatioDataMC"]->GetNbinsX(); i++){
        double R=hName[variable_name+MC_name+"_RatioDataMC"]->GetBinContent(i);
        if(R<0.01) hName[variable_name+MC_name+"_RatioDataMC"]->SetBinContent(i,0);
        if(R<0.01) hName[variable_name+MC_name+"_RatioDataMC"]->SetBinError(i,0);
		
    }
	
	hName[variable_name+MC_name+"_RatioDataMC"]=(TH1F*)hName[variable_name+MC_name+"_RatioDataMC"]->Rebin(N,"",x); 
	TH1F *h=(TH1F*)hName[variable_name+MC_name]->Rebin(N,"h_tmp",x); // temporary rebinned histogram 
		
	if(print)cout << "Nbins Num: "<< hName[variable_name+MC_name+"_RatioDataMC"]->GetNbinsX() << endl;
	if(print)cout << "Nbins Den: " << h->GetNbinsX() << endl;
	
    hName[variable_name+MC_name+"_RatioDataMC"]->GetYaxis()->SetTitle("Data/MC");
    hName[variable_name+MC_name+"_RatioDataMC"]->SetAxisRange(0.0,2.75,"Y");
    hName[variable_name+MC_name+"_RatioDataMC"]->SetMarkerStyle(20);
    hName[variable_name+MC_name+"_RatioDataMC"]->SetMarkerSize(0.5);
    hName[variable_name+MC_name+"_RatioDataMC"]->Divide(h);
    
	// for(int i=1; i<=hName[variable_name+MC_name+"_RatioDataMC"]->GetNbinsX();i++) cout << "R: " << hName[variable_name+MC_name+"_RatioDataMC"]->GetBinContent(i) << endl;
	
	delete h; 
}



void compute_ratio(TString variable_name, TString MC_name){
    bool print=false;
    if(print)cout << "Compute Ratio: " << variable_name << endl;
    clone_histogram(variable_name+"_singleMu",variable_name+MC_name+"_RatioDataMC");
    for (int i=1; i<=hName[variable_name+MC_name+"_RatioDataMC"]->GetNbinsX(); i++){
        double R=hName[variable_name+MC_name+"_RatioDataMC"]->GetBinContent(i);
        if(R<0.01) {
            hName[variable_name+MC_name+"_RatioDataMC"]->SetBinContent(i,0);
            hName[variable_name+MC_name+"_RatioDataMC"]->SetBinError(i,0);
            //hName[variable_name+MC_name]->SetBinContent(i,0);
        }
        hName[variable_name+MC_name]->SetBinError(i,0);
        
    }
    
    hName[variable_name+MC_name+"_RatioDataMC"]->GetYaxis()->SetTitle("Data/MC");
    if(MC_name=="_allExp") hName[variable_name+MC_name+"_RatioDataMC"]->GetYaxis()->SetTitle("Data/Exp");
    hName[variable_name+MC_name+"_RatioDataMC"]->SetAxisRange(0.0,2.75,"Y");
    hName[variable_name+MC_name+"_RatioDataMC"]->SetMarkerStyle(20);
    hName[variable_name+MC_name+"_RatioDataMC"]->SetMarkerSize(0.5);
    hName[variable_name+MC_name+"_RatioDataMC"]->Divide(hName[variable_name+MC_name]);
    
    for (int i=1; i<=hName[variable_name+MC_name+"_RatioDataMC"]->GetNbinsX(); i++){
        if(hName[variable_name+MC_name+"_RatioDataMC"]->GetBinContent(i)<=0){
            hName[variable_name+MC_name+"_RatioDataMC"]->SetBinContent(i,-20);
            cout << "set bin content 0: " << i << endl;
        }
        
        if(variable_name.Contains("EventCategories_1Mu_1El_2btag_SF") && i==1){
            cout << "N: " << hName[variable_name+"_singleMu"]->GetBinContent(i) << " D: " <<hName[variable_name+MC_name]->GetBinContent(i) << endl;
            cout << "Ratio: " <<hName[variable_name+MC_name+"_RatioDataMC"]->GetBinContent(i) << " +/- " << hName[variable_name+MC_name+"_RatioDataMC"]->GetBinError(i) << endl;
        }
    }
    
    
   // for(int i=1; i<=hName[variable_name+MC_name+"_RatioDataMC"]->GetNbinsX();i++) cout << "R: " << hName[variable_name+MC_name+"_RatioDataMC"]->GetBinContent(i) << endl;
}

double getZBI(double sig, double bkg, double bkge){
    double obs=sig+bkg; //obs = signal + bkg
// obs = expected number of _observed_ events (i.e. signal+background)
// back = expected number of background events
// backe = uncertainty on the expected number of background events
    if(bkg<=0) {
        //cout << "bkg must be greater than 0." << endl;
        return 0;
    }
    if(obs==0) return 0;

    double tau = bkg/(bkge*bkge);
    double n_off = tau*bkg;
    double P_BI = TMath::BetaIncomplete(1./(1.+tau), obs, n_off+1);
    double ZBI=sqrt(2.0)*TMath::ErfcInverse(2*P_BI);
    return ZBI;
}

int category_bin(int nJ, int ist){
    int iBin=-1;
    if(ist==0) iBin=nJ;
    if(ist==1) iBin=nJ+nJetmax;
    if(ist==2) iBin=nJ+2*nJetmax;
    return iBin;
}

TString print_unc(int bkg, int bkg_tot, double unc){
    TString out_string="";
    if(bkg>=0){
        for(int ibkg=0; ibkg<bkg_tot; ibkg++){
            if(ibkg!=bkg)out_string=out_string+"\t-";
            else out_string=out_string+Form("\t%.4f",unc);
        }
    }
    else{
        for(int ibkg=0; ibkg<bkg_tot; ibkg++){
            out_string=out_string+Form("\t%.3f",unc);
        }

    }
    
    return out_string;
}

void table(int nbtag,int ist){
    TString variable_name=Form("EventCategories_1Mu_1El_%dbtag_SF",nbtag);
    TString variable_namectrl=Form("EventCategories_1Mu_1El_%dbtag_SF",2);

    double qcd_per=0.5;
    if(ist==2) qcd_per=1.2;
    
    double R[]={0,0,1.03,1.096,1.11,1.23,1.56,1.56,1.56}; //dy correction 0-8 jets

    
    ofstream out;
    out.open(Form("table_EventCounts_%dbtag_%dst.txt",nbtag,ist));
    out << "\\begin{tabular}{|c|c|c|c|c|c|c|c|} \\hline" << endl;
    out << "n-jets & $M_{\\tilde{q}}$=600 GeV & observed & expected & top & Drell-Yan & diboson & Non-prompt \\\\ \\hline " << endl;
    
    
    double ExpTot[4][3];
    double ObsTot[4][3];
    double ExpUTot[4][3];
   // double Sig[4];

    for(int nJ=4; nJ<=nJetmax; nJ++){
        int ibin=category_bin(nJ,ist);
        int ibincontrol=category_bin(nJ,0);
        double exp=0;
        
        exp+=hName[variable_name+"_ttbarExp"]->GetBinContent(ibin);
        exp+=hName[variable_name+"_dy_datadriven"]->GetBinContent(ibin);
        exp+=hName[variable_name+"_diboson"]->GetBinContent(ibin);
        exp+=hName[variable_name+"_QCD"]->GetBinContent(ibin);

        ExpTot[nJ-4][ist]=exp;
        
        double deltaDY=TMath::Abs((hName[variable_name+"_dy_datadriven"]->GetBinContent(ibin)-hName[variable_name+"_dy"]->GetBinContent(ibin))/2);
        double deltaDB=hName[variable_name+"_diboson"]->GetBinContent(ibin)*(R[nJ]-1)/2;

        double deltaQCD=hName[variable_name+"_QCD"]->GetBinContent(ibin)*qcd_per;
        double datactrl=hName[variable_namectrl+"_singleMu"]->GetBinContent(ibincontrol);

        double top=hName[variable_name+"_ttbarExp"]->GetBinContent(ibin); //nominal top
        double topbcP=hName[variable_name+"Pbc_ttbarExp"]->GetBinContent(ibin);//top bc+
        double topbcM=hName[variable_name+"Mbc_ttbarExp"]->GetBinContent(ibin);//top bc-

        double toplP=hName[variable_name+"Pl_ttbarExp"]->GetBinContent(ibin);
        double toplM=hName[variable_name+"Ml_ttbarExp"]->GetBinContent(ibin);
        
        double bP=TMath::Sqrt(TMath::Power(top-topbcP,2)+TMath::Power(top-toplP,2));
        double bM=TMath::Sqrt(TMath::Power(top-topbcM,2)+TMath::Power(top-toplM,2));

        double btagSys=TMath::Max(bP,bM);
        
        double deltatop=hName[variable_name+"_ttbarExp"]->GetBinContent(ibin)*TMath::Sqrt(0.15*0.15+0.1*0.1+TMath::Power(poisson(datactrl)/datactrl,2));
        deltatop=TMath::Sqrt(deltatop*deltatop+btagSys*btagSys);
        
        double expU=TMath::Max(grName[variable_name+"_allExp_systematic"]->GetEYhigh()[ibin-1],grName[variable_name+"_allExp_systematic"]->GetEYlow()[ibin-1]);
        
        ExpUTot[nJ-4][ist]=expU;
        
        if(nJ<nJetmax)out << nJ << " & ";
        else out << "$\\geq$ " << nJ << " & ";

        out << Form("%.1f & ",hName[variable_name+"_stealth_600_300"]->GetBinContent(ibin));
        
        //Sig=hName[variable_name+"_stealth_600_300"]->GetBinContent(ibin);
        
        ObsTot[nJ-4][ist]=hName[variable_name+"_singleMu"]->GetBinContent(ibin);
        
        out << Form("%.0f & ",hName[variable_name+"_singleMu"]->GetBinContent(ibin));
        if(ist==2){
            out << Form("%.2f ",exp);
            out << Form("$\\pm$ %.2f & ",expU);
            
            out << Form("%.2f ",hName[variable_name+"_ttbarExp"]->GetBinContent(ibin));
            out << Form("$\\pm$ %.2f & ",deltatop);
            
            out << Form("%.2f ",hName[variable_name+"_dy_datadriven"]->GetBinContent(ibin));
            out << Form("$\\pm$ %.2f & ",deltaDY);
            
            out << Form("%.2f ",hName[variable_name+"_diboson"]->GetBinContent(ibin));
            out << Form("$\\pm$ %.2f & ",deltaDB);
            
            out << Form("%.2f ",hName[variable_name+"_QCD"]->GetBinContent(ibin));
            out << Form("$\\pm$ %.2f ",deltaQCD) << " \\\\ \\hline" << endl;
        }
        else {
            out << Form("%.1f ",exp);
            out << Form("$\\pm$ %.1f & ",expU);
            
            out << Form("%.1f ",hName[variable_name+"_ttbarExp"]->GetBinContent(ibin));
            out << Form("$\\pm$ %.1f & ",deltatop);
            
            out << Form("%.1f ",hName[variable_name+"_dy_datadriven"]->GetBinContent(ibin));
            out << Form("$\\pm$ %.1f & ",deltaDY);
            
            out << Form("%.1f ",hName[variable_name+"_diboson"]->GetBinContent(ibin));
            out << Form("$\\pm$ %.1f & ",deltaDB);
            
            out << Form("%.1f ",hName[variable_name+"_QCD"]->GetBinContent(ibin));
            out << Form("$\\pm$ %.1f ",deltaQCD) << " \\\\ \\hline" << endl;
        }

    }
    out << "\\end{tabular}" << endl;
    
    out.close();
    
    out.open(Form("table_compressed_EventCounts_%dbtag_%dst.txt",nbtag,ist));
    /*out << "\\begin{tabular}{cccccc}" << endl;
    out << "\\hline\\hline" << endl;
    out << "& & $N_{jets} = 4$ & $N_{jets}$ = 5$ & $N_{jets} = 6$ & $N_{jets} \geq 7$ \\\\ " << endl;
    out << "\\hline" << endl;
    */
    out << " & Observed events & ";
    for(int i=0; i<4; i++) {
        out << Form("%.0f",ObsTot[i][ist]) << " ";
        if(i<3) out << "& ";
    }
    out <<"\\tabularnewline" << endl;
    
    out << " & Total background & ";
    for(int i=0; i<4; i++) {
        //if(ist==2) out << Form("%.2f",ExpTot[i][ist]) << " $\\pm$ " << Form("%.2f",ExpUTot[i][ist]) << " ";
        //else
        out << Form("%.1f",ExpTot[i][ist]) << " $\\pm$ " << Form("%.1f",ExpUTot[i][ist]) << " ";
        if(i<3) out << "& ";
    }
    out << "\\tabularnewline" << endl;

    out << " & Signal & ";
    
    for(int nJ=4; nJ<=7; nJ++){
        double acc=hName["Acceptance_uW_stealth_600_300"]->GetBinContent(category_bin(nJ,ist));
        double Ngen=hName["h_cutflow_table_stealth_600_300"]->GetBinContent(1);
        
        double finiteMC=TMath::Sqrt(acc*(1-acc)/Ngen)/acc;
        double jec=0.05;
        double lumi=0.026;
        
        double sigU=TMath::Sqrt(finiteMC*finiteMC+jec*jec+lumi*lumi);
        
        double sig=hName[variable_name+"_stealth_600_300"]->GetBinContent(category_bin(nJ,ist));
        out << Form("%.1f",sig) << " $\\pm$ " <<  Form("%.1f",sig*sigU) << " ";
        if(nJ<7) out << "& ";
    }
    out << "\\tabularnewline" << endl;
    out.close();
    
    
}

void table_simulation(int nbtag, int ist){
    ofstream out;

    TString variable_name=Form("EventCategories_1Mu_1El_%dbtag_SF",nbtag);
    out.open(Form("table_compressed_EventCounts_%dbtag_%dst.txt",nbtag,ist));
    /*out << "\\begin{tabular}{cccccc}" << endl;
     out << "\\hline\\hline" << endl;
     out << "& & $N_{jets} = 4$ & $N_{jets}$ = 5$ & $N_{jets} = 6$ & $N_{jets} \geq 7$ \\\\ " << endl;
     out << "\\hline" << endl;
     */
    out << " & Observed events & ";
    for(int nJ=4; nJ<=7; nJ++) {
        double obs = hName[variable_name+"_singleMu"]->GetBinContent(category_bin(nJ,ist));
        out << Form("%.0f",obs) << " ";
        if(nJ<7) out << "& ";
    }
    out <<"\\tabularnewline" << endl;
    
    out << " & Total background & ";
    for(int nJ=4; nJ<=7; nJ++) {
        double exp = hName[variable_name+"_allMC"]->GetBinContent(category_bin(nJ,ist));
        double expU = TMath::Max(grName[variable_name+"_allMC_MCsystematic"]->GetEYhigh()[category_bin(nJ,ist)-1],grName[variable_name+"_allMC_MCsystematic"]->GetEYlow()[category_bin(nJ,ist)-1]);

        //if(ist==2) out << Form("%.2f",exp) << " $\\pm$ " << Form("%.2f",expU) << " ";
        out << Form("%.1f",exp) << " $\\pm$ " << Form("%.1f",expU) << " ";
        if(nJ<7) out << "& ";
    }
    out << "\\tabularnewline" << endl;
    out.close();
}

void acceptance_table(){
    ofstream out;
    out.open("Acceptance_Table.tex");

    int ist=0;
    
    std::vector<TString> signal_names;
    signal_names.push_back("_stealth_300_200");
    signal_names.push_back("_stealth_400_200");
    signal_names.push_back("_stealth_500_300");
    signal_names.push_back("_stealth_600_300");
    signal_names.push_back("_stealth_700_400");
    signal_names.push_back("_stealth_800_400");
    signal_names.push_back("_stealth_900_500");
    int mass=300;
    for(int imass=0; imass<signal_names.size(); imass++){
        TString Sig=signal_names.at(imass);
        if(Sig=="_stealth_300_200")ist=0;
        if(Sig=="_stealth_400_200")ist=1;
        if(Sig=="_stealth_500_300")ist=1;
        if(Sig=="_stealth_600_300")ist=2;
        if(Sig=="_stealth_700_400")ist=2;
        if(Sig=="_stealth_800_400")ist=2;
        if(Sig=="_stealth_900_500")ist=2;
        out << mass << " & ";
        for(int nJ=4; nJ<=7; nJ++){
            int ibin=category_bin(nJ, ist);
            double acc=hName["Acceptance_uW"+Sig]->GetBinContent(ibin);
            
            double Ngen=hName["h_cutflow_table"+Sig]->GetBinContent(1);
            
            double finiteMC=TMath::Sqrt(acc*(1-acc)/Ngen)/acc;
            double jec=0.05;
            double lumi=0.026;
            double sigU=TMath::Sqrt(finiteMC*finiteMC+jec*jec+lumi*lumi);
            
            if(nJ<7)out << Form("%.2f",acc*100) << " $\\pm$ " << Form("%.2f",acc*sigU*100) << " & ";
            else out << Form("%.2f",acc*100) << " $\\pm$ " << Form("%.2f",acc*sigU*100) << " ";
        }
        out << "\\\\" << endl;
        mass+=100;
    }
}

void print_loose_events(){
    TString variable_name="EventCategories_1Mu_1El_0btag_loose_noSF";
    cout << "loose working point: " << endl;
    cout << "bkg: ";
    for(int nJ=4; nJ<=7; nJ++){
        int ibin=category_bin(nJ,2);
        cout << hName[variable_name+"_allMC"]->GetBinContent(ibin) << " ";
    }
    
    cout << endl;
    cout << "sig: ";
    for(int nJ=4; nJ<=7; nJ++){
        int ibin=category_bin(nJ,2);
        cout << hName[variable_name+"_stealth_600_300"]->GetBinContent(ibin) << " ";
        
    }
    cout << endl;
    cout << "medium working point: " << endl;

    variable_name="EventCategories_1Mu_1El_0btag_SF";

    cout << "bkg: ";
    for(int nJ=4; nJ<=7; nJ++){
        int ibin=category_bin(nJ,2);
        cout << hName[variable_name+"_allMC"]->GetBinContent(ibin) << " ";
    }
    cout << endl;
    
    cout << "sig: ";
    for(int nJ=4; nJ<=7; nJ++){
        int ibin=category_bin(nJ,2);
        cout << hName[variable_name+"_stealth_600_300"]->GetBinContent(ibin) << " ";
    }
    
    cout << endl;
}

void datacards_new(TString variable_name, TString variable_name_cntrl, TString signal_names){
    cout << "datacards new: " << endl;
    int nJ=1;
	int ist=0;
    
    int stmin=0;
    if(signal_names=="_stealth_300_200")stmin=0;
    if(signal_names=="_stealth_400_200")stmin=1;
    if(signal_names=="_stealth_500_300")stmin=1;
    if(signal_names=="_stealth_600_300")stmin=2;
    if(signal_names=="_stealth_700_400")stmin=2;
    if(signal_names=="_stealth_800_400")stmin=2;
    if(signal_names=="_stealth_900_500")stmin=2;

    TString data_name="_singleMu";
    
    vector<TString> process_name;
    vector<TString> process_title;
    
    process_name.push_back(signal_names);
    process_name.push_back("_ttbarExp");
    process_name.push_back("_dy_datadriven");
    process_name.push_back("_diboson");
    process_name.push_back("_QCD");

    process_title.push_back("sig");
    process_title.push_back("top");
    process_title.push_back("dy");
    process_title.push_back("db");
    process_title.push_back("np");

    
    int nbkg = 4; //top,dy,db,np
    int njetmin=4;
    int nch=1;
    
    int Nnuisance=15;
    double R[]={0,0,1.03,1.096,1.11,1.23,1.56,1.56,1.56}; //dy correction 0-8 jets
    
    
    for(int ich=njetmin; ich<=nJetmax; ich++){
        int ibin=category_bin(ich,stmin);
        int ibin_st0=category_bin(ich,0);
        
        TString datacardName=Form("datacards/datacard_nJets%d",ich)+signal_names+".txt";
        ofstream out;
        out.open(datacardName);
        
        out << Form("jmax %d number of channels ",nch) << endl;
        out << Form("jmax %d number of backgrounds",nbkg) << endl;
        out << Form("kmax %d number of nuisance parameters (sources of systematical uncertainties)",Nnuisance) << endl;
        out << Form("bin\tb%dj\t",ich);

        out << endl;
        //out <<"bin b1 b2 b3 b4 b5" << endl;
        out << " ------------"<< endl;
        
        out << "observation " ;
        double obs=hName[variable_name+data_name]->GetBinContent(ibin);
        out << Form("%.0f\t",obs);
        
        out << endl;
        out << " ------------"<< endl;
        out << "bin\t";
        
        for(int ip=0; ip<process_title.size(); ip++) {
            out << Form("b%dj\t",ich);
        }//loop over process
        
        out << endl;
        
        out << "process\t";
        
        for(int ip=0; ip<process_title.size(); ip++) {
            out << process_title.at(ip)<< "\t";
            //<<Form("b%dj",ich)<< "\t";
        }//loop over process
        
        out << endl;
        
        out << "process\t";
        
        for(int ip=0; ip<process_title.size(); ip++) {
            out << ip << "\t";
        }//loop over process
        
        out << endl;
        
        out << "rate\t";
        
        for(int ip=0; ip<process_name.size(); ip++) {
            //TString VN=variable_name;
            //if(process_name.at(ip)=="_ttbarExp")variable_name="EventCategories_1Mu_1El_0btag_TopCor";
            out << hName[variable_name+process_name.at(ip)]->GetBinContent(ibin) << "\t";
            //if(process_name.at(ip)=="_ttbarExp")variable_name=VN;

        }//loop over process
        
        out << endl;
        out << " ------------"<< endl;
        
        double acc=hName["Acceptance_uW"+signal_names]->GetBinContent(ibin);
        double Ngen=hName["h_cutflow_table"+signal_names]->GetBinContent(1);
        
        double finiteMC=1+TMath::Sqrt(acc*(1-acc)/Ngen)/acc;
        double jec=1.05;
        double lumi=1.026;
        
        double db=(R[ich]-1)/2;
        db=1+TMath::Sqrt(db*db+0.3*0.3);
        
        double topshape_N=hName[variable_name_cntrl+data_name]->GetBinContent(ibin_st0);
        double topshape_alpha=0;
        if(topshape_N>0)topshape_alpha=hName[variable_name+"_ttbarExp"]->GetBinContent(ibin)/topshape_N;
        //if(topshape_N>0)topshape_alpha=hName["EventCategories_1Mu_1El_0btag_TopCor_ttbarExp"]->GetBinContent(ibin)/topshape_N;

        double topNorm[]={1.15,1.15,1.15};//top Normalization in 3 bins
        
        double top=hName[variable_name+"_ttbarExp"]->GetBinContent(ibin);
        //double top=hName["EventCategories_1Mu_1El_0btag_TopCor_ttbarExp"]->GetBinContent(ibin);
        double topU=hName[variable_name+"_ttbarExp_scaleup"]->GetBinContent(ibin);
        double topD=hName[variable_name+"_ttbarExp_scaledown"]->GetBinContent(ibin);

        
        double bcP=hName[variable_name+"Pbc_ttbarExp"]->GetBinContent(ibin);
        double lP=hName[variable_name+"Pl_ttbarExp"]->GetBinContent(ibin);

        double bcM=hName[variable_name+"Mbc_ttbarExp"]->GetBinContent(ibin);
        double lM=hName[variable_name+"Ml_ttbarExp"]->GetBinContent(ibin);

        //cout << "top: " << top << " bcP " << bcP << " lP: " << lP << endl;
        //cout << "nj: " << ich << " ist: " << ist << endl;
        //cout << "top: " << top << " up " << topU << " down " << topD << endl;
        
        double btagsysP=TMath::Sqrt(TMath::Power(bcP-top,2)+TMath::Power(lP-top,2));
        double btagsysM=TMath::Sqrt(TMath::Power(bcM-top,2)+TMath::Power(lM-top,2));

        double btagsys=1+TMath::Max(btagsysP,btagsysP)/top;
        
        double dy_fit=1+(R[ich]-1)/2;
        //Ndy=Ncontrol*alpha -> alpha=Ndy/Ncontrol
   //     cout << "dy finite sample size: " << endl;
        double ndy_cntrl=hName["EventCategories_2Mu_onZ_0btag_SF_singleMu"]->GetBinContent(category_bin(ich,ist));//get the number of dy events in 2mu on Z control
   //     cout << "compute expected number of dy events:" << endl;
        double dy_alpha=hName[variable_name+"_dy_datadriven"]->GetBinContent(category_bin(ich,ist))/ndy_cntrl;
        //on z: 81-101 GeV
        double np[]={1.5,1.5,2.2};
        //cout << "finite MC bkg: " << endl;
        double finite_bkg[3];
        finite_bkg[0]=hName["Acceptance_nEvents_ttbar"]->GetBinContent(category_bin(ich,ist));
        finite_bkg[0]+=hName["Acceptance_nEvents_singleTop"]->GetBinContent(category_bin(ich,ist));
        finite_bkg[1]+=hName["Acceptance_nEvents_dy"]->GetBinContent(category_bin(ich,ist));
        finite_bkg[2]+=hName["Acceptance_nEvents_diboson"]->GetBinContent(category_bin(ich,ist));
        
        for(int i=0; i<3; i++){
            finite_bkg[i]=1+poisson(finite_bkg[i])/finite_bkg[i];
        }
        
        out << Form("sigfiniteMCb%dj\tlnN\t",ich) << print_unc(0,5,finiteMC) << endl;
        out << Form("sigjec\tlnN\t",ich)  << print_unc(0,5,jec) << endl;
        out << Form("siglumi\tlnN\t",ich)  << print_unc(0,5,lumi) << endl;
        out << Form("top_btag\tlnN\t",ich)  << print_unc(1,5,btagsys) << endl;
        out << Form("top_Norm\tlnN\t",ich)  << print_unc(1,5,topNorm[ist]) << endl;
        out << Form("top_shape\tlnN\t",ich)  << print_unc(1,5,1.1) << endl;
        out << Form("topb%dj\tgmN\t ",ich) << topshape_N << "\t"<< print_unc(1,5,topshape_alpha) << endl;;
       
        out << Form("dy_fitb%dj\tlnN\t",ich) << print_unc(2,5,dy_fit) << endl;
       // out << Form("dy_FiniteControlb%dj\t gmN\t",ich) << ndy_cntrl<<"\t" << print_unc(2,5,dy_alpha) << endl;

        out << Form("dbb%dj\tlnN\t",ich) << print_unc(3,5,db) << endl;
        out << Form("npb%dj\tlnN\t",ich) << print_unc(4,5,np[ist]) << endl;
        
        //finite MC sample size
        out << Form("topfiniteMC%dj\tlnN\t",ich) << print_unc(1,5,finite_bkg[0]) << endl;
        out << Form("dyfiniteMC%dj\tlnN\t",ich) << print_unc(2,5,finite_bkg[1]) << endl;
        out << Form("dbfiniteMC%dj\tlnN\t",ich) << print_unc(3,5,finite_bkg[2]) << endl;

        //lepton systematics
        out << Form("muonTrigID\tlnN\t",ich) << print_unc(0,5,1.01) << endl;
        out << Form("EleID\tlnN\t",ich) << print_unc(0,5,1.03) << endl;

        
    }//loop over channels
    
    
   
}

void datacards(TString variable_name, TString data_name, TString signal_names, TString MC_name){
	int nJ=1;
	int ist=0;
	
	double Nexp[5][3];
	double NexpSys[5][3];
	double Nobs[5][3];
	double Nsig[5][3];
    double NsigE[5][3];
    
	cout << "variable: " << variable_name << endl;
	cout << "signal name: " << signal_names << endl;
	
    double Scale=1;
    TString SN=signal_names;
    if(signal_names=="_stealth_800_400")Scale=0.355;
    
    if(signal_names=="_stealth_800_400")signal_names="_stealth_700_400";
  
     
	for(int ibin=1; ibin<=hName[variable_name+data_name]->GetNbinsX(); ibin++){
		
		bool print=true;
		if(nJ<4) print=false; 
		double expSys=TMath::Max(grName[variable_name+MC_name+"_Rsystematic"]->GetEYhigh()[ibin-1],grName[variable_name+MC_name+"_Rsystematic"]->GetEYlow()[ibin-1]); 
		double exp=hName[variable_name+MC_name]->GetBinContent(ibin); 
		double obs=hName[variable_name+data_name]->GetBinContent(ibin); 
		
		TString sigName=variable_name+signal_names;
		double sig=Scale*hName[sigName]->GetBinContent(ibin);
        //double sigE=1+Scale*hName[sigName]->GetBinError(ibin)/sig;
        double acc=hName["Acceptance_uW"+signal_names]->GetBinContent(nJ);
        double Ngen=hName["h_cutflow_table"+signal_names]->GetBinContent(1);
        double sigE = 1+TMath::Sqrt(acc*(1-acc)/Ngen);
        
        
		if(print){
			//cout << "nJ: " << nJ << "  st: " << ist << endl; 
			//cout << " exp: "<< exp<< endl; 
			Nexp[nJ-4][ist]=exp;
            if(exp==0)Nexp[nJ-4][ist]=0.0001;
			NexpSys[nJ-4][ist]=1+expSys;
			Nobs[nJ-4][ist]=obs; 		
			Nsig[nJ-4][ist]=sig;
            NsigE[nJ-4][ist]=sigE;

		}
		//if(print)cout << endl << endl; 
		nJ++; 
		
		if(nJ>nJetmax){
			nJ=1; 
			ist++; 
		}
	}	
	
	for(int ist=0; ist<3; ist++){
		TString datacardName="datacard"+SN+Form("_%dst.txt",ist);
		ofstream out;
		out.open(datacardName); 
		
		out << "jmax 5 number of channels " << endl; 
		out << "jmax 1 number of backgrounds" << endl; 
		out << "kmax 4 number of nuisance parameters (sources of systematical uncertainties)" << endl;
		out <<"bin b1 b2 b3 b4 b5" << endl; 
		out << " ------------"<< endl; 
		
		out << "observation   " ; 
		for(int nJ=0; nJ<5; nJ++){
			out << Nobs[nJ][ist] << " "; 
		}
		out << endl; 
		out << "------------" << endl; 
		out << "bin b1 b1 b2 b2 b3 b3 b4 b4 b5 b5" << endl; 
		out << "process sig bkg sig bkg sig bkg sig bkg sig bkg " << endl; 
		out << "process  0 1  0 1  0 1  0 1  0 1" << endl; 
		out << "rate "; 
		for(int nJ=0; nJ<5; nJ++){
			out << Nsig[nJ][ist] << " " << Nexp[nJ][ist] << " ";  
		}
		out << endl; 
		out << "------------" << endl; 
		out << "bkg   lnN - " ;
		for(int nJ=0; nJ<5; nJ++){
            out << NexpSys[nJ][ist];
            if(nJ<4) out << " - ";
		}
		out << endl;
        
        out << "finiteMC   lnN " ;
		for(int nJ=0; nJ<5; nJ++){
			out << NsigE[nJ][ist] << " - ";
		}
		out << endl;
        
        out << "jec   lnN " ;
		for(int nJ=0; nJ<5; nJ++){
			out  << "1.05" << " - ";
		}
		out << endl;
        
        out << "lumi   lnN " ;
		for(int nJ=0; nJ<5; nJ++){
			out << "1.022" << " - ";
		}
		out << endl;
        
		out.close();
	}
	//cout << endl << endl; 
	
	/*
	 imax 1 number of channels
	 jmax 1 number of backgrounds
	 kmax 1 number of nuisance parameters (sources of systematical uncertainties)
	 bin b1
	 ------------
	 observation 2
	 ------------
	 bin     b1     b1
	 process sig    bkg
	 process 0      1
	 rate   20    2
	 -------------
	 bkg  lnN  -      1.5
	 
	 */

	
	for(int ist=0; ist<3; ist++){
		//cout << "st: "<< ist << endl; 
		for(int nJ=0; nJ<5; nJ++){
			ofstream out; 
			TString datacardName="datacard"+SN+Form("_%dst_%djets.txt",ist,nJ+4);
			out.open(datacardName); 
	
			out << "jmax 1 number of channels " << endl; 
			out << "jmax 1 number of backgrounds" << endl; 
			out << "kmax 1 number of nuisance parameters (sources of systematical uncertainties)" << endl; 
			out <<"bin b1" << endl; 
			out << " ------------"<< endl; 
			
			out << "observation   " ; 
			out << Nobs[nJ][ist] << " " << endl; 
			out << " ------------"<< endl; 
			out << "bin    b1     b1 " << endl; 
			out << "process sig   bkg " << endl; 
			out << "process 0    1 "<< endl; 
			out << "rate "; 
			
			out << Nsig[nJ][ist] << " " << Nexp[nJ][ist] << " ";  
			
			out << endl; 
			out << "-------------" << endl; 
			out << "bkg lnN -    "; 
			
			out << NexpSys[nJ][ist] << " " << endl; 
			out.close();
			
		}
		
	}
}


void model_independent_datacard(TString variable, int nJet, int nbtag){
	
	TString histNameMC=variable+Form("_nJets%d_1Mu_1El_%dbtag_SF_allMC",nJet,nbtag); 
	//TString histNameMC_scaleUp=variable+Form("_nJets%d_1Mu_1El_%dbtag_SF_allMC_scaleup",nJet,nbtag);
	TString histNameMC_sys=variable+Form("_nJets%d_1Mu_1El_%dbtag_SF_allMC_RMCsystematic",nJet,nbtag);

	
	
	TString histNameData=variable+Form("_nJets%d_1Mu_1El_%dbtag_SF_singleMu",nJet,nbtag); 

	cout << "model independent datacard nJets: " << nJet << "  " << nbtag << endl; 
	
	//cout << histNameMC << endl; 
	//cout << histNameData << endl; 
	
	for(int i=1; i<=hName[histNameMC]->GetNbinsX(); i++){
		float st=hName[histNameMC]->GetBinCenter(i); 
		if(st<300) continue; 
		ofstream out; 
		TString datacardName=Form("modelInd/datacard_%.0fst_",st)+variable+Form("_nJets%d_%dbtags.txt",nJet,nbtag);
		
		out.open(datacardName); 
		
		out << "jmax 1 number of channels " << endl; 
		out << "jmax 1 number of backgrounds" << endl; 
		out << "kmax 1 number of nuisance parameters (sources of systematical uncertainties)" << endl; 
		out <<"bin b1" << endl; 
		out << " ------------"<< endl; 
		
		out << "observation   " ; 
		out << hName[histNameData]->GetBinContent(i) << " " << endl; 
		out << " ------------"<< endl; 
		out << "bin    b1     b1 " << endl; 
		out << "process sig   bkg " << endl; 
		out << "process 0    1 "<< endl; 
		out << "rate "; 
		float exp=hName[histNameMC]->GetBinContent(i);
        if(exp==0) exp==0.00001;
		out << 1 << " " << exp << " ";  
		
		out << endl; 
		out << "-------------" << endl; 
		out << "bkg lnN -    "; 
		
		//float scaleU=TMath::Abs(hName[histNameMC_scaleUp]->GetBinContent(i)-exp)/exp;
		//float scaleD=TMath::Abs(hName[histNameMC_scaleDown]->GetBinContent(i)-exp)/exp;
        
        
		float sys=1+TMath::Max(grName[histNameMC_sys]->GetEYhigh()[i-1],grName[histNameMC_sys]->GetEYlow()[i-1]);

		
		out << sys << " " << endl; 
		out.close();
	}

}

void fill_stack_eventcategories(TString variable_name, TString data_name, std::vector<TString> signal_names, std::vector<TString> names,
								std::vector<TString> legend_names, TString MC_name){
    //gStyle->SetErrorX(0);
	bool print=false;
	if(print) cout << variable_name << "  " << data_name << endl; 
    TString canvas_name=variable_name+"_stack_canvas";
    if(expected_bkg)canvas_name=variable_name+"_expectedbkg_stack_canvas";
    int Cx=800;
    int Cy=800;
    if(overlay_stlabels)Cx=900;
    Cx=800;
    CreateCanvas(canvas_name,"",Cx,Cy);
    TString stack_name=variable_name+"_stack";
    CreateStack(stack_name,"");
	
    int Nbkg=names.size();
    /*
    vector<int> colors;
    if(drawQCD)colors.push_back(5);
    else colors.push_back(7);
    colors.push_back(8);
    if(!expected_bkg){
		colors.push_back(5);
		colors.push_back(419);
    }
    else colors.push_back(419);
     names.push_back("_dy");
     names.push_back("_diboson");
     names.push_back("_singleTop");
     names.push_back("_ttbar");
	*/
    //8,30,46,
    vector<int> colors;
    for(int i=0; i<Nbkg; i++){
        if(names.at(i)=="_dy")colors.push_back(593);
        if(names.at(i)=="_dy_datadriven")colors.push_back(602);

        if(names.at(i)=="_ttbar")colors.push_back(625);
        if(names.at(i).Contains("singleTop"))colors.push_back(411);
        if(names.at(i)=="_ttbarExp")colors.push_back(634);
        if(names.at(i).Contains("QCD"))colors.push_back(5);
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
    TString labelTxt;
    if(variable_name.Contains("0btag")) labelTxt="0 b-tag";
    if(variable_name.Contains("1btag")) labelTxt="1 b-tag";
    if(variable_name.Contains("2btag")) labelTxt="#geq 2 b-tags";

    if(variable_name.Contains("1Mu_1El")) labelTxt+=" e,#mu";
    if(variable_name.Contains("2Mu_onZ")) labelTxt+=" #mu#mu on-Z";
    if(variable_name.Contains("2Mu_offZ")) labelTxt+=" #mu#mu off-Z";

    TLatex binlabel(0.6,0.85,labelTxt);
    
    
    TLegend *L = new TLegend(0.6,0.5,0.89,0.8);
	L->SetFillColor(10);
	L->SetLineColor(10);
	L->SetLineWidth(0);
    L->SetTextSize(0.04);
    
    L->AddEntry(hName[variable_name+"_singleMu"], "Data","p");
    
    for(int ibkg=0; ibkg<Nbkg; ibkg++){
        TString histName=variable_name+names.at(ibkg);
        hName[histName]->SetFillColor(colors.at(ibkg));
		if(expected_bkg && histName.Contains("ttbarExp")){
			//cout << "set fill style 3001" << endl; 
			hName[histName]->SetFillStyle(3004); 
		}
        if(expected_bkg && names.at(ibkg).Contains("_dy_datadriven"))hName[histName]->SetFillStyle(3005);
        //if(expected_bkg && names.at(ibkg).Contains("QCD"))hName[histName]->SetFillStyle(3007);

        stackName[stack_name]->Add(hName[histName]);
        
    }
    
    for(int ibkg=Nbkg-1; ibkg>=0; ibkg--){
        TString histName=variable_name+names.at(ibkg);
        L->AddEntry(hName[histName], legend_names.at(ibkg),"F");
    }
    
    
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
    gPad->SetLogy();
    
    hName[variable_name+data_name]->GetYaxis()->SetTitleOffset(0.65);
    hName[variable_name+data_name]->GetYaxis()->SetTitleSize(0.08);
    hName[variable_name+data_name]->GetYaxis()->SetLabelSize(0.05);
    hName[variable_name+data_name]->GetXaxis()->SetLabelSize(0);
    hName[variable_name+data_name]->GetXaxis()->SetTitleSize(0);
    
    hName[variable_name+data_name]->SetMarkerSize(0.75);
    hName[variable_name+data_name]->SetMarkerStyle(20);
    hName[variable_name+data_name]->SetMinimum(0.12);
    hName[variable_name+data_name]->SetMaximum(50000);

    if(print) cout << "draw stacks: " << endl;
	if(print) cout << "notblind: " << notblind << endl; 
	if(notblind){
		if(print) cout << "draw data: " << endl; 
        hName[variable_name+data_name]->DrawCopy("E1");
		if(print) cout << "draw stack: " << endl; 
        stackName[stack_name]->DrawClone("histo same");
        
    }
    else stackName[stack_name]->DrawClone("histo");
	if(notblind) hName[variable_name+data_name]->DrawCopy("E1 same");// draw data on top of MC and MC on top of data, so that data will always be visible
	if(print) cout << "systematics. " << endl; 
    if(plot_systematics){
        L->AddEntry(grName[variable_name+MC_name+"_systematic"],"Systematic unc.","F");
        grName[variable_name+MC_name+"_systematic"]->DrawClone("E2 same");
        hName[variable_name+data_name]->DrawCopy("E1 same");
    }
    
    if(scaleUpDown){
        grName[variable_name+MC_name+"_MCsystematic"]->DrawClone("E2 same");
        hName[variable_name+data_name]->DrawCopy("E1 same");

    }

    
    TLatex txt;
    txt.SetNDC(kTRUE);
    if(overlay_stlabels){
        txt.DrawLatex(0.21,0.73,"S_{T} > 300 GeV");
        txt.DrawLatex(0.4,0.51,"S_{T} > 700 GeV");
        txt.DrawLatex(0.67,0.35,"S_{T} > 1200 GeV");
    }
    if(print) cout << "draw signal: " << endl; 
	if(print) cout << "Nsig: " << signal_names.size() << endl; 
	
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
    CMS_lumi(events_pad,2,10);
    
    
    ratio_pad->cd();
    gPad->SetGridy();
    TString dataMC=variable_name+MC_name+"_RatioDataMC";
	
    hName[dataMC]->GetYaxis()->SetTitleOffset(0.35);
	hName[dataMC]->GetYaxis()->SetTitleSize(0.15);
	hName[dataMC]->GetYaxis()->SetTitleOffset(0.35);
	hName[dataMC]->GetYaxis()->SetLabelSize(0.12);
    hName[dataMC]->GetXaxis()->SetTitleOffset(0.8);
    hName[dataMC]->GetXaxis()->SetTitleSize(0.15);

	hName[dataMC]->SetMarkerSize(0.75);
    if(dataMC.Contains("EventCategories") || dataMC.Contains("met") || dataMC.Contains("npv") || dataMC.Contains("st")) {
        hName[dataMC]->GetXaxis()->LabelsOption("h");
        hName[dataMC]->GetXaxis()->SetLabelSize(0.2);
		if(dataMC.Contains("st") || dataMC.Contains("met")) hName[dataMC]->GetXaxis()->SetLabelSize(0.15);
        hName[dataMC]->GetXaxis()->SetLabelOffset(0.03);
        hName[dataMC]->GetXaxis()->SetTitleOffset(1.2);
    }


    if(print) cout << "draw ratio: " << endl;
    
    if(notblind){
        hName[dataMC]->GetYaxis()->SetNdivisions(3);
        if(ErrorX)gStyle->SetErrorX(0.5);
		hName[dataMC]->DrawCopy("E0");
        

		/*
        if(scaleUpDown){
			TString dataMCUp=variable_name+"_allMC_scaleup_RatioDataMC";
			TString dataMCDown=variable_name+"_allMC_scaledown_RatioDataMC";
			
			hName[dataMCUp]->SetLineColor(kRed);
			hName[dataMCDown]->SetLineColor(kGreen); 
			//hName[dataMCUp]->DrawCopy("histo same");
			//hName[dataMCDown]->DrawCopy("histo same");

		}
         */
	}
	
	
	
    if(plot_systematics){
        cout << "draw systematic: " << endl;
		grName[variable_name+MC_name+"_Rsystematic"]->DrawClone("E2 same");
        if(ErrorX)gStyle->SetErrorX(0.5);
        hName[dataMC]->DrawCopy("E1 same");
       
        
    }
    
    
    if(scaleUpDown) cout << "draw scale up-down. " << endl;
    if(scaleUpDown)MC_name="_allMC";
    //if(scaleUpDown)grName[variable_name+MC_name+"_RMCsystematic"]->SetFillColor(kRed);
    if(scaleUpDown){
        grName[variable_name+MC_name+"_RMCsystematic"]->DrawClone("E2 same");
        hName[dataMC]->DrawCopy("E1 same");
    }
    //reset gStyle
	//gStyle->SetErrorX(0);
}

void draw_header(){
	
	L1.SetNDC(kTRUE);
	L2.SetNDC(kTRUE);
	L1.SetTextSize(0.03);
	L2.SetTextSize(0.03);
	L1.DrawLatex(0.15,0.92, cms_pre);
	L2.DrawLatex(0.43,0.92, lumi);
	
}

void plot_jet_distributions(){
    overlay_stlabels=false;
    TString scales[]={"","_SF"};
    
    int Nscales=sizeof(scales)/sizeof(TString);
    
    std::vector<TString> signal_names;
    //signal_names.push_back("_stealth_300_200");
    signal_names.push_back("_stealth_600_300");
    
    TString variable[]={"nJets_1Mu_1El"};
    int Nvar=sizeof(variable)/sizeof(TString);
    
    vector<TString> names;
    names.push_back("_diboson");
    names.push_back("_dy");
    names.push_back("_singleTop");
    names.push_back("_ttbar");
	
	vector<TString> legend_names;
    legend_names.push_back("Diboson");
	legend_names.push_back("Drell-Yan"); 
	legend_names.push_back("Single t");
	legend_names.push_back("t#bar{t}");
	
	legend_names.push_back("M_{#tilde{q}} = 600 GeV");
	
    for(int ivar=0; ivar<Nvar;ivar++){
        for(int ib=0; ib<=2; ib++){
            for(int i=0; i<Nscales; i++){
                TString variable_name=variable[ivar]+Form("_%dbtag",ib)+scales[i];
                //cout << variable_name << endl;
                
                combine_histograms(variable_name,names,"_allMC");
                compute_ratio(variable_name, "_allMC");
                fill_stack_eventcategories(variable_name,"_singleMu",signal_names,names,legend_names,"_allMC");
            }//scales
            
        }//b-tags
    }//variables
}

void print_fractions(){
    TString tag="_SF";
    TString sample[]={"_singleTop","_dy","_diboson","_ttbar"};
    int Ns=sizeof(sample)/sizeof(TString);
	
	for(int ib=0; ib<=2; ib++){
		for(int is=0; is<Ns;is++){
			TString variable_nameAll=Form("EventCategories_1Mu_1El_%dbtag",ib)+tag+"_allMC";
			TString variable_name=Form("EventCategories_1Mu_1El_%dbtag",ib)+tag+sample[is];
			
			cout << "b-tags: " << ib << " sample: "<< sample[is] << " " << hName[variable_name]->GetBinContent(2)/hName[variable_nameAll]->GetBinContent(2) << endl;
		}
	}
}

void plot_mass(){
    overlay_stlabels=false;
    notblind=false;
    int nJet=6;
    TString variable_name=Form("h_mumu_nJet%d_0btag",nJet);
    
    vector<TString> names;
    names.push_back("_diboson");
    names.push_back("_dy");
    names.push_back("_singleTop");
    names.push_back("_ttbar");
	
	vector<TString> legend_names;
    legend_names.push_back("Diboson");
	legend_names.push_back("Drell-Yan"); 
	legend_names.push_back("Single t");
	legend_names.push_back("t#bar{t}");
	

	legend_names.push_back("M_{#tilde{q}} = 600 GeV");
    
    std::vector<TString> signal_names;
    //signal_names.push_back("_stealth_300_200");
    signal_names.push_back("_stealth_600_300");
    
    combine_histograms(variable_name,names,"_allMC");
    compute_ratio(variable_name, "_allMC");
	
    fill_stack_eventcategories(variable_name,"_singleMu",signal_names,names,legend_names,"_allMC");
    notblind=true;
}

void plot_distributions(){
    bool print=false;
	cout << "plot distributions: "<< endl;
    overlay_stlabels=false;
    TString scales[]={"_SF","_SFPbc","_SFMbc","_SFPl","_SFMl","_TopCor"};
    
    int Nscales=sizeof(scales)/sizeof(TString);
    
    std::vector<TString> signal_names;
  //  signal_names.push_back("_stealth_300_200");
    signal_names.push_back("_stealth_600_300");
   // signal_names.push_back("_stealth_600_300_compressed_ptmu20_pte20");
    
    
    TString variable[]={"EventCategories_2Mu_onZ","EventCategories_2Mu_offZ","st_1Mu_1El", "st23_1Mu_1El", "st_SS_1Mu_1El","st_2Mu_onZ","stInclusive_nJets4_1Mu_1El","stInclusive_nJets5_1Mu_1El","stInclusive_nJets6_1Mu_1El","stInclusive_nJets7_1Mu_1El",
		"stInclusive_Inclusive_nJets4_1Mu_1El","stInclusive_Inclusive_nJets5_1Mu_1El","stInclusive_Inclusive_nJets6_1Mu_1El","stInclusive_Inclusive_nJets7_1Mu_1El",
		"met_1Mu_1El","bjet1pt_1Mu_1El", "bjet2pt_1Mu_1El","ptl1_1Mu_1El","ptl2_1Mu_1El",
        "h_dilepton_mass_2Mu_nJets2_stg300","h_dilepton_mass_2Mu_nJets3_stg300","h_dilepton_mass_2Mu_nJets4_stg300","h_dilepton_mass_2Mu_nJets5_stg300",
    "h_dilepton_mass_2Mu_nJets6_stg300","h_dilepton_mass_2Mu_nJets7_stg300"};
    int Nvar=sizeof(variable)/sizeof(TString);
	
	double RebinX1[]={0,300,500,1000,3000};
	double RebinX2[]={0,100,200,300,400,500};
	double RebinX3[]={0,100,200,300,400,500};
	
	int N1=sizeof(RebinX1)/sizeof(double)-1; 
	int N2=sizeof(RebinX2)/sizeof(double)-1; 
	int N3=sizeof(RebinX3)/sizeof(double)-1; 

    vector<TString> names;
    names.push_back("_diboson");
    names.push_back("_dy");
    names.push_back("_singleTop");
    names.push_back("_ttbar");
    
	vector<TString> namesScaleup;
    namesScaleup.push_back("_diboson");
	namesScaleup.push_back("_dy");
    namesScaleup.push_back("_singleTop");
    namesScaleup.push_back("_ttJets_scaleup");
	
	
	vector<TString> namesScaledown;
    namesScaledown.push_back("_diboson");
	namesScaledown.push_back("_dy");
    namesScaledown.push_back("_singleTop");
    namesScaledown.push_back("_ttJets_scaledown");
	
	vector<TString> legend_names;
    legend_names.push_back("Diboson");
	legend_names.push_back("Drell-Yan"); 
	legend_names.push_back("Single t");
	legend_names.push_back("t#bar{t}");
	
	legend_names.push_back("M_{#tilde{q}} = 600 GeV");
    
    
    for(int ivar=0; ivar<Nvar;ivar++){
        // cout << "variable: " << variable[ivar] << endl;
        for(int ib=0; ib<=2; ib++){
            for(int i=0; i<Nscales; i++){
                TString variable_name=variable[ivar]+Form("_%dbtag",ib)+scales[i];
                cout << "compute systematics: " << variable_name << endl;
                //if(variable[ivar]=="st_1Mu_1El")variable_name=variable[ivar]+Form("_%dbtag",ib)+scales[i];
                if(print)cout << variable_name << endl;
                
                combine_histograms(variable_name,names,"_allMC");
				combine_histograms(variable_name,namesScaledown,"_allMC_scaledown");
				combine_histograms(variable_name,namesScaleup,"_allMC_scaleup");
                
            }
        }
    }
    
    compute_systematic_MC("st_1Mu_1El",names);
    
    cout << "plot" << endl;
    for(int ivar=0; ivar<Nvar;ivar++){
       // cout << "variable: " << variable[ivar] << endl;
        for(int ib=0; ib<=2; ib++){
            for(int i=0; i<Nscales; i++){
                TString variable_name=variable[ivar]+Form("_%dbtag",ib)+scales[i];
                cout << "plot: " << variable_name << endl;
                //if(variable[ivar]=="st_1Mu_1El")variable_name=variable[ivar]+Form("_%dbtag",ib)+scales[i];
                if(print)cout << variable_name << endl;
                bool rebin=false; 
                
                //combine_histograms(variable_name,names,"_allMC");
				//combine_histograms(variable_name,namesScaledown,"_allMC_scaledown");
				//combine_histograms(variable_name,namesScaleup,"_allMC_scaleup");
                
                //cout << "outname: " << variable_name+"_allMC_scaledown" << endl;
                
                if(variable[ivar]=="st_1Mu_1El" && scales[i]=="_SF"){
                    //compute_systematic_MC("st_1Mu_1El",names);
                    //plot_systematics=true;
                    scaleUpDown=true;
                }
				
                //combine_histograms(variable_name,names,"_allMC");

				ErrorX=false;
				/*
				if(variable[ivar]=="st_1Mu_1El"){
					rebin=true;
                    ErrorX=true;
					compute_ratio(variable_name,"_allMC",N1,RebinX1); 
					compute_ratio(variable_name,"_allMC_scaledown",N1,RebinX1); 
					compute_ratio(variable_name,"_allMC_scaleup",N1,RebinX1); 


				}
				if(variable[ivar].Contains("met")){
					rebin=true;
                    ErrorX=true;

					compute_ratio(variable_name,"_allMC",N2,RebinX2);
					compute_ratio(variable_name,"_allMC_scaledown",N2,RebinX2); 
					compute_ratio(variable_name,"_allMC_scaleup",N2,RebinX2); 


				}
				if(variable[ivar].Contains("bjet") || variable[ivar].Contains("ptl")){
					rebin=true;
                    ErrorX=true;

					compute_ratio(variable_name,"_allMC",N3,RebinX3);
					compute_ratio(variable_name,"_allMC_scaledown",N3,RebinX3); 
					compute_ratio(variable_name,"_allMC_scaleup",N3,RebinX3); 

				}
				*/
                if(!rebin) {
					compute_ratio(variable_name, "_allMC");// if not one of the standard rebinnings, just use histogram bins 
					compute_ratio(variable_name, "_allMC_scaleup");// if not one of the standard rebinnings, just use histogram bins 
					compute_ratio(variable_name, "_allMC_scaledown");// if not one of the standard rebinnings, just use histogram bins 

				}
				//cout << "fill stack: " << endl;
				//cout << "scaleupdown: "<< scaleUpDown << endl; 
				//scaleUpDown=false;
                //if()
                fill_stack_eventcategories(variable_name,"_singleMu",signal_names,names,legend_names,"_allMC");
				scaleUpDown=false;
                plot_systematics=false;
            }//scales
            
        }//b-tags
    }//variables
    
    //plot pile-up
    cout << "plot pile-up" << endl;
    scaleUpDown=false;
    TString variable_name="h_before_npv";
    combine_histograms(variable_name,names,"_allMC");
    compute_ratio(variable_name,"_allMC");
    fill_stack_eventcategories(variable_name,"_singleMu",signal_names,names,legend_names,"_allMC");

    variable_name="h_npv";
    combine_histograms(variable_name,names,"_allMC");
    compute_ratio(variable_name,"_allMC");
    fill_stack_eventcategories(variable_name,"_singleMu",signal_names,names,legend_names,"_allMC");
    
    
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
	// cout << "combine histograms: " << variable_name << endl;
	// for(int i=0; i<N; i++) cout << names.at(i) << endl;
	// cout << "clone histogram: " << endl;
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
        //if(keyName.Contains("st"))cout << "histogram: " << keyName << endl;
        if(className=="TH1F" ){
            TH1F *h=(TH1F*)key->ReadObj();
            h->SetStats(kFALSE);
            h->SetBinErrorOption(TH1::kPoisson);
            hName[h->GetName()+file]=h;
            TString NN=h->GetName()+file;
            //if(NN.Contains("EventCategories")==0 && NN.Contains("Acceptance_uW")==0)continue;
            if(NN.Contains("EventCategories")){
                hName[NN]->GetXaxis()->SetTitle("N_{jets}");
                int nJ=1;
                for(int ibin=1; ibin<=hName[NN]->GetNbinsX(); ibin++){
                    if(nJ==1) hName[NN]->SetBinContent(ibin,0);
                    /*
                    TString binLabel=Form("%d",nJ);
                    if(nJ==nJetmax)binLabel=Form("#geq %d",nJ);
                    if(nJ<5) binLabel="";
                    hName[NN]->GetXaxis()->SetBinLabel(ibin,binLabel);
                     */
                    TString BL=hName[NN]->GetXaxis()->GetBinLabel(ibin);
                    if(BL.Contains("7")) BL=" #geq7";
                    hName[NN]->GetXaxis()->SetBinLabel(ibin,BL);
                    
                    nJ++;
                    if(nJ>nJetmax)nJ=1;
                }
            }//eventcat
            
            if(NN.Contains("h_mumu")){
                hName[NN]->SetAxisRange(50,150);
            }//eventcat
            /*
             int nJ=1;
             for(int i=1; i<=hName[Form("EventCategories_1Mu_1El_%dbtag_noSF",nb)]->GetNbinsX();i++){
             TString binLabel=Form("%d",nJ);
             if(nJ==nJetmax)binLabel=Form("#geq %d",nJ);
             if(nJ==1) binLabel="";
             for(int is=0; is<Nscale; is++) hName[Form("EventCategories_1Mu_1El_%dbtag%s",nb,scales[is].c_str())]->GetXaxis()->SetBinLabel(i,binLabel);
             for(int is=0; is<Nscale; is++) hName[Form("EventCategories_SS_1Mu_1El_%dbtag%s",nb,scales[is].c_str())]->GetXaxis()->SetBinLabel(i,binLabel);
             for(int is=0; is<Nscale; is++) hName[Form("EventCategories_2Mu_offZ_%dbtag%s",nb,scales[is].c_str())]->GetXaxis()->SetBinLabel(i,binLabel);
             for(int is=0; is<Nscale; is++) hName[Form("EventCategories_2Mu_onZ_%dbtag%s",nb,scales[is].c_str())]->GetXaxis()->SetBinLabel(i,binLabel);
             
             nJ++;
             if(nJ>nJetmax)nJ=1;
             }
             */
        }
    }
    
}

