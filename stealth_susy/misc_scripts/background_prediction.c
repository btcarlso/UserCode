//
//  background_prediction.c
//  
//
//  Created by Benjamin Carlson on 4/4/14.
//
//
#include "background_prediction.h"

void background_prediction(){
    gROOT->SetBatch();
 	open_file("dy");
	open_file("ttbar");
	open_file("singleTop");
	open_file("diboson");
	open_file("singleMu");
    open_file("allMC");
    
    combine_jetbins("_singleMu");
    combine_jetbins("_singleTop");
    combine_jetbins("_diboson");
    combine_jetbins("_dy");
    combine_jetbins("_ttbar");
    
    subtract_smallbkg();
    int Nbins=91;
    CreateHistogram("transfer_factor_ttbar","","","",Nbins,0,10);
    CreateHistogram("transfer_factorP_ttbar","","","",Nbins,0,10);
    CreateHistogram("transfer_factorM_ttbar","","","",Nbins,0,10);
    CreateHistogram("transfer_factor_bkgsub","","","",Nbins,0,10);
    CreateHistogram("transfer_factorP_bkgsub","","","",Nbins,0,10);
    CreateHistogram("transfer_factorM_bkgsub","","","",Nbins,0,10);

    CreateHistogram("transfer_factor_nb_nj_ttbar","","","",6,0.5,6.5,3,-0.5,2.5);
    CreateHistogram("transfer_factor_nb_nj_bkgsub","","","",6,0.5,6.5,3,-0.5,2.5);

   // for (std::map<TString,TH1F*>::iterator it=hName.begin(); it!=hName.end(); it++) cout << it->first << endl;

    for(int nb_out=0; nb_out<=2; nb_out++){
        for(int nb_in=0; nb_in<=2;nb_in++){
            if(nb_out<nb_in) continue;
            CreateHistogram(Form("transfer_factor_%d-%d_nJ2_ttbar",nb_in,nb_out),"","n-jets",Form("T[%d/%d]",nb_out,nb_in),6,0.5,6.5);
            CreateHistogram(Form("transfer_factor_%d-%d_nJ3_ttbar",nb_in,nb_out),"","n-jets",Form("T[%d/%d]",nb_out,nb_in),6,0.5,6.5);
            
            CreateHistogram(Form("transfer_factorP_%d-%d_nJ2_ttbar",nb_in,nb_out),"","n-jets",Form("T[%d/%d]",nb_out,nb_in),6,0.5,6.5);
            CreateHistogram(Form("transfer_factorP_%d-%d_nJ3_ttbar",nb_in,nb_out),"","n-jets",Form("T[%d/%d]",nb_out,nb_in),6,0.5,6.5);
            
            CreateHistogram(Form("transfer_factorM_%d-%d_nJ2_ttbar",nb_in,nb_out),"","n-jets",Form("T[%d/%d]",nb_out,nb_in),6,0.5,6.5);
            CreateHistogram(Form("transfer_factorM_%d-%d_nJ3_ttbar",nb_in,nb_out),"","n-jets",Form("T[%d/%d]",nb_out,nb_in),6,0.5,6.5);
            
            for(int nJ=2; nJ<=6; nJ++){
                transfer_factor(nb_in,nb_out,nJ,2,"_ttbar");
            }
            for(int nJ=2; nJ<=6; nJ++){
                transfer_factor(nb_in,nb_out,nJ,3,"_ttbar");
            }
        }
    }
    //now fill data
    iBin=1;
    for(int nb_out=0; nb_out<=2; nb_out++){
        for(int nb_in=0; nb_in<=2;nb_in++){
            if(nb_out<nb_in) continue;

            CreateHistogram(Form("transfer_factor_%d-%d_nJ2_bkgsub",nb_in,nb_out),"","n-jets",Form("T[%d/%d]",nb_out,nb_in),6,0.5,6.5);
            CreateHistogram(Form("transfer_factor_%d-%d_nJ3_bkgsub",nb_in,nb_out),"","n-jets",Form("T[%d/%d]",nb_out,nb_in),6,0.5,6.5);
            
            for(int nJ=2; nJ<=6; nJ++){
                transfer_factor(nb_in,nb_out,nJ,2,"_bkgsub");
            }
            for(int nJ=2; nJ<=6; nJ++){
                transfer_factor(nb_in,nb_out,nJ,3,"_bkgsub");
            }
        }
    }
    
    
    cout << "total bins filled: " << iBin << endl;
    CreateCanvas("summary","",800,400);
    CName["summary"]->cd();
    //cout << "bin error: " <<  hName["transfer_factor_bkgsub"]->GetBinError(1) << endl;

    hName["transfer_factor_ttbar"]->Draw("histo");
    hName["transfer_factorP_ttbar"]->SetLineColor(kRed);
    hName["transfer_factorM_ttbar"]->SetLineColor(kGreen);
    hName["transfer_factorP_ttbar"]->Draw("histo same");
    hName["transfer_factorM_ttbar"]->Draw("histo same");
    hName["transfer_factor_bkgsub"]->Draw("E1 same");
  //  cout << "bin error: " <<  hName["transfer_factor_bkgsub"]->GetBinError(1) << endl;
    
    CreateCanvas("ratio","",800,400);
    CName["ratio"]->cd();

    clone_histogram("transfer_factor_bkgsub","T_ratio");
    clone_histogram("transfer_factor_bkgsub","T_ratioP");
    clone_histogram("transfer_factor_bkgsub","T_ratioM");
    
    hName["T_ratio"]->Divide(hName["transfer_factor_ttbar"]);
    hName["T_ratioP"]->Divide(hName["transfer_factorP_ttbar"]);
    hName["T_ratioM"]->Divide(hName["transfer_factorM_ttbar"]);

    hName["T_ratio"]->Draw("E1");
    hName["T_ratioP"]->SetLineColor(kRed);
    hName["T_ratioM"]->SetLineColor(kGreen);
    
    hName["T_ratioP"]->Draw("histo same");
    hName["T_ratioM"]->Draw("histo same");
    
    CreateCanvas("Ratio2D","",700,600);
    CName["Ratio2D"]->cd();
    
    TF2 *fit=new TF2("fit",transfer_cor,2,6,0,2,3);
    fit->SetParameters(1,1.1,1,0.01);
    
    hName2D["transfer_factor_nb_nj_bkgsub"]->Divide( hName2D["transfer_factor_nb_nj_ttbar"]);
    hName2D["transfer_factor_nb_nj_bkgsub"]->Fit(fit,"LR");

    hName2D["transfer_factor_nb_nj_bkgsub"]->Draw("COLZ TEXT");
    fit->Draw("same");
    
    
    CreateCanvas("TransferFactor_nJ2","",900,600);
    //00,11,22,01,02,12
    CName["TransferFactor_nJ2"]->Divide(6,2);
    int iCan=1;
    for(int nb_in=0; nb_in<=2; nb_in++){
        for(int nb_out=0; nb_out<=2; nb_out++){
            if(nb_in!=nb_out && nb_out<nb_in) continue;
            CName["TransferFactor_nJ2"]->cd(iCan);
            TLatex txt(0.2,0.9, Form("%d/%d",nb_out,nb_in));
            txt.SetNDC(kTRUE);
            hName[Form("transfer_factor_%d-%d_nJ2_ttbar",nb_in,nb_out)]->Draw("histo");
            hName[Form("transfer_factorP_%d-%d_nJ2_ttbar",nb_in,nb_out)]->SetLineColor(kRed);
            hName[Form("transfer_factorM_%d-%d_nJ2_ttbar",nb_in,nb_out)]->SetLineColor(kGreen);
            
            hName[Form("transfer_factorP_%d-%d_nJ2_ttbar",nb_in,nb_out)]->Draw("histo same");
            hName[Form("transfer_factorM_%d-%d_nJ2_ttbar",nb_in,nb_out)]->Draw("histo same");

            hName[Form("transfer_factor_%d-%d_nJ2_bkgsub",nb_in,nb_out)]->DrawCopy("E1 same");
            txt.DrawClone();
            clone_histogram(Form("transfer_factor_%d-%d_nJ2_bkgsub",nb_in,nb_out),Form("Ratio_dataMC_%d-%d_nJ2",nb_in, nb_out));
            clone_histogram(Form("transfer_factor_%d-%d_nJ2_bkgsub",nb_in,nb_out),Form("RatioP_dataMC_%d-%d_nJ2",nb_in, nb_out));
            clone_histogram(Form("transfer_factor_%d-%d_nJ2_bkgsub",nb_in,nb_out),Form("RatioM_dataMC_%d-%d_nJ2",nb_in, nb_out));

            hName[Form("RatioP_dataMC_%d-%d_nJ2",nb_in, nb_out)]->SetLineColor(kRed);
            hName[Form("RatioM_dataMC_%d-%d_nJ2",nb_in, nb_out)]->SetLineColor(kGreen);
            
            hName[Form("Ratio_dataMC_%d-%d_nJ2",nb_in, nb_out)]->Divide(hName[Form("transfer_factor_%d-%d_nJ2_ttbar",nb_in,nb_out)]);
            hName[Form("RatioP_dataMC_%d-%d_nJ2",nb_in, nb_out)]->Divide(hName[Form("transfer_factorP_%d-%d_nJ2_ttbar",nb_in,nb_out)]);
            hName[Form("RatioM_dataMC_%d-%d_nJ2",nb_in, nb_out)]->Divide(hName[Form("transfer_factorM_%d-%d_nJ2_ttbar",nb_in,nb_out)]);

            hName[Form("Ratio_dataMC_%d-%d_nJ2",nb_in, nb_out)]->GetYaxis()->SetTitle("Ratio data/MC");
            CName["TransferFactor_nJ2"]->cd(iCan+6);
            
            TF1 *f=new TF1("f",transfer_cor_nJ,2,6,4);
            f->SetParameters(fit->GetParameter(0),fit->GetParameter(1),fit->GetParameter(2),fit->GetParameter(3),nb_out-nb_in);
            
            hName[Form("Ratio_dataMC_%d-%d_nJ2",nb_in, nb_out)]->Draw("E1");
            hName[Form("RatioP_dataMC_%d-%d_nJ2",nb_in, nb_out)]->Draw("histo same");
            hName[Form("RatioM_dataMC_%d-%d_nJ2",nb_in, nb_out)]->Draw("histo same");
            //f->DrawClone("same");
            delete f;
            iCan++;
        }//nbout
    }//nbin
    for(int nb_out=0; nb_out<=2; nb_out++){
        for(int nb_in=0; nb_in<=2;nb_in++){
            if(nb_out<nb_in) continue;
            transfer_factor_st( nb_in,  nb_out,  -1,  -1, "_ttbar");
            transfer_factor_st( nb_in,  nb_out,  -1,  -1, "_bkgsub");
            
            for(int nj_in=2; nj_in<=6; nj_in++){
                for(int nj_out=2; nj_out<=6; nj_out++){
                    if(nj_out<nj_in) continue;
                    transfer_factor_st( nb_in,  nb_out,  nj_in,  nj_out, "_ttbar");
                    transfer_factor_st( nb_in,  nb_out,  nj_in,  nj_out, "_bkgsub");
            }
        }
    }
    }
    
    CreateCanvas("TransferFactor_st","",900,600);
    CName["TransferFactor_st"]->Divide(3,2);
    iCan=1;
    int nj_in=-1;
    int nj_out=-1;
    for(int nb=0; nb<=2; nb++){
        CName["TransferFactor_st"]->cd(iCan);
        hName[Form("TransferFactor_st_%d-%db_%d-%dnJ_ttbar",nb,nb,nj_in,nj_out)]->Draw("histo");
        hName[Form("TransferFactor_st_%d-%db_%d-%dnJ_bkgsub",nb,nb,nj_in,nj_out)]->Draw("E1 same");
        hName[Form("TransferFactorP_st_%d-%db_%d-%dnJ_ttbar",nb,nb,nj_in,nj_out)]->Draw("histo same");
        hName[Form("TransferFactorM_st_%d-%db_%d-%dnJ_ttbar",nb,nb,nj_in,nj_out)]->Draw("histo same");
        
        clone_histogram(Form("TransferFactor_st_%d-%db_%d-%dnJ_bkgsub",nb,nb,nj_in,nj_out),Form("Ratio_DataMC_%d-%d_nJ2_st",nb,nb));
        clone_histogram(Form("TransferFactor_st_%d-%db_%d-%dnJ_bkgsub",nb,nb,nj_in,nj_out),Form("RatioP_DataMC_%d-%d_nJ2_st",nb,nb));
        clone_histogram(Form("TransferFactor_st_%d-%db_%d-%dnJ_bkgsub",nb,nb,nj_in,nj_out),Form("RatioM_DataMC_%d-%d_nJ2_st",nb,nb));

        hName[Form("Ratio_DataMC_%d-%d_nJ2_st",nb,nb)]->Divide(hName[Form("TransferFactor_st_%d-%db_%d-%dnJ_ttbar",nb,nb,nj_in,nj_out)]);
        hName[Form("RatioP_DataMC_%d-%d_nJ2_st",nb,nb)]->Divide(hName[Form("TransferFactorP_st_%d-%db_%d-%dnJ_ttbar",nb,nb,nj_in,nj_out)]);
        hName[Form("RatioM_DataMC_%d-%d_nJ2_st",nb,nb)]->Divide(hName[Form("TransferFactorM_st_%d-%db_%d-%dnJ_ttbar",nb,nb,nj_in,nj_out)]);

        hName[Form("RatioP_DataMC_%d-%d_nJ2_st",nb,nb)]->SetLineColor(kRed);
        hName[Form("RatioM_DataMC_%d-%d_nJ2_st",nb,nb)]->SetLineColor(kGreen);

        CName["TransferFactor_st"]->cd(iCan+3);
        hName[Form("Ratio_DataMC_%d-%d_nJ2_st",nb,nb)]->SetMinimum(0.5);
        hName[Form("Ratio_DataMC_%d-%d_nJ2_st",nb,nb)]->SetMaximum(1.5);
        hName[Form("Ratio_DataMC_%d-%d_nJ2_st",nb,nb)]->Draw("E1");
        hName[Form("RatioP_DataMC_%d-%d_nJ2_st",nb,nb)]->Draw("histo same");
        hName[Form("RatioM_DataMC_%d-%d_nJ2_st",nb,nb)]->Draw("histo same");

        iCan++;
    }
    TFile *output_file = new TFile("transfer_factor_plots.root","RECREATE");
    output_file->cd();
    
	for (std::map<TString,TCanvas*>::iterator it=CName.begin(); it!=CName.end(); it++) {
		CName[it->first]->Write();
    }
    for (std::map<TString,TH1F*>::iterator it=hName.begin(); it!=hName.end(); it++) {
		//hName[it->first]->Write();
    }
    
    output_file->Close();
    
}

double transfer_cor_nJ(double *x, double *par){
    double Nb=par[3];
    double Nj=x[0];
    //return par[0]+1+par[1]*Nb+par[2]*Nj+par[3]*Nj*Nj;
    return TMath::Power(1+par[0],-0.5*Nb)+par[1]*Nj+par[2]*Nj*Nj;
}

double transfer_cor(double *x, double *par){
    double Nb=x[1];
    double Nj=x[0];
   // return par[0]+1+par[1]*Nb+par[2]*Nj+par[3]*Nj*Nj;

    return TMath::Power(1+par[0],-0.5*Nb)+par[1]*Nj+par[2]*Nj*Nj;
}

TString histname(int nbtag, TString file){
    return Form("nJets_1Mu_1El_%dbtag_tight_SF",nbtag)+file;
}

TString histnameP(int nbtag, TString file){
    return Form("nJets_1Mu_1El_%dbtag_tight_SFP",nbtag)+file;
}

TString histnameM(int nbtag, TString file){
    return Form("nJets_1Mu_1El_%dbtag_tight_SFM",nbtag)+file;
}

TString histname(int nbtag,int nJ, TString file){
    if(nJ<0) return Form("st_nJetsInclusive_1Mu_1El_%dbtag_tight_SF_topreweighted",nbtag)+file;
    return Form("st_nJets%d_1Mu_1El_%dbtag_tight_SF_topreweighted",nJ,nbtag)+file;
}

TString histnameP(int nbtag, int nJ,TString file){
    if(nJ<0) return Form("st_nJetsInclusive_1Mu_1El_%dbtag_tight_SFP_topreweighted",nbtag)+file;
    return Form("st_nJets%d_1Mu_1El_%dbtag_tight_SFP_topreweighted",nJ,nbtag)+file;
}

TString histnameM(int nbtag,int nJ, TString file){
    if(nJ<0) return Form("st_nJetsInclusive_1Mu_1El_%dbtag_tight_SFM_topreweighted",nbtag)+file;
    return Form("st_nJets%d_1Mu_1El_%dbtag_tight_SFM_topreweighted",nJ,nbtag)+file;
}


void transfer_factor_st(int nb_in, int nb_out, int nj_in, int nj_out, TString file){
    //cout << "file : " << file << " transfer factor: " << nb_in << " " << nb_out << " nj " << nj_in << " " << nj_out << endl;
    // cout << "file : " << file << endl;
    TString outname=Form("TransferFactor_st_%d-%db_%d-%dnJ",nb_in,nb_out,nj_in,nj_out)+file;
    TString outnameP=Form("TransferFactorP_st_%d-%db_%d-%dnJ",nb_in,nb_out,nj_in,nj_out)+file;
    TString outnameM=Form("TransferFactorM_st_%d-%db_%d-%dnJ",nb_in,nb_out,nj_in,nj_out)+file;

    clone_histogram(histname(nb_in,nj_in,file),outname);
    if(file!="_bkgsub") {
        clone_histogram(histname(nb_in,nj_in,file),outnameP);
        clone_histogram(histname(nb_in,nj_in,file),outnameM);
        hName[outnameP]->SetLineColor(kRed);
        hName[outnameM]->SetLineColor(kGreen);

    }

    int IST0=2;


    for(int ist=2;ist<=hName[histname(nb_in,nj_in,file)]->GetNbinsX();ist++){
    
        double num=hName[histname(nb_in,nj_in,file)]->GetBinContent(ist);
        double numE=hName[histname(nb_in,nj_in,file)]->GetBinError(ist);
    
        double den=hName[histname(nb_out,nj_out,file)]->GetBinContent(IST0);
        double denE=hName[histname(nb_out,nj_out,file)]->GetBinError(IST0);
 
        hName[outname]->SetBinContent(ist,0);
        hName[outname]->SetBinError(ist,0);
        double T=num/den;
        hName[outname]->SetBinContent(ist,T);
        double TE=(T)*TMath::Sqrt(TMath::Power(numE/num,2) + TMath::Power(denE/den,2));
        hName[outname]->SetBinError(ist,TE);

        if(file=="_bkgsub") continue;
        
        
        double numP=hName[histnameP(nb_in,nj_in,file)]->GetBinContent(ist);
        double denP=hName[histnameP(nb_out,nj_out,file)]->GetBinContent(IST0);
        
        double numM=hName[histnameM(nb_in,nj_in,file)]->GetBinContent(ist);
        double denM=hName[histnameM(nb_out,nj_out,file)]->GetBinContent(IST0);
        
        double TP=numP/denP;
        double TM=numM/denM;

        hName[outnameP]->SetBinContent(ist,TP);
        hName[outnameM]->SetBinContent(ist,TM);

    }
    
}

void transfer_factor(int nb_in, int nb_out, int nj_in, int nj_out, TString file){
    //cout << "file : " << file << " transfer factor: " << nb_in << " " << nb_out << " nj " << nj_in << " " << nj_out << endl;
   // cout << "file : " << file << endl;
    
    double num=hName[histname(nb_in,file)]->GetBinContent(nj_in);
    double numE=hName[histname(nb_in,file)]->GetBinError(nj_in);

    double den=hName[histname(nb_out,file)]->GetBinContent(nj_out);
    double denE=hName[histname(nb_out,file)]->GetBinError(nj_out);

    TString label=Form("T%d-%db_%d-%dnJ",nb_in,nb_out, nj_in,nj_out);

    double T=num/den;
    if(file.Contains("_bkgsub")){
        //remove signal region (aka blind)
        if(nb_in ==0 || nb_out==0){
            if(nj_in>3 || nj_out>3) T=0;
        }
        if(nb_in==1 || nb_out==1){
            if(nj_in>4 || nj_out>4) T=0;
        }
      
    }
    double TE=(T)*TMath::Sqrt(TMath::Power(numE/num,2) + TMath::Power(denE/den,2));
    TString Name= Form("transfer_factor_%d-%d_nJ%d",nb_in,nb_out,nj_out)+file;
    hName[Name]->SetBinContent(nj_in,T);
    hName[Name]->SetBinError(nj_in,TE);

    hName["transfer_factor"+file]->SetBinContent(iBin, T);
    hName["transfer_factor"+file]->GetXaxis()->SetBinLabel(iBin,label);
    hName["transfer_factor"+file]->SetBinError(iBin, TE);
    if(nj_in-nj_out>0 && nj_out==2){
       // cout << " deltaNb: " << nb_out-nb_in << " delta NJ: " <<nj_in-nj_out<< endl;
        hName2D["transfer_factor_nb_nj"+file]->SetBinContent(nj_in-nj_out+1, nb_out-nb_in+1,T);
        hName2D["transfer_factor_nb_nj"+file]->SetBinError(nj_in-nj_out+1, nb_out-nb_in+1,TE);
    }
    //cout << "bin content: " << hName["transfer_factor_bkgsub"]->GetBinContent(iBin) << " " << hName["transfer_factor_bkgsub"]->GetBinError(iBin) << endl;

  //  cout << "bin error: " <<    hName["transfer_factor"+file]->GetBinError(iBin);
    if(file.Contains("_bkgsub")){
        iBin++;
        return;
    }
    double numP=hName[histnameP(nb_in,file)]->GetBinContent(nj_in);
    double denP=hName[histnameP(nb_out,file)]->GetBinContent(nj_out);
    
    double numM=hName[histnameM(nb_in,file)]->GetBinContent(nj_in);
    double denM=hName[histnameM(nb_out,file)]->GetBinContent(nj_out);
    
    TString NameP= Form("transfer_factorP_%d-%d_nJ%d",nb_in,nb_out,nj_out)+file;
    TString NameM= Form("transfer_factorM_%d-%d_nJ%d",nb_in,nb_out,nj_out)+file;

    hName[NameP]->SetBinContent(nj_in,numP/denP);
    hName[NameM]->SetBinContent(nj_in,numM/denM);
    
    hName["transfer_factorP"+file]->SetBinContent(iBin, numP/denP);
    hName["transfer_factorP"+file]->GetXaxis()->SetBinLabel(iBin,label);

    hName["transfer_factorM"+file]->SetBinContent(iBin, numM/denM);
    hName["transfer_factorM"+file]->GetXaxis()->SetBinLabel(iBin,label);


    iBin++;
}

void subtract_smallbkg(){
	for(int nbtag=0; nbtag<=2; nbtag++){
        clone_histogram(histname(nbtag,"_singleMu"),histname(nbtag,"_bkgsub"));
        hName[histname(nbtag,"_bkgsub")]->Add(hName[histname(nbtag,"_singleTop")],-1);
        hName[histname(nbtag,"_bkgsub")]->Add(hName[histname(nbtag,"_dy")],-1);
        hName[histname(nbtag,"_bkgsub")]->Add(hName[histname(nbtag,"_diboson")],-1);
        
        clone_histogram(histname(nbtag,-1,"_singleMu"),histname(nbtag,-1,"_bkgsub"));
        hName[histname(nbtag,-1,"_bkgsub")]->Add(hName[histname(nbtag,-1,"_singleTop")],-1);
        hName[histname(nbtag,-1,"_bkgsub")]->Add(hName[histname(nbtag,-1,"_dy")],-1);
        hName[histname(nbtag,-1,"_bkgsub")]->Add(hName[histname(nbtag,-1,"_diboson")],-1);
        
        for(int nJ=2; nJ<=8; nJ++){
            clone_histogram(histname(nbtag,nJ,"_singleMu"),histname(nbtag,nJ,"_bkgsub"));
            hName[histname(nbtag,nJ,"_bkgsub")]->Add(hName[histname(nbtag,nJ,"_singleTop")],-1);
            hName[histname(nbtag,nJ,"_bkgsub")]->Add(hName[histname(nbtag,nJ,"_dy")],-1);
            hName[histname(nbtag,nJ,"_bkgsub")]->Add(hName[histname(nbtag,nJ,"_diboson")],-1);
        }
        
        //cout << "subtracted background: " << histname(nbtag,"_bkgsub") << endl;
    }
    
}

void combine_jetbins(TString file){
    for(int nbtag=0; nbtag<=2; nbtag++){
        clone_histogram(histname(nbtag,2,file),histname(nbtag,-1,file));
        clone_histogram(histnameP(nbtag,2,file),histnameP(nbtag,-1,file));
        clone_histogram(histnameM(nbtag,2,file),histnameM(nbtag,-1,file));

        for(int nJ=3; nJ<=8; nJ++){
            hName[histname(nbtag,-1,file)]->Add(hName[histname(nbtag,nJ,file)]);
            hName[histnameP(nbtag,-1,file)]->Add(hName[histnameP(nbtag,nJ,file)]);
            hName[histnameM(nbtag,-1,file)]->Add(hName[histnameM(nbtag,nJ,file)]);

        }
    }
}

void clone_histogram(TString name1, TString clone_name){
    TH1F *h=(TH1F*)hName[name1]->Clone(clone_name);
    //h->Sumw2();
    hName[h->GetName()]=h;
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
        if(className=="TH1F" && (keyName.Contains("nJets_1Mu_1El") || keyName.Contains("st_nJets") )){
            TH1F *h=(TH1F*)key->ReadObj();
            h->SetStats(kFALSE);
            hName[h->GetName()+file]=h;
        }
    }
    
}
