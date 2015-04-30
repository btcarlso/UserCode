#include "flavor_binned.h"
void flavor_binned(){
    TFile F("output_file_ttFullLept.root","READ");
    for(int itag=0; itag<=2; itag++){
        TH1F *g=(TH1F*)F.Get(Form("nJets_1Mu_1El_%dbtag_tight_SF",itag));
        hName[g->GetName()]=g;
        for(int ib=0; ib<=2; ib++){
            TH1F *h=(TH1F*)F.Get(Form("nJets_1Mu_1El_nbtruth%d_%dbtag_tight_SF",ib,itag));
            if(ib==0) h->SetFillColor(kRed);
            if(ib==1)h->SetFillColor(kBlue);
            if(ib==2)h->SetFillColor(kGreen);
            hName[h->GetName()]=h;
        }
    }
    
    THStack *stack_b0=new THStack();
    THStack *stack_b2=new THStack();
    for(int ib=0; ib<=2; ib++)stack_b0->Add(hName[Form("nJets_1Mu_1El_nbtruth%d_%dbtag_tight_SF",ib,0)]);
    
    for(int ib=0; ib<=2; ib++)stack_b2->Add(hName[Form("nJets_1Mu_1El_nbtruth%d_%dbtag_tight_SF",ib,2)]);

    
    TCanvas *c = new TCanvas();
    c->Divide(2,2);
    c->cd(1);
    gPad->SetLogy();
    stack_b0->DrawClone("histo");
    c->cd(2);
    gPad->SetLogy();
    stack_b2->DrawClone("histo");
    
    for(int ib=0; ib<=2; ib++){
        hName[Form("nJets_1Mu_1El_nbtruth%d_%dbtag_tight_SF",ib,0)]->Divide(hName[Form("nJets_1Mu_1El_%dbtag_tight_SF",0)]);
        c->cd(3);
        hName[Form("nJets_1Mu_1El_nbtruth%d_%dbtag_tight_SF",ib,0)]->SetFillStyle(0);
        hName[Form("nJets_1Mu_1El_nbtruth%d_%dbtag_tight_SF",ib,2)]->SetFillStyle(0);
        if(ib==0){
            hName[Form("nJets_1Mu_1El_nbtruth%d_%dbtag_tight_SF",ib,0)]->SetLineColor(kRed);
            hName[Form("nJets_1Mu_1El_nbtruth%d_%dbtag_tight_SF",ib,2)]->SetLineColor(kRed);
        }
        if(ib==1){
            hName[Form("nJets_1Mu_1El_nbtruth%d_%dbtag_tight_SF",ib,0)]->SetLineColor(kBlue);
            hName[Form("nJets_1Mu_1El_nbtruth%d_%dbtag_tight_SF",ib,2)]->SetLineColor(kBlue);
        }
        
        if(ib==2){
            hName[Form("nJets_1Mu_1El_nbtruth%d_%dbtag_tight_SF",ib,0)]->SetLineColor(kGreen);
            hName[Form("nJets_1Mu_1El_nbtruth%d_%dbtag_tight_SF",ib,2)]->SetLineColor(kGreen);
        }
        
        
        hName[Form("nJets_1Mu_1El_nbtruth%d_%dbtag_tight_SF",ib,0)]->SetMinimum(0);
        hName[Form("nJets_1Mu_1El_nbtruth%d_%dbtag_tight_SF",ib,0)]->SetMaximum(1.1);
        if(ib==0)hName[Form("nJets_1Mu_1El_nbtruth%d_%dbtag_tight_SF",ib,0)]->DrawCopy("histo ");
        else hName[Form("nJets_1Mu_1El_nbtruth%d_%dbtag_tight_SF",ib,0)]->DrawCopy("histo same");

        hName[Form("nJets_1Mu_1El_nbtruth%d_%dbtag_tight_SF",ib,2)]->Divide(hName[Form("nJets_1Mu_1El_%dbtag_tight_SF",2)]);
        c->cd(4);
        hName[Form("nJets_1Mu_1El_nbtruth%d_%dbtag_tight_SF",ib,2)]->SetMinimum(0);
        hName[Form("nJets_1Mu_1El_nbtruth%d_%dbtag_tight_SF",ib,2)]->SetMaximum(1.1);
        if(ib==0)hName[Form("nJets_1Mu_1El_nbtruth%d_%dbtag_tight_SF",ib,2)]->DrawCopy("histo ");
        else hName[Form("nJets_1Mu_1El_nbtruth%d_%dbtag_tight_SF",ib,2)]->DrawCopy("histo same");
    }
    
    
}
