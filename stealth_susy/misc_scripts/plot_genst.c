void plot_genst(){

  Float_t st=-9; 
  Float_t gen_st=-9; 
  Float_t gen_st2=-9; 
  Int_t nJets=-9; 

  TFile *input = new TFile("wJetsToLNu_gen/wJets_gen.root","READ"); 
  TFile *input2= new TFile("w2JetsToLNu_gen/w2Jets_gen.root","READ");
  TTree *tree = (TTree*)input->Get("tree"); 
  TTree *tree2 = (TTree*)input2->Get("tree"); 
  tree->SetBranchAddress("st",&st); 
  tree->SetBranchAddress("jets_n",&nJets); 
  tree->SetBranchAddress("gen_st",&gen_st); 
  tree2->SetBranchAddress("gen_st",&gen_st2); 
Long64_t entries = tree->GetEntries(); 

    TH1F *stg_2 = new TH1F ("stg_2","s_{T} gen 2 jets", 300, 0,3000); 
    TH1F *stg_3 =new TH1F ("stg_3","s_{T} gen 3jets", 300, 0,3000); 
    TH1F *stg_4 =new TH1F ("stg_4","s_{T} gen 4jets", 300, 0,3000);
    TH1F *stg_5 =new TH1F ("stg_5","s_{T} gen 5jets", 300, 0,3000);
    TH1F *stg_6 =new TH1F ("stg_6","s_{T} gen 6jets", 300, 0,3000);

    stg_2->SetLineColor(kRed);
    stg_3->SetLineColor(kRed);
    stg_4->SetLineColor(kRed); 
    stg_5->SetLineColor(kRed);
    stg_6->SetLineColor(kRed);

    TH1F *st_2 =new TH1F ("st_2","s_{T} 2jets", 300, 0,3000);
    TH1F *st_3 =new TH1F ("st_3","s_{T} 3jets", 300, 0,3000);
    TH1F *st_4 =new TH1F ("st_4","s_{T} 4jets", 300, 0,3000);
    TH1F *st_5 =new TH1F ("st_5","s_{T} 5jets", 300, 0,3000);
    TH1F *st_6 =new TH1F ("st_6","s_{T} 6jets", 300, 0,3000);



    TH1F *st_incl4= new TH1F("st_incl4","s_{T}",150,0,3000); 
    TH1F *st_incl2=new TH1F("st_incl2","s_{T}", 150,0,3000); 


    Long64_t entries2=tree2->GetEntries(); 
    for(Long64_t i=0; i<entries2; i++){
      tree2->GetEntry(i);
      st_incl2->Fill(gen_st2); 
    }

for(Long64_t i=0; i<entries; i++){
    tree->GetEntry(i);
    st_incl4->Fill(gen_st); 
    if(nJets==2){
       stg_2->Fill(gen_st); 
       st_2->Fill(st); 
    }

    if(nJets==3){
       stg_3->Fill(gen_st);
       st_3->Fill(st);
    }


    if(nJets==4){
       stg_4->Fill(gen_st);
       st_4->Fill(st);
    }
    if(nJets==5){
      stg_5->Fill(gen_st);
      st_5->Fill(st);
    }

    if(nJets==6){
      stg_6->Fill(gen_st);
      st_6->Fill(st);
    }



}

TCanvas * C =new TCanvas(); 
C->Divide(3,3); 
C->cd(1); 
 gPad->SetLogy();
st_2->Draw();
stg_2->Draw("same"); 

C->cd(2); 
 gPad->SetLogy();
st_3->Draw();
stg_3->Draw("same"); 

C->cd(3); 
 gPad->SetLogy();
st_4->DrawCopy();
stg_4->DrawCopy("same"); 

C->cd(4); 
 gPad->SetLogy();
 st_5->DrawCopy();
 stg_5->DrawCopy("same");
 
 C->cd(5);
 gPad->SetLogy();
 st_6->DrawCopy();
 stg_6->DrawCopy("same"); 
 C->cd(6);

stg_4->Divide(stg_3); 
stg_4->Draw(); 
 C->cd(7);
 stg_5->Divide(stg_3);
 stg_5->Draw();
 C->cd(8);
 stg_6->Divide(stg_3);
 stg_6->Draw();

 TCanvas *C2= new TCanvas();
 st_incl4->Divide(st_incl2);
 st_incl4->Draw(); 

}
