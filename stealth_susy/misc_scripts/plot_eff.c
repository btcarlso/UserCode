void plot_eff(){

  TFile DEN_FILE("singleMu_All/singleMu_All.root","READ");
  TFile NUM_FILE("singleMu_All_HT750/singleMu_All_HT750.root","READ");

  if(DEN_FILE.IsOpen()!=1)cout << "Not Open" << endl; 

  TH1F *num=new TH1F("num","Num",100,0,3000);
  TH1F *den = new TH1F("den","Den",100,0,3000);

  TTree *TNum=(TTree*)NUM_FILE.FindObjectAny("tree");
  TTree *TDen=(TTree*)DEN_FILE.FindObjectAny("tree");

  float stN; 
  TNum->SetBranchAddress("st",&stN); 


  float stD;
  TDen->SetBranchAddress("st",&stD);

  int _entriesN=TNum->GetEntriesFast();
  int _entriesD=TDen->GetEntriesFast();

  cout << "Entries: " << _entriesN << endl;

  for(int i=0; i<_entriesN; i++){
    TNum->GetEntry(i);
    num->Fill(stN);
    }

  for(int i=0; i<_entriesD; i++){
    TDen->GetEntry(i);
    if(den>400)
      den->Fill(stD);
  } 


  cout << "Num: " << num->GetEntries() << " Den: "<< den->GetEntries() << endl; 

  TCanvas * C = new TCanvas();
  C->Divide(2,2);
  C->cd(1);
  num->DrawCopy();
  C->cd(2);
  gPad->SetLogy();
  den->DrawCopy();

  num->Divide(den);
  C->cd(3);
  num->SetStats(kFALSE);
  num->DrawCopy();


}
