void plot_jets(){
  TCanvas *c = new TCanvas(); 
  c->Divide(2,3); 
  for(int nJ=1; nJ<=6; nJ++){
    TProfile *jf_0b=(TProfile*)_file0.Get(Form("jet_flavor_nJets%d_0btag",nJ));
    TProfile *jf_1b=(TProfile*)_file0.Get(Form("jet_flavor_nJets%d_1btag",nJ)); 
    TProfile *jf_2b=(TProfile*)_file0.Get(Form("jet_flavor_nJets%d_2btag",nJ)); 
    cout << "c->cd(nJ)" << nJ << endl; 
    c->cd(nJ); 
    jf_0b->SetLineColor(kBlack);
    jf_1b->SetLineColor(kBlue);
    jf_2b->SetLineColor(kRed);
    jf_0b->SetMinimum(0);
    jf_0b->SetMaximum(1); 
    jf_0b->DrawClone("histo");
    jf_1b->DrawClone("histo same");
    jf_2b->DrawClone("histo same");
  }
}
