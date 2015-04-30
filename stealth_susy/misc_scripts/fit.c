double fcn(double *x, double *par){
  return par[0]*(1+par[1]/x[0]); 
}

void fit(){
  TH1F *stMC=(TH1F*)_file0->Get("st_nJets3_W");
  TH1F *stdata=(TH1F*)_file1->Get("st_nJets3_W");

  stdata->Divide(stMC); 
  TF1 *fcn = new TF1("fcn",fcn,300,3000,2); 
 
  stdata->Fit("fcn","R"); 
  fcn->SetRange(300,3000);
  stdata->Draw();
  fcn->Draw("same");
}
