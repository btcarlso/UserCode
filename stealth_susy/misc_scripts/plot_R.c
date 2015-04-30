void plot_R(){
  TGraphErrors *gr = new TGraphErrors(5); 

  double R[]={0,0,1.02,1.09,1.1,1.22,1.56}; 
  double RE[]={0,0,0.01,0.01,0.02,0.06,0.25}; 
  for(int nj=2; nj<=6; nj++){
    gr->SetPoint(nj-2,nj,R[nj]);
    gr->SetPointError(nj-2,0,RE[nj]); 
}
  
  TCanvas *c = new TCanvas(); 
  c->cd();
  gr->Draw("ap"); 
  gr->GetXaxis()->SetNdivisions(5); 
  gr->GetYaxis()->SetTitle("R"); 
  gr->GetXaxis()->SetTitle("n-jets"); 
  gr->Draw("ap"); 
}
