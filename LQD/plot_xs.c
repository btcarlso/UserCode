void plot_xs(){
  float M[]={300,400,500,600,700,800,900,1000};
  float R[]={0.0285,0.0402,0.1066,0.1846,0.5133,0.6352,0.97,1.99};

  float stopxs[]={1.99608,0.35683, 0.0855847,0.0248009,0.0081141,0.00289588,0.00109501,0.000435488}; 
  int N=sizeof(M)/sizeof(float);

  for(int i=0; i<N; i++)R[i]=stopxs[i]*R[i]; 

  TGraph *gr=new TGraph(N,M,R); 
  TGraph *th=new TGraph(N,M,stopxs); 

  th->SetLineColor(kRed);
  
  gr->Draw("al");
  gr->GetXaxis()->SetTitle("M_{#tilde{t}} (GeV)"); 
  gr->GetYaxis()->SetTitle("#sigma (pb)"); 
  gr->GetYaxis()->SetRangeUser(0.0001,0.1);
  gr->GetXaxis()->SetRangeUser(300,1000);
  gr->Draw("al"); 
  th->Draw("l same"); 
  gPad->SetLogy(); 
  
  TLegend *Leg=new TLegend(0.20,.2,0.4,0.4);
  Leg->SetFillColor(10); 
  Leg->AddEntry(gr,"EXO-12-041","l");
  Leg->AddEntry(th,"#sigma_{#tilde{t}#tilde{t}}","l");
  Leg->Draw();

}
