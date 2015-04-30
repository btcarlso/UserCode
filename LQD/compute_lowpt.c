void compute_lowpt(){

  TFile *fbkg=new TFile("output_file_allMC.root","READ"); 

  for(int M=300; M<=900; M+=100){
    TFile *fsig=new TFile(Form("output_file_RPV_LQD221_M%d.root",M),"READ"); 
    cout << "M: " << M << endl; 
    for(int nJ=5; nJ<=7; nJ++){
      TH1F *hsig=(TH1F*)fsig->Get(Form("h_nJets%d_stmin_lowpt",nJ));
      TH1F *hbkg=(TH1F*)fbkg->Get(Form("h_nJets%d_stmin_lowpt",nJ));
      std::vector<double> sb; 
      std::vector<double> sig; 
      std::vector<double> bkg; 
      std::vector<double> stmin; 
      for(int ist=1; ist<=hsig->GetNbinsX(); ist++){
	if(hsig->GetBinContent(ist)>0 && hbkg->GetBinContent(ist)>0 && hsig->GetBinCenter(ist)>300){
	  sig.push_back(hsig->GetBinContent(ist));
	  bkg.push_back(hbkg->GetBinContent(ist)); 
	  stmin.push_back(hsig->GetBinCenter(ist)); 
	  sb.push_back(hsig->GetBinContent(ist)/TMath::Sqrt(hbkg->GetBinContent(ist))); 
	  }
	
      }
      int iMax = TMath::LocMax(sb.size(),&sb[0]);
      //cout << "M: " << M << " nJ: " << "stmin: " << stmin.at(iMax) <<  " sig: " << sig.at(iMax) << "  bkg " << bkg.at(iMax) << endl;
      cout << sig.at(iMax) << " "  << bkg.at(iMax) << " " ; 
    }
    cout << endl; 
  }

}
