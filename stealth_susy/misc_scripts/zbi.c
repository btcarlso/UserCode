double zbi(double sig, double bkg, double bkge){
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
