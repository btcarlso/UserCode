
void pois(double n){

  const double alpha = 1 - 0.6827;
  double l = 0.5*TMath::ChisquareQuantile(alpha/2,2*n);
  double h = (n==0) ? ( 0.5*TMath::ChisquareQuantile(1-alpha,2*(n+1)) ) : ( 0.5*TMath::ChisquareQuantile(1-alpha/2,2*(n+1)) );
   
   cout << "n: " << n << " +/- " << h << " " << l << endl; 
   if(n>0)cout << "n: " << n << " +/- [%] " << 100*(h-n)/n << " " << 100*(n-l)/n << endl;
}
