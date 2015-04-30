
double poisson(double N, double CL,double step){

  //  cout << "N: " << N << " CL: " << CL << " step: "<< step << endl; 
  double rhs=CL/2+TMath::Gamma(N,N);
  double rhs_=TMath::Gamma(N,N)-CL/2;
  double epsilon=1; 
  double mu_p=N;
  double mu_m=N;

  while(epsilon>0.01){
    epsilon=TMath::Abs(rhs-TMath::Gamma(N,mu_p));
    cout << "mu_p: " << mu_p << " epsilon: " << epsilon << endl; 
    mu_p+=step;
  }
  epsilon=1;
  
  while(epsilon>0.01 && mu_m>0){
    epsilon=TMath::Abs(rhs_-TMath::Gamma(N,mu_m)); 
     cout << "mu_m: " << mu_m << " epsilon: " << epsilon << endl; 
    mu_m-=step;
  }
  cout << "EYhi: " << mu_p-N << " EYlo: " << N-mu_m << endl; 
  

  return mu_p-static_cast<double>(N);
}
