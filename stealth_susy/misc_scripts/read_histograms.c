void read_histograms(){
  TFile *_file0 = new TFile("output_file_allMC.root","READ"); 
  TIter nextkey(_file0->GetListOfKeys()); 
  TKey *key;
 
  while((key=(TKey*)nextkey())){
    TString className=key->ReadObj()->ClassName(); 
    //    cout << "className: " << className << endl; 
    cout << "Name: " << key->GetName()<<endl; 
    if(className=="TH1F"){
      TH1F *h=(TH1F*)key->ReadObj(); 
      cout << "Histogram: " << h->GetName() << endl; 
    }
  }

}
