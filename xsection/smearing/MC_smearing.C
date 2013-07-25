
#include "MC_smearing.h"

void MC_smearing(){
	if(MC->IsOpen()!=kTRUE){
		cout << "MC failed to open." << endl; 
		return; 
	}
	
	data = new TFile(data_file.c_str(), "READ"); 
	cout << "Opening Data File: " <<  data_file << endl; 
	
	if(data->IsOpen()!=kTRUE){
		cout << "data failed to open." << endl; 
		return; 
	}
	cout << "Output Dir: " << output_dir << endl; 
	cout << "Output File Name: " << FN << endl; 
	
	if(outfile->IsOpen()!=kTRUE){
		cout << "outfile failed to open." << endl; 
		return; 
	}
	

	cout <<  "bins: " << endl << endl;
	for(int i=0; i<=fNpt2; i++)
		cout << fPTbin2[i] << " ";
	cout << endl; 
	
	cout << "Mass Windows: " << endl;
	cout << "Y1S: "<< y1m[0] << "-" << y1m[1] << endl; 
	cout << "Y2S: "<< y2m[0] << "-" << y2m[1] << endl; 
	cout << "Y3S: "<< y3m[0] << "-" << y3m[1] << endl; 

	
	TStopwatch ts;
	ts.Start();
	
    // resetFiles();
     for(int ups=1; ups<=3; ups++){
		// if(ups!=1) continue;
		 peak(ups);
	 }//loop over ups type	

  ts.Stop(); 
  cout << "Real Time: " << ts.RealTime() << "s" << endl; 
  cout << "CPU Time: " << ts.CpuTime() << "s" << endl; 
  cout << "Closing output file. " << endl; 
  outfile->Close();
	
	
}

int get_memory_usage_kb(long* vmrss_kb, long* vmsize_kb)
{
	/* Get the the current process' status file from the proc filesystem */
	FILE* procfile = fopen("/proc/self/status", "r");
	
	const long to_read = 8192;
	char buffer[to_read];
	int read = fread(buffer, sizeof(char), to_read, procfile);
	fclose(procfile);
	
	short found_vmrss = 0;
	short found_vmsize = 0;
	char* search_result;
	
	/* Look through proc status contents line by line */
	char delims[] = "\n";
	char* line = strtok(buffer, delims);
	
	while (line != NULL && (found_vmrss == 0 || found_vmsize == 0) )
    {
		search_result = strstr(line, "VmRSS:");
		if (search_result != NULL)
		{
			sscanf(line, "%*s %ld", vmrss_kb);
			found_vmrss = 1;
		}
		
		search_result = strstr(line, "VmSize:");
		if (search_result != NULL)
		{
			sscanf(line, "%*s %ld", vmsize_kb);
			found_vmsize = 1;
		}
		
		line = strtok(NULL, delims);
    }
	
	return (found_vmrss == 1 && found_vmsize == 1) ? 0 : 1;
}

double get_dm(TH1D *dm_distribution){
	//Function randomly samples a dm distribution, and returns a dm value in units of GeV 

	double dm=999; 
	
	TH1D *dm_clone = (TH1D*)dm_distribution->Clone("dm_clone"); 
	dm_clone->SetAxisRange(40,dm_max); 
	
	while (dm*1000>dm_max) {
		dm=dm_clone->GetRandom()/1000.;
	}

	delete dm_clone;
	
	return dm; 
	
}

double get_corr(TH2D *dm_m,TH2D *ref, int peak,double pm){
	
	TH1D *massplot=dm_m->ProjectionX("massplot",0,dm_m->GetNbinsY());
	int iYbin1=0; 
	int iYbin2=0; 
	if(peak==1){
		iYbin1=massplot->FindBin(9.26); 
		iYbin2=massplot->FindBin(pm);
	}
	if(peak==2){
		iYbin1=massplot->FindBin(9.95); 
		iYbin2=massplot->FindBin(pm);
	}
	if(peak==3){
		iYbin1=massplot->FindBin(10.15); 
		iYbin2=massplot->FindBin(pm);
	}
	
	TH1D *dm_distribution=dm_m->ProjectionY("dm_distribution",iYbin1,iYbin2);
	TH1D *dm_distribution_ref=ref->ProjectionY("dm_distribution_ref",iYbin1,iYbin2);

	double shift=(dm_distribution->GetMean()-dm_distribution_ref->GetMean())/1000;

	delete dm_distribution;
	delete dm_distribution_ref;
	delete massplot; 

	return shift;
	
	
}

void peak(int ups){
	string mass_func_name=Form("mass_func_%dS",ups);
	mass_function=(TF1*)MC->FindObjectAny(mass_func_name.c_str());
	
	long vmrss; 
	long vmsize; 
	
	for(int ipt=0; ipt<fNpt2; ++ipt){
			//if(fPTbin2[ipt]!=30) continue;
			get_memory_usage_kb(&vmrss,&vmsize);
			cout << "ipt: " << ipt << " iy: 0" << " Memory: " << vmrss << endl; 
			PDF(ups, 0,ipt);
	}//pt loop
	
}

double shape_m(double *x, double *par){
	
	double f=0; 
	double neventd=static_cast<double> (N_event);
	
	for(int i=0; i<N_event; i++){
		int im=2*i;//index for mean
		int is=2*i+1;//index for sigma
		double meanprime=par[im]; // +x[1]
		double sigmaprime=par[is];
		double norm=1./( TMath::Sqrt(2*TMath::Pi())*sigmaprime); // Normalization
		f+=norm*TMath::Gaus(x[0],meanprime,sigmaprime);// Gaussian at point x, with position mean, with width sigma, normalized by 1/2pi
		
	}
	
	return f/neventd; 
}

double PDF_shape(double *x, double *par){
	
	double f=0; 
	double neventd=static_cast<double> (N_event);
	
	for(int i=0; i<N_event; i++){
		int im=2*i;//index for mean
		int is=2*i+1;//index for sigma
		double meanprime=par[im]; // +x[1]
		double sigmaprime=par[is]*x[1];
		double norm=1./( TMath::Sqrt(2*TMath::Pi())*sigmaprime); // Normalization
		f+=norm*TMath::Gaus(x[0],meanprime,sigmaprime);// Gaussian at point x, with position mean, with width sigma, normalized by 1/2pi
	
	}
	
	return f/neventd; 
}

void PDF(int peak, int iy, int ipt){
	
	bool print=0; 
	if(print==1) cout << "Generate PDF. " << endl; 
	
	//string file_name=output_dir+"smeared_nocorr.root"; 
	

	
	string LShist_name="INIT";
	int ref_bin=fNpt2-3; 
	string dm_m_name="INIT"; 
	string dm_m_name_ref="INIT";
	
	LShist_name=LS_namew(peak,iy,ipt); 
	dm_m_name=dm_m_hist_name_weighted(iy, ipt,mode); 
	//dm_m_name=dm_m_hist_name(iy,ipt); 
	//dm_m_name="dm_m";
	
	
	string genmass_filled_name=gen_name(peak,iy,ipt); 
	string dm_used_name=dm_name(peak,iy,ipt); 
	
	string bin_name=BN(iy,ipt); 

	if(print) cout << "bin:" << bin_name << endl; 
	if(print) cout << "Before Creating Histograms. " << endl; 
	
	double dm_scale_max;
	double dm_scale_min;
	if(fPTbin2[ipt]==add_params){
		dm_scale_max=1.0+dm_scale_width_last[0];
		dm_scale_min=1.0-dm_scale_width_last[1];
	}
	else {
		dm_scale_min=1.0-dm_scale_width; 
		dm_scale_max=1.0+dm_scale_width; 
	}

	
	TF2 *LS = new TF2(LShist_name.c_str(), PDF_shape, 8.7,11.2,dm_scale_min,dm_scale_max,N_event*2); 
	//LS->SetNpx(1000); 
	TH1D *mass_filled= new TH1D(genmass_filled_name.c_str(),genmass_filled_name.c_str(),1250,8.7,11.2);
 

	if(print) cout << "After Creating Histograms. " << endl; 
	
	double l=0;
	double h=0; 
	double R=0;
	double pm;
	if(peak==1){
		l=9.0;//was 9
		h=9.455;
		pm=9.4603; // PDG value for peak mass 
		R=7.5; //7.5***actual value
	}
	
	if(peak==2){
		l=9.6; 
		h=10.05;
		pm=10.02326;
		R=6.11;//6.11
	}
	if(peak==3){
		l=9.9;
		h=10.355;
		pm=10.3552;
		R=6.68;//6.68

	}
	

	if(print) cout << dm_m_name << endl; 
	TH2D *dm_m=(TH2D*)data->FindObjectAny(dm_m_name.c_str());
	TH1D *dm_used = new TH1D(dm_used_name.c_str(),dm_used_name.c_str(),dm_m->GetNbinsY(),0,1000);

	//TH2D *dm_m_ref=dm_m_ref=(TH2D*)data->FindObjectAny(dm_m_name_ref.c_str());
	
	double mean_corr=0.0; 
	/*
	if(ipt>ref_bin){	
		mean_corr=get_corr(dm_m,dm_m_ref,peak,pm); 
		dm_m=dm_m_ref; 
	}
	 */
	
	int Ntail=static_cast<int>(static_cast<double>(N_event)/(R+1) );
	if(print==1) cout << "Number of tail events: " << Ntail << " Total Number of events: " << N_event << " Ratio: " << static_cast<double>(N_event-Ntail)/static_cast<double>(Ntail) << endl; 

	outfile->cd();

	TH1D *massPDFdata=dm_m->ProjectionX("massPDFdata",0,dm_m->GetNbinsY());
	int iYbin1=0; 
	int iYbin2=0; 
	
	if(peak==1){
		iYbin1=massPDFdata->FindBin(y1m[0]); 
		//iYbin2=massPDFdata->FindBin(pm);
		iYbin2=massPDFdata->FindBin(y1m[1]);
	}
	if(peak==2){
		iYbin1=massPDFdata->FindBin(y2m[0]); 
		//iYbin2=massPDFdata->FindBin(pm);
		iYbin2=massPDFdata->FindBin(y2m[1]);
	}
	if(peak==3){
		iYbin1=massPDFdata->FindBin(y3m[0]); 
		//iYbin2=massPDFdata->FindBin(pm);
		iYbin2=massPDFdata->FindBin(y3m[1]);
	}
	
	TH1D *dm_distribution=dm_m->ProjectionY("dm_distribution",iYbin1,iYbin2);
	
	if(print) cout << "dm distribution mean: " << dm_distribution->GetMean() << endl; 
	
	TF1 *gaus = new TF1("gaus","gaus", 8.7,11.2); 
	
	TH1D *LS1D=new TH1D("LS1D","massPDF", 2500, 8.7,11.2); 
	
	TH1D *mass_tail= new TH1D("mass_tail","mass_tail",1250,8.7,11.2);
	
	double x1;
	double x2; 
	mass_function->GetRange(x1,x2); 
	mass_function->SetRange(8.7,pm-0.002); 
	mass_tail->Add(mass_function); 
	mass_function->SetRange(x1,x2); 
	
	int N1=N_event/1000+5;
	int N2=Ntail/1000+5;
	
	cout <<  "N1: " << N1 << " N2: " << N2 << endl;
	
	TGraph *gr_ks = new TGraph(N1); 
	gr_ks->SetName("gr_ks"); 
	gr_ks->GetXaxis()->SetTitle("Number of dm samplings"); 
	gr_ks->GetYaxis()->SetTitle("KS"); 

	
	TGraph *gr_area = new TGraph(N2);
	gr_area->GetXaxis()->SetTitle("Number of FSR tail sampling"); 
	gr_area->GetXaxis()->SetTitle("Ratio of High to Low Integral"); 
	gr_area->SetName("gr_area"); 
	
	TGraph *gr_tail_ks=new TGraph(N2); 
	
	gr_tail_ks->GetXaxis()->SetTitle("Number of FSR tail sampling"); 
	gr_tail_ks->GetXaxis()->SetTitle("KS test for tail"); 
	gr_tail_ks->SetName("gr_tail_ks"); 
	
	int ip1=1; 
	int ip2=1; 
	//Loop over events 
	for(int i=0; i<N_event; i++){	
		int im=2*i; 
		int is=2*i+1; 

		 //Random sampling method 
	
	  //double dm=dm_distribution->GetRandom(30,130)/1000.;
	  double dm=-9;
	  double m=-9;
	
	//Get mass from histograms 	
	//mass_function->SetRange(8.7,9.4);	
	if(i<Ntail)m=mass_function->GetRandom(l,pm-0.002); //should be pm-0.002
	else m=pm;

	dm=get_dm(dm_distribution)+mean_corr; //+mean_corr..+4*mean_corr

	if(print==1) cout << "dm: " << dm << endl; 
	if(print==1) cout << "Event Number: " << i << endl; 
	if(print==1) cout << "Good Event found at dm: "<< dm << endl; 
	if(print==1) cout << "m "<< m << endl; 
		
     mass_filled->Fill(m);
     dm_used->Fill(dm*1000);
    if(print==1) cout << "mass and dm used written to histogram." << endl; 	

		
		LS->SetParameter(im,m);
		LS->SetParameter(is,dm); 
		gaus->SetParameter(0,1/(TMath::Sqrt(2*TMath::Pi())*dm));
		gaus->SetParameter(1,m);
		gaus->SetParameter(2,dm); 
		LS1D->Add(gaus);
		//GF->SetNpy(1000); 
		//GF->SetNpz(1000); 
		
	if(print) cout << LShist_name << endl; 
	/*
	if((i<1000 && i%200==0)||(i>=1000 && i%1000==0)){
			cout << "Event: "<< i << endl;
			TH1D *tmp = (TH1D*)LS1D->Clone("tmp"); 
			tmp->Scale(1./tmp->Integral()); 
			if(print) cout << "mass shape Integral: " << tmp->Integral() << endl; 
			double KS=dm_distribution->KolmogorovTest(dm_used);
			double A=tmp->Integral(tmp->FindBin(9.56),tmp->FindBin(9.65))/tmp->Integral(tmp->FindBin(9.1),tmp->FindBin(9.345));
			if(print) cout << "Kolmogorov Test: " <<  KS << endl; 
			if(print)cout << "Graph point: " << ip1 << endl; 
			gr_ks->SetPoint(ip1,i,KS); 
			ip1++; 
			if(i<Ntail){
				TH1D *tmp_tail = (TH1D*)mass_filled->Clone(); 
				tmp_tail->SetAxisRange(8.7,pm-0.002); 
				double KStail=mass_tail->KolmogorovTest(tmp_tail); 
				
				if(print)cout << "High Int/Low Int: " << A << endl; 
				if(print)cout << "Tail KS: " << KStail << endl; 
				gr_area->SetPoint(ip2,i,A); 
				gr_tail_ks->SetPoint(ip2,i,KStail); 
				if(print)cout << "Graph Point: " << ip2 << endl; 
				ip2++;
				delete tmp_tail; 
				}
			delete tmp; 
		}
	*/
	}//end event loop 
	
	
	double BW=mass_filled->GetBinWidth(1);
	double peak_events= mass_filled->GetBinContent(mass_filled->FindBin(pm));
	double tail_events=mass_filled->Integral(mass_filled->FindBin(l),mass_filled->FindBin(pm-BW));
	
	if(print==1) cout << "Actual Peak: " << peak_events << " Actual tail: " << tail_events << " Actual peak to tail: " << peak_events/tail_events << endl << endl; 

	cout << "Writing Line shape with " << Npx << " points " << endl; 
	outfile->cd();
	LS->SetNpx(Npx);
	LS->SetNpy(Npy);
	
	
	if(ipt==0 && peak==1){
		LS1D->Write();
		gr_area->Write();
		gr_ks->Write(); 
		gr_tail_ks->Write(); 
	}
	
	delete LS1D; 
	delete gaus; 
	
	delete gr_area;
	delete gr_ks;
	delete gr_tail_ks;
	LS->Write();
	delete LS;
	mass_filled->Write();
	delete mass_filled; 
	dm_used->Write();
	delete dm_used;
	delete mass_tail;
	
	delete dm_distribution;
	delete massPDFdata;
	/*
	if(ipt>=ref_bin)
		delete dm_m; 
	else {
		delete dm_m;
		delete dm_m_ref;
	}
	 */
	delete dm_m;

	
}//end of PDF function
