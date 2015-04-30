#include "analyze_trees.h"

void analyze_trees(int jobnumber){
   expo_cor = new TF1("expo_cor", "expo",0,1000);
    expo_cor->SetParameters(0.156,-0.00137);
    setup_files(jobnumber);
    TString dir=gSystem->pwd();
    load_btagEff();
    load_btag_sys();
	bookHisto();

	open_files();
	
    gSystem->cd(dir);
	//divide_BW();
	writeHisto();
	gSystem->Exit(1,1);
}


void setup_files(int jobnumber){
    
	if(jobnumber==-2){
        //test job, not standard
		sample_list.push_back("ttFullLept");
		test=true;
		isMC=true;
		output_file_name="test";
		selections="_Eletrigger";
        //_ptmu3_pte3
        //_compressed_ptmu3_pte3
		//_stealth_600_300_ptmu3_pte3 Nominal
        //stealth_600_300_singlino200_ptmu3_pte3
        //stealth_600_300_singlino225_ptmu3_pte3
	}
    
    if(jobnumber==-1){
        //test job, not standard
		sample_list.push_back("singleEle");
		test=true;
		isMC=false;
		output_file_name="singleEle";
		selections="_trigger";
        //_ptmu3_pte3
        //_compressed_ptmu3_pte3
		//_stealth_600_300_ptmu3_pte3 Nominal
        //stealth_600_300_singlino200_ptmu3_pte3
        //stealth_600_300_singlino225_ptmu3_pte3
	}
    
	int JN=0;
    
	if(jobnumber==JN) {
		sample_list.push_back("singleMu");
		isMC=false;
		output_file_name="singleMu";
		selections="_trigger";
	}
    JN++;
    
    if(jobnumber==JN) {
		sample_list.push_back("ttFullLept");
        for(int nJ=1; nJ<=4; nJ++) {
			TString tag=Form("dy%dJets",nJ);
			sample_list.push_back(tag);
		}
        
        sample_list.push_back("WWJetsTo2L2Nu");
		sample_list.push_back("WZJetsTo2L2Q");
		sample_list.push_back("WZJetsTo3LNu");
		
		sample_list.push_back("ZZJetsTo2L2Q");
		sample_list.push_back("ZZJetsTo2L2Nu");
		sample_list.push_back("ZZJetsTo4L");
        
        sample_list.push_back("TBar_t");
		sample_list.push_back("TBar_s");
		sample_list.push_back("TBar_tW");
		
		sample_list.push_back("T_t");
		sample_list.push_back("T_s");
		sample_list.push_back("T_tW");
        
		isMC=true;
		output_file_name="allMC";
		selections="_trigger";
	}
    JN++;
    
	if(jobnumber==JN) {
		for(int nJ=1; nJ<=4; nJ++) {
			TString tag=Form("dy%dJets",nJ);
			sample_list.push_back(tag);
		}
		isMC=true;
		output_file_name="dy";
		selections="_trigger";
        
	}
    JN++;
    
	   
    //3
	if(jobnumber==JN) {
		sample_list.push_back("ttSemiLept");
        
		isMC=true;
		output_file_name="ttSemiLept";
		selections="_trigger";
        
	}
    JN++;
    
	if(jobnumber==JN) {
		sample_list.push_back("ttFullLept");
		
		isMC=true;
		output_file_name="ttFullLept";
		selections="_trigger";
	}
    JN++;
    
	if(jobnumber==JN){
		sample_list.push_back("TBar_t");
		sample_list.push_back("TBar_s");
		sample_list.push_back("TBar_tW");
		
		sample_list.push_back("T_t");
		sample_list.push_back("T_s");
		sample_list.push_back("T_tW");
		isMC=true;
		output_file_name = "singleTop";
		selections="_trigger";
        
	}
    JN++;
	//
	if(jobnumber==JN){
		sample_list.push_back("WWJetsTo2L2Nu");
		sample_list.push_back("WZJetsTo2L2Q");
		sample_list.push_back("WZJetsTo3LNu");
		
		sample_list.push_back("ZZJetsTo2L2Q");
		sample_list.push_back("ZZJetsTo2L2Nu");
		sample_list.push_back("ZZJetsTo4L");
		
		isMC=true;
		QCD=false;
		output_file_name = "diboson";
		selections="_trigger";
		
	}
    JN++;
    //19
    if(jobnumber==JN){
		sample_list.push_back("TTZJets");
        sample_list.push_back("TTWJets");

		
		isMC=true;
		QCD=false;
		output_file_name = "ttRare";
		selections="_trigger";
		
	}
    JN++;
    
    if(jobnumber==JN){
		sample_list.push_back("w1Jets");
        sample_list.push_back("w2Jets");
		sample_list.push_back("w3Jets");
		sample_list.push_back("w4Jets");
        
		
		isMC=true;
		QCD=false;
		output_file_name = "wJets";
		selections="_trigger";
		
	}
    JN++;
    
    
    
	//RPV 9
	if(jobnumber==JN) {
		sample_list.push_back("RPV_LQD221_M300");
		
		isMC=true;
		output_file_name="RPV_LQD221_M300";
		selections="_trigger";
		
	}
    cout << " RPV 300" << JN << endl;
    JN++;
	
	if(jobnumber==JN) {
		sample_list.push_back("RPV_LQD221_M400");
		
		isMC=true;
		output_file_name="RPV_LQD221_M400";
		selections="_trigger";
	}
	JN++;
	
	if(jobnumber==JN) {
		sample_list.push_back("RPV_LQD221_M500");
		
		isMC=true;
		output_file_name="RPV_LQD221_M500";
		selections="_trigger";
	}
    JN++;
    
	if(jobnumber==JN) {
		sample_list.push_back("RPV_LQD221_M600");
		
		isMC=true;
		output_file_name="RPV_LQD221_M600";
		selections="_trigger";
	}
    JN++;
	
	if(jobnumber==JN) {
		sample_list.push_back("RPV_LQD221_M700");
		
		isMC=true;
		output_file_name="RPV_LQD221_M700";
		selections="_trigger";
	}
    JN++;
	
	if(jobnumber==JN) {
		sample_list.push_back("RPV_LQD221_M800");
		
		isMC=true;
		output_file_name="RPV_LQD221_M800";
		selections="_trigger";
		
	}
    JN++;
	
	
	if(jobnumber==JN) {
		sample_list.push_back("RPV_LQD221_M900");
		
		isMC=true;
		output_file_name="RPV_LQD221_M900";
		selections="_trigger";
	}
    JN++;
	
	if(jobnumber==JN) {
		sample_list.push_back("RPV_LQD221_M1000");
		
		isMC=true;
		output_file_name="RPV_LQD221_M1000";
		selections="_trigger";
		
	}
	JN++;
    
	if(jobnumber==JN){
		sample_list.push_back("UDD300");
		isMC=true;
		QCD=false;
		output_file_name="UDD300";
		selections="_trigger";
	}
	JN++;
    
    if(jobnumber==JN) {
		sample_list.push_back("RPV_LQD221_M300");
		
		isMC=true;
		output_file_name="RPV_LQD221_M300_jecP";
		selections="_trigger_jecP";
		
	}
    cout << " RPV 300 JEC+" << JN << endl;
    JN++;
    
    if(jobnumber==JN) {
		sample_list.push_back("RPV_LQD221_M300");
		
		isMC=true;
		output_file_name="RPV_LQD221_M300_jecM";
		selections="_trigger_jecM";
		
	}
    cout << " RPV 300 JEC+" << JN << endl;
    JN++;
    
    if(jobnumber==JN) {
		sample_list.push_back("RPV_LQD221_M400");
		
		isMC=true;
		output_file_name="RPV_LQD221_M400_jecP";
		selections="_trigger_jecP";
		
	}
    cout << " RPV 300 JEC+" << JN << endl;
    JN++;
    
    if(jobnumber==JN) {
		sample_list.push_back("RPV_LQD221_M400");
		
		isMC=true;
		output_file_name="RPV_LQD221_M400_jecM";
		selections="_trigger_jecM";
		
	}
    JN++;
    
    if(jobnumber==JN) {
		sample_list.push_back("RPV_LQD221_M900");
		
		isMC=true;
		output_file_name="RPV_LQD221_M900_jecP";
		selections="_trigger_jecP";
	}
    JN++;
    
    if(jobnumber==JN) {
		sample_list.push_back("RPV_LQD221_M900");
		
		isMC=true;
		output_file_name="RPV_LQD221_M900_jecM";
		selections="_trigger_jecM";
	}
    JN++;
    
   	
}

void open_files(){
	prepareNormalization();
    
	//TString sample_list[]="wJetsht150","wJetsht200","wJetsht250","wJetsht300","wJetsht400",
	//TString sample_list[]={"RPV500"};
	
	int N_files=sample_list.size();
	cout << "N Files: " << N_files << endl;
    
	weight=1;

	cout << "output file: " << output_file_name << endl;
	TString output_FILE_NAME="output_file_"+output_file_name+".root";
	output_file=new TFile(output_FILE_NAME,"RECREATE");
	cout << "Analyzing samples: " << endl;
	for (int iF=0; iF<N_files; iF++){
		cout << sample_list.at(iF) << " " << fileName[sample_list.at(iF)] << endl;
	}
	cout << endl;
	
	
	for (int iF=0; iF<N_files; iF++) {
		cout << sample_list[iF] << endl;
		if(isMC) weight=xs[sample_list[iF]]*19700;
		TString FN=dir+fileName[sample_list.at(iF)]+selections;
        if(fileName[sample_list.at(iF)].Contains("tt")) reweight_top = true;
        else reweight_top=false;
		analyze_file(FN);
        
	}
	
}

void prepareNormalization(){
   std::cout <<"prepare for normalization .... " << std::endl;
   
   // create maps of filename, xsec and generated events
   xs["singleEle"]=1.0;
   xs["singleMu"]=1.0;

   fileName["singleMu"]="singleMu";
   fileName["singleEle"]="single";
    
   fileName["singleEleA"]="singleElectron_A";
   fileName["singleEleB"]="singleElectron_B";
   fileName["singleEleC"]="singleElectron_C";
   fileName["singleEleD"]="singleElectron_D";

//
   for(int nJ=1; nJ<=4; nJ++){
      fileName[Form("w%dJets",nJ)]=Form("w%djets",nJ);
      fileName[Form("dy%dJets",nJ)]=Form("dy%dJetsToLL",nJ);
   }
   xs["w1Jets"]=6662;
   xs["w2Jets"]=2159;
   xs["w3Jets"]=640;
   xs["w4Jets"]=264;
  
    //
   xs["dy1Jets"]=666;
   xs["dy2Jets"]=215;
   xs["dy3Jets"]=61;
   xs["dy4Jets"]=27;
   
//
   xs["ttSemiLept"]=107.6;
   xs["ttFullLept"]=25.8;
   
   fileName["ttSemiLept"]="ttJetsSemiLept";
   fileName["ttFullLept"]="ttJetsFullLept";

//   
   fileName["TBar_t"]="TBar_t";
   fileName["TBar_s"]="TBar_s";
   fileName["TBar_tW"]="TBar_tW";
   fileName["T_t"]="T_t";
   fileName["T_s"]="T_s";
   fileName["T_tW"]="T_tW";
   
   xs["TBar_t"]=30.7;
   xs["TBar_s"]=1.76;
   xs["TBar_tW"]=11.1;
   xs["T_t"]=56.4;
   xs["T_s"]=3.79;
   xs["T_tW"]=11.1;

//
   fileName["WWJetsTo2L2Nu"]="WWJetsTo2L2Nu";
   fileName["WZJetsTo2L2Q"]="WZJetsTo2L2Q";
   fileName["WZJetsTo3LNu"]="WZJetsTo3LNu";
   fileName["WZJetsTo2QLNu"]="WZJetsTo2QLNu";
   fileName["ZZJetsTo2L2Q"]="ZZJetsTo2L2Q";
   fileName["ZZJetsTo2L2Nu"]="ZZJetsTo2L2Nu";
   fileName["ZZJetsTo4L"]="ZZJetsTo4L";
   
   xs["WWJetsTo2L2Nu"]=7.3; //4.7
   xs["WZJetsTo2L2Q"]=5.995; //1.755;
   xs["WZJetsTo3LNu"]=1.057; //0.8674;
   xs["WZJetsTo2QLNu"]=3.1;
   xs["ZZJetsTo2L2Q"]=0.91;
   xs["ZZJetsTo2L2Nu"]=0.32; //0.28;
   xs["ZZJetsTo4L"]=0.1296;

//
   fileName["TTTT"]    ="TTTT";
   fileName["TTWJets"] ="TTWJets";
   fileName["TTZJets"] ="TTZJets";
   fileName["TTWWJets"]="TTWWJets";
   
   xs["TTTT"]     = 0.000716;
   xs["TTWJets"]  = 0.232; 
   xs["TTZJets"]  = 0.208;
   xs["TTWWJets"] = 0.002037;

//
   xs["QCD"]=1.0;
    
   fileName["QCD_A"]="singleElectronA_QCD";
   fileName["QCD_B"]="singleElectronB_QCD";   
   fileName["QCD_C"]="singleElectronC_QCD"; 
   fileName["QCD_D"]="singleElectronD_QCD"; 
//
   xs["stealthT2_300_200_wwToLNu"]=19.8283;
   fileName["stealthT2_300_200_wwToLNu"]="stealthT2_300_200_wwToLNu";

//
   xs["stealthT2_500_400_wwToLNu"]=0.847051;
   fileName["stealthT2_500_400_wwToLNu"]="stealthT2_500_400_wwToLNu";
 
//
    
    
    xs["RPV_LQD221_M200"]=18.5245;
    fileName["RPV_LQD221_M200"]="RPV_LQD221_M200";
    
    //
    xs["RPV_LQD221_M300"]=1.99608;
    fileName["RPV_LQD221_M300"]="RPV_LQD221_M300";
    
    //
    xs["RPV_LQD221_M400"]=0.35683;;
    fileName["RPV_LQD221_M400"]="RPV_LQD221_M400";
    
    //
    xs["RPV_LQD221_M500"]=0.0855847;
    fileName["RPV_LQD221_M500"]="RPV_LQD221_M500";
    
    //
    xs["RPV_LQD221_M600"]=0.0248009;
    fileName["RPV_LQD221_M600"]="RPV_LQD221_M600";
    
    //
    xs["RPV_LQD221_M700"]=0.0081141;
    fileName["RPV_LQD221_M700"]="RPV_LQD221_M700";
    
    //
    xs["RPV_LQD221_M800"]=0.00289588;
    fileName["RPV_LQD221_M800"]="RPV_LQD221_M800";
    
    //
    xs["RPV_LQD221_M900"]=0.00109501;
    fileName["RPV_LQD221_M900"]="RPV_LQD221_M900";
    
    //
    xs["RPV_LQD221_M1000"]=0.000435488;
    fileName["RPV_LQD221_M1000"]="RPV_LQD221_M1000";

//LQD121
    
   xs["RPV_LQD121_M200"]=18.5245;
   fileName["RPV_LQD121_M200"]="RPV_LQD121_M200";

//
   xs["RPV_LQD121_M300"]=1.99608;
   fileName["RPV_LQD121_M300"]="RPV_LQD121_M300";
 
//
   xs["RPV_LQD121_M400"]=0.35683;;
   fileName["RPV_LQD121_M400"]="RPV_LQD121_M400";

//
   xs["RPV_LQD121_M500"]=0.0855847;
   fileName["RPV_LQD121_M500"]="RPV_LQD121_M500";

//
   xs["RPV_LQD121_M600"]=0.0248009;
   fileName["RPV_LQD121_M600"]="RPV_LQD121_M600";
 
//
   xs["RPV_LQD121_M700"]=0.0081141;
   fileName["RPV_LQD121_M700"]="RPV_LQD121_M700";
 
//
   xs["RPV_LQD121_M800"]=0.00289588;
   fileName["RPV_LQD121_M800"]="RPV_LQD121_M800";

//
   xs["RPV_LQD121_M900"]=0.00109501;
   fileName["RPV_LQD121_M900"]="RPV_LQD121_M900";

//
   xs["RPV_LQD121_M1000"]=0.000435488;
   fileName["RPV_LQD121_M1000"]="RPV_LQD121_M1000";
//
   xs["wJets"]=37509;
   fileName["wJets"]="wjets";

//
   xs["dyJets"]=3504;
   fileName["dyJets"]="dyJetsToLL";

//
   xs["ttJetsMatchDown"]=245.8;
   fileName["ttJetsMatchDown"]="ttJetsMatchDown";

// 
   xs["ttJetsMatchUp"]=245.8;
   fileName["ttJetsMatchUp"]="ttJetsMatchUp";
 
// 
   xs["ttJetsScaleDown"]=245.8;
   fileName["ttJetsScaleDown"]="ttJetsScaleDown";
 
// 
   xs["ttJetsScaleUp"]=245.8;
   fileName["ttJetsScaleUp"]="ttJetsScaleUp";

// 
}

int list_files(TString dirname, TString ext){
    int N=0;
    int maxEvt=100000;
    TString startdir=gSystem->pwd();
    
    TSystemDirectory dir(dirname, dirname);
    TList *files = dir.GetListOfFiles();
    if (files) {
        TSystemFile *file;
        TString fname;
        TIter next(files);
        while ((file=(TSystemFile*)next())) {
            fname = file->GetName();
            if (!file->IsDirectory() && fname.EndsWith(ext)) {
                cout << fname.Data() << " " ;
                TFile *fTmp=new TFile(dirname+"/"+fname,"READ");
                if(fTmp->ReadKeys()!=0){
                    TH1F *h=(TH1F*)fTmp->Get("h_Nevents");
                    N+=h->GetBinContent(1);
                    cout <<h->GetBinContent(1) << endl;
                    delete h;
                }
                //if(N>=maxEvt) return N;
                delete fTmp;
            }
        }
    }
    gSystem->cd(startdir);
    return N;
}

void analyze_file(TString folder){

    TChain *tree = new TChain("tree");
    TString FN=folder+"/hist_analysis*.root";
    cout << "FN: " << FN << endl;
    int Nfiles=tree->Add(FN);
    int NGEN=list_files(folder, ".root");
    NG=NGEN/filter_eff;
    hName["Nevents"]->SetBinContent(1,NGEN/filter_eff);
    cout << "Actual number processed: " << NGEN << endl;
    
    if(isMC) weight=weight/(NGEN/filter_eff);
    
    cout << "actual weight used: "<< weight << endl;

    initialize_tree(tree);

    event_loop(tree);
    delete tree;

}

void initialize_tree(TChain *tree){
	entries=tree->GetEntries();
    
//	tree->SetBranchAddress("eventNo", &eventNo);
//   tree->SetBranchAddress("runNo", &runNo);
//   tree->SetBranchAddress("lumiNo", &lumiNo);
    
//	tree->SetBranchAddress("jets_n",&nJets);
	
//	tree->SetBranchAddress("st",&st);
    
	tree->SetBranchAddress("vertices_n",&nPv);
	tree->SetBranchAddress("met_et",&met);
	tree->SetBranchAddress("met_phi",&met_phi);
    
	tree->SetBranchAddress("puWeight_nom",&puWt_nom);
//	tree->SetBranchAddress("puWeight_up",&puWt_up);
//	tree->SetBranchAddress("puWeight_down",&puWt_down);
	
	tree->SetBranchAddress("electron_px",&electron_px, &b_electron_px);
	tree->SetBranchAddress("electron_py",&electron_py, &b_electron_py);
	tree->SetBranchAddress("electron_pz",&electron_pz, &b_electron_pz);
	tree->SetBranchAddress("electron_e",&electron_e, &b_electron_e);
	tree->SetBranchAddress("electron_charge",&electron_charge, &b_electron_charge);
	
	//loose electrons
	tree->SetBranchAddress("loose_electron_px",&loose_electron_px, &b_loose_electron_px);
	tree->SetBranchAddress("loose_electron_py",&loose_electron_py, &b_loose_electron_py);
	tree->SetBranchAddress("loose_electron_pz",&loose_electron_pz, &b_loose_electron_pz);
	tree->SetBranchAddress("loose_electron_e",&loose_electron_e, &b_loose_electron_e);
	tree->SetBranchAddress("loose_electron_charge",&loose_electron_charge, &b_loose_electron_charge);
	
    
	//get tight muons
	tree->SetBranchAddress("muon_px",&muon_px, &b_muon_px);
	tree->SetBranchAddress("muon_py",&muon_py, &b_muon_py);
	tree->SetBranchAddress("muon_pz",&muon_pz, &b_muon_pz);
	tree->SetBranchAddress("muon_e",&muon_e, &b_muon_e);
	tree->SetBranchAddress("muon_charge",&muon_charge,&b_muon_charge);
    if(iso){
        tree->SetBranchAddress("muon_noIso_px",&muon_noIso_px, &b_muon_noIso_px);
        tree->SetBranchAddress("muon_noIso_py",&muon_noIso_py, &b_muon_noIso_py);
        tree->SetBranchAddress("muon_noIso_pz",&muon_noIso_pz, &b_muon_noIso_pz);
        tree->SetBranchAddress("muon_noIso_e",&muon_noIso_e, &b_muon_noIso_e);
        tree->SetBranchAddress("muon_noIso_PFIso04",&muon_noIso_PFIso04, &b_muon_noIso_PFIso04);
        tree->SetBranchAddress("muon_noIso_PFIso03",&muon_noIso_PFIso03, &b_muon_noIso_PFIso03);
        tree->SetBranchAddress("muon_noIso_charge",&muon_noIso_charge,&b_muon_noIso_charge);
    }
	//get loose muons
	
	tree->SetBranchAddress("loose_muon_px",&loose_muon_px, &b_loose_muon_px);
	tree->SetBranchAddress("loose_muon_py",&loose_muon_py, &b_loose_muon_py);
	tree->SetBranchAddress("loose_muon_pz",&loose_muon_pz, &b_loose_muon_pz);
	tree->SetBranchAddress("loose_muon_e",&loose_muon_e, &b_loose_muon_e);
	tree->SetBranchAddress("loose_muon_charge",&loose_muon_charge,&b_loose_muon_charge);
		
	tree->SetBranchAddress("jet_px",&jet_px, &b_jet_px);
	tree->SetBranchAddress("jet_py",&jet_py, &b_jet_py);
	tree->SetBranchAddress("jet_pz",&jet_pz, &b_jet_pz);
	tree->SetBranchAddress("jet_e",&jet_e, &b_jet_e);
	tree->SetBranchAddress("jet_unc",&jet_unc, &b_jet_unc);
//	tree->SetBranchAddress("jet_bTagL",&jet_bTagL, &b_jet_bTagL);
	tree->SetBranchAddress("jet_bTagM",&jet_bTagM, &b_jet_bTagM);
//	tree->SetBranchAddress("jet_bTagT",&jet_bTagT, &b_jet_bTagT);
	tree->SetBranchAddress("jet_algFlavor",&jet_algFlavor, &b_jet_algFlavor);
	tree->SetBranchAddress("jet_phyFlavor",&jet_phyFlavor, &b_jet_phyFlavor);
    
	/*
	tree->SetBranchAddress("w_px",&w_px, &b_w_px);
	tree->SetBranchAddress("w_py",&w_py, &b_w_py);
	tree->SetBranchAddress("w_pz",&w_pz, &b_w_pz);
	tree->SetBranchAddress("w_e",&w_e, &b_w_e);
	
	tree->SetBranchAddress("z_px",&z_px, &b_z_px);
	tree->SetBranchAddress("z_py",&z_py, &b_z_py);
	tree->SetBranchAddress("z_pz",&z_pz, &b_z_pz);
	tree->SetBranchAddress("z_e",&z_e, &b_z_e);
	*/
	tree->SetBranchAddress("top_px",&top_px, &b_top_px);
	tree->SetBranchAddress("top_py",&top_py, &b_top_py);
	tree->SetBranchAddress("top_pz",&top_pz, &b_top_pz);
	tree->SetBranchAddress("top_e",&top_e, &b_top_e);
    
	
	
	
}
void event_loop(TChain *tree){
	int maxEvt=10000;
	//entries=maxEvt;
    
    double time_evt=0;
	double W=weight;
	cout << "Entries: " << entries << endl;
	for (int iEvt=0; iEvt<entries; iEvt++){
		if(iEvt%1000000==0)cout << iEvt << endl;
	
		tree->GetEntry(iEvt);
        
		//count_muons();
		//count_electrons();
        //cout << "fill_jets: " << endl;
        fill_jets();
        //cout << "fill leptons: " << endl;
        loop_leptons();
        //cout << "calc st: " << endl;
        calcSt();
        //cout << "calc mass: " << endl;
        calDilepMass();
        //cout << "calc muJ: " << endl;
        calcMuJ();
		weight=W;
		general_plots(); 

		//cout << "weight before eff: " << weight << endl;
        fill_isolation_efficiency();
		if(preselection_selection()==false) continue;
        // pileup_reweight();
        //print_event();
        //cout << "fill LQ: " << endl;
        fill_LQ();
        //cout << "fill OPT: " << endl;
        fill_optimization();
        fill_analysis();
	}
    
    
}

//-------------------------------------------------------------------------------

void print_event(){
    
    cout << "hT: " << hT << " sT: " << sT << " met: " << met << endl;
    cout << "nElecrons: " << nElectrons << " nMuons " << nMuons << endl;
    cout << "dimuon: " << dimuon_enriched << " " << dimuon_mass << endl;
    cout << "e,mu: " << tt_enriched << " " << endl;
    cout << "mu_j: " << m_muj << endl;
    cout << "dielectron: " << dielectron_enriched << " " << dielectron_mass << endl;
    
    for(int i=0; i<muonCollection.size();i++){
        cout << "muon pt: " << muonCollection.at(i).Perp() << endl;
    }
    
    for(int i=0; i<jetCollection.size();i++){
        cout << "jet pt: " << jetCollection.at(i).Perp() << endl;
    }
    
    cout << "nJets: " << nJets << "  " << nBtags << endl << endl;

}

void fill_jets(){
    TLorentzVector pJ;
    hT=0;
    
    int nB=0;
    jetCollection.clear();
    
    B_weight=1;
    nJets=0;
    
    
    for(unsigned int ijet=0; ijet<jet_px->size(); ijet++){
        pJ.SetPxPyPzE(jet_px->at(ijet), jet_py->at(ijet), jet_pz->at(ijet),jet_e->at(ijet));
        if(jetCollection.size()==0 && pJ.Perp()<jetPt1_cut) continue;
        if(jetCollection.size()==1 && pJ.Perp()<jetPt2_cut) continue;
        if(jetCollection.size()>=2 && pJ.Perp()<jetPt_cut) continue;
        
        if(jet_bTagM->at(ijet))nB++;
        //b-tag SF
        if(isMC){
            double SF=getSF(pJ, TMath::Abs(jet_algFlavor->at(ijet)), "CSVM","");
            double eff=read_btag_efficiency(pJ,TMath::Abs(jet_algFlavor->at(ijet)) ,"M");
            
            if(jet_bTagM->at(ijet)==1) {
                B_weight*=SF*eff/eff;
            }
            else {
                B_weight*=(1-SF*eff)/(1-eff);
            }
            
        }//MC     

        nJets++;
        jetCollection.push_back(pJ);
        //if(ijet==0 || ijet==1)hT+=pJ.Perp(); only for LQ analysis
        hT+=pJ.Perp();
    }
    nBtags=nB;
}

void general_plots(){
	float MuW=1;
    float EleW=1;
    
    if(isMC){
        top_cor();
        MuW= muonSF(muonCollection);
        EleW=electronSF(electronCollection);
    }
	float EvtW=MuW*EleW*TTbar_corr;
	
	if(jetCollection.size()>=1){
		//cout << " 1 jet" << endl; 
		if(jetCollection.at(0).Perp()>50){
			//cout << " pt threshold" << endl; 

			if(electronCollection.size()>=1){
				//cout << " 1 electron" << endl; 
				//cout << "st: " << sT << endl; 
				if(met>150){
					hName["h_st_1El_met150"]->Fill(sT,EvtW); 
				}
				else {
					hName["h_st_1El_met0"]->Fill(sT,EvtW); 
				}
				
			}
		}
	}
	
	if(m_muj>0)hName["h_muj_all"]->Fill(m_muj,EvtW);
	if(m_ej>0)hName["h_ej"]->Fill(m_ej,EvtW);
	
}

void fill_isolation_efficiency(){
    //if(dimuon_enriched==0) return;
    //if(SS==1) return;
    if(met>100)return;
    if(nElectrons!=1) return;
    if(sT<stMin) return;
    
    if(nBtags>=1){
        if(nMuons_R03==1){
            hName["nJets_Iso_total_1btag"]->Fill(nJets,1);
        }
        if(nMuons_R04==1){
            hName["nJets_Iso_pass_1btag"]->Fill(nJets,1);
        }
    }
    
    if(nBtags>=2){
        if(nMuons_R03==1){
            hName["nJets_Iso_total"]->Fill(nJets,1);
            hName["h_met_Iso_total"]->Fill(met,1);
            hName["h_muonpt_Iso_total"]->Fill(muonCollection_noIso.at(0).Perp(),1);
        }
        if(nMuons_R04==1){
            hName["nJets_Iso_pass"]->Fill(nJets,1);
            hName["h_met_Iso_pass"]->Fill(met,1);
            hName["h_muonpt_Iso_pass"]->Fill(muonCollection.at(0).Perp(),1);
        }
    }
}

void fill_analysis(){
    
    float MuW=1;
    float EleW=1;
    
    if(isMC){
        top_cor();
        MuW= muonSF(muonCollection);
        EleW=electronSF(electronCollection);
    }
    
    float EvtW=MuW*EleW*B_weight*TTbar_corr;
    
    float st_=sT;
    if(st_>1600)st_=1601;
    
    int nb=nBtags;
    if(nb>=2)nb=2;
    
	if(SS) return;
	if(sT>875 && nJets>=5 && nBtags>=1 && dimuon_enriched && dimuon_mass>130) {
		for(int imet=1; imet<=hName["h_metmax"]->GetNbinsX(); imet++){
            double binCenter=hName["h_metmax"]->GetBinCenter(imet);
            if(met<binCenter)hName["h_metmax"]->Fill(binCenter,weight*EvtW);
        }
	}
	
	if(sT>875 && nJets>=5 && nBtags>=1 && dimuon_enriched && dimuon_mass>105 && met<100) {
		for(int imet=1; imet<=hName["h_dimuomassmin"]->GetNbinsX(); imet++){
            double binCenter=hName["h_dimuomassmin"]->GetBinCenter(imet);
            if(dimuon_mass>binCenter)hName["h_dimuomassmin"]->Fill(binCenter,weight*EvtW);
        }
	}
	
    if(met>100) return;
    if(dimuon_enriched && dimuon_mass>50 && dimuon_mass<130)hName[Form("h_nJets_%dbtags_lowM_noST",nb)]->Fill(nJets,weight*EvtW);
    if(sT<stMin) return;


    
    if(dimuon_enriched && dimuon_mass>50 && dimuon_mass<130)hName[Form("h_nJets_%dbtags_lowM",nb)]->Fill(nJets,weight*EvtW);
    if(tt_enriched)hName[Form("h_nJets_%dbtags_MuE",nb)]->Fill(nJets,weight*EvtW);
    if(tt_enriched)hName["h_nJets_MuE"]->Fill(nJets,weight*EvtW); 
    if(dimuon_enriched)hName[Form("h_nJets_%dbtags",nb)]->Fill(nJets,weight*EvtW);

    
    if(nBtags>=1  && tt_enriched){
        float st2_=sT;
        if(sT>800) st2_=801.;
        hName[Form("h_nJets%d_MuE_st",nJets)]->Fill(st2_,weight*EvtW);
    }
    
    if(nBtags>=1 && dimuon_enriched && dimuon_mass>130){
        //if(nJets==5) print_event();
        hName[Form("h_nJets%d_st",nJets)]->Fill(st_,weight*EvtW);
        hName["h_nJets"]->Fill(nJets,weight*EvtW);
    }
    
}


void top_cor(){
	TTbar_corr=1;
    
	if(!isMC) return;
    if(!reweight_top) return;
	for(int i=0; i<top_px->size(); i++){
		TLorentzVector top;
		top.SetPxPyPzE(top_px->at(i), top_py->at(i), top_pz->at(i), top_e->at(i));
		double	 pt_top=top.Perp();
		if(pt_top>400)pt_top=400;
		TTbar_corr=TTbar_corr*expo_cor->Eval(pt_top);
	}
    TTbar_corr=TMath::Sqrt(TTbar_corr);

}

void fill_optimization(){
    float st_=sT;
    if(st_>1500)st_=1499;
   // cout << "optimization: " << endl;
    //if(dimuon_enriched==0) return;
    if(dimuon_enriched && dimuon_mass>130)hName[Form("h_met_nJets%d",nJets)]->Fill(met,weight);
    if(nJets>=5 && !Z_enriched && sT>300)hName["h_met"]->Fill(met);
    if(met>100)return;
    
    if(nBtags>=1 && !Z_enriched && SS==1){
        for(int ist=1; ist<=hName[Form("h_nJets%d_stminSS",nJets)]->GetNbinsX(); ist++){
            double binCenter=hName[Form("h_nJets%d_stminSS",nJets)]->GetBinCenter(ist);
            if(sT>binCenter)hName[Form("h_nJets%d_stminSS",nJets)]->Fill(binCenter,weight);
        }
        
    }
    
    if(nJets>=5 && nBtags>=1 && !Z_enriched && sT>300){
        float MuW= muonSF(muonCollection);
        float EleW=electronSF(electronCollection);
        
        if(dimuon_enriched)hName["nJets_muSF"]->Fill(5,MuW);
        if(dielectron_enriched)hName["nJets_eleSF"]->Fill(5,EleW);
        if(dielectron_enriched)hName["nJets_noSF"]->Fill(5,1);
        hName["nJets_EleMuSF"]->Fill(5,EleW/MuW);

    }
    
    if(nBtags>=1 && sT>300 && tt_enriched){
        float st2_=sT;
        if(sT>800) st2_=799.;
        hName[Form("h_nJets%d_MuE_st",nJets)]->Fill(st2_,weight);
    }

    if(nBtags>=1 && sT>300 && tt_enriched){
        double pt=muonCollection.at(0).Perp();
        if(SS==0 && nJets>=5 && pt>150)hName["h_misQID"]->Fill("OS, N_{jets}#geq5, p_{T} > 150",weight);
        if(SS==0 && nJets<5 && pt>150)hName["h_misQID"]->Fill("OS, N_{jets}<5, p_{T} > 150",weight);

        if(SS==1 && nJets>=5 && pt>150)hName["h_misQID"]->Fill("SS, N_{jets}#geq5, p_{T} > 150",weight);
        if(SS==1 && nJets<5 && pt>150)hName["h_misQID"]->Fill("SS, N_{jets}<5, p_{T} > 150",weight);
        
        if(SS==0 && nJets>=5 && pt<=150)hName["h_misQID"]->Fill("OS, N_{jets}#geq5, p_{T} #leq 150",weight);
        if(SS==0 && nJets<5 && pt<=150)hName["h_misQID"]->Fill("OS, N_{jets}<5, p_{T} #leq 150",weight);
        
        if(SS==1 && nJets>=5 && pt<=150)hName["h_misQID"]->Fill("SS, N_{jets}#geq5, p_{T} #leq 150",weight);
        if(SS==1 && nJets<5 && pt<=150)hName["h_misQID"]->Fill("SS, N_{jets}<5, p_{T} #leq 150",weight);
        
    }
    float MuW= muonSF(muonCollection);
    if(nJets>=5 && nBtags>=1 && sT>300){
        //float MuW= muonSF(muonCollection);
        float EleW=electronSF(muonCollection);
        hName["h_mumu_muSF"]->Fill(dimuon_mass,MuW);
        hName["h_mumu_noSF"]->Fill(dimuon_mass,1);
        hName["h_mumu_eleSF"]->Fill(dimuon_mass,EleW);

    }
    //print_event();
    if(SS==1) return;
    
    if(nJets>=5 && dimuon_enriched){
        if(nBtags==0)hName["h_mumu_0btag"]->Fill(dimuon_mass,weight);
        if(nBtags>=1)hName["h_mumu_1btag"]->Fill(dimuon_mass,weight);

    }
    
    //if(nBtags>=1 && dimuon_mass>105)hName[Form("h_nJets%d_st",nJets)]->Fill(sT,weight);
    if(nBtags>=1 && dimuon_enriched && !Z_enriched){
   
        hName[Form("h_nJets%d_st_10GeV",nJets)]->Fill(sT,weight);
        if(dimuon_mass>130)hName[Form("h_nJets%d_st",nJets)]->Fill(st_,weight);
        hName[Form("h_nJets%d_st_uW",nJets)]->Fill(sT);
        for(int ist=1; ist<=hName[Form("h_nJets%d_stmin",nJets)]->GetNbinsX(); ist++){
            double binCenter=hName[Form("h_nJets%d_stmin",nJets)]->GetBinCenter(ist)-hName[Form("h_nJets%d_stmin",nJets)]->GetBinWidth(ist)/2;
            if(sT>binCenter)hName[Form("h_nJets%d_stmin_nominal",nJets)]->Fill(binCenter,weight*MuW*B_weight);
            if(dimuon_mass>130){
                if(sT>binCenter)hName[Form("h_nJets%d_stmin",nJets)]->Fill(binCenter,weight*MuW*B_weight);
                if(sT>binCenter)hName[Form("h_nJets%d_stmin_uW",nJets)]->Fill(binCenter,1);
                if(lept_pt_binlow && sT>binCenter)hName[Form("h_nJets%d_stmin_lowpt",nJets)]->Fill(binCenter,weight);
                if(lept_pt_binhigh && sT>binCenter)hName[Form("h_nJets%d_stmin_highpt",nJets)]->Fill(binCenter,weight);
            }
        }
        
    }
    if(nJets<5) return;

    //fill_optimization_hist(sT, dimuon_mass, m_muj, weight);

}

void fill_LQ(){
    //cout << "fill LQ: " << endl;
    float MLQ[]={300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000};
    float st_bins[]={380, 460, 540, 615, 685, 755, 820, 880, 935, 990, 1040, 1090, 1135, 1175, 1210};
    float mass_bins[]={100, 115, 125, 140, 150, 165, 175, 185, 195, 205, 215, 220, 230, 235, 245};
    float mj[]={115, 115, 120, 135, 155, 180, 210, 250, 295, 345, 400, 465, 535, 610, 690};
    int N=sizeof(st_bins)/sizeof(float);
    
    //cout << m_muj << " " << weight << endl;
    
   
    if(dimuon_enriched)hName["h_mumu"]->Fill(dimuon_mass,weight);
    if(dimuon_enriched)hName["h_muj"]->Fill(m_muj,weight);
    if(dimuon_enriched)hName["h_st"]->Fill(sT,weight);

    
    for(int i=0; i<N; i++){
        if(dimuon_enriched && (sT>st_bins[i]) && (dimuon_mass>mass_bins[i]) && (m_muj>mj[i]) ) {
          hName["h_mumujj_LQ"]->Fill(MLQ[i]-1,weight);
          if(MLQ[i]==900)hName["nJets_LQ900"]->Fill(nJets,weight);

        }
    }
    
}

void loop_leptons(){
    //no cuts implemented yet for leptons. Do at ntuple level
    TLorentzVector pL;
    lepPt=0;
    muonCollection.clear();
    electronCollection.clear();
    muonCollection_noIso.clear();
    muonQ.clear();
    nMuons_R03=0;
    nMuons_R04=0;
    
    lept_pt_binlow=0;
    lept_pt_binhigh=0;
    
    ElectronQ.clear();
    for(unsigned int i=0; i<electron_px->size(); i++){
        pL.SetPxPyPzE(electron_px->at(i),electron_py->at(i),electron_pz->at(i),electron_e->at(i));
        if(pL.Perp()<lepton_Pt_cut) continue;
        electronCollection.push_back(pL);
        lepPt+=pL.Perp();
        ElectronQ.push_back(electron_charge->at(i));

    }
    
    for(unsigned int i=0; i<muon_px->size(); i++){
        pL.SetPxPyPzE(muon_px->at(i),muon_py->at(i),muon_pz->at(i),muon_e->at(i));
        double pt=pL.Perp();
		double Energy=pL.E(); 
        double zeta=-0.05;
        double zeta_low=-0.002;
        double eta=TMath::Abs(pL.Eta()); 
        //zeta*=2;
        //zeta_low*=2;
        
        //if(pt>200)pL.SetPerp((1+zeta*pt/1000)*pt);
        //if(pt<200)pL.SetPerp((1+zeta_low)*pt);
		//if(eta>1.44)pL.SetE((1.-0.004)*Energy); 
		//else pL.SetE((1.-0.008)*Energy);
        if(pL.Perp()<lepton_Pt_cut) continue;
        muonCollection.push_back(pL);
        if(!iso)lepPt+=pL.Perp();
        muonQ.push_back(muon_charge->at(i));
    }
    nMuons=muonCollection.size();
    
    if(iso){
        
    for(unsigned int i=0; i<muon_noIso_px->size(); i++){
        pL.SetPxPyPzE(muon_noIso_px->at(i),muon_noIso_py->at(i),muon_noIso_pz->at(i),muon_noIso_e->at(i));
        double pt=pL.Perp();
        if(pL.Perp()<lepton_Pt_cut) continue;
        // continue;
        if(muon_noIso_PFIso03->at(i) < 0.12){
            muonCollection_noIso.push_back(pL);
            lepPt+=pL.Perp();
            nMuons_R03++;
        }
        if(muon_noIso_PFIso04->at(i) < 0.12){
            nMuons_R04++;
        }
        
    }
         /*
        
        for(unsigned int i=0; i<loose_muon_px->size(); i++){
            pL.SetPxPyPzE(loose_muon_px->at(i),loose_muon_py->at(i),loose_muon_pz->at(i),loose_muon_e->at(i));
            double pt=pL.Perp();
            if(pL.Perp()<lepton_Pt_cut) continue;
            muonCollection_noIso.push_back(pL);
            lepPt+=pL.Perp();
            nMuons_noIso++;
        }*/
    }
    nMuons=muonCollection.size();
    
    if(muonCollection.size()==2){
        //cout << " pt1: " << muonCollection.at(0).Perp() << " " << muonCollection.at(1).Perp() << endl;
        hName2D["h_pt1_pt2"]->Fill(muonCollection.at(0).Perp(),muonCollection.at(1).Perp());
        if(muonCollection.at(0).Perp()>lepton_Pt_cut && muonCollection.at(1).Perp()>lepton_Pt_cut )lept_pt_binhigh=true;
        if(muonCollection.at(0).Perp()<lepton_Pt_cut || muonCollection.at(1).Perp()<lepton_Pt_cut )lept_pt_binlow=true;
        if(lept_pt_binlow)hName2D["h_pt1_pt2_low"]->Fill(muonCollection.at(0).Perp(),muonCollection.at(1).Perp());
        if(lept_pt_binhigh)hName2D["h_pt1_pt2_high"]->Fill(muonCollection.at(0).Perp(),muonCollection.at(1).Perp());
        
    }
    
    nElectrons=electronCollection.size();
}

void calcSt(){
	//cout << "ht: " << hT << " " << leptPt << endl; 
    sT=hT+lepPt;
}

bool preselection_selection(){
    //if(nLooseElectrons>=1) return false; // no additional loose electrons
    //if(nLooseMuons>=1) return false;
    if(nJets<2) return false;
    if(nJets>=nJetmax) nJets=nJetmax;
    //if(sT<stMin) return false;
    
    //3 cases: ee,mumu,emu
    
    dimuon_enriched=false;
    dielectron_enriched=false;
    tt_enriched=false;
    Z_enriched=false;
    
    SS=false;
    
    //require exactly 2 muons
	if(nElectrons==0 && nMuons==2 && dimuon_mass>min_diLeptonMass) dimuon_enriched=true;
    //require exactly 2 electrons
    if(nMuons==0 && nElectrons==2 && dielectron_mass>min_diLeptonMass) dielectron_enriched=true;
    //require 1 mu,e
	if(nMuons==1 && nElectrons==1) tt_enriched=true;

	if(dimuon_enriched && muonQ.at(0)*muonQ.at(1)==1) SS=true;
    if(dielectron_enriched && ElectronQ.at(0)*ElectronQ.at(1)==1) SS=true;
	if(tt_enriched && muonQ.at(0)*ElectronQ.at(0)==1) SS=true;
    
	if(dimuon_enriched && SS==false){
		if(fabs(dimuon_mass-Zmass)<ZVeto) Z_enriched=true;
	}
    if(dielectron_enriched && SS==false){
		if(fabs(dielectron_mass-Zmass)<ZVeto) Z_enriched=true;
	}
    
    if(dimuon_enriched==0 && dielectron_enriched==0 && tt_enriched==0) return false; //not valid event
    
}

bool count_electrons(){
	TLorentzVector pLoose;
	TLorentzVector pTight;
	
	nLooseElectrons=0;
	nElectrons=electron_px->size();
	
	int nMatched=0;
	std::vector<int> matched;
	
	for(int iloose=0; iloose<loose_electron_px->size();iloose++){
		pLoose.SetPxPyPzE(loose_electron_px->at(iloose),loose_electron_py->at(iloose),loose_electron_pz->at(iloose),loose_electron_e->at(iloose));
		for(int itight=0; itight<electron_px->size();itight++){
			pTight.SetPxPyPzE(electron_px->at(itight),electron_py->at(itight),electron_pz->at(itight),electron_e->at(itight));
			bool skip=false;
			for(int ii=0; ii<matched.size();ii++){
				if(itight==matched[ii]) skip=true;
			}
			if(skip) continue;
			
			if(pLoose.DeltaR(pTight)<0.01){
				matched.push_back(itight);
				nMatched++;
			}
		}
	}
	nLooseElectrons=loose_electron_px->size()-nMatched;
}

bool count_muons(){
	TLorentzVector pLoose;
	TLorentzVector pTight;
	
	nLooseMuons=0;
	nMuons=muon_px->size();
	
	int nMatched=0;
	std::vector<int> matched;
	
	for(int iloose=0; iloose<loose_muon_px->size();iloose++){
		pLoose.SetPxPyPzE(loose_muon_px->at(iloose),loose_muon_py->at(iloose),loose_muon_pz->at(iloose),loose_muon_e->at(iloose));
		for(int itight=0; itight<muon_px->size();itight++){
			pTight.SetPxPyPzE(muon_px->at(itight),muon_py->at(itight),muon_pz->at(itight),muon_e->at(itight));
			bool skip=false;
			for(int ii=0; ii<matched.size();ii++){
				if(itight==matched[ii]) skip=true;
			}
			if(skip) continue;
			
			if(pLoose.DeltaR(pTight)<0.01){
				matched.push_back(itight);
				nMatched++;
			}
		}
	}
	nLooseMuons=loose_muon_px->size()-nMatched;
     
    
    
}

void calcMuJ(){
    bool print=false;
    //2 config
    //config 1
    //LQ mu1,j1
    //LQ~ mu2,j2
    
    //config 2
    //LQ mu1,j2
    //LQ~ mu2,j1
    m_muj=0;
    float MLQ1[]={0,0};//two comb for LQ mass
    float MLQ2[]={0,0};
    
    for(int imu=0; imu<muonCollection.size(); imu++){
        TLorentzVector pTmp(0,0,0,0);
        pTmp+=muonCollection.at(imu);

        for(int ijet=0; ijet<jetCollection.size(); ijet++){
            if(ijet>1) continue;
            pTmp+=jetCollection.at(ijet);
            if(print)cout << "imu: " << imu  << " ijet " << ijet <<" " << pTmp.M() << endl;
            if(imu==0 && ijet==0)MLQ1[0]=pTmp.M();
            if(imu==1 && ijet==1)MLQ1[1]=pTmp.M();

            if(imu==0 && ijet==1)MLQ2[0]=pTmp.M();
            if(imu==1 && ijet==0)MLQ2[1]=pTmp.M();
        }
    }
    
    if(print)cout << "MLQ1: " << MLQ1[0] << " MLQ1~: " << MLQ1[1] << endl;
    if(print)cout << "MLQ2: " << MLQ2[0] << " MLQ2~: " << MLQ2[1] << endl;

    float delta1=fabs(MLQ1[0]-MLQ1[1]);
    float delta2=fabs(MLQ2[0]-MLQ2[1]);

    if(print)cout << "delta1: " << delta1 << " delta2: " << delta2 << endl;
    if(TMath::Min(delta1,delta2)==delta1) {
        if(print)cout << "config 1: " << endl;
        if(print)cout << "MLQ: " << TMath::Min(MLQ1[0],MLQ1[1]) << endl;
        m_muj=TMath::Min(MLQ1[0],MLQ1[1]);
    }
    else {
        if(print)cout << "config 2: " << endl;
        if(print)cout << "MLQ: " << TMath::Min(MLQ2[1],MLQ2[1]) << endl;
        m_muj=TMath::Min(MLQ2[1],MLQ2[1]);
    }
    
	//electrons
	m_ej=0;
	MLQ1[0]=0; 
	MLQ1[1]=0; 
	MLQ2[0]=0; 
	MLQ2[1]=0; 
    
    for(int imu=0; imu<electronCollection.size(); imu++){
        TLorentzVector pTmp(0,0,0,0);
        pTmp+=electronCollection.at(imu);
		
        for(int ijet=0; ijet<jetCollection.size(); ijet++){
            if(ijet>1) continue;
            pTmp+=jetCollection.at(ijet);
            if(print)cout << "imu: " << imu  << " ijet " << ijet <<" " << pTmp.M() << endl;
            if(imu==0 && ijet==0)MLQ1[0]=pTmp.M();
            if(imu==1 && ijet==1)MLQ1[1]=pTmp.M();
			
            if(imu==0 && ijet==1)MLQ2[0]=pTmp.M();
            if(imu==1 && ijet==0)MLQ2[1]=pTmp.M();
        }
    }
    
	delta1=fabs(MLQ1[0]-MLQ1[1]);
	delta2=fabs(MLQ2[0]-MLQ2[1]);

    if(TMath::Min(delta1,delta2)==delta1) {
        m_ej=TMath::Min(MLQ1[0],MLQ1[1]);
    }
    else {
        m_ej=TMath::Min(MLQ2[1],MLQ2[1]);
    }
	
    if(print)cout << endl;
}

void calDilepMass(){
    dimuon_mass=0;
    dielectron_mass=0;
    if(muonCollection.size()==2){
        TLorentzVector pM1=muonCollection.at(0)+muonCollection.at(1);
        dimuon_mass=pM1.M();
        
    }
    
    if(electron_px->size()==2){
        TLorentzVector pM1;
        TLorentzVector pM2;
        pM1.SetPxPyPzE(electron_px->at(0),electron_py->at(0),electron_pz->at(0),electron_e->at(0));
        pM2.SetPxPyPzE(electron_px->at(1),electron_py->at(1),electron_pz->at(1),electron_e->at(1));
        pM1+=pM2;
        dielectron_mass=pM1.M();
    }

}

float muonSF(std::vector<TLorentzVector> muP4){
   //if(muP4.size() > 1){cout << muP4.size() << endl;}
   Double_t muPt,muEta; Double_t SF1=1.,SF2=1.;
   for(unsigned int m=0; m<muP4.size(); m++){
      muPt = muP4[m].Perp();
      muEta = muP4[m].Eta();
      if(m==0){
         SF1 = calcMuIdSF(muEta, muPt);
      }
      if(m==1){
         SF2 = calcMuIdSF(muEta, muPt);
      }
   }
   return SF1*SF2;
}

float  calcMuIdSF(Double_t eta, Double_t pt){
   float abseta = abs(eta);
   float SF = 1.0;
   if (abseta >= 0.0 && abseta < 0.90){
      if (pt > 10. && pt < 20.){ SF = 0.9327;}
      else if (pt > 20.  && pt < 25.){ SF = 0.9639;}
      else if (pt > 25.  && pt < 30.){ SF = 0.9910;}
      else if (pt > 30.  && pt < 35.){ SF = 0.9881;}
      else if (pt > 35.  && pt < 40.){ SF = 0.9876;}
      else if (pt > 40.  && pt < 50.){ SF = 0.9866;}
      else if (pt >= 50. && pt < 60.){ SF = 0.9866;}
      else if (pt >= 60. && pt < 90.){ SF = 0.9895;}
      else if (pt >= 90. && pt < 140.){ SF = 1.0049;}
      else if (pt >= 140.&& pt < 300.){ SF = 1.0283;}
      else if (pt >= 300.&& pt < 500.){ SF = 1.0198;}
      else {SF = 1.0198;}
   }
   else if (abseta >= 0.90 && abseta < 1.2){
      if (pt > 10. && pt < 20.){ SF = 0.9393;}
      else if (pt > 20.  && pt < 25.){ SF = 0.9759;}
      else if (pt > 25.  && pt < 30.){ SF = 0.9950;}
      else if (pt > 30.  && pt < 35.){ SF = 0.9892;}
      else if (pt > 35.  && pt < 40.){ SF = 0.9895;}
      else if (pt >= 40. && pt < 50.){ SF = 0.9878;}
      else if (pt >= 50. && pt < 60.){ SF = 0.9900;}
      else if (pt >= 60. && pt < 90.){ SF = 0.9855;}
      else if (pt >= 90. && pt < 140.){ SF = 1.0115;}
      else if (pt >= 140. && pt < 300.){ SF = 0.9525;}
      else if (pt >= 300. && pt < 500.){ SF = 1.0078;}
      else {SF = 1.0078;}
   }                 
   else if (abseta >= 1.2 and abseta < 2.1){
      if (pt > 10. && pt < 20.){ SF = 0.9922;}
      else if (pt > 20.  && pt < 25.){ SF = 0.9989;}
      else if (pt > 25.  && pt < 30.){ SF = 1.0070;}
      else if (pt > 30.  && pt < 35.){ SF = 1.0044;}
      else if (pt > 35.  && pt < 40.){ SF = 1.0003;}
      else if (pt >= 40.  &&  pt < 50.){ SF = 0.9995;}
      else if (pt >= 50.  && pt < 60.){ SF = 0.9988;}
      else if (pt >= 60.  && pt < 90.){ SF = 0.9946;}
      else if (pt >= 90.  && pt < 140.){ SF = 1.0191;}
      else if (pt >= 140. && pt < 300.){ SF = 1.0164;}
      else if (pt >= 300. && pt < 500.){ SF = 0.6173;}
      else { SF = 0.6173;}
   }
   return SF;
}

float electronSF( std::vector<TLorentzVector> eleP4){
   //if(eleP4.size() > 1){cout << eleP4.size() << endl;}
   Double_t elePt,eleEta; Double_t SF1=1.,SF2=1.;
   for(unsigned int e=0; e<eleP4.size(); e++){
      elePt = eleP4[e].Perp();
      eleEta = eleP4[e].Eta();
      if(e==0){
         SF1 = calcEleIdSF(eleEta, elePt) * calcEleTrigSF(eleEta, elePt);
      }
      if(e==1){
         SF2 = calcEleIdSF(eleEta, elePt);// * calcEleTrigSF(eleEta, elePt);//since event only passes single ele trigger
      }
   }
   //if(electron_px->size()== 2){
   //cout << "SF1 = " << SF1 << ", SF2 = " << SF2 << endl;
   //}
   return SF1*SF2;
   
}

float calcEleTrigSF(Double_t eta, Double_t pt){
   float abseta = abs(eta);
   float SF = 1.0;

   if(abseta >= 0.0 && abseta < 0.80) {
      if (pt >= 20. && pt < 30.){SF = 0.695;}
      else if (pt >= 30. && pt < 40.){SF = 0.984;}
      else if (pt >= 40. && pt < 50.){SF = 0.999;}  
      else if (pt >  50. && pt < 200.){SF = 0.999;}
      else {SF = 0.999;}
   }
   else if (abseta >= 0.80 and abseta < 1.48){
      if (pt >= 20. && pt < 30.){SF = 0.462;}
      else if (pt >= 30. && pt < 40.){SF = 0.967;}
      else if (pt >= 40. && pt < 50.){SF = 0.988;}  
      else if (pt >  50. && pt < 200.){SF = 0.988;}
      else {SF = 0.988;}
   }
   else if (abseta >= 1.48 and abseta <  2.50){
      if (pt >= 20. && pt < 30.){SF = 0.804;}
      else if (pt >= 30. && pt < 40.){SF = 0.991;}
      else if (pt >= 40. && pt < 50.){SF = 1.018;}  
      else if (pt >  50. && pt < 200.){SF = 0.977;}
      else {SF = 0.977;}
   }
   return SF;
}

float calcEleIdSF(Double_t eta, Double_t pt){
   float abseta = abs(eta);
   float SF = 1.0;

   if(abseta >= 0.0 && abseta < 0.80) {
      if (pt >= 20. && pt < 30.){SF = 0.972;}
      else if (pt >= 30. && pt < 40.){SF = 0.950;}
      else if (pt >= 40. && pt < 50.){SF = 0.966;}  
      else if (pt >  50. && pt < 200.){SF = 0.961;}
      else {SF = 0.961;}
   }
   else if (abseta >= 0.80 and abseta < 1.48){
      if (pt >= 20. && pt < 30.){SF = 0.928;}
      else if (pt >= 30. && pt < 40.){SF =  0.957;}
      else if (pt >= 40. && pt < 50.){SF =  0.961;}  
      else if (pt >  50. && pt < 200.){SF = 0.963;}
      else {SF = 0.963;}
   }
   else if (abseta >= 1.48 and abseta <  2.50){
      if (pt >= 20. && pt < 30.){SF = 0.834;}
      else if (pt >= 30. && pt < 40.){SF =  0.922;}
      else if (pt >= 40. && pt < 50.){SF =  0.941;}  
      else if (pt >  50. && pt < 200.){SF = 0.971;}
      else {SF = 0.971;}
   }
   return SF;
}



void load_btagEff(){
    //Need to have std::map<TString,TEfficiency*> effName defined in header
   cout << "loading efficiency map  ..." << endl;
    TString fname=output_file_name;
    fname.ReplaceAll("_jecP","");
    fname.ReplaceAll("_jecM","");
    
    fname.ReplaceAll("dy1Jets","dy");
    fname.ReplaceAll("dy2Jets","dy");
    fname.ReplaceAll("dy3Jets","dy");
    fname.ReplaceAll("dy4Jets","dy");
    
    fname.ReplaceAll("_LQD221_M","");
    
    
    if(fname=="ttZ") fname="ttbar";
    fname="ttbar";

    TString file_name="/eos/uscms/store/user/btcarlso/efficiencies/efficiency_file_"+fname+".root";

   fBtagEff = new TFile(file_name,"READ");
   
   TH1F *h=(TH1F*)fBtagEff->FindObjectAny("TEfficiency_names");
   
   for(int i=1; i<=h->GetNbinsX();i++){
      TString name=h->GetXaxis()->GetBinLabel(i);
      if(name=="")continue;
      if(name=="TH1F_names" || name=="TH2F_names" || name=="TProfile_names" || name=="TEfficiency_names") continue;
      //cout << "loading file: " << name << endl;
      TEfficiency *pEff=(TEfficiency*)fBtagEff->FindObjectAny(name);
      //cout << "name--------->" << pEff->GetName() << endl;
      effName[pEff->GetName()]=pEff;
   }
   
}

double getSF(TLorentzVector pJi, int algFlv, TString Atagger, TString sys){
    
   // TLorentzVector pJi;
   // pJi.SetPxPyPzE(jet_px->at(ijet),jet_py->at(ijet),jet_pz->at(ijet),jet_e->at(ijet)); // 4-vector for jet ijet
    
    double x=pJi.Perp();

    int pt_bin=hName["btag_systematic_medium"]->FindBin(x);
    if(pt_bin>hName["btag_systematic_medium"]->GetNbinsX())pt_bin=hName["btag_systematic_medium"]->GetNbinsX();

    double SF=0;
    double SFsys=0;
   // int flv=TMath::Abs(jet_algFlavor->at(ijet));
    if(algFlv==4 || algFlv==5){
        if(Atagger=="CSVM") SF=(0.938887+(0.00017124*x))+(-2.76366e-07*(x*x));
        if(Atagger=="CSVT") SF=(0.927563+(1.55479e-05*x))+(-1.90666e-07*(x*x));
        if(Atagger=="CSVM" && sys=="min")SFsys=SF-hName["btag_systematic_medium"]->GetBinContent(pt_bin);
        if(Atagger=="CSVM" && sys=="max")SFsys=SF+hName["btag_systematic_medium"]->GetBinContent(pt_bin);

        if(Atagger=="CSVT" && sys=="min")SFsys=SF-hName["btag_systematic_tight"]->GetBinContent(pt_bin);
        if(Atagger=="CSVT" && sys=="max")SFsys=SF+hName["btag_systematic_tight"]->GetBinContent(pt_bin);

    }
    else{
        if(Atagger=="CSVM")SF=SF_light(Atagger,"mean",pJi);
        if(Atagger=="CSVT")SF=SF_light(Atagger,"mean",pJi);

        if(Atagger=="CSVM" && sys=="min")SFsys=SF-SF_light(Atagger,"min",pJi);
        if(Atagger=="CSVM" && sys=="max")SFsys=SF+SF_light(Atagger,"max",pJi);

        if(Atagger=="CSVT" && sys=="min")SFsys=SF-SF_light(Atagger,"min",pJi);
        if(Atagger=="CSVT" && sys=="max")SFsys=SF+SF_light(Atagger,"max",pJi);
    }

    if(sys=="") return SF;
    else if(sys=="min" || sys=="max")return SFsys;
    else return 1.0;
}

double SF_light(TString Atagger, TString mode, TLorentzVector jetP){
   //stopwatch["SF_light"].Start(kFALSE);
   double eta=TMath::Abs(jetP.Eta()); 
   double x=jetP.Perp();
   if(x>850)x=850; 
   double SF=0; 
   if(Atagger=="CSVM" && eta>=0 && eta<=0.8 ){
      double mean=((1.07541+(0.00231827*x))+(-4.74249e-06*(x*x)))+(2.70862e-09*(x*(x*x)));
      if(mode=="mean")SF=((1.07541+(0.00231827*x))+(-4.74249e-06*(x*x)))+(2.70862e-09*(x*(x*x))); 
      if(mode=="min")SF=mean-((0.964527+(0.00149055*x))+(-2.78338e-06*(x*x)))+(1.51771e-09*(x*(x*x)));
      if(mode=="max")SF=((1.18638+(0.00314148*x))+(-6.68993e-06*(x*x)))+(3.89288e-09*(x*(x*x)))-mean;
   }
   
   if(Atagger=="CSVM" && eta>0.8 && eta<=1.6 ){
      double mean=((1.05613+(0.00114031*x))+(-2.56066e-06*(x*x)))+(1.67792e-09*(x*(x*x)));
      if(mode=="mean")SF=((1.05613+(0.00114031*x))+(-2.56066e-06*(x*x)))+(1.67792e-09*(x*(x*x)));
      if(mode=="min")SF=mean-((0.946051+(0.000759584*x))+(-1.52491e-06*(x*x)))+(9.65822e-10*(x*(x*x)));
      if(mode=="max")SF=((1.16624+(0.00151884*x))+(-3.59041e-06*(x*x)))+(2.38681e-09*(x*(x*x)))-mean;
   }
   
   
   if(Atagger=="CSVM" && eta>1.6 && eta<=2.4 ){
      double mean=((1.05625+(0.000487231*x))+(-2.22792e-06*(x*x)))+(1.70262e-09*(x*(x*x)));
      if(mode=="mean")SF=((1.05625+(0.000487231*x))+(-2.22792e-06*(x*x)))+(1.70262e-09*(x*(x*x)));
      if(mode=="min")SF=mean-((0.956736+(0.000280197*x))+(-1.42739e-06*(x*x)))+(1.0085e-09*(x*(x*x)));
      if(mode=="max")SF=((1.15575+(0.000693344*x))+(-3.02661e-06*(x*x)))+(2.39752e-09*(x*(x*x)))-mean; 
   }
   if(Atagger=="CSVT"){
      double mean=((1.00462+(0.00325971*x))+(-7.79184e-06*(x*x)))+(5.22506e-09*(x*(x*x)));
      if(mode=="mean")SF=((1.00462+(0.00325971*x))+(-7.79184e-06*(x*x)))+(5.22506e-09*(x*(x*x)));
      if(mode=="min")SF=mean-((0.845757+(0.00186422*x))+(-4.6133e-06*(x*x)))+(3.21723e-09*(x*(x*x)));
      if(mode=="max")SF=((1.16361+(0.00464695*x))+(-1.09467e-05*(x*x)))+(7.21896e-09*(x*(x*x)))-mean;
   }
   
   if((mode=="min" || mode=="max") && x>850) SF=SF*2; //double systematic for high pt
   //stopwatch["SF_light"].Stop();
   
   return SF;
}

double read_btag_efficiency(TLorentzVector pJ, int flavor,string Tagger){
        int ieta=0;
        flavor=TMath::Abs(flavor);
        if(flavor<4) flavor=3; //udgs light flavor
        //if(flavor==21) flavor=3; //udgs light flavor
        double etaBin[]={0,0.8,1.6,2.4};
        
        double eta=TMath::Abs(pJ.Eta()); 
         
        for(int i=0; i<3; i++){
                if(eta>etaBin[i] && eta<=etaBin[i+1]) ieta=i;
        }
        
        TH1F *h=(TH1F*)effName[Form("btag%s_eff_eta%d_flavor%d",Tagger.c_str(),ieta,flavor)]->GetPassedHistogram();
        int ptbin=h->FindBin(pJ.Perp());
        if(ptbin > h->GetNbinsX()) ptbin=h->GetNbinsX();
        
        double eff=effName[Form("btag%s_eff_eta%d_flavor%d",Tagger.c_str(),ieta,flavor)]->GetEfficiency(ptbin);
        return eff;

}

void load_btag_sys(){
   float ptbin_b[] = {20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600,800};
   int Nb=sizeof(ptbin_b)/sizeof(double)-1;
   
   createHistogram("btag_systematic_medium","","jet p_{T} (GeV)","",Nb,ptbin_b); 
   createHistogram("btag_systematic_tight","","jet p_{T} (GeV)","",Nb,ptbin_b); 
   
   
   double SFb_error_M[] = {
      0.0415707,
      0.0204209,
      0.0223227,
      0.0206655,
      0.0199325,
      0.0174121,
      0.0202332,
      0.0182446,
      0.0159777,
      0.0218531,
      0.0204688,
      0.0265191,
      0.0313175,
      0.0415417,
      0.0740446,
      0.0596716 };
   
   double SFb_error_T[] = {
      0.0515703,
      0.0264008,
      0.0272757,
      0.0275565,
      0.0248745,
      0.0218456,
      0.0253845,
      0.0239588,
      0.0271791,
      0.0273912,
      0.0379822,
      0.0411624,
      0.0786307,
      0.0866832,
      0.0942053,
      0.102403 };
   
   for(int i=1; i<=hName["btag_systematic_medium"]->GetNbinsX();i++){
      hName["btag_systematic_medium"]->SetBinContent(i,SFb_error_M[i-1]); 
      hName["btag_systematic_tight"]->SetBinContent(i,SFb_error_T[i-1]); 
      
   }
}

void createHistogram(const char* name, const char* title,
                     const char* xTitle, const char* yTitle,
                     Int_t       nBinsX, Double_t    xLow, Double_t xUp)
{
   TH1F* h = new TH1F(name, title, nBinsX, xLow, xUp);
   h->GetXaxis()->SetTitle(xTitle);
   h->GetYaxis()->SetTitle(yTitle);
   h->Sumw2();
   hName[name] = h;
}

void createHistogram(const char* name, const char* title,
                     const char* xTitle, const char* yTitle,
                     Int_t       nBinsX, const Float_t* xBins)
{
   TH1F* h = new TH1F(name, title, nBinsX, xBins);
   h->GetXaxis()->SetTitle(xTitle);
   h->GetYaxis()->SetTitle(yTitle);
   h->Sumw2();
   hName[name] = h;
}

void createHistogram(const char* name,   const char* title,
                     const char* xTitle, const char* yTitle,
                     Int_t nBinsX, Double_t xLow, Double_t xUp,
                     Int_t nBinsY,Double_t yLow, Double_t yUp)
{
   TH2F* h = new TH2F(name, title, nBinsX, xLow, xUp, nBinsY, yLow, yUp);
   h->GetXaxis()->SetTitle(xTitle);
   h->GetYaxis()->SetTitle(yTitle);
   h->Sumw2();
   hName2D[name] = h;
}

void createProfile(const char* name,   const char* title,
                   const char* xTitle, const char* yTitle,
                   Int_t       nBinsX, Double_t    xLow, Double_t xUp)
{
        TProfile* pr = new TProfile(name, title, nBinsX, xLow, xUp);
        pr->GetXaxis()->SetTitle(xTitle);
        pr->GetYaxis()->SetTitle(yTitle);
        prName[name] = pr;
}

void bookHisto(){
  //create histograms
    //loop over jets & st
    for(int nJ=2; nJ<=nJetmax; nJ++){
        //plot st distribution
        for(int ist=100; ist<=1200; ist+=100){
            //create mass plots
        }
    }
    
    for(int ist=100; ist<=1200; ist+=100){
        //create n-jet distributions
    }
    createHistogram("h_misQID","","","Events",10,0,10);
    createHistogram("h_pt1_pt2","","p_{T}(1) (GeV)","p_{T}(2) (GeV)",50,0,500,50,0,500);
    createHistogram("h_pt1_pt2_low","","p_{T}(1) (GeV)","p_{T}(2) (GeV)",50,0,500,50,0,500);
    createHistogram("h_pt1_pt2_high","","p_{T}(1) (GeV)","p_{T}(2) (GeV)",50,0,500,50,0,500);

    float MLQ[]={300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000};
    int Nx=sizeof(MLQ)/sizeof(float)-1;
    createHistogram("Nevents","","","",10,0,10);
    createHistogram("h_nJets","","N_{jets}","Events",nJetmax,0.5,nJetmax+0.5);
    createHistogram("h_nJets_MuE","","N_{jets}","Events",nJetmax,0.5,nJetmax+0.5);

	createHistogram("h_st_1El_met150","","S_{T} (GeV)","Events",30,0,3000); 
	createHistogram("h_st_1El_met0","","S_{T} (GeV)","Events",30,0,3000); 
	createHistogram("h_muj_all","","M(#mu,j) (GeV)","Events",100,0,1000); 
	createHistogram("h_ej","","M(e,j) (GeV)","Events",100,0,1000); 

	
    createHistogram("nJets_LQ900","","S_{T} (GeV)","Events",nJetmax,0.5,nJetmax+0.5);

    createHistogram("nJets_Iso_pass","","N_{jets}","Events",nJetmax,0.5,nJetmax+0.5);
    createHistogram("nJets_Iso_total","","N_{jets}","Events",nJetmax,0.5,nJetmax+0.5);
    
    createHistogram("nJets_Iso_pass_1btag","","N_{jets}","Events",nJetmax,0.5,nJetmax+0.5);
    createHistogram("nJets_Iso_total_1btag","","N_{jets}","Events",nJetmax,0.5,nJetmax+0.5);
    
    createHistogram("h_muonpt_Iso_pass","","p_{T}(#mu)","Events",50,0,500);
    createHistogram("h_muonpt_Iso_total","","p_{T}(#mu)","Events",50,0,500);

    createHistogram("h_met_Iso_pass","","MET","Events",50,0,500);
    createHistogram("h_met_Iso_total","","MET","Events",50,0,500);
    
    createHistogram("nJets_muSF","","N_{jets}","Events",nJetmax,0.5,nJetmax+0.5);
    createHistogram("nJets_eleSF","","N_{jets}","Events",nJetmax,0.5,nJetmax+0.5);
    createHistogram("nJets_noSF","","N_{jets}","Events",nJetmax,0.5,nJetmax+0.5);
    createHistogram("nJets_EleMuSF","","N_{jets}","Events",nJetmax,0.5,nJetmax+0.5);

    createHistogram("h_mumu_muSF","","M_{#mu#mu} (GeV)","Events",50,50,130);
    createHistogram("h_mumu_noSF","","M_{#mu#mu} (GeV)","Events",50,50,130);
    createHistogram("h_mumu_eleSF","","M_{#mu#mu} (GeV)","Events",50,50,130);

    
    createHistogram("h_met","","MET (GeV)","Events",100, 0, 500);
    for(int nJ=2; nJ<=nJetmax; nJ++) createHistogram(Form("h_met_nJets%d",nJ),"","MET (GeV)","Events",100, 0, 500);
	createHistogram("h_metmax","","MET (GeV)","Events",100, 0, 500);
	createHistogram("h_dimuomassmin","","M(ll) (GeV)","Events",100, 0, 500);

	
    createHistogram("h_st","","S_{T} (GeV)","Events",50,0,2000);
    
    float stbins[]={0,100,200,300,400,500,600,800,1000,1200,1400,1600,2000};
    float stbins2[]={0,100,200,300,400,500,600,800,1000};

    int N=sizeof(stbins)/sizeof(float)-1;
    int N2=sizeof(stbins2)/sizeof(float)-1;

    for(int nJ=2; nJ<=nJetmax; nJ++){
        createHistogram(Form("h_nJets%d_stmin_uW",nJ),"","S_{T} (GeV)","Events",120,0,3000);
        createHistogram(Form("h_nJets%d_stmin",nJ),"","S_{T} (GeV)","Events",120,0,3000);
        createHistogram(Form("h_nJets%d_stmin_nominal",nJ),"","S_{T} (GeV)","Events",120,0,3000);

        createHistogram(Form("h_nJets%d_stmin_lowpt",nJ),"","S_{T} (GeV)","Events",120,0,3000);
        createHistogram(Form("h_nJets%d_stmin_highpt",nJ),"","S_{T} (GeV)","Events",120,0,3000);

        createHistogram(Form("h_nJets%d_stminSS",nJ),"","S_{T} (GeV)","Events",120,0,3000);
        
        createHistogram(Form("h_nJets%d_st_10GeV",nJ),"","S_{T} (GeV)","Events",300,0,3000);
        createHistogram(Form("h_nJets%d_st_uW",nJ),"","S_{T} (GeV)","Events",50,0,3000);
    }
    
    for(int nJ=2; nJ<=nJetmax; nJ++){
        createHistogram(Form("h_nJets%d_st",nJ),"","S_{T} (GeV)","Events",N,stbins);
        createHistogram(Form("h_nJets%d_MuE_st",nJ),"","S_{T} (GeV)","Events",N2,stbins2);
    }

    for(int nB=0; nB<=2; nB++){
        createHistogram(Form("h_nJets_%dbtags",nB),"","N_{jets}","Events",nJetmax,0.5,nJetmax+0.5);
        createHistogram(Form("h_nJets_%dbtags_lowM",nB),"","N_{jets}","Events",nJetmax,0.5,nJetmax+0.5);
        createHistogram(Form("h_nJets_%dbtags_lowM_noST",nB),"","N_{jets}","Events",nJetmax,0.5,nJetmax+0.5);

        createHistogram(Form("h_nJets_%dbtags_MuE",nB),"","N_{jets}","Events",nJetmax,0.5,nJetmax+0.5);
    }
    
    
    
    createHistogram("h_mumu_0btag","0 b-tag","M_{#mu#mu} (GeV)","Events",250,0,1000);
    createHistogram("h_mumu_1btag","1 b-tag","M_{#mu#mu} (GeV)","Events",250,0,1000);
    

    createHistogram("h_mumu","","M_{#mu#mu} (GeV)","Events",50,0,1000);

    createHistogram("h_muj","","M_{#muj} (GeV)","Events",50,0,1000);

    createHistogram("h_mumujj_LQ","","M_{LQ} (GeV)","Events",Nx,MLQ);
    
    createHistogram("h_optimization_nJetsg5","","","",1000,0,1000);
    createHistogram("h_optimization_nJetsg6","","","",1000,0,1000);
    createHistogram("h_optimization_nJetsg7","","","",1000,0,1000);

    
    createHistogram("h_optimization_nJetsg5_g1b","","","",1000,0,1000);
    createHistogram("h_optimization_nJetsg5_g2b","","","",1000,0,1000);
    
    createHistogram("h_optimization_nJetsg6_g1b","","","",1000,0,1000);
    createHistogram("h_optimization_nJetsg6_g2b","","","",1000,0,1000);
    
    
    createHistogram("h_optimization_nJetsg7_g1b","","","",1000,0,1000);
    createHistogram("h_optimization_nJetsg7_g2b","","","",1000,0,1000);
    
    createHistogram("h_optimization_nJets5_g1b","","","",1000,0,1000);
    createHistogram("h_optimization_nJets6_g1b","","","",1000,0,1000);

    
    fill_optimization_hist(StMax+10, massMax+10, m_mujMax+10, 0);
    
}
void fill_optimization_hist(float st_, float mass_, float muj_, float w_ ){
   // cout << "st: " << st_ << " mass " << mass_ <<  " muj_ " << muj_ << endl;
    if(st_>=StMax)st_=StMax;
    if(mass_>=massMax)mass_=massMax;
    if(Z_enriched) mass_=0;
    if(muj_>=m_mujMax)muj_=m_mujMax;
    
    for(float ist=300; ist<=st_; ist+=100){
        for(float imass=55; imass<=mass_; imass+=25){
            for(float muj=0; muj<=muj_; muj+=100){
                TString name=Form("St>%.0f,M>%.0f,Mmuj>%.0f",ist,imass,muj);
                //cout << name << endl;
                //TString name=Form("St>%.0f,M>%.0f",ist,imass);
                if(w_<=0){
                    hName["h_optimization_nJetsg5"]->Fill(name,w_);
                    hName["h_optimization_nJetsg5_g1b"]->Fill(name,w_);
                    hName["h_optimization_nJetsg5_g2b"]->Fill(name,w_);
                    
                    hName["h_optimization_nJetsg6_g1b"]->Fill(name,w_);
                    hName["h_optimization_nJetsg6_g2b"]->Fill(name,w_);
                    
                    hName["h_optimization_nJetsg7_g1b"]->Fill(name,w_);
                    hName["h_optimization_nJetsg7_g2b"]->Fill(name,w_);
                    
                    hName["h_optimization_nJetsg6"]->Fill(name,w_);
                    hName["h_optimization_nJetsg7"]->Fill(name,w_);
                    
                    hName["h_optimization_nJets5_g1b"]->Fill(name,w_);
                    hName["h_optimization_nJets6_g1b"]->Fill(name,w_);


                }
                else{
                 
                    if(nBtags>=1){
                        hName["h_optimization_nJetsg5_g1b"]->Fill(name,w_); //nj>=5
                        if(nJets==5)hName["h_optimization_nJets5_g1b"]->Fill(name,w_);
                        if(nJets==6)hName["h_optimization_nJets6_g1b"]->Fill(name,w_);
                        
                        if(nJets>=6)hName["h_optimization_nJetsg6_g1b"]->Fill(name,w_);
                        if(nJets>=7)hName["h_optimization_nJetsg7_g1b"]->Fill(name,w_);
                    }
                    if(nBtags>=2){
                        hName["h_optimization_nJetsg5_g2b"]->Fill(name,w_);
                        if(nJets>=6)hName["h_optimization_nJetsg6_g2b"]->Fill(name,w_);
                        if(nJets>=7)hName["h_optimization_nJetsg7_g2b"]->Fill(name,w_);
                    }
                    hName["h_optimization_nJetsg5"]->Fill(name,w_);//no b-tag requirement
                    if(nJets>=6)hName["h_optimization_nJetsg6"]->Fill(name,w_);
                    if(nJets>=7)hName["h_optimization_nJetsg7"]->Fill(name,w_);

                }
            }
        }
    }
}
void writeHisto(){
   output_file->cd();
   for (std::map<TString,TH1F*>::iterator it=hName.begin(); it!=hName.end(); it++) {
      hName[it->first]->Write();
   }
   for (std::map<TString,TH2F*>::iterator it=hName2D.begin(); it!=hName2D.end(); it++){
      hName2D[it->first]->Write();
   }
   
   output_file->Close();
   delete output_file; 
};


