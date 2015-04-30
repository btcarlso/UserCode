#include "analyze_trees.h"

void analyze_trees(int jobnumber){
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
    
	if(jobnumber==-1){
        //test job, not standard
		sample_list.push_back("ttFullLept");
		test=true;
		isMC=true;
		output_file_name="test";
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
		for(int nJ=1; nJ<=4; nJ++) {
			TString tag=Form("dy%dJets",nJ);
			sample_list.push_back(tag);
		}
		isMC=true;
		output_file_name="dy";
		selections="_trigger";
        
	}
    JN++;
    
	//dy systematics
	if(jobnumber==JN) {
		sample_list.push_back("dyJets_scaleup");
		
		isMC=true;
		output_file_name="dyJets_scaleup";
		selections="_trigger";
		
	}
	JN++;
    
	if(jobnumber==JN) {
		sample_list.push_back("dyJets_scaledown");
		
		isMC=true;
		output_file_name="dyJets_scaledown";
		selections="_trigger";
		
	}
	JN++;
    
	if(jobnumber==JN) {
		sample_list.push_back("dyJets_matchingup");
		
		isMC=true;
		output_file_name="dyJets_matchingup";
		selections="_trigger";
		
	}
    JN++;
	
	if(jobnumber==JN) {
		sample_list.push_back("dyJets_matchingdown");
		
		isMC=true;
		output_file_name="dyJets_matchingdown";
		selections="_trigger";
		
	}
	JN++;
    
	if(jobnumber==JN) {
		sample_list.push_back("dyJets_scaleup");
		
		isMC=true;
		output_file_name="dyJets_scaleup";
		selections="_trigger";
		
	}
    JN++;
    
    //7
	//ttbar semi+fully leptonic
	if(jobnumber==JN) {
		sample_list.push_back("ttSemiLept");
		sample_list.push_back("ttFullLept");
		//sample_list.push_back("ttG");
		
		isMC=true;
		output_file_name="ttbar";
		selections="_trigger";
        
	}
    JN++;
	//8
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
    
	if(jobnumber==JN) {
		sample_list.push_back("ttJets_scaleup");
		
		isMC=true;
		test=true;
		output_file_name="ttJets_scaleup";
		selections="_trigger";
		
	}
    JN++;
	
	if(jobnumber==JN) {
		sample_list.push_back("ttJets_scaledown");
		
		isMC=true;
		output_file_name="ttJets_scaledown";
		selections="_trigger";
	}
	JN++;
    
	if(jobnumber==JN) {
		sample_list.push_back("ttJets_matchingup");
		
		isMC=true;
		output_file_name="ttJets_matchingup";
		selections="_trigger";
	}
	JN++;
    
	if(jobnumber==JN) {
		sample_list.push_back("ttJets_matchingdown");
		
		isMC=true;
		output_file_name="ttJets_matchingdown";
		selections="_trigger";
	}
	JN++;
    
    if(jobnumber==JN) {
		sample_list.push_back("ttJetsisMCaNLO");
		
		isMC=true;
		output_file_name="ttaMCNLO";
		selections="_trigger";
	}
	JN++;
    
    if(jobnumber==JN) {
		sample_list.push_back("ttJetsHadronic");
		
		isMC=true;
		output_file_name="ttJetsHadronic";
		selections="_trigger";
	}
	JN++;
	
	if(jobnumber==JN) {
		sample_list.push_back("QCD_30-50");
		sample_list.push_back("QCD_50-80");
		sample_list.push_back("QCD_80-120");
		sample_list.push_back("QCD_120-170");
		sample_list.push_back("QCD_170-300");
		sample_list.push_back("QCD_300-470");
		sample_list.push_back("QCD_470-600");
		sample_list.push_back("QCD_600-800");
		sample_list.push_back("QCD_800-1000");
        
		isMC=true;
		output_file_name="QCD";
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
	//18
	if(jobnumber==JN){
		sample_list.push_back("WWJetsTo2L2Nu");
		sample_list.push_back("WZJetsTo2L2Q");
		sample_list.push_back("WZJetsTo3LNu");
		sample_list.push_back("WZJetsTo2QLNu");
		
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
        
		
		isMC=true;
		QCD=false;
		output_file_name = "ttZ";
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
    
	//stealth 21
	if(jobnumber==JN) {
		sample_list.push_back("stealth_300_200");
		
		isMC=true;
        _fastSim=true;
        filter_eff=0.543;
		output_file_name="stealth_300_200";
		selections="_trigger";
		
	}
	JN++;
    
	if(jobnumber==JN) {
		sample_list.push_back("stealth_400_200");
		
		isMC=true;
        _fastSim=true;
        filter_eff=0.543;
		output_file_name="stealth_400_200";
		selections="_trigger";
	}
	JN++;
    
	if(jobnumber==JN) {
		sample_list.push_back("stealth_500_300");
		
		isMC=true;
        _fastSim=true;
        filter_eff=0.543;
		output_file_name="stealth_500_300";
		selections="_trigger";
		
	}
	JN++;
	
	if(jobnumber==JN) {
		sample_list.push_back("stealth_600_300");
		
		isMC=true;
        _fastSim=true;
        filter_eff=0.543;
		output_file_name="stealth_600_300";
		selections="_trigger";
		
	}
    JN++;
    
	if(jobnumber==JN) {
		sample_list.push_back("stealth_700_400");
		
		isMC=true;
        _fastSim=true;
        filter_eff=0.543;
		output_file_name="stealth_700_400";
		selections="_trigger";
		
	}
    JN++;
    
    if(jobnumber==JN) {
		sample_list.push_back("stealth_800_400");
		
		isMC=true;
        _fastSim=true;
        filter_eff=0.543;
		output_file_name="stealth_800_400";
		selections="_trigger";
		
	}
    JN++;
    
    if(jobnumber==JN) {
		sample_list.push_back("stealth_900_500");
		
		isMC=true;
        _fastSim=true;
        filter_eff=0.543;
		output_file_name="stealth_900_500";
		selections="_trigger";
		
	}
    JN++;
    
    
    //stealth 28
	if(jobnumber==JN) {
		sample_list.push_back("stealth_300_200");
		
		isMC=true;
        _fastSim=true;
        filter_eff=0.543;
		output_file_name="stealth_300_200_jecP";
		selections="_trigger_jecP";
		
	}
	JN++;
    
	if(jobnumber==JN) {
		sample_list.push_back("stealth_400_200");
		
		isMC=true;
        _fastSim=true;
        filter_eff=0.543;
		output_file_name="stealth_400_200_jecP";
		selections="_trigger_jecP";
	}
	JN++;
    
	if(jobnumber==JN) {
		sample_list.push_back("stealth_500_300");
		
		isMC=true;
        _fastSim=true;
        filter_eff=0.543;
		output_file_name="stealth_500_300_jecP";
		selections="_trigger_jecP";
		
	}
	JN++;
	
	if(jobnumber==JN) {
		sample_list.push_back("stealth_600_300");
		
		isMC=true;
        _fastSim=true;
        filter_eff=0.543;
		output_file_name="stealth_600_300_jecP";
		selections="_trigger_jecP";
		
	}
    JN++;
    
	if(jobnumber==JN) {
		sample_list.push_back("stealth_700_400");
		
		isMC=true;
        _fastSim=true;
        filter_eff=0.543;
		output_file_name="stealth_700_400_jecP";
		selections="_trigger_jecP";
		
	}
    JN++;
    
    if(jobnumber==JN) {
		sample_list.push_back("stealth_800_400");
		
		isMC=true;
        _fastSim=true;
        filter_eff=0.543;
		output_file_name="stealth_800_400_jecP";
		selections="_trigger_jecP";
		
	}
    JN++;
    
    if(jobnumber==JN) {
		sample_list.push_back("stealth_900_500");
		
		isMC=true;
        _fastSim=true;
        filter_eff=0.543;
		output_file_name="stealth_900_500_jecP";
		selections="_trigger_jecP";
		
	}
    JN++;
    //35
    
	if(jobnumber==JN) {
		sample_list.push_back("stealth_300_200");
		
		isMC=true;
        _fastSim=true;
        filter_eff=0.543;
		output_file_name="stealth_300_200_jecM";
		selections="_trigger_jecM";
		
	}
	JN++;
    
	if(jobnumber==JN) {
		sample_list.push_back("stealth_400_200");
		
		isMC=true;
        _fastSim=true;
        filter_eff=0.543;
		output_file_name="stealth_400_200_jecM";
		selections="_trigger_jecM";
	}
	JN++;
    
	if(jobnumber==JN) {
		sample_list.push_back("stealth_500_300");
		
		isMC=true;
        _fastSim=true;
        filter_eff=0.543;
		output_file_name="stealth_500_300_jecM";
		selections="_trigger_jecM";
		
	}
	JN++;
	
	if(jobnumber==JN) {
		sample_list.push_back("stealth_600_300");
		
		isMC=true;
        _fastSim=true;
        filter_eff=0.543;
		output_file_name="stealth_600_300_jecM";
		selections="_trigger_jecM";
		
	}
    JN++;
    
	if(jobnumber==JN) {
		sample_list.push_back("stealth_700_400");
		
		isMC=true;
        _fastSim=true;
        filter_eff=0.543;
		output_file_name="stealth_700_400_jecM";
		selections="_trigger_jecM";
		
	}
    JN++;
    
    if(jobnumber==JN) {
		sample_list.push_back("stealth_800_400");
		
		isMC=true;
        _fastSim=true;
        filter_eff=0.543;
		output_file_name="stealth_800_400_jecM";
		selections="_trigger_jecM";
		
	}
    JN++;
    
    if(jobnumber==JN) {
		sample_list.push_back("stealth_900_500");
		
		isMC=true;
        _fastSim=true;
        filter_eff=0.543;
		output_file_name="stealth_900_500_jecM";
		selections="_trigger_jecM";
		
	}
    JN++;
    
    
    if(jobnumber==JN) {
        TString tag=Form("dy%dJets",1);
        sample_list.push_back(tag);
        isMC=true;
		output_file_name="dy1Jets";
		selections="_trigger";
        
	}
    JN++;
    
    if(jobnumber==JN) {
        TString tag=Form("dy%dJets",2);
        sample_list.push_back(tag);
        isMC=true;
		output_file_name="dy2Jets";
		selections="_trigger";
        
	}
    JN++;
    
    if(jobnumber==JN) {
        TString tag=Form("dy%dJets",3);
        sample_list.push_back(tag);
        isMC=true;
		output_file_name="dy3Jets";
		selections="_trigger";
        
	}
    JN++;
    
    if(jobnumber==JN) {
        TString tag=Form("dy%dJets",4);
        sample_list.push_back(tag);
        isMC=true;
		output_file_name="dy4Jets";
		selections="_trigger";
        
	}
    JN++;
    
    if(jobnumber==JN) {
        TString tag="dyJets";
        sample_list.push_back(tag);
        isMC=true;
		output_file_name="dyJets_fullSim";
		selections="_noTrigger";
        
	}
    JN++;
    
    if(jobnumber==JN) {
        TString tag="dyJets";
        sample_list.push_back(tag);
        isMC=true;
		output_file_name="dyJets_fastSim";
		selections="_fastSim_noTrigger";
        
	}
    JN++;
    
	//RPV
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
		for(int nJ=1; nJ<=4; nJ++) {
			TString tag=Form("w%dJets",nJ);
			sample_list.push_back(tag);
		}
		isMC=true;
		output_file_name="wJets";
		selections="_trigger";
	}
    JN++;
	
	if(jobnumber==JN) {
		for(int nJ=1; nJ<=4; nJ++) {
			TString tag=Form("w%dJets",nJ);
			sample_list.push_back(tag);
		}
		for(int nJ=1; nJ<=4; nJ++) {
			TString tag=Form("dy%dJets",nJ);
			sample_list.push_back(tag);
		}
		sample_list.push_back("ttSemiLept");
		sample_list.push_back("ttFullLept");
		
		sample_list.push_back("QCD_30-50");
		sample_list.push_back("QCD_50-80");
		sample_list.push_back("QCD_80-120");
		sample_list.push_back("QCD_120-170");
		sample_list.push_back("QCD_170-300");
		sample_list.push_back("QCD_300-470");
		sample_list.push_back("QCD_470-600");
		sample_list.push_back("QCD_600-800");
		sample_list.push_back("QCD_800-1000");
		
		
		sample_list.push_back("TBar_t");
		sample_list.push_back("TBar_s");
		sample_list.push_back("TBar_tW");
		
		sample_list.push_back("T_t");
		sample_list.push_back("T_s");
		sample_list.push_back("T_tW");
		
		sample_list.push_back("WWJetsTo2L2Nu");
		sample_list.push_back("WZJetsTo2L2Q");
		sample_list.push_back("WZJetsTo3LNu");
		sample_list.push_back("WZJetsTo2QLNu");
		
		sample_list.push_back("ZZJetsTo2L2Q");
		sample_list.push_back("ZZJetsTo2L2Nu");
		sample_list.push_back("ZZJetsTo4L");
		
		//sample_list.push_back("ttG");
		
		isMC=true; 
		output_file_name="allMC"; 
		selections="_trigger";
		
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
        
		analyze_file(FN);
        
	}
	
}

void prepareNormalization(){
   std::cout <<"prepare for normalization .... " << std::endl;
   
   // create maps of filename, xsec and generated events
   xs["singleEle"]=1.0;
   xs["singleMu"]=1.0;

   fileName["singleMu"]="singleMu";
    
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
	
	tree->SetBranchAddress("top_px",&top_px, &b_top_px);
	tree->SetBranchAddress("top_py",&top_py, &b_top_py);
	tree->SetBranchAddress("top_pz",&top_pz, &b_top_pz);
	tree->SetBranchAddress("top_e",&top_e, &b_top_e);
     */
	
	
	
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
        
		count_muons();
		count_electrons();
        fill_jets();
        loop_leptons();
        calcSt();
        calDilepMass();
        calcMuJ();
		weight=W;

		//cout << "weight before eff: " << weight << endl;
        
		if(preselection_selection()==false) continue;
        // pileup_reweight();
        //print_event();
        fill_LQ();
        
	}
    
    
}

//-------------------------------------------------------------------------------

void print_event(){
    
    cout << "hT: " << hT << " sT: " << sT << endl;
    cout << "nElecrons: " << nElectrons << " nMuons " << nMuons << endl;
    cout << "dimuon: " << dimuon_enriched << " " << dimuon_mass << endl;
    cout << "e,mu: " << tt_enriched << " " << endl;
    cout << "mu_j: " << m_muj << endl;
    cout << "dielectron: " << dielectron_enriched << " " << dielectron_mass << endl;
    cout << "nJets: " << nJets << endl;

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
        if(ijet==0 && pJ.Perp()<jetPt1_cut) continue;
        if(ijet==1 && pJ.Perp()<jetPt2_cut) continue;
        if(ijet>=2 && pJ.Perp()<jetPt_cut) continue;
        
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
        if(ijet==0 || ijet==1)hT+=pJ.Perp();
    }
    
}

void fill_LQ(){
    //cout << "fill LQ: " << endl;
    float MLQ[]={300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000};
    float st_bins[]={380, 460, 540, 615, 685, 755, 820, 880, 935, 990, 1040, 1090, 1135, 1175, 1210};
    float mass_bins[]={100, 115, 125, 140, 150, 165, 175, 185, 195, 205, 215, 220, 230, 235, 245};
    float mj[]={115, 115, 120, 135, 155, 180, 210, 250, 295, 345, 400, 465, 535, 610, 690};
    int N=sizeof(st_bins)/sizeof(float);
    
    //cout << m_muj << " " << weight << endl;
    
    if(dimuon_enriched)hName["nJets"]->Fill(nJets,weight);
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
    for(unsigned int i=0; i<electron_px->size(); i++){
        pL.SetPxPyPzE(electron_px->at(i),electron_py->at(i),electron_pz->at(i),electron_e->at(i));
        electronCollection.push_back(pL);
        lepPt+=pL.Perp();
    }
    
    for(unsigned int i=0; i<muon_px->size(); i++){
        pL.SetPxPyPzE(muon_px->at(i),muon_py->at(i),muon_pz->at(i),muon_e->at(i));
        muonCollection.push_back(pL);
        lepPt+=pL.Perp();
    }
}

void calcSt(){
    sT=hT+lepPt;
}

bool preselection_selection(){
    if(nLooseElectrons>=1) return false; // no additional loose electrons
    if(nLooseMuons>=1) return false;
    if(nJets<2) return false;
    if(nJets>=nJetmax) nJets=nJetmax;
    if(sT<stMin) return false;
    
    //3 cases: ee,mumu,emu
    
    dimuon_enriched=false;
    dielectron_enriched=false;
    tt_enriched=false;
    Z_enriched=true;
    
    SS=false;
    
    //require exactly 2 muons
	if(nElectrons==0 && nMuons==2 && dimuon_mass>min_diLeptonMass) dimuon_enriched=true;
    //require exactly 2 electrons
    if(nMuons==0 && nElectrons==2 && dielectron_mass>min_diLeptonMass) dielectron_enriched=true;
    //require 1 mu,e
	if(nMuons==1 && nElectrons==1) tt_enriched=true;
	if(dimuon_enriched && muon_charge->at(0)*muon_charge->at(1)==1) SS=true;
    if(dielectron_enriched && muon_charge->at(0)*muon_charge->at(1)==1) SS=true;
	if(tt_enriched && muon_charge->at(0)*electron_charge->at(0)==1) SS=true;
    
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
    
    if(print)cout << endl;
}

void calDilepMass(){
    if(muon_px->size()==2){
        TLorentzVector pM1;
        TLorentzVector pM2;
        pM1.SetPxPyPzE(muon_px->at(0),muon_py->at(0),muon_pz->at(0),muon_e->at(0));
        pM2.SetPxPyPzE(muon_px->at(1),muon_py->at(1),muon_pz->at(1),muon_e->at(1));
        pM1+=pM2;
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
    float MLQ[]={300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000};
    int Nx=sizeof(MLQ)/sizeof(float)-1;
    
    createHistogram("nJets","","S_{T} (GeV)","Events",nJetmax,0.5,nJetmax+0.5);
    createHistogram("nJets_LQ900","","S_{T} (GeV)","Events",nJetmax,0.5,nJetmax+0.5);

    
    createHistogram("h_st","","S_{T} (GeV)","Events",50,0,2000);

    createHistogram("h_mumu","","M_{#mu#mu} (GeV)","Events",50,0,1000);

    createHistogram("h_muj","","M_{#muj} (GeV)","Events",50,0,1000);

    createHistogram("h_mumujj_LQ","","M_{LQ} (GeV)","Events",Nx,MLQ);
    
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


