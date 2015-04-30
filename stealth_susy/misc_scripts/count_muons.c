/*
 *  count_muons.c
 *  
 *
 *  Created by Benjamin Carlson on 9/24/13.
 *  Copyright 2013 Carnegie Mellon University. All rights reserved.
 *
 */

#include "count_muons.h"

void count_muons(){
	initialize_xs();
	open_file("ttSemiLept","ttJetsSemiLept.root"); 
	open_file("ttFullLept","ttJetsFullLept.root"); 
	
	open_file("wJetsHT150","wJets_HT_150To200.root"); 
	open_file("wJetsHT200","wJets_HT_200To250.root"); 
	open_file("wJetsHT250","wJets_HT_250To300.root"); 
	open_file("wJetsHT300","wJets_HT_300To400.root"); 
	open_file("wJetsHT400","wJets_HT_400Toinf.root"); 
	
	open_file("w1Jets","w1JetsToLNu.root"); 
	open_file("w2Jets","w2JetsToLNu.root"); 
	open_file("w3Jets","w3JetsToLNu.root"); 
	open_file("w4Jets","w4JetsToLNu.root"); 
	
	open_file("dy1Jets","dy1JetsToLL.root"); 
	open_file("dy2Jets","dy2JetsToLL.root"); 
	open_file("dy3Jets","dy3JetsToLL.root"); 
	open_file("dy4Jets","dy4JetsToLL.root");
	
	open_file("singleMu","singleMu.root"); 

	for (std::map<TString,TFile*>::iterator it=FName.begin(); it!=FName.end(); it++) {
		cout << "File: " << it->first << endl; 
		FILENAME=it->first;
		analyze_file(it->first); 
	}
	print();
	

}

void initialize_xs(){
	xs["w1Jets"]=6662;
	xs["w2Jets"]=2159;
	xs["w3Jets"]=640;
	xs["w4Jets"]=264;
	
	Ngen["w1Jets"]=23141598;
	Ngen["w2Jets"]=34044921;
	Ngen["w3Jets"]=15539503;
	Ngen["w4Jets"]=13382803;
	
	xs["wJetsHT150"]=290.7;
	xs["wJetsHT200"]=111.37;
	xs["wJetsHT250"]=59.24;
	xs["wJetsHT300"]=47.25;
	xs["wJetsHT400"]=31.12;
	
	Ngen["wJetsHT150"]=21686209;
	Ngen["wJetsHT200"]=10039771;
	Ngen["wJetsHT250"]=6575572;
	Ngen["wJetsHT300"]=6840509;
	Ngen["wJetsHT400"]=6619654;
	
	xs["dy1Jets"]=666;
	xs["dy2Jets"]=215;
	xs["dy3Jets"]=61;
	xs["dy4Jets"]=27;
	
	Ngen["dy1Jets"]=24045248;
	Ngen["dy2Jets"]=21852156;
	Ngen["dy3Jets"]=11015445;
	Ngen["dy4Jets"]=6402827;
	
	xs["ttSemiLept"]=107.6;
	xs["ttFullLept"]=25.6;

	Ngen["ttSemiLept"]=25364818;
	Ngen["ttFullLept"]=12119013;
	
	xs["singleMu"]=1.0;
	Ngen["singleMu"]=1.0; 
	
}

void print(){
	
	int N=4;
	ofstream output("table1.tex"); 
	output <<"\\begin{tabular}{"; 
	for(int i=0; i<N; i++) output << "c"; 
	output << "}\\hline" << endl; 
	output << "Actual Events & Expected Events \\\\ \\hline"<< endl; 
	output << "sample & weight & 0 $\\mu$ & 1$\\mu$ & 2$\\mu$ & $\\ge$2$\mu$ "; 
	output <<" 0 $\\mu$ & 1$\\mu$ & 2$\\mu$ & $\\ge$2$\mu$ \\\\" << endl;
	for (std::map<TString,TFile*>::iterator it=FName.begin(); it!=FName.end(); it++) {
		output << it->first << " & ";
		double w=xs[it->first]*19600/Ngen[it->first];
		output << w << " & ";
		output<< Mu0[it->first] << " & ";
		output << Mu1[it->first] << " & ";
		output << Mu2[it->first] << " & ";
		output << Mug2[it->first] << " & "; 
		
		output<< w*Mu0[it->first] << " & ";
		output << w*Mu1[it->first] << " & ";
		output << w*Mu2[it->first] << " & ";
		output << w*Mug2[it->first] << " & ";
		
		output << "\\\\" << endl;
	}
	output << "\\hline" << endl << "\\end{tabular}" << endl; 
	output.close(); 
	
}

void analyze_file(TString file){
	TTree *tree = (TTree*)FName[file]->Get("tree");
	cout << "Intitialize Tree: " << endl; 
	initialize_tree(tree);
	event_loop(tree); 
	
}
void event_loop(TTree *tree){
	TStopwatch t; 
	t.Start(kFALSE); 
	
	cout << "Entries: " << entries << endl; 
	for (int iEvt=0; iEvt<entries; iEvt++){
		tree->GetEntry(iEvt);
		if(nJets<4 ) continue; 
		if(st>1000 && nMuons==0) Mu0[FILENAME]+=1.0; 
		if(st>1000 && nMuons==1) Mu1[FILENAME]+=1.0; 
		if(st>1000 && nMuons==2) Mu2[FILENAME]+=1.0; 
		if(st>1000 && nMuons>2) Mug2[FILENAME]+=1.0; 


	}
}

void open_file(TString name, TString full_name){
	TString dir="/eos/uscms/store/user/btcarlso/MC_trees/trigger/new/";
	if(name=="singleMu") dir="/eos/uscms/store/user/btcarlso/data_trees/";
	full_name=dir+full_name; 
	TFile *F = new TFile(full_name,"READ"); 
	FName[name]=F; 
}

void initialize_tree(TTree *tree){
	entries=tree->GetEntries();
	
	tree->SetBranchAddress("jets_n",&nJets); 
	tree->SetBranchAddress("st",&st); 
	tree->SetBranchAddress("muons_n",&nMuons); 
}