#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <math.h>

#include "TROOT.h"
#include "TRint.h"
#include "TChain.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TString.h"
#include "TString.h"
#include "TRandom.h"
#include "TUnixSystem.h"

#include "xsReader.hh"
#include "treeReaderXS.hh"

using namespace std;


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %% Usage: ./runXSReaders -f test.root
// %%        ./runXSReaders -c chains/bg-test -D root
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int main(int argc, char *argv[]) {

  int processID = gSystem->GetPid();
  cout << "Running under process ID  " << processID << endl;

  string progName  = argv[0]; 
  string writeName, fileName;
  int file(0);
  int dirspec(0);
  int nevents(-1), start(-1);
  bool skim_analysis; 
  int randomSeed(processID);

  // -- Some defaults
  string dirBase("./");               // this could point to "/home/ursl/data/root/."
  string dirName("."); dirspec = 0;   // and this to, e.g. "bmm", "bee", "bem", ...
  string cutFile("tree.defaults.cuts");

  string treeName("T1");
  string evtClassName("TAna01Event");

  string readerName("xsReader"); 
  TString histfile("");
  TString skim_folder("");

  // -- command line arguments
  for (int i = 0; i < argc; i++){
    if (!strcmp(argv[i],"-h")) {
	std::cout << "List of arguments:" << std::endl;
	std::cout << "-c filename   chain definition file" << std::endl;
	std::cout << "-C filename   file with cuts" << std::endl;
	std::cout << "-D path       where to put the output" << std::endl;
	std::cout << "-f filename   single file instead of chain" << std::endl;
	std::cout << "-n integer    number of events to run on" << std::endl;
	std::cout << "-r class name which tree reader class to run" << std::endl;
	std::cout << "-s number     seed for random number generator" << std::endl;
	std::cout << "-o filename   set output file" << std::endl;
	std::cout << "-h            prints this message and exits" << std::endl;
	std::cout << "-a skim_analysis  analysis or skim" << std::endl;
	std::cout << "-l skimmed folder" << std::endl; 
	return 0;
    }
    if (!strcmp(argv[i],"-c"))  {fileName   = string(argv[++i]); file = 0; }     // file with chain definition
    if (!strcmp(argv[i],"-C"))  {cutFile    = string(argv[++i]);           }     // file with cuts
    if (!strcmp(argv[i],"-D"))  {dirName    = string(argv[++i]);  dirspec = 1; } // where to put the output
    if (!strcmp(argv[i],"-f"))  {fileName   = string(argv[++i]); file = 1; }     // single file instead of chain
    if (!strcmp(argv[i],"-n"))  {nevents    = atoi(argv[++i]); }                 // number of events to run 
    if (!strcmp(argv[i],"-r"))  {readerName = string(argv[++i]); }               // which tree reader class to run
    if (!strcmp(argv[i],"-s"))  {randomSeed = atoi(argv[++i]); }                 // set seed for random gen.
    if (!strcmp(argv[i],"-o"))  {histfile   = TString(argv[++i]); }              // set output file
	if (!strcmp(argv[i],"-a")) {skim_analysis = atoi(argv[++i]);}                    // run in analysis or skim mode. If skim_anal = 1 run in analysis mode. 
	  //If 0 run in skim mode
	if(!strcmp(argv[i],"-l"))   {skim_folder = TString(argv[++i]);}
  }

  // -- Prepare histfilename variation with (part of) cut file name
  TString fn(cutFile);
  fn.ReplaceAll("cuts/", "");
  fn.ReplaceAll(".cuts", "");
  fn.ReplaceAll("tree", "");
  
  // -- Determine filename for output histograms and 'final' small/reduced tree
  TString meta = fileName;
  if(histfile == "") {
    TString  barefile(fileName), chainFile, meta;
    if (file == 0) {
      // -- input from chain
      if (barefile.Contains("chains/")) {
	barefile.ReplaceAll("chains/", "");
	histfile = barefile + "." + fn + ".root";
	if (dirspec) {
	  if (dirName[0] == '/') {
	    histfile = dirName + "/" + histfile;
	  } else {
	    histfile = dirBase + "/" + dirName + "/" + histfile;
	  }
	}
      } else {
	histfile =  barefile + "." + fn + ".root";
	if (dirspec) {
	  if (dirName[0] == '/') {
	    histfile = dirName + "/" + histfile;
	  } else {
	    histfile = dirBase + "/" + dirName + "/" + histfile;
	  }
	}
      }
    // -- The following lines strip everything from the string up to and including the last '/'
    int fl = barefile.Last('/');
    TString bla(barefile);
    bla.Replace(0, fl+1, ' '); bla.Strip(TString::kLeading, ' ');  bla.Remove(0,1);
    histfile =  bla + "." + fn + ".root";
    if (dirspec) {
      histfile = dirBase + "/" + dirName + "/" + histfile;
    }
  }  else if (file == 1) {
    // -- single file input
    // -- The following lines strip everything from the string up to and including the last '/'
    int fl = barefile.Last('/');
    TString bla(barefile);
    bla.Replace(0, fl+1, ' '); bla.Strip(TString::kLeading, ' ');  bla.Remove(0,1);
    histfile =  bla;
    histfile.ReplaceAll(".root", "");
    histfile +=  "." + fn + ".root";
    if (dirspec) {
	if (dirName[0] == '/') {
	  histfile = dirName + "/" + histfile;
	} else {
	  histfile = dirBase + "/" + dirName + "/" + histfile;
	}
      }
    }
  }

  cout << "Opening " << histfile.Data() << " for output histograms" << endl;
  cout << "Opening " << fileName.c_str() << " for input" << endl;


  // -- Set up chain
  cout << "treeName: " << treeName << endl; 
  TChain *chain = new TChain(TString(treeName));
  cout << "Chaining ... " << treeName << endl;
  char pName[2000]; 
  int nentries; 
  if (file == 0) {
    // -- non-trivial chain input
    ifstream is(meta);  
    while(meta.ReadLine(is) && (!meta.IsNull())){ 
      nentries = -1;
      if (meta.Data()[0] == '#') continue; 
      sscanf(meta.Data(), "%s %d", pName, &nentries); 
      if (nentries > -1) {
        cout << pName << " -> " << nentries << " entries" << endl; 
        chain->Add(pName, nentries); 
      } else {
        cout << meta << endl; 
        chain->Add(meta.Data()); 
      }
    }
    is.close();
  }
  else if (file == 1) {
    // -- single file input
    cout << fileName << endl;
    chain->Add(TString(fileName));
  }

  // -- Now instantiate the tree-analysis class object, initialize, and run it ...
  treeReaderXS *a = NULL;
  if (readerName == "xsReader") a = new xsReader(chain, TString(evtClassName));
  else {
    cout << "please provide a class name to instantiate" << endl;
  }
  
	if(skim_analysis==1)
		cout << endl << endl << endl << "IN ANALYSIS MODE!" << endl << endl << endl;
	if(skim_analysis==0)
		cout << endl << endl << endl << "IN SKIMMING MODE!" << endl << endl << endl;
	
	TString skimmed_region1=skim_folder + "/r1.root";
	TString skimmed_region2=skim_folder + "/r2.root";
	TString skimmed_region3=skim_folder + "/r3.root";
	TString discarded = skim_folder + "/discarded.root";	
	
	
  if (a){
  if(skim_analysis==1)
	  a->openHistFile(histfile); 
  else{
      a->openHistFile(histfile);
      a->openCopyFile1(skimmed_region1);
      a->openCopyFile2(skimmed_region2);
      a->openCopyFile3(skimmed_region3);
      a->openDiscardedFile(discarded);
	  
  }

  a->readCuts(cutFile, 1);
  a->bookHist();
  if(skim_analysis==1)
    //  a->bookHist();
  a->startAnalysis(); 
  a->loop(nevents, start,skim_analysis);
  a->print_numbers(nevents,skim_analysis);
  if(skim_analysis==1){
	a->normalize();
	a->closeHistFile(); 
  }
  else {
	  a->closeCopyFile1();
	  a->closeCopyFile2();
	  a->closeCopyFile3();
	  a->closeDiscardedFile();
  } 	  
  } else 
    cerr << "Readerclass '" << readerName << "' not found" << endl;
  
  delete a; // so we can dump some information in the destructor
  
  return 0; 
  
}
