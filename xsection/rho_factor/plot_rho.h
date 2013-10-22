
#include "TObject.h"
#include "TEnv.h"
#include <TROOT.h>
#include <TSystem.h>
#include "TClass.h"
#include "TStopwatch.h"

//formatting 
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TGraphAsymmErrors.h"
#include "TGraph2D.h"
#include "TPaveText.h"
#include "TMultiGraph.h"
#include <TFitResult.h>
#include <THStack.h>
#include <TLegend.h>
#include <TAttMarker.h>
#include <TProfile.h>
#include <TGaxis.h> 
#include <TAxis.h>
#include <TEfficiency.h>
#include <TMarker.h>
#include <TPave.h>
#include <TPaveStats.h>

//random stuff
#include <TText.h>
#include <TError.h>

//io stuff
#include <TTree.h>
#include <TDirectory.h>
#include "TString.h"
#include "TFile.h"

//Math stuff
#include <TMath.h>
#include <TF1.h>
#include <TF2.h>
#include <TF3.h>
#include <TRandom3.h>
#include "TVector3.h"
#include "TLorentzVector.h"
#include <TMinuit.h>
#include "TVirtualFitter.h"


//TFile *rho_MC = new TFile("output_canvas.root","READ");
//TH1F *rho_ptMC =(TH1F*)rho_MC->FindObjectAny("rho_pt_TH1"); 

double MC_shape(double *x, double *par);
void CreateCanvas(TString Name,TString Title, int x, int y );
void eff(TH1F *passed, TH1F *total, float m1,float m2, float &Np, float &NpE, float &Nt, float &NtE, bool constant, int can_);
void plot_rho();

TString binning;
std::map<TString, TCanvas*> CName;           // Map for histograms
