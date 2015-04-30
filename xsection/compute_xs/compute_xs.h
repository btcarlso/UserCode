/*
 *  create_plots.h
 *  
 *
 *  Created by Benjamin Carlson on 8/9/13.
 *  Copyright 2013 Carnegie Mellon University. All rights reserved.
 *
 */

#include "/uscms/home/btcarlso/stealth_code/root_headers.h"
#include "/uscms/home/btcarlso/code/bins_final.card"
#include "/uscms/home/btcarlso/code/names.C"
#include "/uscms/home/btcarlso/code/smearing_parameters.c"

std::map<TString, TH1F*> hName;
std::map<TString, TProfile*> prName;
std::map<TString, TH2F*> hName2D;
std::map<TString, TCanvas*> CName; 
std::map<TString, TGraphAsymmErrors*> grName; 
std::map<TString, TF1*> f1Name; //1D TF1 name

TString x_label="p_{T} [GeV]";
//TString xs_y="#scale[0.7]{d#sigma/dp_{T} #times Br(#mu^{+}#mu^{#font[122]{\55}}) [pb/GeV]}";
TString xs_y="d#sigma/dp_{T} #times Br(#mu^{+}#mu^{#font[122]{\55}}) [pb/GeV]";
TString Ratio_y3="Corrected Ratio R_{31}";
//"#sigma(3S)#timesBr(#mu^{+}#mu^{-})/#sigma(1S)#timesBr(#mu^{+}#mu^{-})";
TString Ratio_y2="Corrected Ratio R_{21}";
//"#sigma(2S)#timesBr(#mu^{+}#mu^{-})/#sigma(1S)#timesBr(#mu^{+}#mu^{-})";

TString cms_pre = "CMS Preliminary";
TString cms_sim = "CMS Simulation";
TString lumi = "#sqrt{s} = 8 TeV, L = 19.7 fb^{-1}";
TLatex L1;
TLatex L2;
TLatex L1sim;

//stuff for lumi text

TString cmsText     = "CMS";
float cmsTextFont   = 61;  // default is helvetic-bold

bool writeExtraText = false;
TString extraText   = "Preliminary";
float extraTextFont = 52;  // default is helvetica-italics

// text sizes and text offsets with respect to the top frame
// in unit of the top margin size
float lumiTextSize     = 0.6;
float lumiTextOffset   = 0.2;
float cmsTextSize      = 0.75;
float cmsTextOffset    = 0.1;  // only used in outOfFrame version

float relPosX    = 0.045;
float relPosY    = 0.035;
float relExtraDY = 1.2;

// ratio of "CMS" and extra text size
float extraOverCmsTextSize  = 0.76;

TString lumi_13TeV = "20.1 fb^{-1}";
TString lumi_8TeV  = "19.7 fb^{-1}";
TString lumi_7TeV  = "4.9 fb^{-1}";

bool drawLogo      = false;


void compute_xs();
void create_plots();
void read_NLO_xs(int ups, int iy);
void readHisto();
void draw_muon_pt();
void writeHisto();
void draw_zeta();
void draw_header();
void CMS_lumi( TPad* pad, int iPeriod, int iPosX ); 
void draw_headersim();
void draw_dR();
void draw_zeta_m();
void draw_zeta_profile();
void draw_pt();
void yield_table();
void draw_efficiency();
void efficiency_closure();
void draw_m();
void draw_zeta_eta();
void atlas_cmsratio();
void efficiency_table(int iy);

void acceptance_ratio(int iy, int ups, string mode);
void acceptance_ratio_plot();
void load_acceptance();
void load_sg(int iy);
void Rebin(TString name);
void rebin_acceptance();
void weighted_acceptance(int iy, int ups);
void dsigma_pt(int iy, int ups);
void uncertainty(int iy, int ups); 
void uncertaintyR(int iy, int num);
void RN1(int iy, int num);
void draw_unc(int iy);
void draw_unc_ratio(int iy);

int set_precision(double x);
int set_precisionE(double x);
void xs_table(int iy);
void xs_tablehepdata(int iy);
void ratio_table(int iy); 
void ratio_tablehepdata(int iy);

void get_atlas_ratios();
void draw_xs();
void rap_dep2010();
void rapidity_dep_atlas();
void table_rapidity_ratio();
void draw_rapidity_ratio();
void draw_xs_fits(int iy,int ups);
void draw_xs_fits_wNLO(int iy,int ups); 
void compare_ratios(int num);

void plot_atlas_pt(int ups);
void draw_xs_atlas();

void set_CMS2010(TGraphAsymmErrors *gr);
void set_CMS2011(TGraphAsymmErrors *gr,string MODE, int ups);

void bin_center(int iy, int ups);

void plot_orderDist();
void fit_nuisance(int iy);
void xs_pull(int iy, int ups);
void draw_ratio(int num);
TString divide_graphs(TString num, TString den);
void plot_cms2010(int ups);
void scale_graph(TGraphAsymmErrors *gr, double SF);

void xs_fit(int iy, int ups);
void xs_fit_po_study(int iy, int ups);
void acceptance_table(int iy, int ups);
void acceptance_table_smearing(int iy, int ups);
void acceptance_table_summary(int iy, int ups);
void acceptance_table_summaryhepdata(int iy, int ups);
double power_law(double *x, double*par);
double power_law_variablepo(double *x, double*par);
double ratio_function(double *x, double *par);
void fit_table(int iy);

void acceptance_ratio(int iy);
void draw_acc_ratio(int iy);
void combine_rapidity_bins(int ups);


void CreateCanvas(TString Name,TString Title, int x, int y );
void CreateHistogram(const char* name,   const char* title,
					 const char* xTitle, const char* yTitle,
					 Int_t       nBinsX, Double_t    xLow, Double_t xUp);
void CreateHistogram(const char* name,   const char* title,
					 const char* xTitle, const char* yTitle,
					 Int_t       nBinsX, const Float_t* xBins);
int Cx=800,Cy=600; 
TFile *hist_file; 
TFile *yields_file; 
TFile *canvas_file;

TFile *eff1S; 
TFile *eff2S; 
TFile *eff3S; 




//TString cms_sim = "CMS Simulation"; 
TString lumi_header = "#sqrt{s} = 7 TeV, L = 4.9 fb^{-1}"; 
TLatex L1_header; 
TLatex L2_header; 
TLatex L1sim_header; 