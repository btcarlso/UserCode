/*
 *  compute_xs.h
 *  
 *
 *  Created by Benjamin Carlson on 1/23/13.
 *  Copyright 2013 Carnegie Mellon University. All rights reserved.
 *
 */
//C Headers
#include <iostream> 
#include <iomanip>
#include <fstream>
#include <math.h>
#include <vector>


//Headers for running locally 
#include "/Users/carlsonbt1/results/code/root_headers.h"
#include "/Users/carlsonbt1/results/code/bins_use.card" // also bins_new.card
#include "/Users/carlsonbt1/results/code/opt_params.C"
#include "/Users/carlsonbt1/results/code/names.C"

//open universal files 

TFile *data = new TFile("/Users/carlsonbt1/results/data_files/data_fullmass.root","READ"); //data_feb11_newbins.root
TFile *lineshape_file = new TFile("/Users/carlsonbt1/results/data_files/smearing_newbins_feb11.root","READ"); //smearing_newbins_feb11.root
TFile *acceptance_file1S = new TFile("/Users/carlsonbt1/results/data_files/accept1S.root","READ"); 
TFile *acceptance_file2S = new TFile("/Users/carlsonbt1/results/data_files/accept2S.root","READ"); 
TFile *acceptance_file3S = new TFile("/Users/carlsonbt1/results/data_files/accept3S.root","READ"); 

TFile *output_files = new TFile("/Users/carlsonbt1/results/compute_xs/plots.root","RECREATE"); 
//TFile *yields_file = new TFile("/Users/carlsonbt1/results/public_frozen/output_feb11.root","READ");

double lumi2011=4.9;
double label_text_size=0.035; 
double label_spacing=0.04;
double label_x=0.37; 
double label_y=0.85; 
TString cms_pre="CMS Preliminary";
TString lumi_string="pp  #sqrt{s} = 7 TeV, L_{int} = " + TString::Format("%.1f fb^{-1}",lumi2011); 

TLatex prelimText(label_x, label_y, cms_pre);
TLatex prelimTextSim(0.3,0.5,"CMS Preliminary Simulation"); 
TLatex lumiText(label_x, label_y-label_spacing, lumi_string);

TString lumi_exc="luminosity unc. (2.2%) excluded";


TLatex prelimTextL(0.28, 0.3, cms_pre);

TLatex lumiTextL(0.28, 0.3-label_spacing, lumi_string);



string x_label="p_{T}(#mu#mu) [GeV]";
string xs_y="#scale[0.6]{d#sigma/dp_{T}#timesBr(#mu^{+}#mu^{-}) [pb GeV^{-1}]}";



double PDF_shape(double *x, double *par);						 

TF2 *LS1_shape; //1S
TF2 *LS2_shape; //2S
TF2 *LS3_shape; //3S

TH2D *dm_m_data; 
TH2D *dm_m_data_unweighted;
TH1D *mdata; 


TH1D *S1_pt = new TH1D("S1_pt", "Fit Parameter S1; P_{T} [GeV]", fNpt2,fPTbin2); 
TH1D *S2_pt = new TH1D("S2_pt", "Fit Parameter S2; P_{T} [GeV]", fNpt2,fPTbin2); 
TH1D *S3_pt = new TH1D("S3_pt", "Fit Parameter S3; P_{T} [GeV]", fNpt2,fPTbin2); 

TH1D *chi2_pt = new TH1D("chi2_pt", "#chi^{2}/NDOF vs P_{T}; P_{T} [GeV]" , fNpt2, fPTbin2); 

TH1D *deltaM_pt = new TH1D("deltaM_pt", "; P_{T}(#mu#mu) GeV; #deltaM MeV", fNpt2,fPTbin2);
TH1D *dm_scale_pt = new TH1D("dm_scale_pt", ";P_{T}(#mu#mu) GeV;c_{w}", fNpt2,fPTbin2);

TH1D *res_pt_P1S = new TH1D("res_pt_P1S","Signal Resolution + vs. P_{T} GeV; P_{T} GeV; #sigma_{+}", fNpt2,fPTbin2); 
TH1D *res_pt_M1S = new TH1D("res_pt_M1S","Signal Resolution - vs. P_{T} GeV; P_{T} GeV; #sigma_{-}", fNpt2,fPTbin2); 
/*
TH1D *yield1S = new TH1D("yield1S","Yield #Upsilon(1S) vs. P_{T} GeV;P_{T} GeV", fNpt2, fPTbin2); 
TH1D *yield2S = new TH1D("yield2S","Yield #Upsilon(2S) vs. P_{T} GeV;P_{T} GeV", fNpt2, fPTbin2); 
TH1D *yield3S = new TH1D("yield3S","Yield #Upsilon(3S) vs. P_{T} GeV;P_{T} GeV", fNpt2, fPTbin2);
*/
TH1D *bin_width = new TH1D("bin_width","Bin Width for each P_{T} bin", fNpt2, fPTbin2); 


TH1D *Y1_bkg_4 = new TH1D("Y1_bkg_4", Form("Integral(%.2f-%.2f) for 4th order",y1m[0],y1m[1]),fNpt2,fPTbin2);
TH1D *Y1_bkg_3 = new TH1D("Y1_bkg_3", Form("Integral(%.2f-%.2f) for 3rd order",y1m[0],y1m[1]),fNpt2,fPTbin2);
TH1D *Y1_bkg_2 = new TH1D("Y1_bkg_2", Form("Integral(%.2f-%.2f) for 2nd order",y1m[0],y1m[1]),fNpt2,fPTbin2);
TH1D *Y1_bkg_1 = new TH1D("Y1_bkg_1", Form("Integral(%.2f-%.2f) for 1st order",y1m[0],y1m[1]),fNpt2,fPTbin2);
TH1D *Y1_bkg_0 = new TH1D("Y1_bkg_0", Form("Integral(%.2f-%.2f) for constant ",y1m[0],y1m[1]),fNpt2,fPTbin2);


TH1D *Y1_bkg_chi2_4 = new TH1D("Y1_bkg_chi2_4", Form("#chi^{2}/NDOF for 4th order",y1m[0],y1m[1]),fNpt2,fPTbin2);
TH1D *Y1_bkg_chi2_3 = new TH1D("Y1_bkg_chi2_3", Form("#chi^{2}/NDOF for 3rd order",y1m[0],y1m[1]),fNpt2,fPTbin2);
TH1D *Y1_bkg_chi2_2 = new TH1D("Y1_bkg_chi2_2", Form("#chi^{2}/NDOF for 2nd order",y1m[0],y1m[1]),fNpt2,fPTbin2);
TH1D *Y1_bkg_chi2_1 = new TH1D("Y1_bkg_chi2_1", Form("#chi^{2}/NDOF for 1st order",y1m[0],y1m[1]),fNpt2,fPTbin2);
TH1D *Y1_bkg_chi2_0 = new TH1D("Y1_bkg_chi2_0", Form("#chi^{2}/NDOF for constant",y1m[0],y1m[1]),fNpt2,fPTbin2);

TH1D *eff=new TH1D("eff","1/#epsilon(#mu#mu)",fNpt2,fPTbin2); 

//TH1D *r12 = new TH1D("r12","r12",fNpt2, fPTbin2); 
//TH1D *r13 = new TH1D("r13","r13",fNpt2, fPTbin2); 


void CreateCanvas(string Name,string Title, int x, int y );
void book_canvas();
TLegend *error_legend = new TLegend(0.1,0.1,0.9,0.9,"","NDC"); 

std::map<TString, TCanvas*> CName;           // Map for histograms

void initialize();
void compute_xs();

void compute_yields();
int get_mem(long* vmrss_kb, long* vmsize_kb);
double BKG_cheb(double *x, double *par);
double sig_bkg(double *x, double *par);
double sig_bkg2(double *x, double *par);
double signal_shape(double *x, double *par);
double ratio_function(double *x, double *par);

void scale_graph(TGraphAsymmErrors *gr, double SF);
void clear_points(TGraphAsymmErrors *gr, double Xmin, double Xmax);
void set_titles(TMultiGraph *GR);
void set_CMS2010(TGraphAsymmErrors *gr,string MODE);
void set_CMS2011(TGraphAsymmErrors *gr,string MODE, int ups);
void set_ATLAS(TGraphAsymmErrors *gr);

void estimate_yield(double *y);
double fit_bkg_N(TH1D *mdata,TF1 *BKG_fixed, int ipt);
void fit_mass_range(TF1 *FF, int ipt, double *Y, double *YE,double *R, double *RE, double &chi2, int &NDOF, string fit_method);
void set_bkg_order(int n,TF1 *LS_bkg);
void compute_mass_sys(int iy, int ipt, double &chi2_min, int &NDOF_min);
void fill_toy_yield(int ups,int iy);
void fit_sys_error(int iy, int ipt);
void fit(int iy, int ipt);
double compute_resid(int iy, int ipt,TH1D *h, TF1 *fit_shape,int &NDOF);
void compute_chi2(TH1D *h, TF1 *F, double a, double b, double &chi2, int &NDOF);
void make_fit_plot(int iy, int ipt, string mode);
void draw_background_errors();

void acceptance(int ups);
void acceptance_2D(int ups);
void acceptance_summary_plot();
void set_acceptance_info(TH1D *h,int ups);
void reweight_acceptance(int ups, int iy);
void yield_GeV(int iy, string mode);
void total_systematic(int iy, int ups);
void total_ratio_uncertainty(int num);
void stat_error(int iy, int ups);
void plot_total_error(int iy,int ups);
void compare_bkg_systematics(int iy, int ups);
void set_name_reco(TH1D *h);
void set_name_all(TH1D *h);
						 
void compare_xs(int ups);
void overlay_xs(string experiment);
void compare_cms2010(int ups);
void plot_atlas(int ups);
void plot_cms2010(int ups);
void compare_atlas(int ups);
void fill_yields(int iy, string mode);
void make_xs_stat_graph(int iy,int ups);
void xs_fit(int ups);
void xs_fit_hist(TF1 *xs_fit, int ups);
void divide_accept(int iy,int ups,string mode);
void print_tables(int iy, string mode);
void compute_ratio_function(int num,int den);
void fill_ratio_yield(int iy, int num);
void plot_ratios(int num, int den, string PP,bool show_fits);

void check_background_model(int iy, int ipt);
double integral_error(TF1 *FF, double a, double b);
double Cheb_int(double a,double b, int order);

double ups_errors(TF1 *f,TMatrixDSym cor, int ups);

void get_dm_m(int iy, int ipt, string mode);
void get_LS(int iy, int ipt, string mode);
double get_bin_center(int ipt, int ups);
void make_mdata(int iy, int ipt);
double adjust_uncertainties(TH1D *h, TH1D *hun_w);
void write_histos();
void setBW();
double expo_pol2(double *x, double *par);
double power_law(double *x, double*par);
double xs_full_func(double *x, double *par);
double resolutionP(TF1 *signal, int ups );
double resolutionM(TF1 *signal, int ups );
void make_sample_LS();
