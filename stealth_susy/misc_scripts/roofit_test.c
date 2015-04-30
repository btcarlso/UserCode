/*
 *  roofit_test.c
 *  
 *
 *  Created by Benjamin Carlson on 8/13/13.
 *  Copyright 2013 Carnegie Mellon University. All rights reserved.
 *
 */

//#include "roofit_test.h"

#include <iostream> 
#include <iomanip>
#include <fstream>
#include <math.h>
#include <vector>

#include <TMath.h>
#include <TF1.h>

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooFormulaVar.h"
#include "RooProdPdf.h"
#include "RooEfficiency.h"
#include "RooPolynomial.h"
#include "RooChi2Var.h"
#include "RooMinuit.h"
#include "RooChebychev.h"
#include "RooCategory.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooExtendPdf.h"
#include "RooArgList.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooRealSumPdf.h"

#ifndef __CINT__
#include "RooCFunction1Binding.h" 
#include "RooCFunction2Binding.h" 
#include "RooCFunction3Binding.h"
#endif
#include "RooTFnBinding.h" 

using namespace RooFit ;

double fit_weighted(double wgt, bool error_method){
	
	
	RooRealVar x1("x1","x1",8,11);
	RooRealVar sigma1("sigma1","sigma1",0.1,0.05,0.2) ;
	RooRealVar mean1("mean1","mean1",9.46,8,11) ;
	RooRealVar n1("yield1","number of events",1000,0,10000000) ;
	RooGaussian gauss1("gauss1","gauss1",x1,mean1,sigma1) ;
	RooExtendPdf egauss1("egauss1","extended gaussian PDF", gauss1,n1); 
	
	RooRealVar w("w","w",wgt); // weight variable 
	
	RooDataSet* data1 = egauss1.generate(x1,100);
	data1->addColumn(w);
	RooDataSet *wdata = new RooDataSet("wdata","wdata",*data1->get(),Import(*data1), WeightVar("w")); 
	
	RooFitResult* fitres(0); 
	if(error_method==0) fitres= egauss1.fitTo(*wdata,Extended(kTRUE),Save(),SumW2Error(kFALSE),PrintEvalErrors(-1),PrintLevel(-1),Verbose(kFALSE)) ;
	if(error_method==1) fitres= egauss1.fitTo(*wdata,Extended(kTRUE),Save(),SumW2Error(kTRUE),PrintEvalErrors(-1),PrintLevel(-1),Verbose(kFALSE)) ;

	
	double yield     = ((RooRealVar*)fitres->floatParsFinal().find("yield1"))->getVal();
	double yield_err = ((RooRealVar*)fitres->floatParsFinal().find("yield1"))->getError();
	
	//cout << "extended maximum likelihood: " << endl; 
	//printf("yield: %6.0f +/- %4.0f \t rel. error: %4.2f%%\n",yield, yield_err, yield_err/yield*100);
	
	return yield_err*100/yield;
	
}

double sqrt_fcn(double *x, double *par){
	return par[0]*TMath::Sqrt(1/x[0]); 
}

void roofit_test(){

	
	TGraph *weighed=new TGraph(20); 
	TGraph *weighed_true=new TGraph(20); 

	weighed->SetTitle("Weighted events: SumW2Error(kFALSE)"); 
	weighed_true->SetTitle("Weighted Events: SumW2Error(kTRUE)"); 
	
	double frac_error=0; 
	int iP=0; 
	for(double w=1.0; w<=100; w+=5){
		frac_error=fit_weighted(w,0);
		weighed->SetPoint(iP,w,frac_error);
		cout << "fractional error SumW2Error(kFALSE): " << frac_error << endl; 
		frac_error=fit_weighted(w,1);
		cout << "fractional error SumW2Error(kTRUE): " << frac_error << endl; 

		weighed_true->SetPoint(iP,w,frac_error);

		iP++;
	}
	TCanvas *C = new TCanvas();
	C->Divide(1,2); 
	C->cd(1); 
	TF1 *sqrtfcn = new TF1("sqrtfcn",sqrt_fcn,1,100,1); 
	weighed->Fit("sqrtfcn","R"); 
	weighed->Draw("ap");
	weighed->GetXaxis()->SetTitle("weight"); 
	weighed->GetYaxis()->SetTitle("fractional uncertainty"); 
	weighed->Draw("ap"); 


	C->cd(2); 
	weighed_true->Draw("ap"); 
	

	
}


