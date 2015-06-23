//
//  plotter.c
//  
//
//  Created by Benjamin Carlson on 8/26/14.
//
//


#include "plotter.h"

void plotter(){
    cout << "plotter. " << endl;
    map_all();
    ifstream input;
    input.open("CMS_Cool_Test_140821_1043.txt");
    cout << "opened file. " << endl;
    std::vector<TString> label;
    std::vector<double> value;
    std::vector<TString> date;
    std::vector<TString> time;
    std::vector<TString> notes;
    std::vector<TString> flow;

    int L=0;
    int N=26+17;
    while(!input.eof()){
        
        double x;
        string txt;
        string line;
        std::vector<double> value_tmp;
        
        for(int i=0; i<N; i++){
            if(L==0){
                input >>txt;
                label.push_back(txt);
            }
         
            if(L>=1){
                if(i==0){
                    input >> txt;
                    date.push_back(txt);
                }
                if(i==1){
                    input >> txt;
                    time.push_back(txt);
                }
                if(i>1){
                    input >> x;
                    value.push_back(x);
                }
                if(i==N-1){
                    input >> txt;
                    notes.push_back(txt);
                    input >> txt;
                    flow.push_back(txt);
                }
            }

        }//loop over i
        
        getline(input,line);
        //cout << line << endl;
        L++;
    }
    cout << "Number of lines: " << L << endl;
    cout << "size of label: " << label.size() << endl;

    cout << "size of value: " << value.size() << endl;
    /*
    int hI=0;
    
    for(int imeas=0; imeas<L-2;imeas++){
        cout << "measurement: " << imeas <<  " ";
        TDatime DatTime=StringsToTime(date.at(imeas),time.at(imeas));
        cout << DatTime.GetHour() << " ";

        cout << label.at(hI+2+Nsens) << " " << value.at(heaterIndex(imeas,hI)) << endl;
    }
    */
    CreateCanvas("TimeSummmary", 800,600);
    CreateGraph("RT101");
    int sI=9;
    double T0=0;
    for(int imeas=0; imeas<L-2;imeas++){
       // cout << "measurement: " << imeas <<  " ";
        TDatime DatTime=StringsToTime(date.at(imeas),time.at(imeas));
        if(imeas==0)T0=DatTime.GetHour()*60+DatTime.GetMinute();
        cout << DatTime.GetHour() << " ";
        
        cout << label.at(sI+2) << " ";
        cout << value.at(heaterIndex(imeas,0)) << " ";
        cout << value.at(sensIndex(imeas,sI)) << " ";
        cout << notes.at(imeas) << " " << flow.at(imeas) << endl;
        grList.at(0)->SetPoint(imeas,DatTime.GetHour()*60+DatTime.GetMinute()-T0,value.at(sensIndex(imeas,sI)));
    }
    CName["TimeSummmary"]->cd();
    grList.at(0)->SetMarkerStyle(8);
    grList.at(0)->SetMarkerSize(1.5);
    grList.at(0)->SetMarkerColor(kRed);
    grList.at(0)->Draw("ap");
    grList.at(0)->GetYaxis()->SetRangeUser(-29.5,-28);
    grList.at(0)->GetYaxis()->SetTitle("C");
    grList.at(0)->GetXaxis()->SetTitle("min");

    CreateGraph2D("TempPlot4");
    TH2F *h2d=new TH2F("h2d",";X (cm); Y (cm)",15,0,12*2.54,30,0,15*2.54);
    

    int imeas=11;
    for(int isens=0; isens<Nsens;isens++){
        // cout << "measurement: " << imeas <<  " ";
        cout << "heat: " <<value.at(heaterIndex(imeas,0)) << endl;
        //cout << label.at(isens+2) << " " << value.at(sensIndex(imeas,isens)) << " ";
        cout << notes.at(imeas) << endl;
        h2d->Fill(pos[label.at(isens+2)]->X(),pos[label.at(isens+2)]->Y(),value.at(sensIndex(imeas,isens))  );
        //grName2D["TempPlot4"]->SetPoint(isens,pos[label.at(isens+2)]->X(),pos[label.at(isens+2)]->Y(),value.at(sensIndex(imeas,isens)));
    }
    CreateCanvas("TempPlot2D",600,800);
    CName["TempPlot2D"]->cd();
    CName["TempPlot2D"]->SetTopMargin(0.1);
    CName["TempPlot2D"]->SetRightMargin(0.25);
    int colors[]={7,417,3,633,2,802,0};
    double min= h2d->GetBinContent(h2d->GetMinimumBin());
    double max=-27.9;
    double step=(max-min)/5;
    
    cout << "min: " << min << " max " << max << " " << min+5*step << endl;
    
    double levels[]={min,min+step,min+2*step,min+3*step,min+4*step,max,0};
    h2d->SetContour(sizeof(levels)/sizeof(double),levels);
    h2d->SetMaximum(max);
    h2d->SetMinimum(min-0.1);
    //h2d->GetZaxis()->SetTitle("#circ C");
    gStyle->SetPalette(sizeof(colors)/sizeof(int),colors);
    h2d->SetStats(kFALSE);
    h2d->Draw("COLZ TEXT");
}

TDatime StringsToTime(TString date, TString time){
    TObjArray *Date_I=(TObjArray*)date.Tokenize("/");
    TObjArray *Time_I=(TObjArray*)time.Tokenize(":");
    int hr=0;
    int min=0;
    int sec=0;
    
    int yr=0;
    int mo=0;
    int day=0;
    
    for(int i=0; i<3; i++){
        TObjString *dt=(TObjString*)Date_I->At(i);
        TObjString *st=(TObjString*)Time_I->At(i);
        TString tmp=st->GetString();
        TString tmp_dat=dt->GetString();
        if(i==0)hr=tmp.Atof();
        if(i==1)min=tmp.Atof();
        if(i==2)sec=tmp.Atof();
        
        if(i==2)yr=tmp_dat.Atof();
        if(i==0)mo=tmp_dat.Atof();
        if(i==1)day=tmp_dat.Atof();
        
        delete st;
        delete dt;
    }
    yr+=2000;
    TDatime DatTime(yr,mo,day,hr,min,sec);
    return DatTime;
}

void CreateGraph(TString name){
    int N=100;
    TGraph *gr = new TGraph(N);
    gr->SetName(name);
    grList.push_back(gr);
}

void CreateGraph2D(TString name){
    int N=100;
    TGraph2D *gr = new TGraph2D();
    gr->SetName(name);
    grName2D[name]=gr;
}

int index(int imeas, int index){
    int N=41;
    return imeas*N+index;
}
int heaterIndex(int imeas, int iHeat){
    int N=41;
    return imeas*N+2+Nsens+iHeat;
}

int sensIndex(int imeas, int iSens){
    int N=41;
    return imeas*N+1+iSens;
}

void CreateCanvas(TString name, int Cx,int Cy){
    TCanvas *c = new TCanvas(name,"",Cx,Cy);
    CName[name]=c;
}

void map_all(){
    map_pos("RTD101",3,0.25);
    map_pos("RTD102",8.625,0.25);
    
    map_pos("RTD103",3,2.3750);
    map_pos("RTD104",6,2.3750);
    map_pos("RTD105",8.625,2.3750);
    
    map_pos("RTD106",1.5,5.75);
    map_pos("RTD107",3,5.75);
    map_pos("RTD108",4.5,5.75);
    map_pos("RTD109",6,5.75);
    map_pos("RTD110",7.5,5.75);
    map_pos("RTD111",9.0,5.75);
    map_pos("RTD112",10.5,5.75);

    map_pos("RTD113",1.5,9.25);
    map_pos("RTD114",3,9.25);
    map_pos("RTD115",4.5,9.25);
    map_pos("RTD116",6,9.25);
    map_pos("RTD117",7.5,9.25);
    map_pos("RTD118",9.0,9.25);
    map_pos("RTD119",10.5,9.25);
    
    map_pos("RTD120",3,12.6250);
    map_pos("RTD121",6,12.6250);
    map_pos("RTD122",8.6250,12.6250);
    
    map_pos("RTD123",3,14.25);
    map_pos("RTD124",8.6250,14.25);

    
}

void map_pos(TString name, double x, double y){
    x*=2.54;//convert in -> cm
    y*=2.54;
    TVector2 *p=new TVector2(x,y);
    pos[name]=p;
}


