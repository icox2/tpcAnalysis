#include "TSystem.h"
#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TProfile.h"
#include "TChain.h"
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <ProcessorRootStruc.hpp>

double QDCcalculator(vector<double> trace, unsigned int lowBound, unsigned int highBound);

void tracePlotterTPC(){
  TFile *_file0 = TFile::Open("yso_vault05_bigEvent_DD.root");
  TTree *GS = (TTree*)_file0->Get("PixTree");
  TTreeReader singe;
  cout<<"readermade"<<endl;
  singe.SetTree( GS );
  //TTreeReaderValue<std::vector<unsigned int>> trace = {singe,"GS_Trace"};
  //TTreeReaderValue<Double_t> tdiff = {singe,"GS_BetaGammaTDiff"};
  //TTreeReaderValue<std::string> type = {singe,"GS_Type"};
  TTreeReaderArray<processor_struct::ROOTDEV> rd = {singe, "root_dev_vec_"};  //gives vector of stucture 
  TCanvas* c1 = new TCanvas("c1","c1",1200,600);
  TCanvas* c2 = new TCanvas("c2","c2",1200,600);
  c1->SetLogz();
  c2->SetLogz();
  TH2D* dynTraceLow=new TH2D("dynTraceLow","dynTraces Low Gain",5000.,0.,5000.,5000,0.,5000.);
  TH2D* dynTraceHigh=new TH2D("dynTraceHigh","dynTraces High Gain",250.,0.,250.,5000,0.,5000.);
  /*TH2D* xaTrace=new TH2D("xaTrace","xaTraces",250.,0.,250.,5000,0.,5000.);
  TH2D* xbTrace=new TH2D("xbTrace","xbTraces",250.,0.,250.,5000,0.,5000.);
  TH2D* yaTrace=new TH2D("yaTrace","yaTraces",250.,0.,250.,5000,0.,5000.);
  TH2D* ybTrace=new TH2D("ybTrace","ybTraces",250.,0.,250.,5000,0.,5000.); */
  TH2D* positionLow=new TH2D("positionLow","Positions Low Gain",250.,0.,25.,250.,0.,25.);
  TH2D* positionHigh=new TH2D("positionHigh","Positions High Gain",250.,0.,25.,250.,0.,25.);
  TH2D* siTrace=new TH2D("siTrace","siTraces",2500.,0.,2500.,2500,0.,2500.);
  TH1D* timeDiff=new TH1D("timeDiff","timeDifferences",1000.,-500.,500.);
  TH1D* siEnergy=new TH1D("siEnergy","siEnergies",1800.,0.,18000.);
  TH1D* dynEnergy=new TH1D("dynEnergy","dynEnergies",4000.,0.,40000.);
  TH2D* qdcRatio=new TH2D("qdcRatio","ShortQDC vs LongQDC",1000.,0.,10000.,1000.,0.,10000.);

  int eventNum=0;

  while(singe.Next()){
    //cout<<"while Loop"<<endl;
    /*std::vector<unsigned int> traceV = *trace;
    std::string dtype = *type;
    double bgTDiff = *tdiff;
    */	
    double xal=0, xbl=0, yal=0, ybl=0;
    double xposl=0, yposl=0;  //low gain
    double xah=0, xbh=0, yah=0, ybh=0;
    double xposh=0, yposh=0;  //high gain 
    double siTime=0, dynTime=0;
    double dynodeElow = 0, siE = 0;

    for(auto itC = rd.begin(); itC!=rd.end();itC++){
	std::vector<unsigned> trace = itC->trace;
	double energy = itC->energy;
	double time = itC->time;
	std::string dtype = itC->subtype.Data();
	std::string dgroup = itC->group.Data();
	int channel = itC->chanNum;  //used to try and gate on channel number rather than group/subtype

    if (dtype == "dynode_low"){
      for (unsigned int it=0,itend = trace.size();it<itend;it++){
	dynTraceLow->Fill(it,trace.at(it));
	}
	dynTime = time;
	dynEnergy->Fill(energy);
	dynodeElow = energy;
      }
     else if (channel == 6){
      for (unsigned int it=0,itend = trace.size();it<itend;it++){
	siTrace->Fill(it,trace.at(it));
	}
	siTime = time;
	siEnergy->Fill(energy);
	siE = energy;
      }
    else if (dtype == "anode_low" && dgroup == "xa"){
      /*for (unsigned int it=0,itend = trace.size();it<itend;it++){
          xaTrace->Fill(it,trace.at(it));
      }*/
	xal = energy;
      }
    else if (dtype == "anode_low" && dgroup == "xb"){
      /*for (unsigned int it=0,itend = trace.size();it<itend;it++){
          xbTrace->Fill(it,trace.at(it));
      }*/
	xbl = energy;
      }
    else if (dtype == "anode_low" && dgroup == "ya"){
      /*for (unsigned int it=0,itend = trace.size();it<itend;it++){
          yaTrace->Fill(it,trace.at(it));
      }*/
	yal = energy;
      }
    else if (dtype == "anode_low" && dgroup == "yb"){
      /*for (unsigned int it=0,itend = trace.size();it<itend;it++){
          ybTrace->Fill(it,trace.at(it));
      }*/
	ybl = energy;
      }
    else if(dtype == "dynode_high"){
        vector<double> dynodeTrace;
        double qdcShort = 0.0, qdcLong = 0.0;
        unsigned int lowBoundShort = 5;
        unsigned int highBoundShort = 7 + lowBoundShort;
        unsigned int lowBoundLong = 1 + highBoundShort;
        unsigned int highBoundLong = 40 + lowBoundLong;
      for (unsigned int it=0,itend = trace.size();it<itend;it++){
          dynTraceHigh->Fill(it,trace.at(it));
          dynodeTrace.push_back(trace.at(it));
      }
        qdcShort = QDCcalculator(dynodeTrace, lowBoundShort, highBoundShort);
        qdcLong = QDCcalculator(dynodeTrace, lowBoundLong, highBoundLong);
        qdcRatio->Fill(qdcShort,qdcLong);
      dynTime = time;
	  dynEnergy->Fill(energy);
	  dynodeElow = energy;
      
    }
    /*else if (dtype == "anode_high" && dgroup == "xa"){
      for (unsigned int it=0,itend = trace.size();it<itend;it++){
	//xaTrace->Fill(it,trace.at(it));
	}
	xah = energy;
      }
    else if (dtype == "anode_high" && dgroup == "xb"){
      for (unsigned int it=0,itend = trace.size();it<itend;it++){
	//xbTrace->Fill(it,trace.at(it));
	}
	xbh = energy;
      }
    else if (dtype == "anode_high" && dgroup == "ya"){
      for (unsigned int it=0,itend = trace.size();it<itend;it++){
	//yaTrace->Fill(it,trace.at(it));
	}
	yah = energy;
      }
    else if (dtype == "anode_high" && dgroup == "yb"){
      for (unsigned int it=0,itend = trace.size();it<itend;it++){
	//ybTrace->Fill(it,trace.at(it));
	}
	ybh = energy;
      } */
    }

   if(xah > 0 && xbh > 0 && yah > 0 && ybh > 0){
	//xpos = 25 * ((xa + xb) - (ya +yb)) / (xa + xb + ya + yb);
	//ypos = 25 * ((xa + yb) - (xb +ya)) / (xa + xb + ya + yb);
        xposh = 25 * (ybh + xah) / (xah + xbh + yah + ybh);
        yposh = 25 * (xah + xbh) / (xah + xbh + yah + ybh);
	positionHigh->Fill(xposh,yposh);
       }
   if(xal > 0 && xbl > 0 && yal > 0 && ybl > 0){
	//xpos = 25 * ((xa + xb) - (ya +yb)) / (xa + xb + ya + yb);
	//ypos = 25 * ((xa + yb) - (xb +ya)) / (xa + xb + ya + yb);
        xposl = 25 * (ybl + xal) / (xal + xbl + yal + ybl);
        yposl = 25 * (xal + xbl) / (xal + xbl + yal + ybl);
	positionLow->Fill(xposl,yposl);
	//siEnergy->Fill(siE);
       }
   timeDiff->Fill(siTime - dynTime);

   eventNum++;

  }

  c1->cd();
  positionLow->Draw("colz");
  
  c2->cd();
  siEnergy->Draw();
}

//Function for calculating the QDC

double QDCcalculator(vector<double> trace, unsigned int lowBound, unsigned int highBound){
    //first subtract baseline
    baseline = trace.begin();
    for (int j=0; j<trace.size();j++){
        trace[j] = trace[j] - baseline;
    }
    
    //calculate integral after baseline
    double integral = 0.0;
    for (int i=lowBound; i<highBound;i++){
        integral += 0.5 * (double(trace[i-1] + trace[i]))
    }
    return integral;
}
