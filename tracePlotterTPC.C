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
int maxCalculator(vector<double> trace);

void tracePlotterTPC(){
  //TFile *_file0 = TFile::Open("yso_vault05_bigEvent_DD.root");
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
  TH2D* dynTraceLow=new TH2D("dynTraceLow","dynTraces Low Gain",250.,0.,250.,5000,0.,5000.);
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
  TH2D* qdcRatio=new TH2D("qdcRatio","Long Integral vs. Short Integral",300.,0.,6000.,700.,-2000.,12000.);

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
    double dynodeElow = 0, dynodeEhigh=0, siE = 0;
    double qdcShortL = 0.0, qdcLongL = 0.0;
    double qdcShortH = 0.0, qdcLongH = 0.0;

    for(auto itC = rd.begin(); itC!=rd.end();itC++){
	std::vector<unsigned> trace = itC->trace;
	double energy = itC->energy;
	double time = itC->time;
	std::string dtype = itC->subtype.Data();
	std::string dgroup = itC->group.Data();
	int channel = itC->chanNum;  //used to try and gate on channel number rather than group/subtype

    if (dtype == "dynode_low"){
       vector<double> dynodeTrace;
      for (unsigned int it=0,itend = trace.size();it<itend;it++){
          dynTraceLow->Fill(it,trace.at(it));
          dynodeTrace.push_back(trace.at(it));
	  }
       	int maxLoc = maxCalculator(dynodeTrace);
        unsigned int lowBoundShort = maxLoc - 28;
        unsigned int highBoundShort = maxLoc + 20;
        unsigned int lowBoundLong = 1 + highBoundShort;
        unsigned int highBoundLong = 250;
        qdcShortL = QDCcalculator(dynodeTrace, lowBoundShort, highBoundShort);
        qdcLongL = QDCcalculator(dynodeTrace, lowBoundLong, highBoundLong);
        dynTime = time;
        dynodeElow = energy;
    }
     else if (channel == 6){
      for (unsigned int it=0,itend = trace.size();it<itend;it++){
          siTrace->Fill(it,trace.at(it));
      }
      siTime = time;
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
       for (unsigned int it=0,itend = trace.size();it<itend;it++){
           dynTraceHigh->Fill(it,trace.at(it));
           dynodeTrace.push_back(trace.at(it));
       }
        int maxLoc = maxCalculator(dynodeTrace);
        unsigned int lowBoundShort = maxLoc - 28;
        unsigned int highBoundShort = maxLoc + 20;
        unsigned int lowBoundLong = 1 + highBoundShort;
        unsigned int highBoundLong = 250;
        qdcShortH = QDCcalculator(dynodeTrace, lowBoundShort, highBoundShort);
        qdcLongH = QDCcalculator(dynodeTrace, lowBoundLong, highBoundLong);
        dynTime = time;
        dynodeEhigh = energy;
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

      //Section for filling plots
      
      siEnergy->Fill(siE);
      dynEnergy->Fill(dynodeEhigh); //can be changed to low
      qdcRatio->Fill(qdcShortH,qdcLongH);  //for the H for high gain L for low gain plotting
      
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

  c2->cd();
  qdcRatio->Draw("colz");
  
  c1->cd();
  dynTraceHigh->Draw("colz");
}

//Function for calculating the QDC

double QDCcalculator(vector<double> trace, unsigned int lowBound, unsigned int highBound){
    //first subtract baseline
    double baseSum = 0;
    for(int j=0;j<20;j++){
	baseSum += trace[j];
    }
    double  baseline = baseSum/20;
    for (int j=0; j<trace.size();j++){
        trace[j] = trace[j] - baseline;
    }
    /*if(highBound>trace.size()){
	cout<<"Upper bound is greater than trace length"<<endl;
	return 0;
    } */
    //calculate integral after baseline
    double integral = 0.0;
    for (int i=lowBound; i<highBound;i++){
        integral += 0.5 * (double(trace[i-1] + trace[i]));
    }
    return integral;
}

//Function for calculating the location of the maximum value in a trace

int maxCalculator(vector<double> trace){
	int maxLoc = 0;
	for(int i=1;i<trace.size();i++){
		    if(trace[i]>trace[maxLoc]){maxLoc = i;}
	}
	return maxLoc;
}
