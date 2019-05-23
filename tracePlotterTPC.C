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

double QDCcalculator(vector<unsigned int> trace, unsigned int lowBound, unsigned int highBound);
pair <int, int> boundsCalc(vector<unsigned int> trace, int maxPos);
pair <double, double> baselineCalc(vector<unsigned int> trace);

void tracePlotterTPC(){
  TFile *_file0 = TFile::Open("yso_csClose_DD.root");
  TTree *GS = (TTree*)_file0->Get("PixTree");
  TTreeReader singe;
  cout<<"readermade"<<endl;
  singe.SetTree( GS );
  //TTreeReaderValue<std::vector<unsigned int>> trace = {singe,"GS_Trace"};
  //TTreeReaderValue<Double_t> tdiff = {singe,"GS_BetaGammaTDiff"};
  //TTreeReaderValue<std::string> type = {singe,"GS_Type"};
  TTreeReaderArray<processor_struct::ROOTDEV> rd = {singe, "root_dev_vec_"};  //gives vector of stucture 
  TH2D* timeDiff=new TH2D("timeDiff","timeDifferences",600.,-1.5e-6,1.5e-6,4000.,0.,4000.);
  TH1D* siEnergy=new TH1D("siEnergy","siEnergies",1800.,0.,18000.);
  TH1D* siGate=new TH1D("siGate","Si gated by dynode",4000.,0.,4000.);
  TH1D* dynEnergy=new TH1D("dynEnergy","dynEnergies",8000.,-1000.,36000.);
  TH1D* dynGate=new TH1D("dynGate","Dynode gated by Si",8000.,-1000.,36000.);
   
  double xal=0.0, xbl=0.0, yal=0.0, ybl=0.0;
  double xposl=0.0, yposl=0.0;  //low gain
  double xah=0.0, xbh=0., yah=0., ybh=0.;
  double xposh=0., yposh=0.;  //high gain
  double dynTime=0.;
  double siTime = 0.0;
  double dynodeElow = 0., dynodeEhigh=0., siE = 0.;
  double qdcShortA = 0.0, qdcLongA = 0.0;
  double qdcShortH = 0.0, qdcLongH = 0.0;
  int maxLoc = 0, iDyn = 0, ixa=0,ixb=0,iya=0,iyb=0, minimumLoc=0;
  pair<int,int> Bound;
  bool goodTrace = false;
  vector<unsigned int> dynodeTrace, siliconTrace, xaTrace, xbTrace, yaTrace, ybTrace, anodeSum;
  int eventNum=0, iSi=0, undershoot=0;    
  double traceMax =0.0, stddev = 0.0, traceMin = 0.0;
  double qdcSxa = 0.0;
  double qdcSxb = 0.0;
  double qdcSya = 0.0;
  double qdcSyb = 0.0;
  double qdcLxa = 0.0;
  double qdcLxb = 0.0;
  double qdcLya = 0.0;
  double qdcLyb = 0.0;
	
  TFile *newFile = new TFile("newFile.root","RECREATE");
  TTree *newTree = new TTree("newTree","tree filled with traces and energies etc.");
  
  newTree->Branch("eventNum", &eventNum);
  newTree->Branch("dynodeTrace", &dynodeTrace);
  newTree->Branch("dynodeEnergy", &dynodeEhigh);
  newTree->Branch("siEnergy", &siE);
  newTree->Branch("siTrace", &siliconTrace);
  newTree->Branch("shortQDC", &qdcShortA);
  newTree->Branch("longQDC",&qdcLongA);
  newTree->Branch("shortQDCdyn", &qdcShortH);
  newTree->Branch("longQDCdyn",&qdcLongH);
  newTree->Branch("xpos", &xposl);
  newTree->Branch("ypos", &yposl);
  newTree->Branch("iDyn", &iDyn);
  newTree->Branch("ixa", &ixa);
  newTree->Branch("ixb", &ixb);
  newTree->Branch("iya", &iya);
  newTree->Branch("iyb", &iyb);
  newTree->Branch("iSi", &iSi);
  newTree->Branch("siTime", &siTime);
  newTree->Branch("dynTime", &dynTime);
  newTree->Branch("tmax", &traceMax);
  newTree->Branch("tmin", &traceMin);
  newTree->Branch("stddev", &stddev);
  newTree->Branch("undershootcount", &undershoot);
//    newTree->Branch("xaTrace", &xaTrace);
//    newTree->Branch("xbTrace", &xbTrace);
//    newTree->Branch("yaTrace", &yaTrace);
//    newTree->Branch("ybTrace", &ybTrace);

std::vector<unsigned> *trace;
  while(singe.Next()){
    //cout<<"while Loop"<<endl;
    /*std::vector<unsigned int> traceV = *trace;
    std::string dtype = *type;
    double bgTDiff = *tdiff;
    */
	//clearing the previous values from last event
	iDyn = 0;
	ixa = 0;
	ixb = 0;
	iya = 0;
	iyb = 0;
	iSi = 0; undershoot=0;
	qdcShortA = -999.; qdcLongA = -999.;
	qdcShortH = -999.; qdcLongH = -999.;
	siE = -999.;
	dynodeEhigh = -999.;
	xposl = 0.0;
	yposl = 0.0;
    xal=0.0; xbl=0.0; yal=0.0; ybl=0.0;
    stddev = 0.0; traceMax = 0.0; traceMin = 0.0;
	qdcSxa = 0.0; qdcSxb = 0.0; qdcSya = 0.0; qdcSyb = 0.0;
	qdcLxa = 0.0; qdcLxb = 0.0; qdcLya = 0.0; qdcLyb = 0.0;
	
	if (eventNum%50000==0) cout<<eventNum<<endl;

	goodTrace = false;      	
	dynodeTrace.clear();
	siliconTrace.clear();
    xaTrace.clear();
    xbTrace.clear();
    yaTrace.clear();
    ybTrace.clear();
    anodeSum.clear();
    siTime=0.0;
 
     for(auto itC = rd.begin(); itC!=rd.end();itC++){
	trace = &(itC->trace);
	double energy = itC->energy;
	double time = itC->time;
	std::string dtype = itC->subtype.Data();
	std::string dgroup = itC->group.Data();
	int channel = itC->chanNum;  //used to try and gate on channel number rather than group/subtype



/*    if (dtype == "dynode_low" && false){
	dynodeTraceLow = *trace;
        unsigned int lowBoundShort = maxLoc - 28;
        unsigned int highBoundShort = maxLoc + 20;
        unsigned int lowBoundLong = 1 + highBoundShort;
        unsigned int highBoundLong = 250;
        qdcShortL = QDCcalculator(dynodeTrace, lowBoundShort, highBoundShort);
        qdcLongL = QDCcalculator(dynodeTrace, lowBoundLong, highBoundLong);
        dynTime = time;
        dynodeElow = energy;
    }*/
      if (channel == 6){
	iSi++;
	siliconTrace = *trace;
        siTime = time;
        siE = energy;
	siEnergy->Fill(energy);
      }
    else if (dtype == "anode_low" && dgroup == "xa"){
      ixa++;
      xaTrace = *trace;
      xal = energy;
      maxLoc = std::distance(xaTrace.begin(),std::max_element(xaTrace.begin(),xaTrace.end()));
      Bound = boundsCalc(xaTrace, maxLoc);
      qdcSxa = QDCcalculator(xaTrace, Bound.first, Bound.second);
      qdcLxa = QDCcalculator(xaTrace, Bound.second+1, 250);
      }
    else if (dtype == "anode_low" && dgroup == "xb"){
      ixb++;
      xbTrace = *trace;
      xbl = energy;
      maxLoc = std::distance(xbTrace.begin(),std::max_element(xbTrace.begin(),xbTrace.end()));
      Bound = boundsCalc(xbTrace, maxLoc);
      qdcSxb = QDCcalculator(xbTrace, Bound.first, Bound.second);
      qdcLxb = QDCcalculator(xbTrace, Bound.second+1, 250);
      }
    else if (dtype == "anode_low" && dgroup == "ya"){
	iya++;
	yaTrace = *trace;
	yal = energy;
	maxLoc = std::distance(yaTrace.begin(),std::max_element(yaTrace.begin(),yaTrace.end()));
      Bound = boundsCalc(yaTrace, maxLoc);
      qdcSya = QDCcalculator(yaTrace, Bound.first, Bound.second);
      qdcLya = QDCcalculator(yaTrace, Bound.second+1, 250); 
      }
    else if (dtype == "anode_low" && dgroup == "yb"){
	iyb++;
	ybTrace = *trace;
	ybl = energy;
      maxLoc = std::distance(ybTrace.begin(),std::max_element(ybTrace.begin(),ybTrace.end()));
      Bound = boundsCalc(ybTrace, maxLoc);
      qdcSyb = QDCcalculator(ybTrace, Bound.first, Bound.second);
      qdcLyb = QDCcalculator(ybTrace, Bound.second+1, 250);
      }
    else if(dtype == "dynode_high"){
       iDyn++;
       if(iDyn==1){
	dynodeTrace = *trace;
        maxLoc = std::distance(dynodeTrace.begin(),std::max_element(dynodeTrace.begin(),dynodeTrace.end()));
	minimumLoc = std::distance(dynodeTrace.begin(),std::min_element(dynodeTrace.begin(),dynodeTrace.end()));
        Bound = boundsCalc(dynodeTrace, maxLoc);
	traceMax = dynodeTrace[maxLoc];
	traceMin = dynodeTrace[minimumLoc];
        stddev = baselineCalc(dynodeTrace).second;
	if(traceMin<baselineCalc(dynodeTrace).first-6*stddev){undershoot++;}
        unsigned int lowBoundShort = Bound.first;
        unsigned int highBoundShort = Bound.second;
        unsigned int lowBoundLong = 1 + highBoundShort;
        unsigned int highBoundLong = 250;
        dynTime = time;
        dynodeEhigh = energy;
	dynEnergy->Fill(energy);
        if(traceMax<4090){   //avoiding saturated traces
            qdcShortH = QDCcalculator(dynodeTrace, lowBoundShort, highBoundShort);
	    qdcLongH = QDCcalculator(dynodeTrace, lowBoundLong, highBoundLong);
            goodTrace = true;
	 }
       
     }
    }
    /*else if (dtype == "anode_high" && dgroup == "xa"){
	xaTrace = *trace;      
	xah = energy;
      }
    else if (dtype == "anode_high" && dgroup == "xb"){
      	xbTrace = *trace;
	xbh = energy;
      }
    else if (dtype == "anode_high" && dgroup == "ya"){
	yaTrace = *trace;      
	yah = energy;
      }
    else if (dtype == "anode_high" && dgroup == "yb"){
     	ybTrace = *trace;
	ybh = energy;
      } */
      //if(qdcLongH>0.5*qdcShortH && qdcLongH<2*qdcShortH) 



    }

      //Section for filling plots
    if(ixa>0 && ixb>0 && iya>0 && iyb>0){
        qdcShortH = 0.0;
	qdcLongH = 0.0;
        qdcShortA = qdcSxa + qdcSxb + qdcSya + qdcSyb;
        qdcLongA = qdcLxa + qdcLxb + qdcLya + qdcLyb;
   }

    if(iSi>0 && iDyn>0){
      siGate->Fill(siE);
      dynGate->Fill(dynodeEhigh);
      timeDiff->Fill(siTime-dynTime,siE);
    }
    
   if(xah > 0 && xbh > 0 && yah > 0 && ybh > 0){
	//xpos = 25 * ((xa + xb) - (ya +yb)) / (xa + xb + ya + yb);
	//ypos = 25 * ((xa + yb) - (xb +ya)) / (xa + xb + ya + yb);
        xposh = 25 * (ybh + xah) / (xah + xbh + yah + ybh);
        yposh = 25 * (xah + xbh) / (xah + xbh + yah + ybh);
       }
   if(xal > 0 && xbl > 0 && yal > 0 && ybl > 0){
	//xpos = 25 * ((xa + xb) - (ya +yb)) / (xa + xb + ya + yb);
	//ypos = 25 * ((xa + yb) - (xb +ya)) / (xa + xb + ya + yb);
        xposl = 25 * (ybl + xal) / (xal + xbl + yal + ybl);
        yposl = 25 * (xal + xbl) / (xal + xbl + yal + ybl);
       }

   newTree->Fill();
	eventNum++;

  }

   
    newTree->Write();
	timeDiff->Write();
	siEnergy->Write();
	dynEnergy->Write();
	siGate->Write();
	dynGate->Write();
}

//Function for calculating the QDC

double QDCcalculator(vector<unsigned int> trace, unsigned int lowBound, unsigned int highBound){

    double  baseline = baselineCalc(trace).first;
   //calculate integral after baseline
    double integral = 0.0;
    for (int i=lowBound; i<highBound;i++){
        integral += 0.5 * (double(abs(trace[i-1] - baseline) + abs(trace[i] - baseline)));
    }
    return integral;
}

//Function for calculating upper and lower bounds for short due to pulse height
pair <int, int> boundsCalc(vector<unsigned int> trace, int maxPos){
    double maxVal = trace[maxPos];
    return std::make_pair (maxPos-(maxVal/100+1),maxPos+(maxVal/50+1));
    /*if(maxVal<600){
        return std::make_pair (maxPos-5,maxPos+10);
    }
    else{
        return std::make_pair (maxPos-10, maxPos+15);
    }*/
}

pair <double, double> baselineCalc(vector<unsigned int> trace){
    //calculating the baseline
    double baseSum = 0.;
    for(int j=0;j<20;j++){
        baseSum += trace[j];
    }
    double  baseline = baseSum/20.;
    
    //calculating the standard dev
    double stddev = 0.0;
    for(int j=0;j<20;j++){
        stddev += pow(trace[j]-baseline,2);
    }
    stddev = sqrt(stddev/20);
    
    return std::make_pair (baseline, stddev);
}
