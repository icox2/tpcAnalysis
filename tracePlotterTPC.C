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
  TFile *_file0 = TFile::Open("cs_test_1000V_DD.root");
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
  double qdcShortAl = 0.0, qdcLongAl = 0.0;
  double qdcShortAh = 0.0, qdcLongAh = 0.0;
  double qdcShortH = 0.0, qdcLongH = 0.0;
  double anodeRatio[4];
  int maxLoc = 0, iDyn = 0, ixa=0,ixb=0,iya=0,iyb=0, minimumLoc=0;
  pair<int,int> Bound;
  bool goodTrace = false;
  vector<unsigned int> dynodeTrace, siliconTrace, xaTracel, xbTracel, yaTracel, ybTracel, anodeSum;
  vector<unsigned int> xaTraceh, xbTraceh, yaTraceh, ybTraceh;
  int eventNum=0, iSi=0, undershoot=0;    
  double traceMax =0.0, stddev = 0.0, traceMin = 0.0;
  double qdcS[8];//first 4 are low gain, second 4 are high gain
  double qdcL[8];
  unsigned int tmax[8];
  double qdcSum[4];
  double scale = 1.;
	
  TFile *newFile = new TFile("newFile.root","RECREATE");
  TTree *newTree = new TTree("newTree","tree filled with traces and energies etc.");
  
  newTree->Branch("eventNum", &eventNum);
  newTree->Branch("dynodeTrace", &dynodeTrace);
  newTree->Branch("dynodeEnergy", &dynodeEhigh);
  newTree->Branch("siEnergy", &siE);
  newTree->Branch("siTrace", &siliconTrace);
  newTree->Branch("shortQDCl", &qdcShortAl);
  newTree->Branch("longQDCl",&qdcLongAl);
  newTree->Branch("shortQDCh", &qdcShortAh);
  newTree->Branch("longQDCh",&qdcLongAh);
  newTree->Branch("shortQDCdyn", &qdcShortH);
  newTree->Branch("longQDCdyn",&qdcLongH);
  newTree->Branch("xposl", &xposl);
  newTree->Branch("yposl", &yposl);
  newTree->Branch("xposh", &xposh);
  newTree->Branch("yposh", &yposh);
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
  newTree->Branch("anodeRatio", anodeRatio);
  newTree->Branch("qdcSum", qdcSum);
  newTree->Branch("undershootcount", &undershoot);
  newTree->Branch("xaTrace", &xaTracel);
  newTree->Branch("xbTrace", &xbTracel);
  newTree->Branch("yaTrace", &yaTracel);
  newTree->Branch("ybTrace", &ybTracel);

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
	qdcShortAl = -999.; qdcLongAl = -999.;
	qdcShortH = -999.; qdcLongH = -999.;
	siE = -999.;
	dynodeEhigh = -999.;
	xposl = 0.0;
	yposl = 0.0;
    xal=0.0; xbl=0.0; yal=0.0; ybl=0.0;
    stddev = 0.0; traceMax = 0.0; traceMin = 0.0;
	for(int i=0;i<4;i++){
		anodeRatio[i]=-999.;
		scale=1.;
		qdcS[2*i] = -999.;
		qdcL[2*i] = -999.;
		qdcS[2*i+1] = -999.;
		qdcL[2*i+1] = -999.;
		tmax[2*i]=0;
		tmax[2*i+1]=0;
	}
	
	if (eventNum%50000==0) cout<<eventNum<<endl;

	goodTrace = false;      	
	dynodeTrace.clear();
	siliconTrace.clear();
    xaTracel.clear();
    xbTracel.clear();
    yaTracel.clear();
    ybTracel.clear();
    xaTraceh.clear();
    xbTraceh.clear();
    yaTraceh.clear();
    ybTraceh.clear();
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
        unsigned int highBoundLong = 300;
        qdcShortL = QDCcalculator(dynodeTrace, lowBoundShort, highBoundShort);
        qdcLongL = QDCcalculator(dynodeTrace, lowBoundLong, highBoundLong);
        dynTime = time;
        dynodeElow = energy;
    }
      if (channel == 6){
	iSi++;
	siliconTrace = *trace;
        siTime = time;
        siE = energy;
	siEnergy->Fill(energy);
      }
    else*/ if (dtype == "anode_low" && dgroup == "xa"){
      ixa++;
      xaTracel = *trace;
      xal = energy;
      maxLoc = std::distance(xaTracel.begin(),std::max_element(xaTracel.begin(),xaTracel.end()));
      tmax[0] = *std::max_element(xaTracel.begin(),xaTracel.end());
      Bound = boundsCalc(xaTracel, maxLoc);
      if(maxLoc>450) continue;
      qdcS[0] = QDCcalculator(xaTracel, Bound.first, Bound.second);
      qdcL[0] = QDCcalculator(xaTracel, Bound.second+1, 500);
      }
    else if (dtype == "anode_low" && dgroup == "xb"){
      ixb++;
      xbTracel = *trace;
      xbl = energy;
      maxLoc = std::distance(xbTracel.begin(),std::max_element(xbTracel.begin(),xbTracel.end()));
      tmax[1] = *std::max_element(xbTracel.begin(),xbTracel.end());
      Bound = boundsCalc(xbTracel, maxLoc);
      if(maxLoc>450) continue;
      qdcS[1] = QDCcalculator(xbTracel, Bound.first, Bound.second);
      qdcL[1] = QDCcalculator(xbTracel, Bound.second+1, 500);
      }
    else if (dtype == "anode_low" && dgroup == "ya"){
	iya++;
	yaTracel = *trace;
	yal = energy;
	maxLoc = std::distance(yaTracel.begin(),std::max_element(yaTracel.begin(),yaTracel.end()));
      tmax[2] = *std::max_element(yaTracel.begin(),yaTracel.end());
      Bound = boundsCalc(yaTracel, maxLoc);
      if(maxLoc>450) continue;
      qdcS[2] = QDCcalculator(yaTracel, Bound.first, Bound.second);
      qdcL[2] = QDCcalculator(yaTracel, Bound.second+1, 500); 
      }
    else if (dtype == "anode_low" && dgroup == "yb"){
	iyb++;
	ybTracel = *trace;
	ybl = energy;
      maxLoc = std::distance(ybTracel.begin(),std::max_element(ybTracel.begin(),ybTracel.end()));
      tmax[3] = *std::max_element(ybTracel.begin(),ybTracel.end());
      Bound = boundsCalc(ybTracel, maxLoc);
      if(maxLoc>450) continue;
      qdcS[3] = QDCcalculator(ybTracel, Bound.first, Bound.second);
      qdcL[3] = QDCcalculator(ybTracel, Bound.second+1, 500);
      }
    else if(dtype == "dynode_high"){
       iDyn++;
       if(iDyn==1){
	dynodeTrace = *trace;
        maxLoc = std::distance(dynodeTrace.begin(),std::max_element(dynodeTrace.begin(),dynodeTrace.end()));
	minimumLoc = std::distance(dynodeTrace.begin(),std::min_element(dynodeTrace.begin(),dynodeTrace.end()));
        Bound = boundsCalc(dynodeTrace, maxLoc);
	traceMax = dynodeTrace[maxLoc];
	if(traceMax>16380) maxLoc += 5;
	traceMin = dynodeTrace[minimumLoc];
        stddev = baselineCalc(dynodeTrace).second;
	if(traceMin<baselineCalc(dynodeTrace).first-3*stddev){undershoot++;}
        unsigned int lowBoundShort = Bound.first;
        unsigned int highBoundShort = Bound.second;
        unsigned int lowBoundLong = 1 + highBoundShort;
        unsigned int highBoundLong = 300;
        if(maxLoc>250){/*cout<<"Broken Bounds"<<endl;*/ continue;}
	dynTime = time;
        dynodeEhigh = energy;
	dynEnergy->Fill(energy);
        //if(traceMax<4090){   //avoiding saturated traces
            qdcShortH = QDCcalculator(dynodeTrace, lowBoundShort, highBoundShort);
	    qdcLongH = QDCcalculator(dynodeTrace, lowBoundLong, highBoundLong);
            goodTrace = true;
    	    if(qdcShortH<1 || qdcLongH<1){cout<<"QDC < 1"<<endl;}
	 //}
       
     }
    }
    else if (dtype == "anode_high" && dgroup == "xa"){
	xaTraceh = *trace;      
	xah = energy;
      maxLoc = std::distance(xaTraceh.begin(),std::max_element(xaTraceh.begin(),xaTraceh.end()));
      tmax[4] = *std::max_element(xaTraceh.begin(),xaTraceh.end());
      Bound = boundsCalc(xaTraceh, maxLoc);
      if(maxLoc>450) continue;
      qdcS[4] = QDCcalculator(xaTraceh, Bound.first, Bound.second);
      qdcL[4] = QDCcalculator(xaTraceh, Bound.second+1, 500);
      }
    else if (dtype == "anode_high" && dgroup == "xb"){
      	xbTraceh = *trace;
	xbh = energy;
      maxLoc = std::distance(xbTraceh.begin(),std::max_element(xbTraceh.begin(),xbTraceh.end()));
      tmax[5] = *std::max_element(xbTraceh.begin(),xbTraceh.end());
      Bound = boundsCalc(xbTraceh, maxLoc);
      if(maxLoc>450) continue;
      qdcS[5] = QDCcalculator(xbTraceh, Bound.first, Bound.second);
      qdcL[5] = QDCcalculator(xbTraceh, Bound.second+1, 500);
      }
    else if (dtype == "anode_high" && dgroup == "ya"){
	yaTraceh = *trace;      
	yah = energy;
      maxLoc = std::distance(yaTraceh.begin(),std::max_element(yaTraceh.begin(),yaTraceh.end()));
      tmax[6] = *std::max_element(yaTraceh.begin(),yaTraceh.end());
      Bound = boundsCalc(yaTraceh, maxLoc);
      if(maxLoc>450) continue;
      qdcS[6] = QDCcalculator(yaTraceh, Bound.first, Bound.second);
      qdcL[6] = QDCcalculator(yaTraceh, Bound.second+1, 500);
      }
    else if (dtype == "anode_high" && dgroup == "yb"){
     	ybTraceh = *trace;
	ybh = energy;
      maxLoc = std::distance(ybTraceh.begin(),std::max_element(ybTraceh.begin(),ybTraceh.end()));
      tmax[7] = *std::max_element(ybTraceh.begin(),ybTraceh.end());
      Bound = boundsCalc(ybTraceh, maxLoc);
      if(maxLoc>450) continue;
      qdcS[7] = QDCcalculator(ybTraceh, Bound.first, Bound.second);
      qdcL[7] = QDCcalculator(ybTraceh, Bound.second+1, 500);
      } 


    }

      //Section for filling plots
    if(ixa>0 && ixb>0 && iya>0 && iyb>0){
        qdcShortAl = 0.0;
	qdcLongAl = 0.0;
        qdcShortAh = 0.0;
	qdcLongAh = 0.0;
	for(int i=0;i<4;i++){
	  scale = qdcS[0]/qdcS[i];
	  qdcSum[i] = qdcS[i] + qdcL[i];
	  qdcShortAl += scale*qdcS[i];
	  scale = qdcL[0]/qdcL[i];
	  qdcLongAl += scale*qdcL[i];
	  anodeRatio[i] = qdcS[i]/qdcSum[i];
	  scale = qdcS[4]/qdcS[i+4];
	  qdcShortAh += scale*qdcS[i+4];
	  scale = qdcL[4]/qdcL[i+4];
	  qdcLongAh += scale*qdcL[i+4];
	}
   }

    if(iSi>0 && iDyn>0){
      siGate->Fill(siE);
      dynGate->Fill(dynodeEhigh);
      timeDiff->Fill(siTime-dynTime,siE);
    }
    
   if(qdcS[0]+qdcL[0] > 0 && qdcS[1]+qdcL[1] > 0 && qdcS[2]+qdcL[2] > 0 && qdcS[3]+qdcL[3] > 0){
	//xpos = 25 * ((xa + xb) - (ya +yb)) / (xa + xb + ya + yb);
	//ypos = 25 * ((xa + yb) - (xb +ya)) / (xa + xb + ya + yb);
        xposl = 25 * (qdcS[3]+qdcL[3] + qdcS[0]+qdcL[0]) / (qdcS[0]+qdcL[0] + qdcS[1]+qdcL[1] + qdcS[2]+qdcL[2] + qdcS[3]+qdcL[3]);
        yposl = 25 * (qdcS[0]+qdcL[0] + qdcS[1]+qdcL[1]) / (qdcS[0]+qdcL[0] + qdcS[1]+qdcL[1] + qdcS[2]+qdcL[2] + qdcS[3]+qdcL[3]);
       }

   if(qdcS[4]+qdcL[4] > 0 && qdcS[5]+qdcL[5] > 0 && qdcS[6]+qdcL[6] > 0 && qdcS[7]+qdcL[7] > 0){
        xposh = 25 * (qdcS[7]+qdcL[7] + qdcS[4]+qdcL[4]) / (qdcS[4]+qdcL[4] + qdcS[5]+qdcL[5] + qdcS[6]+qdcL[6] + qdcS[7]+qdcL[7]);
        yposh = 25 * (qdcS[4]+qdcL[4] + qdcS[5]+qdcL[5]) / (qdcS[4]+qdcL[4] + qdcS[5]+qdcL[5] + qdcS[6]+qdcL[6] + qdcS[7]+qdcL[7]);
       }

  if(iDyn>0 || (ixa>0 && ixb>0 && iya>0 && iyb>0)) 
   newTree->Fill();
	eventNum++;

  }
    if(iDyn>0 || ixa>0 || ixb>0 || iya>0 || iyb>0 || iSi>0)
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
    return std::make_pair (maxPos-25,maxPos+25);
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
