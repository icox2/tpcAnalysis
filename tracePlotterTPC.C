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
pair <int, int> maxCalculator(vector<double> trace);
pair <int, int> boundsCalc(vector<double> trace, int maxPos);
pair <double, double> baselineCalc(vector<double> trace);

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
//  TCanvas* c1 = new TCanvas("c1","c1",1200,600);
//  TCanvas* c2 = new TCanvas("c2","c2",1200,600);
//  c1->SetLogz();
//  c2->SetLogz();
  //TH2D* dynTraceLow=new TH2D("dynTraceLow","dynTraces Low Gain",250.,0.,250.,5000,0.,5000.);
  TH2D* dynTraceHigh=new TH2D("dynTraceHigh","dynTraces High Gain",250.,0.,250.,5000,0.,5000.);
  /*TH2D* xaTrace=new TH2D("xaTrace","xaTraces",250.,0.,250.,5000,0.,5000.);
  TH2D* xbTrace=new TH2D("xbTrace","xbTraces",250.,0.,250.,5000,0.,5000.);
  TH2D* yaTrace=new TH2D("yaTrace","yaTraces",250.,0.,250.,5000,0.,5000.);
  TH2D* ybTrace=new TH2D("ybTrace","ybTraces",250.,0.,250.,5000,0.,5000.); */
  TH2D* positionLow=new TH2D("positionLow","Positions Low Gain",250.,0.,25.,250.,0.,25.);
  //TH2D* positionHigh=new TH2D("positionHigh","Positions High Gain",250.,0.,25.,250.,0.,25.);
  //TH2D* siTrace=new TH2D("siTrace","siTraces",2500.,0.,2500.,2500,0.,2500.);
  //TH1D* timeDiff=new TH1D("timeDiff","timeDifferences",1000.,-500.,500.);
  //TH1D* siEnergy=new TH1D("siEnergy","siEnergies",1800.,0.,18000.);
  //TH1D* dynEnergy=new TH1D("dynEnergy","dynEnergies",4000.,0.,40000.);
  TH2D* qdcRatio=new TH2D("qdcRatio","Short/Total vs Total",1000.,0.,10000.,100.,0.,1.);
  //TH1D* traceRatio=new TH1D("traceRatio","trace ratio",100.,0.,1.);
    
    double xal=0.0, xbl=0.0, yal=0.0, ybl=0.0;
    double xposl=0.0, yposl=0.0;  //low gain
    double xah=0.0, xbh=0., yah=0., ybh=0.;
    double xposh=0., yposh=0.;  //high gain
    double siTime=0., dynTime=0.;
    double dynodeElow = 0., dynodeEhigh=0., siE = 0.;
    double qdcShortL = 0.0, qdcLongL = 0.0;
    double qdcShortH = 0.0, qdcLongH = 0.0;
    int maxLoc = 0, iDyn = 0, ixa=0,ixb=0,iya=0,iyb=0;
    pair<int,int> Bound;
    bool goodTrace = false;
    vector<double> dynodeTrace, siliconTrace, xaTrace, xbTrace, yaTrace, ybTrace, anodeSum;
    int eventNum=0, iSi=0, undershoot=0;    
    double traceMax =0.0, stddev = 0.0;

    TFile *newFile = new TFile("newFile.root","RECREATE");
    TTree *newTree = new TTree("newTree","tree filled with traces and energies etc.");
  
    newTree->Branch("eventNum", &eventNum);
    newTree->Branch("dynodeTrace", &dynodeTrace);
    newTree->Branch("dynodeEnergy", &dynodeEhigh);
    newTree->Branch("siEnergy", &siE);
    newTree->Branch("siTrace", &siliconTrace);
    newTree->Branch("shortQDC", &qdcShortH);
    newTree->Branch("longQDC",&qdcLongH);
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
    newTree->Branch("stddev", &stddev);
    newTree->Branch("undershootcount", &undershoot);
    newTree->Branch("xaTrace", &xaTrace)
    newTree->Branch("xbTrace", &xbTrace)
    newTree->Branch("yaTrace", &yaTrace)
    newTree->Branch("ybTrace", &ybTrace)
    newTree->Branch("anodeSum", &anodeSum)

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
	qdcShortL = -999.; qdcLongL = -999;
	qdcShortH = -999; qdcLongH = -999;
	siE = -999;
	dynodeEhigh = -999;
	xposl = 0.0;
	yposl = 0.0;
    xal=0.0; xbl=0.0; yal=0.0; ybl=0.0;
    stddev = 0.0;
	
	if (eventNum%50000==0) cout<<eventNum<<endl;
//        if (eventNum==50000) break;
	goodTrace = false;      	
	dynodeTrace.clear();
	siliconTrace.clear();
    xaTrace.clear();
    xbTrace.clear();
    yaTrace.clear();
    ybTrace.clear();
    anodeSum.clear();
 
     for(auto itC = rd.begin(); itC!=rd.end();itC++){
	trace = &(itC->trace);
	double energy = itC->energy;
	double time = itC->time;
	std::string dtype = itC->subtype.Data();
	std::string dgroup = itC->group.Data();
	int channel = itC->chanNum;  //used to try and gate on channel number rather than group/subtype



/*    if (dtype == "dynode_low" && false){
      for (unsigned int it=0,itend = trace.size();it<itend;it++){
          dynTraceLow->Fill(it,trace.at(it));
          //dynodeTrace.push_back(trace.at(it));
	  }
       	//maxLoc = maxCalculator(dynodeTrace);
        unsigned int lowBoundShort = maxLoc - 28;
        unsigned int highBoundShort = maxLoc + 20;
        unsigned int lowBoundLong = 1 + highBoundShort;
        unsigned int highBoundLong = 250;
        qdcShortL = QDCcalculator(dynodeTrace, lowBoundShort, highBoundShort);
        qdcLongL = QDCcalculator(dynodeTrace, lowBoundLong, highBoundLong);
        dynTime = time;
        dynodeElow = energy;
       qdcRatio->Fill(qdcShortL,qdcLongL);  //for the H for high gain L for low gain plotting
	traceRatio->Fill(qdcShortL/(qdcShortL+qdcLongL));    
    }*/
      if (channel == 6){
	iSi++;
      //if(iSi==1){
      for (unsigned int it=0,itend = trace->size();it<itend;it++){
          siTrace->Fill(it,trace->at(it));
          siliconTrace.push_back(trace->at(it));
      }
      siTime = time;
	  siE = energy;
      siEnergy->Fill(siE);
      }
      //}
    else if (dtype == "anode_low" && dgroup == "xa"){
      ixa++;
      for (unsigned int it=0,itend = trace.size();it<itend;it++){
          xaTrace.push_back(trace.at(it));
      }
	xal = energy;
      }
    else if (dtype == "anode_low" && dgroup == "xb"){
	ixb++;
      for (unsigned int it=0,itend = trace.size();it<itend;it++){
          xbTrace.push_back(it,trace.at(it));
      }
	xbl = energy;
      }
    else if (dtype == "anode_low" && dgroup == "ya"){
	iya++;
      for (unsigned int it=0,itend = trace.size();it<itend;it++){
          yaTrace.push_back(it,trace.at(it));
      }
	yal = energy;
      }
    else if (dtype == "anode_low" && dgroup == "yb"){
	iyb++;
      for (unsigned int it=0,itend = trace.size();it<itend;it++){
          ybTrace.push_back(it,trace.at(it));
      }
	ybl = energy;
      }
    else if(dtype == "dynode_high"){
       iDyn++;
       if(iDyn==1){
	for (unsigned int it=0,itend = trace->size();it<itend;it++){
           dynTraceHigh->Fill(it,trace->at(it));
           dynodeTrace.push_back(trace->at(it));
        }
        maxLoc = maxCalculator(dynodeTrace).first;
        Bound = boundsCalc(dynodeTrace, maxLoc);
	traceMax = dynodeTrace[maxLoc];
        stddev = baselineCalc(dynodeTrace).second;
	if(maxCalculator(dynodeTrace).second<baselineCalc(dynodeTrace).first-6*stddev){undershoot++;}
        unsigned int lowBoundShort = Bound.first;
        unsigned int highBoundShort = Bound.second;
        unsigned int lowBoundLong = 1 + highBoundShort;
        unsigned int highBoundLong = 250;
        dynTime = time;
        dynodeEhigh = energy;
         if(traceMax<4090){
            qdcShortH = QDCcalculator(dynodeTrace, lowBoundShort, highBoundShort);
        	qdcLongH = QDCcalculator(dynodeTrace, lowBoundLong, highBoundLong);
        	traceRatio->Fill(qdcShortH/(qdcShortH+qdcLongH));    
                goodTrace = true;
		//cout << "Long: " << qdcLongH << ", ShortH: " << qdcShortH << endl;
                //newTree->Fill();
	 }
       
      dynEnergy->Fill(dynodeEhigh); //can be changed to low
       qdcRatio->Fill((qdcLongH+qdcShortH),qdcShortH/(qdcLongH+qdcShortH));  //for the H for high gain L for low gain plotting
     }
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
      //if(qdcLongH>0.5*qdcShortH && qdcLongH<2*qdcShortH) 



    }

      //Section for filling plots
    if(ixa>0 && ixb>0 && iya>0 && iyb>0){
        for(int i=0;i<250;i++){
            anodeSum.push_back(xaTrace[i]);
            anodeSum[i] += xbTrace[i];
            anodeSum[i] += yaTrace[i];
            anodeSum[i] += ybTrace[i];
        }
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
  //if(goodTrace){ newTree->Fill();}
   newTree->Fill();
	eventNum++;

  }

//  c2->cd();
//  qdcRatio->Draw("colz");
//  
//  c1->cd();
//  dynTraceHigh->Draw("colz");
    
    newTree->Write();
	qdcRatio->Write();
	dynTraceHigh->Write();
	positionLow->Write();
 
 TCut XY = "ixa>0 && ixb>0 && iya>0 && iyb>0";
 XY.Write();

}

//Function for calculating the QDC

double QDCcalculator(vector<double> trace, unsigned int lowBound, unsigned int highBound){

    double  baseline = baselineCalc(trace).first;
   //calculate integral after baseline
    double integral = 0.0;
    for (int i=lowBound; i<highBound;i++){
        integral += 0.5 * (double(abs(trace[i-1] - baseline) + abs(trace[i] - baseline)));
    }
    return integral;
}

//Function for calculating the location of the maximum value in a trace

pair<int, int> maxCalculator(vector<double> trace){
	int mLoc = 0;
	int minLoc=0;
	for(int i=1;i<trace.size();i++){
		    if(trace[i]>trace[mLoc]){mLoc = i;}
		if(trace[i]<trace[minLoc]){minLoc = i;}
	}
	return std::make_pair (mLoc,minLoc);
}

//Function for calculating upper and lower bounds for short due to pulse height
pair <int, int> boundsCalc(vector<double> trace, int maxPos){
    double maxVal = trace[maxPos];
    return std::make_pair (maxPos-(maxVal/100+1),maxPos+(maxVal/50+1));
    /*if(maxVal<600){
        return std::make_pair (maxPos-5,maxPos+10);
    }
    else{
        return std::make_pair (maxPos-10, maxPos+15);
    }*/
}

pair <double, double> baselineCalc(vector<double> trace){
    //calculating the baseline
    double baseSum = 0.;
    for(int j=0;j<20;j++){
        baseSum += trace[j];
    }
    double  baseline = baseSum/20.;
    
    //calculating the standard dev
    double stddev = 0.0;
    for(int j=0;j<20;j++){
        stddev += pow(j-baseline,2);
    }
    stddev = sqrt(stddev/20);
    
    return std::make_pair (baseline, stddev);
}
