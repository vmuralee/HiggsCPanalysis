
#include <iostream>
#include <vector>
#include <map>
#include <iomanip>


#include "TPad.h"
#include "TROOT.h"
#include "TColor.h"
#include "TEfficiency.h"
#include "TMath.h"

using namespace std;

void stack_plot(
		int nBins  =  20,
		float xmin =   40,
		float xmax = 250,
		TString Cuts ="&&pt_plus>20&&pt_neg>20&&Prompt_pT>50",
		double lumi = 41900
		)
{
  gROOT->ProcessLine(".L tdrStyle.C");
  setTDRStyle();
  gStyle->SetOptStat(0);
  TH1D * hist[6];
  double x_sec[5]={5765,5765,5765,870.31,61526};
  TString cuts[6];
  
  cuts[0]="os>0.5&&(gen_match_1==5&&gen_match_2==5)"+Cuts;
  cuts[1]="!(isZL || isZJ)"+Cuts;
  cuts[2]="os>0.5&&isZJ>0.5"+Cuts;
  cuts[3]="os>0.5"+Cuts;
  cuts[4]="os>0.5"+Cuts;
  cuts[5]="os>0.5"+Cuts;
  TFile* file_data = new TFile("DATA_Tau.root");
  TTree * tree_data = (TTree*)file_data->Get("TauCheck");
  TH1D* h_data = new TH1D("h_data","",nBins,xmin,xmax);
  tree_data->Draw("m_vis>>h_data","os>0.5&&pt_plus>20&&pt_neg>20&&Prompt_pT>50");
  TString file_names[6] = {"DYJetsToLL_M-50_TuneCP5.root","DYJetsToLL_M-50_TuneCP5.root","DYJetsToLL_M-50_TuneCP5.root","TTTo2L2Nu_TuneCP5.root","WJetsToLNu_TuneCP5.root"};
  
  for(int i =0; i<5;i++){
    TFile* file = new TFile(file_names[i]);
    TH1D * nWeightedEvents = (TH1D*)file->Get("inputEventsH"); //inputEventsH
    TTree * tree = (TTree*)file->Get("TauCheck");
    float eta_plus;float eta_neg;float pt_tmp=-10000;
    float pt_plus;float pt_neg;
    tree->SetBranchAddress("eta_plus",&eta_plus);
    tree->SetBranchAddress("eta_neg",&eta_neg);
    tree->SetBranchAddress("pt_plus",&pt_plus);
    tree->SetBranchAddress("pt_neg",&pt_neg);
    TString histName   = file_names[i];
    hist[i] = new TH1D(histName,"",nBins,xmin,xmax);
    double norm = x_sec[i]*lumi/nWeightedEvents->GetSumOfWeights();
    tree->Draw("m_vis>>"+histName,cuts[i]);
    hist[i]->Scale(norm);
  }
  //hist[0]->SetLineColor(kRed);
  hist[0]->SetFillColor(kOrange);
  hist[1]->SetFillColor(38);
  hist[2]->SetFillColor(3);
  hist[3]->SetFillColor(40);
  hist[4]->SetFillColor(28);
  TLegend *legend = new TLegend(0.40, 0.30, 1, 0.1);
  legend -> SetTextFont(42);
  //legend -> AddEntry(histData, "Observed", "ple");
  legend -> AddEntry(hist[0],"Z#rightarrow #tau#tau","f");
  legend -> AddEntry(hist[1],"Z#rightarrow #mu#mu/ee","f");
  legend -> AddEntry(hist[2],"Z+jets","f");
  legend -> AddEntry(hist[3],"t#bar{t}","f");
  legend -> AddEntry(hist[4],"W+jets","f");

 
  h_data->GetXaxis()->SetTitle("visible di-#tau mass");
  
  // Add all bkg contributions to one stack plot
  THStack *stack = new THStack("Background","");
  stack->Add(hist[4]);
  stack->Add(hist[3]);
  stack->Add(hist[2]);
  stack->Add(hist[1]);
  stack->Add(hist[0]);
  stack->Draw("HIST");
  h_data->Draw("epsame");
  legend->Draw();
}
