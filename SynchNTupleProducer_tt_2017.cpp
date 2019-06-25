#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>

#include "TFile.h" 
#include "TH1.h" 
#include "TH2.h"
#include "TGraph.h"
#include "TTree.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRFIOFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TChain.h"
#include "TMath.h"

#include "TLorentzVector.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "TRandom.h"

#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"
#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h"
#include "DesyTauAnalyses/NTupleMaker/interface/json.h"
#include "DesyTauAnalyses/NTupleMaker/interface/PileUp.h"
#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"
#include "DesyTauAnalyses/NTupleMaker/interface/functions.h"
#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"
#include "HTT-utilities/RecoilCorrections/interface/MEtSys.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Jets.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
#include "RooWorkspace.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LV;
struct DiTauInfo 
{ 
  DiTauInfo(){}; 
  int diTauCharge_; 
  double sumPt_; 
  double sumIso_;
  int index1_;
  int index2_;
}; 

struct SortDiTauPairs 
{ 
  bool operator() (const DiTauInfo t1, const DiTauInfo t2) 
  { 
   
    if(t1.sumIso_ > t2.sumIso_ ) return true;
    if(t1.sumIso_ < t2.sumIso_ ) return false;
   
    return (t1.sumPt_ > t2.sumPt_);  
  } 
}; 
int main(int argc, char * argv[]) {

  using namespace std;
  //HLT trigger 
  vector<string> filterDiTau;
  // **** configuration
  Config cfg(argv[1]);
  filterDiTau.push_back(cfg.get<string>("filterDiTau1"));
  filterDiTau.push_back(cfg.get<string>("filterDiTau2"));
  filterDiTau.push_back(cfg.get<string>("filterDiTau3"));

  const bool isData = cfg.get<bool>("IsData");
  const bool isDY   = cfg.get<bool>("IsDY");
  const bool isW    = cfg.get<bool>("IsW");
  const bool isTOP = cfg.get<bool>("IsTOP");
  
  const bool applyGoodRunSelection = cfg.get<bool>("ApplyGoodRunSelection");
  const string jsonFile = cfg.get<string>("JSON");
  const string DataPUFile = cfg.get<string>("DataPUFile");
  const string MCPUFile = cfg.get<string>("MCPUFile");
  const string pileUpforMC = cfg.get<string>("pileUpforMC");

  const bool applyInclusiveSelection = cfg.get<bool>("ApplyInclusiveSelection");
  const bool computeDitauMass        = cfg.get<bool>("ComputeDitauMass");
  const string zMassPtWeightsFileName   = cfg.get<string>("ZMassPtWeightsFileName");
  TString ZMassPtWeightsFileName(zMassPtWeightsFileName);
  
  const string zMassPtWeightsHistName   = cfg.get<string>("ZMassPtWeightsHistName");
  TString ZMassPtWeightsHistName(zMassPtWeightsHistName);

  // kinematic cuts on taus
  const float ptTauCut        = cfg.get<float>("ptTauCut");
  const float etaTauCut       = cfg.get<float>("etaTauCut");
  const float dzTauCut        = cfg.get<float>("dzTauCut");
  const bool ApplyTauId = cfg.get<bool>("ApplyTauId");
  const string TauLegXTriggerFile = cfg.get<string>("TauLegXTriggerFile");

  const float dRleptonsCut = cfg.get<float>("dRleptonsCut");
  const float deltaRTrigMatch = cfg.get<float>("DRTrigMatch");
  // jets
  const string bTagDiscriminator = cfg.get<string>("BTagDiscriminator");
  const float jetEtaCut = cfg.get<float>("JetEtaCut");
  const float jetPtLowCut = cfg.get<float>("JetPtLowCut");
  const float jetPtHighCut = cfg.get<float>("JetPtHighCut");
  const float dRJetLeptonCut = cfg.get<float>("dRJetLeptonCut");
  const float bJetEtaCut = cfg.get<float>("bJetEtaCut");
  const float btagCut = cfg.get<float>("btagCut");
  const bool applyJetPfId = cfg.get<bool>("ApplyJetPfId");
  const bool applyJetPuId = cfg.get<bool>("ApplyJetPuId");
  const string recoilFileName   = cfg.get<string>("RecoilFileName");
  TString RecoilFileName(recoilFileName);
  // veto electrons
  const float ptVetoElectronCut   = cfg.get<float>("ptVetoElectronCut");
  const float etaVetoElectronCut  = cfg.get<float>("etaVetoElectronCut");
  const float dxyVetoElectronCut  = cfg.get<float>("dxyVetoElectronCut");
  const float dzVetoElectronCut   = cfg.get<float>("dzVetoElectronCut");
  const float isoVetoElectronCut  = cfg.get<float>("isoVetoElectronCut");
  const bool applyVetoElectronId = cfg.get<bool>("ApplyVetoElectronId");

  // veto muons
  const float ptVetoMuonCut   = cfg.get<float>("ptVetoMuonCut");
  const float etaVetoMuonCut  = cfg.get<float>("etaVetoMuonCut");
  const float dxyVetoMuonCut   = cfg.get<float>("dxyVetoMuonCut");
  const float dzVetoMuonCut   = cfg.get<float>("dzVetoMuonCut");
  const float isoVetoMuonCut   = cfg.get<float>("isoVetoMuonCut");
  const bool applyVetoMuonId = cfg.get<bool>("ApplyVetoMuonId");
  const bool isIsoR03 = cfg.get<bool>("IsIsoR03");
  // 
  const float ptDilepMuonCut = cfg.get<float>("ptDilepMuonCut");
  const float etaDilepMuonCut = cfg.get<float>("etaDilepMuonCut");
  const float dxyDilepMuonCut = cfg.get<float>("dxyDilepMuonCut");
  const float dzDilepMuonCut = cfg.get<float>("dzDilepMuonCut");
  const float isoDilepMuonCut = cfg.get<float>("isoDilepMuonCut");
  const float dRDilepVetoCut = cfg.get<float>("dRDilepVetoCut");

  std::string rootFileName(argv[2]);
  std::ifstream fileList(argv[2]);
  std::ifstream fileList0(argv[2]);
  std::string ntupleName("makeroottree/AC1B");
  std::string initNtupleName("initroottree/AC1B");
  TString BTagDiscriminator(bTagDiscriminator);
  TString TStrName(rootFileName);
  std::cout <<TStrName <<std::endl; 
  // output fileName with histograms
  TFile * file = new TFile(TStrName+TString(".root"),"recreate");
  file->cd("");

  TH1D * inputEventsH = new TH1D("inputEventsH","",1,-0.5,0.5);
  TH1D * histWeightsH = new TH1D("histWeightsH","",1,-0.5,0.5);
  TH1D * ZEE = new TH1D("ZEE","",1,-0.5,0.5);
  TH1D * ZMM = new TH1D("ZMM","",1,-0.5,0.5);
  TH1D * ZTT = new TH1D("ZTT","",1,-0.5,0.5);
  TH1D * ALL = new TH1D("ALL","",1,-0.5,0.5);
  TH1D * GAM = new TH1D("GAM","",1,-0.5,0.5);

  TTree * tree = new TTree("TauCheck","TauCheck");
  // Declaration of leaf types
  Int_t           run;
  Int_t           lumi;
  Int_t           evt;
  Int_t           npv;
  Int_t           npu;
  Float_t rho;
  Float_t         mcweight;
  Float_t         puweight;
  Float_t         trigweight_1;
  Float_t         trigweight_2;
  Float_t         idweight_1;
  Float_t         idweight_2;
  Float_t         isoweight_1;
  Float_t         isoweight_2;
  Float_t         effweight;
  Float_t         fakeweight;
  Float_t         embeddedWeight;
  Float_t         signalWeight;
  Float_t weight;
  Float_t zptmassweight;
  Float_t         nuPx;
  Float_t         nuPy;
  Float_t nuPz;
  Float_t         lepPx;
  Float_t         lepPy;
  Float_t lepPz;
  Float_t         bosonPx;
  Float_t         bosonPy;
  Float_t         bosonPz;
  Float_t         bosonPt;
  Float_t bosonMass;
  Float_t         m_vis;
  Float_t         tauPair_m_vis;
  Bool_t isZLL;
  Bool_t isZMM;
  Bool_t isZEE;
  Bool_t isZTT;
  Bool_t isZL;
  Bool_t isZJ;
  Bool_t          extraelec_veto;
  Bool_t extramuon_veto;
  Bool_t os;

  Double_t Prompt_pT;
  Float_t         againstMuonTight3_2;
  Float_t againstElectronVTightMVA6_2;
  Float_t         againstMuonTight3_1;
  Float_t againstElectronVTightMVA6_1;
  Float_t         pt_tau1;
  Float_t         phi_tau1;
  Float_t         eta_tau1;
  Float_t         m_tau1;

  Float_t         pt_plus;
  Float_t         phi_plus;
  Float_t         eta_plus;
  Float_t         m_plus;
  Int_t           q_plus;
  Float_t         iso_plus;
  Float_t         mva_plus;
  Float_t         d0_plus;
  Float_t         dZ_plus;

  Float_t mt_plus;
  Float_t         pt_neg;
  Float_t         phi_neg;
  Float_t         eta_neg;
  Float_t         m_neg;
  Int_t           q_neg;
  Float_t         iso_neg;
  Float_t         mva_neg;
  Float_t         d0_neg;
  Float_t         dZ_neg;

  Float_t         pt_tau2;
  Float_t         phi_tau2;
  Float_t         eta_tau2;
  Float_t         m_tau2;

  Float_t mt_neg;

  Int_t           njets;
  Int_t           njetspt20;
  Float_t         jpt_1;
  Float_t         jeta_1;
  Float_t         jphi_1;
  Float_t         jptraw_1;
  Float_t         jptunc_1;
  Float_t         jmva_1;
  Float_t         jlrm_1;
  Int_t           jctm_1;
  Float_t         jpt_2;
  Float_t         jeta_2;
  Float_t         jphi_2;
  Float_t         jptraw_2;
  Float_t         jptunc_2;
  Float_t         jmva_2;
  Float_t         jlrm_2;
  Int_t           jctm_2;
  Float_t         mjj;
  Float_t         jdeta;
  Int_t njetingap;

  Float_t         met;
  Float_t         metphi;
  Float_t         metcov00;
  Float_t         metcov01;
  Float_t         metcov10;
  Float_t         metcov11;
  Float_t         genmet;
  Float_t genmetphi;
  Int_t gen_match_1; 
  Int_t gen_match_2; 
  // Declaration of branches
  tree->Branch("run", &run, "run/I");
  tree->Branch("lumi", &lumi, "lumi/I");
  tree->Branch("evt", &evt, "evt/I");
  tree->Branch("npv", &npv, "npv/I");
  tree->Branch("npu", &npu, "npu/I");
  tree->Branch("rho", &rho, "rho/F");
  tree->Branch("mcweight", &mcweight, "mcweight/F");
  tree->Branch("puweight", &puweight, "puweight/F");
  tree->Branch("trigweight_1", &trigweight_1, "trigweight_1/F");
  tree->Branch("trigweight_2", &trigweight_2, "trigweight_2/F");
  tree->Branch("idweight_1", &idweight_1, "idweight_1/F");
  tree->Branch("idweight_2", &idweight_2, "idweight_2/F");
  tree->Branch("isoweight_1", &isoweight_1, "isoweight_1/F");
  tree->Branch("isoweight_2", &isoweight_2, "isoweight_2/F");
  tree->Branch("effweight", &effweight, "effweight/F");
  tree->Branch("fakeweight", &fakeweight, "fakeweight/F");
  tree->Branch("embeddedWeight", &embeddedWeight, "embeddedWeight/F");
  tree->Branch("signalWeight", &signalWeight, "signalWeight/F");
  tree->Branch("weight", &weight, "weight/F");
  tree->Branch("zptmassweight",&zptmassweight,"zptmassweight/F");
  tree->Branch("nuPx",&nuPx,"nuPx/F");
  tree->Branch("nuPy",&nuPy,"nuPy/F");
  tree->Branch("nuPz",&nuPz,"nuPz/F");
  tree->Branch("m_vis",&m_vis,"m_vis/F");

  tree->Branch("tauPair_m_vis",&tauPair_m_vis,"tauPair_m_vis/F");

  tree->Branch("lepPx",&lepPx,"lepPx/F");
  tree->Branch("lepPy",&lepPy,"lepPy/F");
  tree->Branch("lepPz",&lepPz,"lepPz/F");

  tree->Branch("bosonPx",&bosonPx,"bosonPx/F");
  tree->Branch("bosonPy",&bosonPy,"bosonPy/F");
  tree->Branch("bosonPz",&bosonPz,"bosonPz/F");
  tree->Branch("bosonPt",&bosonPt,"bosonPt/F");
  tree->Branch("bosonMass",&bosonMass,"bosonMass/F");
  tree->Branch("isZLL",&isZLL,"isZLL/O");
  tree->Branch("isZEE",&isZEE,"isZEE/O");
  tree->Branch("isZMM",&isZMM,"isZMM/O");
  tree->Branch("isZTT",&isZTT,"isZTT/O");
  tree->Branch("isZL",&isZL,"isZL/O");
  tree->Branch("isZJ",&isZJ,"isZL/O");
  tree->Branch("extraelec_veto", &extraelec_veto, "extraelec_veto/O");
  tree->Branch("extramuon_veto", &extramuon_veto, "extramuon_veto/O");
  tree->Branch("os", &os, "os/O");
  tree->Branch("Prompt_pT", &Prompt_pT, "Prompt_pT/D");
  tree->Branch("againstMuonTight3_2", &againstMuonTight3_2, "againstMuonTight3_2/F");
  //tree->Branch("againstElectronVLooseMVA6_2",&againstElectronVLooseMVA6_2,"againstElectronVLooseMVA6_2/F");
  tree->Branch("againstElectronVTightMVA6_2",&againstElectronVTightMVA6_2,"againstElectronVTightMVA6_2/F");
  tree->Branch("againstMuonTight3_1", &againstMuonTight3_1, "againstMuonTight3_1/F");
  //tree->Branch("againstElectronVLooseMVA6_1",&againstElectronVLooseMVA6_1,"againstElectronVLooseMVA6_1/F");
  tree->Branch("againstElectronVTightMVA6_1",&againstElectronVTightMVA6_1,"againstElectronVTightMVA6_1/F");
  tree->Branch("pt_plus", &pt_plus, "pt_plus/F");
  tree->Branch("phi_plus", &phi_plus, "phi_plus/F");
  tree->Branch("eta_plus", &eta_plus, "eta_plus/F");
  tree->Branch("m_plus", &m_plus, "m_plus/F");

  tree->Branch("pt_tau1", &pt_tau1, "pt_tau1/F");
  tree->Branch("phi_tau1", &phi_tau1, "phi_tau1/F");
  tree->Branch("eta_tau1", &eta_tau1, "eta_tau1/F");
  tree->Branch("m_tau1", &m_tau1, "m_tau1/F");

  tree->Branch("q_plus", &q_plus, "q_plus/I");
  tree->Branch("iso_plus", &iso_plus, "iso_plus/F");
  tree->Branch("mva_plus", &mva_plus, "mva_plus/F");
  tree->Branch("d0_plus", &d0_plus, "d0_plus/F");
  tree->Branch("dZ_plus", &dZ_plus, "dZ_plus/F");
  tree->Branch("mt_plus", &mt_plus, "mt_plus/F");

  tree->Branch("pt_tau2", &pt_tau2, "pt_tau2/F");
  tree->Branch("phi_tau2", &phi_tau2, "phi_tau2/F");
  tree->Branch("eta_tau2", &eta_tau2, "eta_tau2/F");
  tree->Branch("m_tau2", &m_tau2, "m_tau2/F");

  tree->Branch("pt_neg", &pt_neg, "pt_neg/F");
  tree->Branch("phi_neg", &phi_neg, "phi_neg/F");
  tree->Branch("eta_neg", &eta_neg, "eta_neg/F");
  tree->Branch("m_neg", &m_neg, "m_neg/F");
  tree->Branch("q_neg", &q_neg, "q_neg/I");
  tree->Branch("iso_neg", &iso_neg, "iso_neg/F");
  tree->Branch("mva_neg", &mva_neg, "mva_neg/F");
  tree->Branch("d0_neg", &d0_neg, "d0_neg/F");
  tree->Branch("dZ_neg", &dZ_neg, "dZ_neg/F");
  tree->Branch("mt_neg", &mt_neg, "mt_neg/F");

  tree->Branch("njets", &njets, "njets/I");
  tree->Branch("njetspt20", &njetspt20, "njetspt20/I");

  tree->Branch("jpt_1", &jpt_1, "jpt_1/F");
  tree->Branch("jeta_1", &jeta_1, "jeta_1/F");
  tree->Branch("jphi_1", &jphi_1, "jphi_1/F");
  tree->Branch("jptraw_1", &jptraw_1, "jptraw_1/F");
  tree->Branch("jptunc_1", &jptunc_1, "jptunc_1/F");
  tree->Branch("jmva_1", &jmva_1, "jmva_1/F");
  tree->Branch("jlrm_1", &jlrm_1, "jlrm_1/F");
  tree->Branch("jctm_1", &jctm_1, "jctm_1/I");

  tree->Branch("jpt_2", &jpt_2, "jpt_2/F");
  tree->Branch("jeta_2", &jeta_2, "jeta_2/F");
  tree->Branch("jphi_2", &jphi_2, "jphi_2/F");
  tree->Branch("jptraw_2", &jptraw_2, "jptraw_2/F");
  tree->Branch("jptunc_2", &jptunc_2, "jptunc_2/F");
  tree->Branch("jmva_2", &jmva_2, "jlrm_2/F");
  tree->Branch("jlrm_2", &jlrm_2, "jlrm_2/F");
  tree->Branch("jctm_2", &jctm_2, "jctm_2/I");

  tree->Branch("mjj", &mjj, "mjj/F");
  tree->Branch("jdeta", &jdeta, "jdeta/F");
  tree->Branch("njetingap", &njetingap, "njetingap/I");
 

  tree->Branch("met", &met, "met/F");
  tree->Branch("metphi", &metphi, "metphi/F");
  tree->Branch("metcov00", &metcov00, "metcov00/F");
  tree->Branch("metcov01", &metcov01, "metcov01/F");
  tree->Branch("metcov10", &metcov10, "metcov10/F");
  tree->Branch("metcov11", &metcov11, "metcov11/F");

  tree->Branch("genmet",&genmet,"genmet/F");
  tree->Branch("genmetphi",&genmetphi,"genmetphi/F");
  tree->Branch("gen_match_1",&gen_match_1,"gen_match_1/I");
  tree->Branch("gen_match_2",&gen_match_2,"gen_match_2/I");

 // end of Sych tree
  int nTotalFiles = 0;
  std::string dummy;
  // count number of files --->
  while (fileList0 >> dummy) nTotalFiles++;
  std::vector<Period> periods;
  string cmsswBase = (getenv ("CMSSW_BASE"));
  string fullPathToJsonFile = cmsswBase + "/src/DesyTauAnalyses/NTupleMaker/test/json/" + jsonFile;
  if (isData) {
    std::fstream inputFileStream(fullPathToJsonFile.c_str(), std::ios::in);
    if (inputFileStream.fail()) {
      std::cout << "Error: cannot find json file " << fullPathToJsonFile << std::endl;
      std::cout << "please check" << std::endl;
      std::cout << "quitting program" << std::endl;
      exit(-1);
    }
    for(std::string s; std::getline(inputFileStream, s); )
      {
	periods.push_back(Period());
	std::stringstream ss(s);
	ss >> periods.back();
      }
  }
  
  // Official PU reweighting
  PileUp * PUofficial = new PileUp();
  TFile * filePUdistribution_data = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/"+TString(DataPUFile),"read");
  TFile * filePUdistribution_MC = new TFile (TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/"+TString(MCPUFile),"read");
  TH1D * PU_data = (TH1D *)filePUdistribution_data->Get("pileup");
  if(isTOP){
    TH1D * PU_mc = (TH1D *)filePUdistribution_MC->Get("TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8_pileup");
    PUofficial->set_h_MC(PU_mc);  
  }
  else if(isDY){
    TH1D * PU_mc = (TH1D *)filePUdistribution_MC->Get("DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_pileup");
    PUofficial->set_h_MC(PU_mc);
  }
  else if(isW){
    TH1D * PU_mc = (TH1D *)filePUdistribution_MC->Get("WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_pileup");
    PUofficial->set_h_MC(PU_mc);
  }
  PUofficial->set_h_data(PU_data);
  
  // Tau Scale factors
  TString filename_XTrigTauLegSF = TString(cmsswBase)+"/src/"+TString(TauLegXTriggerFile);
  TFile f_XTrigTauLegSF(filename_XTrigTauLegSF);
  RooWorkspace *w_XTrigTauLegSF = (RooWorkspace*)f_XTrigTauLegSF.Get("w");
  // Recoil corrector and met systematics
  RecoilCorrector recoilMetCorrector(RecoilFileName);
  // SV fit mass
  edm::FileInPath inputFileName_visPtResolution("TauAnalysis/SVfitStandalone/data/svFitVisMassAndPtResolutionPDF.root");
  TH1::AddDirectory(false);  
  TFile * inputFile_visPtResolution = new TFile(inputFileName_visPtResolution.fullPath().data());
  //Z pt mass reweighting
  TFile * fileZMassPtWeights = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/"+ZMassPtWeightsFileName);
  if (fileZMassPtWeights->IsZombie()) {
    std::cout << "File " << TString(cmsswBase) << "/src/DesyTauAnalyses/NTupleMaker/data/" << ZMassPtWeightsFileName << "  does not exist!" << std::endl;
    exit(-1);
  }
  TH2D * histZMassPtWeights = (TH2D*)fileZMassPtWeights->Get(ZMassPtWeightsHistName);
  if (histZMassPtWeights==NULL) {
    std::cout << "histogram " << ZMassPtWeightsHistName << " is not found in file " << TString(cmsswBase) << "/src/DesyTauAnalyses/NTupleMaker/data/" << ZMassPtWeightsFileName
	      << std::endl;
    exit(-1);
  }
  // BTag scale factors
  BTagCalibration calib("csvv2", cmsswBase+"/src/DesyTauAnalyses/NTupleMaker/data/CSVv2.csv");
  BTagCalibrationReader reader(BTagEntry::OP_MEDIUM,  // operating point
			       "central");            // systematics type
  reader.load(calib,               // calibration instance
	      BTagEntry::FLAV_B,    // btag flavour
	      "mujets");            // measurement type
  
  reader.load(calib,               // calibration instance
	      BTagEntry::FLAV_C,    // btag flavour
	      "mujets");             // measurement type

  reader.load(calib,               // calibration instance
	      BTagEntry::FLAV_UDSG,    // btag flavour
	      "incl"); // measurement type

  TFile * fileTagging = new TFile(TString(cmsswBase)+TString("/src/DesyTauAnalyses/NTupleMaker/data/tagging_efficiencies.root"));
  TH1F * tagEff_B = (TH1F*)fileTagging->Get("btag_eff_b");
  TH1F * tagEff_C = (TH1F*)fileTagging->Get("btag_eff_c");
  TH1F * tagEff_Light = (TH1F*)fileTagging->Get("btag_eff_oth");
  TRandom3 rand;

  float MaxBJetPt = 670.;
  float MaxLJetPt = 1000.;
  float MinLJetPt = 20.;
  float MinBJetPt = 30.;

  int nEvents = 0;
  int selEvents = 0;
  int nFiles = 0;
  // Looping over RooT files 
  for (int iF=0; iF<nTotalFiles; ++iF) {

    std::string filen;
    fileList >> filen;

    std::cout << "file " << iF+1 << " out of " << nTotalFiles << " filename : " << filen << std::endl;
    TFile * file_ = TFile::Open(TString(filen));
    TTree * _tree = NULL;
    _tree = (TTree*)file_->Get(TString(ntupleName));
    if (_tree==NULL) continue;
    TH1D * histoInputEvents = NULL;
   
    histoInputEvents = (TH1D*)file_->Get("makeroottree/nEvents");
    
    if (histoInputEvents==NULL) continue;
    
    int NE = int(histoInputEvents->GetEntries());
    for (int iE=0;iE<NE;++iE)
      inputEventsH->Fill(0.);

    AC1B analysisTree(_tree);
    Long64_t numberOfEntries = analysisTree.GetEntries();
    std::cout << "      number of entries in Tree = " << numberOfEntries << std::endl;
    for (Long64_t iEntry=0; iEntry<numberOfEntries; iEntry++) { 
      analysisTree.GetEntry(iEntry);
      nEvents++;
      
      if (nEvents%10000==0) 
	cout << "      processed " << nEvents << " events" << endl; 
      if (isData && applyGoodRunSelection) {
	bool lumi = false;
        int n=analysisTree.event_run;
        int lum = analysisTree.event_luminosityblock;

	std::string num = std::to_string(n);
	std::string lnum = std::to_string(lum);
	for(const auto& a : periods)
	  {
	    if ( num.c_str() ==  a.name ) {
              for(auto b = a.ranges.begin(); b != std::prev(a.ranges.end()); ++b) {
                if (lum  >= b->lower && lum <= b->bigger ) lumi = true;
              }
              auto last = std::prev(a.ranges.end());
	      if ( (lum >=last->lower && lum <= last->bigger )) lumi=true;
	    }
	  }
	if (!lumi) continue;	
      }
      vector<int> nDiTauTrig(filterDiTau.size(),-1);
      vector<bool> checkFilterDiTauTrig(filterDiTau.size(), false);
      unsigned int nfilters = analysisTree.run_hltfilters->size();
      for(unsigned int i=0; i<nfilters; ++i){
	TString HLTFilter(analysisTree.run_hltfilters->at(i));
	for(unsigned int i_trig=0; i_trig<filterDiTau.size(); i_trig++){
	  if (HLTFilter==filterDiTau.at(i_trig)){ nDiTauTrig.at(i_trig) = i; checkFilterDiTauTrig.at(i_trig) = true;}
	}
      }
      isZLL = false;
      isZEE = false;
      isZMM = false;
      isZTT = false;
      bool isZfound = false;
      bool isWfound = false;
      bool isHfound = false;

      nuPx = 0;
      nuPy = 0;
      nuPz = 0;

      lepPx = 0;
      lepPy = 0;
      lepPz = 0;

      bosonPx = 0;
      bosonPy = 0;
      bosonPz = 0;

      std::vector<TLorentzVector> promptTausFirstCopy; promptTausFirstCopy.clear();
      std::vector<TLorentzVector> promptTausLastCopy;  promptTausLastCopy.clear();
      std::vector<TLorentzVector> promptElectrons; promptElectrons.clear();
      std::vector<TLorentzVector> promptMuons; promptMuons.clear();
      std::vector<TLorentzVector> promptNeutrinos; promptNeutrinos.clear();
      std::vector<TLorentzVector> nonpromptNeutrinos; nonpromptNeutrinos.clear();
      TLorentzVector promptTausLV; promptTausLV.SetXYZT(0,0,0,0);
      TLorentzVector promptVisTausLV; promptVisTausLV.SetXYZT(0,0,0,0);
      TLorentzVector zBosonLV; zBosonLV.SetXYZT(0,0,0,0);
      TLorentzVector wBosonLV; wBosonLV.SetXYZT(0,0,0,0);
      TLorentzVector hBosonLV; hBosonLV.SetXYZT(0,0,0,0);
      TLorentzVector promptElectronsLV; promptElectronsLV.SetXYZT(0,0,0,0);
      TLorentzVector promptMuonsLV; promptMuonsLV.SetXYZT(0,0,0,0);
      TLorentzVector promptNeutrinosLV;  promptNeutrinosLV.SetXYZT(0,0,0,0);
      TLorentzVector nonpromptNeutrinosLV;  nonpromptNeutrinosLV.SetXYZT(0,0,0,0);
      TLorentzVector wDecayProductsLV; wDecayProductsLV.SetXYZT(0,0,0,0);
      TLorentzVector fullVLV; fullVLV.SetXYZT(0,0,0,0);
      TLorentzVector visVLV; visVLV.SetXYZT(0,0,0,0);

      TLorentzVector genBosonLV; genBosonLV.SetXYZT(0,0,0,0);
      TLorentzVector genVisBosonLV; genVisBosonLV.SetXYZT(0,0,0,0);
      if (!isData) {
	for (unsigned int igentau=0; igentau<analysisTree.gentau_count; ++igentau) {
	  TLorentzVector tauLV; tauLV.SetXYZT(analysisTree.gentau_px[igentau],
					      analysisTree.gentau_py[igentau],
					      analysisTree.gentau_pz[igentau],
					      analysisTree.gentau_e[igentau]);
	  TLorentzVector tauVisLV; tauVisLV.SetXYZT(analysisTree.gentau_visible_px[igentau],
						    analysisTree.gentau_visible_py[igentau],
						    analysisTree.gentau_visible_pz[igentau],
						    analysisTree.gentau_visible_e[igentau]);
	  if (analysisTree.gentau_isPrompt[igentau]&&analysisTree.gentau_isFirstCopy[igentau]) {
	    promptTausFirstCopy.push_back(tauLV);
	    promptTausLV += tauLV;
	    wDecayProductsLV += tauLV;
	  }
	  if (analysisTree.gentau_isPrompt[igentau]&&analysisTree.gentau_isLastCopy[igentau]) {
	    promptTausLastCopy.push_back(tauVisLV);
	    promptVisTausLV += tauVisLV;
	  }

	  for (unsigned int igen=0; igen < analysisTree.genparticles_count; ++igen) {
	    TLorentzVector genLV; genLV.SetXYZT(analysisTree.genparticles_px[igen],
						analysisTree.genparticles_py[igen],
						analysisTree.genparticles_pz[igen],
						analysisTree.genparticles_e[igen]);
	    if (analysisTree.genparticles_pdgid[igen]==23) { 
	      isZfound = true;
	      zBosonLV = genLV;
	    }
	    if (abs(analysisTree.genparticles_pdgid[igen])==24) { 
	      isWfound = true;
	      wBosonLV = genLV;
	    }
	    if (abs(analysisTree.genparticles_pdgid[igen])==25||
		abs(analysisTree.genparticles_pdgid[igen])==35||
		abs(analysisTree.genparticles_pdgid[igen])==36) {
	      isHfound = true;
	      hBosonLV = genLV;
	    }

	    bool fromHardProcessFinalState = analysisTree.genparticles_fromHardProcess[igen]&&analysisTree.genparticles_status[igen]==1;
	    bool isMuon = false;
	    bool isElectron = false;
	    bool isNeutrino = false;
	    bool isDirectHardProcessTauDecayProduct = analysisTree.genparticles_isDirectHardProcessTauDecayProduct[igen];

	      
	    if (fabs(analysisTree.genparticles_pdgid[igen])==11) { 
	      isElectron = true;
	      if (analysisTree.genparticles_isPrompt[igen]&&analysisTree.genparticles_status[igen]==1) {
		promptElectrons.push_back(genLV);
		promptElectronsLV += genLV;
		wDecayProductsLV += genLV;
	      }
	    }
	      
	    if (fabs(analysisTree.genparticles_pdgid[igen])==13) { 
	      isMuon = true;
	      if (analysisTree.genparticles_isPrompt[igen]&&analysisTree.genparticles_status[igen]==1) {
		promptMuons.push_back(genLV);
		promptMuonsLV += genLV;
		wDecayProductsLV += genLV;
	      }
	    }
	    if (fabs(analysisTree.genparticles_pdgid[igen])==12||
		fabs(analysisTree.genparticles_pdgid[igen])==14||
		fabs(analysisTree.genparticles_pdgid[igen])==16)  {
	      isNeutrino = true;
	      if (analysisTree.genparticles_isPrompt[igen]&&
		  !analysisTree.genparticles_isDirectHardProcessTauDecayProduct[igen]&&
		  analysisTree.genparticles_status[igen]==1) {
		promptNeutrinos.push_back(genLV);
		promptNeutrinosLV += genLV;
		wDecayProductsLV += genLV;
	      }
	      if (analysisTree.genparticles_isDirectHardProcessTauDecayProduct[igen]&&
		  analysisTree.genparticles_status[igen]==1) {
		nonpromptNeutrinos.push_back(genLV);
		nonpromptNeutrinosLV += genLV;
	      }
	    }

	    bool isBoson = (fromHardProcessFinalState && (isMuon || isElectron || isNeutrino)) || isDirectHardProcessTauDecayProduct;
	    bool isVisibleBoson = (fromHardProcessFinalState && (isMuon || isElectron)) || (isDirectHardProcessTauDecayProduct && !isNeutrino);
	      
	    if (isBoson)
	      genBosonLV += genLV;
	    if (isVisibleBoson)
	      genVisBosonLV += genLV;
	    
	  }
	}
	
	
      }//isData
      if (isDY) {
	if (promptTausFirstCopy.size()==2) {
	  bosonPx = promptTausLV.Px(); bosonPy = promptTausLV.Py(); bosonPz = promptTausLV.Pz();
	  lepPx = promptVisTausLV.Px(); lepPy = promptVisTausLV.Py(); lepPz = promptVisTausLV.Pz();
	}
	else if (promptMuons.size()==2) {
	  bosonPx = promptMuonsLV.Px(); bosonPy = promptMuonsLV.Py(); bosonPz = promptMuonsLV.Pz();
	  lepPx = promptMuonsLV.Px(); lepPy = promptMuonsLV.Py(); lepPz = promptMuonsLV.Pz();
	}
	else {
	  bosonPx = promptElectronsLV.Px(); bosonPy = promptElectronsLV.Py(); bosonPz = promptElectronsLV.Pz();
	  lepPx = promptElectronsLV.Px(); lepPy = promptElectronsLV.Py(); lepPz = promptElectronsLV.Pz();
	}
      }
      else if (isW) {
	if (isWfound) {
	  bosonPx = wBosonLV.Px(); bosonPy = wBosonLV.Py(); bosonPz = wBosonLV.Pz();
	}
	else {
	  bosonPx = wDecayProductsLV.Px(); bosonPy = wDecayProductsLV.Py(); bosonPz = wDecayProductsLV.Pz();
	}
	if (promptTausLastCopy.size()==1) { 
	  lepPx = promptVisTausLV.Px(); lepPy = promptVisTausLV.Py(); lepPz = promptVisTausLV.Pz();
	}
	else if (promptMuons.size()==1) { 
	  lepPx = promptMuonsLV.Px(); lepPy = promptMuonsLV.Py(); lepPz = promptMuonsLV.Pz();
	}
	else { 
	  lepPx = promptElectronsLV.Px(); lepPy = promptElectronsLV.Py(); lepPz = promptElectronsLV.Pz();
	}
      }
      else {
	TLorentzVector bosonLV = promptTausLV + promptMuonsLV + promptElectronsLV + promptNeutrinosLV;
	bosonPx = bosonLV.Px(); bosonPy = bosonLV.Py(); bosonPz = bosonLV.Pz();
	TLorentzVector lepLV = promptVisTausLV + promptMuonsLV + promptElectronsLV;
	lepPx = lepLV.Px(); lepPy = lepLV.Py(); lepPz = lepLV.Pz();
      }

      bosonPx = genBosonLV.Px();
      bosonPy = genBosonLV.Py();
      bosonPz = genBosonLV.Pz();
      bosonPt = genBosonLV.Pt();
      bosonMass = genBosonLV.M();
      
      lepPx = genVisBosonLV.Px();
      lepPy = genVisBosonLV.Py();
      lepPz = genVisBosonLV.Pz();
      
      isZLL = isZEE || isZMM;
      ALL->Fill(0.);

      zptmassweight = 1;
      if (isDY) { // applying Z pt mass weights
	if (bosonMass>50.0) {
	  float bosonMassX = bosonMass;
	  float bosonPtX = bosonPt;
	  if (bosonMassX>1000.) bosonMassX = 1000.;
	  if (bosonPtX<1.)      bosonPtX = 1.;
	  if (bosonPtX>1000.)   bosonPtX = 1000.;
	  zptmassweight = histZMassPtWeights->GetBinContent(histZMassPtWeights->GetXaxis()->FindBin(bosonMassX),
							    histZMassPtWeights->GetYaxis()->FindBin(bosonPtX));
	}
      }

	
      if (isZEE) ZEE->Fill(0.);
      if (isZMM) ZMM->Fill(0.);
      if (isZTT) ZTT->Fill(0.);
      if (!isZfound) GAM->Fill(0.);
      run = int(analysisTree.event_run);
      lumi = int(analysisTree.event_luminosityblock);
      evt = int(analysisTree.event_nr);

      // weights
      mcweight = analysisTree.genweight;
      puweight = 1;
      trigweight_1 = 1;
      trigweight_2 = 1;
      idweight_1 = 1;
      idweight_2 = 1;
      isoweight_1 = 1;
      isoweight_2 = 1;
      effweight = 1;
      fakeweight = 1;
      embeddedWeight = 1;
      signalWeight = 1;
      weight = 1;
         
      npv = analysisTree.primvertex_count;
      npu = analysisTree.numtruepileupinteractions;
      rho = analysisTree.rho;
      if (!isData)
	puweight = float(PUofficial->get_PUweight(double(analysisTree.numtruepileupinteractions)));
      unsigned int nBTagDiscriminant = 0;
      for (unsigned int iBTag=0; iBTag < analysisTree.run_btagdiscriminators->size(); ++iBTag) {
	TString discr(analysisTree.run_btagdiscriminators->at(iBTag));
	if (discr.Contains(BTagDiscriminator))
	  nBTagDiscriminant = iBTag;
      }
      std::vector<DiTauInfo> sortDiTauInfos;  sortDiTauInfos.clear();
      // tau selection

      vector<int> taus; taus.clear();
      for (unsigned int it = 0; it<analysisTree.tau_count; ++it) {
	bool validDecayMode = analysisTree.tau_decayMode[it]<=4 || analysisTree.tau_decayMode[it]>=10; // excluding 2prongs 
	if (!validDecayMode) continue;
	if (analysisTree.tau_decayModeFinding[it]<=0.5) continue;
	if (analysisTree.tau_byVLooseIsolationMVArun2017v2DBoldDMwLT2017[it] < 0.5) continue;//*applied tight criteria after selecting index pair
	//if(analysisTree.tau_againstMuonTight3[it]<0.5)continue;
	//if(analysisTree.tau_againstElectronVTightMVA6[it]<0.5)continue;
	if (analysisTree.tau_pt[it]<ptTauCut) continue;
	if (fabs(analysisTree.tau_eta[it])>etaTauCut) continue;
	if (fabs(analysisTree.tau_leadchargedhadrcand_dz[it])>dzTauCut) continue;
	if (fabs(analysisTree.tau_charge[it])<0.5||
	    fabs(analysisTree.tau_charge[it])>1.5) continue;
	taus.push_back(it);
      }
      if (taus.size()<2) continue;
      int tauIndex_1 = -1;
      int tauIndex_2 = -1;
      bool isTauTrig;
      bool isTau2Trig;
      float isoTauMin = 50;
      float ptTau = 0;
      float pt_max=0;
      vector<int>tau_first;tau_first.clear();
      vector<int>tau_second;tau_second.clear();
      for (unsigned int it=0; it<taus.size(); ++it) {
     
	for (unsigned int itn=it+1; itn<taus.size(); ++itn) {
	 
	  unsigned int tIndex1 = taus.at(it);
	  unsigned int tIndex2 = taus.at(itn);
	  
	  float dR = deltaR(analysisTree.tau_eta[tIndex1],analysisTree.tau_phi[tIndex1],
			    analysisTree.tau_eta[tIndex2],analysisTree.tau_phi[tIndex2]);
	  
	  if (dR<dRleptonsCut) continue;
	  //make isTrigger to false
	  isTauTrig = false;
	  isTau2Trig = false;
	  for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	    float dRtrigTau1 = deltaR(analysisTree.tau_eta[tIndex1], analysisTree.tau_phi[tIndex1], 
				      analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	    float dRtrigTau2 = deltaR(analysisTree.tau_eta[tIndex2], analysisTree.tau_phi[tIndex2], 
				      analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	  
	   
	    if(dRtrigTau1<deltaRTrigMatch){
	      
	      for(unsigned int i_trig = 0; i_trig<filterDiTau.size(); i_trig++){
		if (nDiTauTrig.at(i_trig) == -1) continue;
		if (analysisTree.trigobject_filters[iT][nDiTauTrig.at(i_trig)]) isTauTrig = true;
	      }
	    }
	    if(dRtrigTau2<deltaRTrigMatch){
	      for(unsigned int i_trig = 0; i_trig<filterDiTau.size(); i_trig++){
		//cout<<nDiTauTrig.at(i_trig)<<endl;
		if (nDiTauTrig.at(i_trig) == -1) continue;
		
		if (analysisTree.trigobject_filters[iT][nDiTauTrig.at(i_trig)]) isTau2Trig = true;
		
	      }
	    }

	  }//trigger obj match
   	  bool trigMatch =isTauTrig && isTau2Trig;
	  if(!trigMatch) continue;
	  
	  tau_first.push_back(tIndex1);
	  tau_second.push_back(tIndex2);
	  float sumPt = analysisTree.tau_pt[tIndex1]+analysisTree.tau_pt[tIndex2];
	  int pairCharge = analysisTree.tau_charge[tIndex1]*analysisTree.tau_charge[tIndex2];
	  float isoTau1 = analysisTree.tau_byIsolationMVArun2v1DBoldDMwLTraw[tIndex1];
	  float isoTau2 = analysisTree.tau_byIsolationMVArun2v1DBoldDMwLTraw[tIndex2];
	  float sumIso = isoTau1+isoTau2;
	  DiTauInfo sortDiTauInfo; 
	  sortDiTauInfo.index1_ = tIndex1;
	  sortDiTauInfo.index2_ = tIndex2;
	  sortDiTauInfo.sumPt_ = sumPt; 
	  sortDiTauInfo.sumIso_ = sumIso;
	  sortDiTauInfo.diTauCharge_ = pairCharge; 
	  sortDiTauInfos.push_back(sortDiTauInfo);
	  
	}//tau-
      }//tau+
      //Selecting isolated minimum tau pair
      std::sort(sortDiTauInfos.begin(), sortDiTauInfos.end(), SortDiTauPairs()); 
      int diTauCounter = -1;  
      for(std::vector<DiTauInfo>::iterator iter = sortDiTauInfos.begin(); iter != sortDiTauInfos.end() ; iter++){ 
	if(diTauCounter >= 0) continue;
    
	tauIndex_1 = iter->index1_;
	tauIndex_2 = iter->index2_;
	LV temp_Leg1_(analysisTree.tau_px[tauIndex_1], analysisTree.tau_py[tauIndex_1], analysisTree.tau_pz[tauIndex_1], analysisTree.tau_e[tauIndex_1]);
	LV temp_Leg2_(analysisTree.tau_px[tauIndex_2], analysisTree.tau_py[tauIndex_2], analysisTree.tau_pz[tauIndex_2], analysisTree.tau_e[tauIndex_2]);
	LV Leg1P4_, Leg2P4_;
	if(analysisTree.tau_byIsolationMVArun2v1DBoldDMwLTraw[tauIndex_1] > analysisTree.tau_byIsolationMVArun2v1DBoldDMwLTraw[tauIndex_2]){
	  Leg1P4_ = temp_Leg1_; 
	  Leg2P4_ = temp_Leg2_;
	}
	else if(analysisTree.tau_byIsolationMVArun2v1DBoldDMwLTraw[tauIndex_1] < analysisTree.tau_byIsolationMVArun2v1DBoldDMwLTraw[tauIndex_2]){
	  Leg1P4_ = temp_Leg2_;tauIndex_1 = iter->index2_;
	  Leg2P4_ = temp_Leg1_;tauIndex_2 = iter->index1_;
	}
	else if(temp_Leg1_.pt() > temp_Leg2_.pt()){
	  Leg1P4_ = temp_Leg1_;
	  Leg2P4_ = temp_Leg2_;
	}
	else {
	  Leg1P4_ = temp_Leg2_; tauIndex_1 = iter->index2_;
	  Leg2P4_ = temp_Leg1_; tauIndex_2 = iter->index1_;
	}
	++diTauCounter;
      }
      //Taus are selected
      
      if(tauIndex_1 >-1 && tauIndex_2 >-1) {
	//cout<<tauIndex_1<<"  "<<tauIndex_2<<endl;
	os=analysisTree.tau_charge[tauIndex_1]*analysisTree.tau_charge[tauIndex_2]<0;
	
	// looking for extra electron
	bool foundExtraElectron = false;
	for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {
	  if (analysisTree.electron_pt[ie]<ptVetoElectronCut) continue;
	  if (fabs(analysisTree.electron_eta[ie])>etaVetoElectronCut) continue;
	  if (fabs(analysisTree.electron_dxy[ie])>dxyVetoElectronCut) continue;
	  if (fabs(analysisTree.electron_dz[ie])>dzVetoElectronCut) continue;
	
	  bool electronMvaId = analysisTree.electron_mva_wp90_nontrig_Spring15_v1[ie];
	  if (!electronMvaId&&applyVetoElectronId) continue;
	  if (!analysisTree.electron_pass_conversion[ie]&&applyVetoElectronId) continue;
	  if (analysisTree.electron_nmissinginnerhits[ie]>1&&applyVetoElectronId) continue;
	  float neutralHadIsoEle = analysisTree.electron_neutralHadIso[ie];
	  float photonIsoEle = analysisTree.electron_photonIso[ie];
	  float chargedHadIsoEle = analysisTree.electron_chargedHadIso[ie];
	  float puIsoEle = analysisTree.electron_puIso[ie];
	  if (isIsoR03) {
	    neutralHadIsoEle = analysisTree.electron_r03_sumNeutralHadronEt[ie];
	    photonIsoEle = analysisTree.electron_r03_sumPhotonEt[ie];
	    chargedHadIsoEle = analysisTree.electron_r03_sumChargedHadronPt[ie];
	    puIsoEle = analysisTree.electron_r03_sumPUPt[ie];
	  }
	  float neutralIsoEle = neutralHadIsoEle + photonIsoEle - 0.5*puIsoEle;
	  neutralIsoEle = TMath::Max(float(0),neutralIsoEle); 
	  float absIsoEle =  chargedHadIsoEle + neutralIsoEle;
	  float relIsoEle = absIsoEle/analysisTree.electron_pt[ie];
	  if (relIsoEle>isoVetoElectronCut) continue;
	  foundExtraElectron = true;
	}
	// looking for extra muon's (dimuon veto)
	bool foundExtraMuon = false;
	vector<int> mu_dimuons; mu_dimuons.clear(); 
	for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {

	  float neutralHadIsoMu = analysisTree.muon_neutralHadIso[im];
	  float photonIsoMu = analysisTree.muon_photonIso[im];
	  float chargedHadIsoMu = analysisTree.muon_chargedHadIso[im];
	  float puIsoMu = analysisTree.muon_puIso[im];
	  if (isIsoR03) {
	    neutralHadIsoMu = analysisTree.muon_r03_sumNeutralHadronEt[im];
	    photonIsoMu = analysisTree.muon_r03_sumPhotonEt[im];
	    chargedHadIsoMu = analysisTree.muon_r03_sumChargedHadronPt[im];
	    puIsoMu = analysisTree.muon_r03_sumPUPt[im];
	  }
	  float neutralIsoMu = neutralHadIsoMu + photonIsoMu - 0.5*puIsoMu;
	  neutralIsoMu = TMath::Max(float(0),neutralIsoMu); 
	  float absIsoMu = chargedHadIsoMu + neutralIsoMu;
	  float relIsoMu = absIsoMu/analysisTree.muon_pt[im];

	  if (analysisTree.muon_pt[im]>ptDilepMuonCut&&
	      fabs(analysisTree.muon_eta[im])<etaDilepMuonCut&&
	          analysisTree.muon_isGlobal[im]&&
	          analysisTree.muon_isTracker[im]&&
	          analysisTree.muon_isPF[im]&&
	      fabs(analysisTree.muon_dxy[im])<dxyDilepMuonCut&&
	      fabs(analysisTree.muon_dz[im])<dzDilepMuonCut&&
	      relIsoMu<isoDilepMuonCut)
	    mu_dimuons.push_back(im);

	  
	  if (analysisTree.muon_pt[im]<ptVetoMuonCut) continue;
	  if (fabs(analysisTree.muon_eta[im])>etaVetoMuonCut) continue;
	  if (fabs(analysisTree.muon_dxy[im])>dxyVetoMuonCut) continue;
	  if (fabs(analysisTree.muon_dz[im])>dzVetoMuonCut) continue;
	  if (applyVetoMuonId && !analysisTree.muon_isMedium[im]) continue;
	  if (relIsoMu>isoVetoMuonCut) continue;
	  foundExtraMuon = true;
	}
	extraelec_veto = foundExtraElectron;
	extramuon_veto = foundExtraMuon;
	if (extraelec_veto) continue;
	if (extramuon_veto) continue;
	bool tauPass1 = 
	  analysisTree.tau_againstElectronVLooseMVA6[tauIndex_1]>0.5 &&
	  analysisTree.tau_againstMuonTight3[tauIndex_1]>0.5 &&
	  analysisTree.tau_byTightIsolationMVArun2v1DBoldDMwLT[tauIndex_1] > 0.5;
	bool tauPass2 = 
	  analysisTree.tau_againstElectronVLooseMVA6[tauIndex_2]>0.5 &&
	  analysisTree.tau_againstMuonTight3[tauIndex_2]>0.5 &&
	  analysisTree.tau_byTightIsolationMVArun2v1DBoldDMwLT[tauIndex_2] > 0.5;
	bool tauPass=tauPass1*tauPass2;
	if(!tauPass)continue;

	TLorentzVector tauLV_1; tauLV_1.SetXYZT(analysisTree.tau_px[tauIndex_1],
					    analysisTree.tau_py[tauIndex_1],
					    analysisTree.tau_pz[tauIndex_1],
					    analysisTree.tau_e[tauIndex_1]);
	TLorentzVector tauLV_2; tauLV_2.SetXYZT(analysisTree.tau_px[tauIndex_2],
					    analysisTree.tau_py[tauIndex_2],
					    analysisTree.tau_pz[tauIndex_2],
					    analysisTree.tau_e[tauIndex_2]);
	TLorentzVector diTauLV = tauLV_1 + tauLV_2;
	m_vis=diTauLV.M(); 
	// filling tau+ variables
	pt_plus = analysisTree.tau_pt[tauIndex_1];
	eta_plus = analysisTree.tau_eta[tauIndex_1];
	phi_plus = analysisTree.tau_phi[tauIndex_1];
	q_plus = -100;
	if (analysisTree.tau_charge[tauIndex_1]>0)
	  q_plus = 1;
	//mva_plus = -9999;
	d0_plus = analysisTree.tau_leadchargedhadrcand_dxy[tauIndex_1];
	dZ_plus = analysisTree.tau_leadchargedhadrcand_dz[tauIndex_1];
	iso_plus = analysisTree.tau_byTightIsolationMVArun2v1DBoldDMwLT[tauIndex_1];
	m_plus = analysisTree.tau_mass[tauIndex_1];
	//filling tau- variables
	pt_neg = analysisTree.tau_pt[tauIndex_2];
	eta_neg = analysisTree.tau_eta[tauIndex_2];
	phi_neg = analysisTree.tau_phi[tauIndex_2];
	q_neg = -100;
	if (analysisTree.tau_charge[tauIndex_2]<0)
	  q_neg = 1;
	//mva_neg = -9999;
	d0_neg = analysisTree.tau_leadchargedhadrcand_dxy[tauIndex_2];
	dZ_neg = analysisTree.tau_leadchargedhadrcand_dz[tauIndex_2];
	iso_neg = analysisTree.tau_byTightIsolationMVArun2v1DBoldDMwLT[tauIndex_2];
	m_neg = analysisTree.tau_mass[tauIndex_2];

      }
      
      againstMuonTight3_2 = analysisTree.tau_againstMuonTight3[tauIndex_2];
      againstElectronVTightMVA6_2 = analysisTree.tau_againstElectronVTightMVA6[tauIndex_2];
      againstMuonTight3_1 = analysisTree.tau_againstMuonTight3[tauIndex_1];
      againstElectronVTightMVA6_1 = analysisTree.tau_againstElectronVTightMVA6[tauIndex_1];
      // counting jets
      vector<unsigned int> jets; jets.clear();
      vector<unsigned int> jetshad; jetshad.clear();
      vector<unsigned int> jetspt20; jetspt20.clear();
      vector<unsigned int> bjets; bjets.clear();
      vector<unsigned int> hadjets; hadjets.clear();

      int indexLeadingJet = -1;
      float ptLeadingJet = -1;

      int indexSubLeadingJet = -1;
      float ptSubLeadingJet = -1;
      
      int indexLeadingBJet = -1;
      float ptLeadingBJet = -1;
      for (unsigned int jet=0; jet<analysisTree.pfjet_count; ++jet) {
	float jetEta = analysisTree.pfjet_eta[jet];
	float absJetEta = fabs(analysisTree.pfjet_eta[jet]);
	if (absJetEta>jetEtaCut) continue;
	float jetPt = analysisTree.pfjet_pt[jet];
	if (jetPt<jetPtLowCut) continue;
	// jetId
	bool isPFJetId = looseJetiD(analysisTree,int(jet));
	if (!isPFJetId) continue;
	float dR1 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
			   eta_plus,phi_plus);
	if (dR1<dRJetLeptonCut) continue;

	float dR2 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
                           eta_neg,phi_neg);
        if (dR2<dRJetLeptonCut) continue;

	jetspt20.push_back(jet);
	if (absJetEta<bJetEtaCut) { // jet within b-tagging acceptance
	  bool tagged = analysisTree.pfjet_btag[jet][nBTagDiscriminant]>btagCut; // b-jet

	  if (!isData) {
	    int flavor = abs(analysisTree.pfjet_flavour[jet]);
	        
	    double jet_scalefactor = 1;
	    double JetPtForBTag = jetPt;
	    double tageff = 1;
	        
	    if (flavor==5) {
	      if (JetPtForBTag>MaxBJetPt) JetPtForBTag = MaxBJetPt - 0.1;
	      if (JetPtForBTag<MinBJetPt) JetPtForBTag = MinBJetPt + 0.1;
	      jet_scalefactor = reader.eval_auto_bounds("central", BTagEntry::FLAV_B, absJetEta, JetPtForBTag);
	      tageff = tagEff_B->Interpolate(JetPtForBTag,absJetEta);
	    }
	    else if (flavor==4) {
	      if (JetPtForBTag>MaxBJetPt) JetPtForBTag = MaxBJetPt - 0.1;
	      if (JetPtForBTag<MinBJetPt) JetPtForBTag = MinBJetPt + 0.1;
	      jet_scalefactor = reader.eval_auto_bounds("central", BTagEntry::FLAV_C, absJetEta, JetPtForBTag);
	      tageff = tagEff_C->Interpolate(JetPtForBTag,absJetEta);
	    }
	    else {
	      if (JetPtForBTag>MaxLJetPt) JetPtForBTag = MaxLJetPt - 0.1;
	      if (JetPtForBTag<MinLJetPt) JetPtForBTag = MinLJetPt + 0.1;
	      jet_scalefactor = reader.eval_auto_bounds("central", BTagEntry::FLAV_UDSG, absJetEta, JetPtForBTag);
	      tageff = tagEff_Light->Interpolate(JetPtForBTag,absJetEta);
	    }
	        
	    if (tageff<1e-5)      tageff = 1e-5;
	    if (tageff>0.99999)   tageff = 0.99999;
	    rand.SetSeed((int)((jetEta+5)*100000));
	    double rannum = rand.Rndm();
	        
	    if (jet_scalefactor<1 && tagged) { // downgrade
	      double fraction = 1-jet_scalefactor;
	      if (rannum<fraction) {
		tagged = false;
		//std::cout << "downgrading " << std::endl;
	      }
	    }
	    if (jet_scalefactor>1 && !tagged) { // upgrade
	      double fraction = (jet_scalefactor-1.0)/(1.0/tageff-1.0);
	      if (rannum<fraction) { 
		tagged = true;
		//std::cout << "upgrading " << std::endl;
	      }
	    }
	  }
	}
	if (jetPt>jetPtHighCut)
	  jets.push_back(jet);
	if (indexLeadingJet>=0) {
	  if (jetPt<ptLeadingJet&&jetPt>ptSubLeadingJet) {
	    indexSubLeadingJet = jet;
	    ptSubLeadingJet = jetPt;
	  }
	}
	
	if (jetPt>ptLeadingJet) {
	  indexSubLeadingJet = indexLeadingJet;
	  ptSubLeadingJet = ptLeadingJet;
	  indexLeadingJet = jet;
	  ptLeadingJet = jetPt;
	}

      }//jet loop

      njets = jets.size();

      jpt_1 = -9999;
      jeta_1 = -9999;
      jphi_1 = -9999;
      jptraw_1 = -9999;
      jmva_1 = -9999;
      if (indexLeadingJet>=0) {
	jpt_1 = analysisTree.pfjet_pt[indexLeadingJet];
	jeta_1 = analysisTree.pfjet_eta[indexLeadingJet];
	jphi_1 = analysisTree.pfjet_phi[indexLeadingJet];
	jptraw_1 = analysisTree.pfjet_pt[indexLeadingJet]*analysisTree.pfjet_energycorr[indexLeadingJet];
	jmva_1 = -1; //analysisTree.pfjet_pu_jet_full_mva[indexLeadingJet];
      }
      jpt_2 = -9999;
      jeta_2 = -9999;
      jphi_2 = -9999;
      jptraw_2 = -9999;
      jmva_2 = -9999;
      if (indexLeadingJet>=0) {
	jpt_2 = analysisTree.pfjet_pt[indexLeadingJet];
	jeta_2 = analysisTree.pfjet_eta[indexLeadingJet];
	jphi_2 = analysisTree.pfjet_phi[indexLeadingJet];
	jptraw_2 = analysisTree.pfjet_pt[indexLeadingJet]*analysisTree.pfjet_energycorr[indexLeadingJet];
	jmva_2 = -1; //analysisTree.pfjet_pu_jet_full_mva[indexLeadingJet];
      }
      mjj =  -9999;
      jdeta = -9999;
      njetingap = 0;
      if (indexLeadingJet>=0 && indexSubLeadingJet>=0) {

	TLorentzVector jet1; jet1.SetPxPyPzE(analysisTree.pfjet_px[indexLeadingJet],
					     analysisTree.pfjet_py[indexLeadingJet],
					     analysisTree.pfjet_pz[indexLeadingJet],
					     analysisTree.pfjet_e[indexLeadingJet]);

	TLorentzVector jet2; jet2.SetPxPyPzE(analysisTree.pfjet_px[indexSubLeadingJet],
					     analysisTree.pfjet_py[indexSubLeadingJet],
					     analysisTree.pfjet_pz[indexSubLeadingJet],
					     analysisTree.pfjet_e[indexSubLeadingJet]);

	mjj = (jet1+jet2).M();
	jdeta = abs(analysisTree.pfjet_eta[indexLeadingJet]-
		    analysisTree.pfjet_eta[indexSubLeadingJet]);
 
	float etamax = analysisTree.pfjet_eta[indexLeadingJet];
	float etamin = analysisTree.pfjet_eta[indexSubLeadingJet];
	if (etamax<etamin) {
	  float tmp = etamax;
	  etamax = etamin;
	  etamin = tmp;
	}
	for (unsigned int jet=0; jet<jetspt20.size(); ++jet) {
	  int index = jetspt20.at(jet);
	  float etaX = analysisTree.pfjet_eta[index];
	  if (index!=indexLeadingJet&&index!=indexSubLeadingJet&&etaX>etamin&&etaX<etamax) 
	    njetingap++;
	}
      }
      // met
      float pfmet_ex = analysisTree.pfmetcorr_ex;
      float pfmet_ey = analysisTree.pfmetcorr_ey;
      met = TMath::Sqrt(pfmet_ex*pfmet_ex + pfmet_ey*pfmet_ey);
      metphi = TMath::ATan2(pfmet_ey,pfmet_ex);

      // genmet
      float genmet_x = analysisTree.genmet_ex;
      float genmet_y = analysisTree.genmet_ey;
      genmet = TMath::Sqrt(genmet_x*genmet_x+genmet_y*genmet_y);
      genmetphi = TMath::ATan2(genmet_y,genmet_x);

      metcov00 = analysisTree.pfmetcorr_sigxx;
      metcov01 = analysisTree.pfmetcorr_sigxy;
      metcov10 = analysisTree.pfmetcorr_sigyx;
      metcov11 = analysisTree.pfmetcorr_sigyy;
      float pfmet_corr_ex = pfmet_ex;
      float pfmet_corr_ey = pfmet_ey;

      int njetsforrecoil = njets;
      if (isW) njetsforrecoil = njets + 1;
      if ((isDY||isW)&&!isData)
	recoilMetCorrector.CorrectByMeanResolution(pfmet_ex,pfmet_ey,bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,pfmet_corr_ex,pfmet_corr_ey);

      pfmet_ex = pfmet_corr_ex;
      pfmet_ey = pfmet_corr_ey;
      met = TMath::Sqrt(pfmet_ex*pfmet_ex+pfmet_ey*pfmet_ey);
      metphi = TMath::ATan2(pfmet_ey,pfmet_ex);

      float dPhiMETTau_1 = dPhiFrom2P(analysisTree.tau_px[tauIndex_1],analysisTree.tau_py[tauIndex_1],
				   pfmet_ex,pfmet_ey);
      mt_plus = TMath::Sqrt(2*met*pt_plus*(1-TMath::Cos(dPhiMETTau_1)));

      float dPhiMETTau_2 = dPhiFrom2P(analysisTree.tau_px[tauIndex_2],analysisTree.tau_py[tauIndex_2],
				    pfmet_ex,pfmet_ey);
      mt_neg = TMath::Sqrt(2*met*pt_neg*(1-TMath::Cos(dPhiMETTau_2)));
      //QCD supresser parameter
      double per_px=(analysisTree.tau_px[tauIndex_1]+analysisTree.tau_px[tauIndex_2]+pfmet_ex);
      double per_py=(analysisTree.tau_py[tauIndex_1]+analysisTree.tau_py[tauIndex_2]+pfmet_ey);
      Prompt_pT=sqrt(per_px*per_px+per_py*per_py);

      gen_match_1 = 6;
      gen_match_2 = 6;

      isZTT = false;
      isZL  = false;
      isZJ  = false;

      float minDR_1 = 0.2;
      float minDR_2 = 0.2;
      if (!isData) {
	for (unsigned int igen=0; igen < analysisTree.genparticles_count; ++igen) {
	  float ptGen = PtoPt(analysisTree.genparticles_px[igen],
			      analysisTree.genparticles_py[igen]);
	  bool type1 = abs(analysisTree.genparticles_pdgid[igen])==11 && analysisTree.genparticles_isPrompt[igen] && ptGen>8;
	  bool type2 = abs(analysisTree.genparticles_pdgid[igen])==13 && analysisTree.genparticles_isPrompt[igen] && ptGen>8;
	  bool type3 = abs(analysisTree.genparticles_pdgid[igen])==11 && analysisTree.genparticles_isDirectPromptTauDecayProduct[igen] && ptGen>8;
	  bool type4 = abs(analysisTree.genparticles_pdgid[igen])==13 && analysisTree.genparticles_isDirectPromptTauDecayProduct[igen] && ptGen>8;
	  bool isAnyType = type1 || type2 || type3 || type4;
	  if (isAnyType) {
	    float etaGen = PtoEta(analysisTree.genparticles_px[igen],
				  analysisTree.genparticles_py[igen],
				  analysisTree.genparticles_pz[igen]);
	    float phiGen = PtoPhi(analysisTree.genparticles_px[igen],
				  analysisTree.genparticles_py[igen]);
	    float deltaR_1 = deltaR(analysisTree.tau_eta[tauIndex_1],analysisTree.tau_phi[tauIndex_1],//bug fix
				    etaGen,phiGen);
	    if (deltaR_1<minDR_1) {
	      minDR_1 = deltaR_1;
	      if (type1) gen_match_1 = 1;
	      else if (type2) gen_match_1 = 2;
	      else if (type3) gen_match_1 = 3;
	      else if (type4) gen_match_1 = 4;
	    }
	        
	    float deltaR_2 = deltaR(analysisTree.tau_eta[tauIndex_2],analysisTree.tau_phi[tauIndex_2],
				    etaGen,phiGen);
	    if (deltaR_2<minDR_2) {
	      minDR_2 = deltaR_2;
	      if (type1) gen_match_2 = 1;
	      else if (type2) gen_match_2 = 2;
	      else if (type3) gen_match_2 = 3;
	      else if (type4) gen_match_2 = 4;
	    }
	  }
	}
	for (unsigned int igen=0; igen < analysisTree.gentau_count; ++igen) {
	  if (analysisTree.gentau_visibleNoLep_pt[igen]>15.) {
	    float deltaR_1 = deltaR(analysisTree.tau_eta[tauIndex_1],analysisTree.tau_phi[tauIndex_1],
				    analysisTree.gentau_visibleNoLep_eta[igen],analysisTree.gentau_visibleNoLep_phi[igen]);
	    if (deltaR_1<minDR_1) {
	      minDR_1 = deltaR_1;
	      gen_match_1 = 5;
	    }
	    float deltaR_2 = deltaR(analysisTree.tau_eta[tauIndex_2],analysisTree.tau_phi[tauIndex_2],
				    analysisTree.gentau_visibleNoLep_eta[igen],analysisTree.gentau_visibleNoLep_phi[igen]);
	    if (deltaR_2<minDR_2) {
	      minDR_2 = deltaR_2;
	      gen_match_2 = 5;
	    }
	  }
	}
	if (gen_match_2==5 && gen_match_1==5) isZTT = true;
	else if (gen_match_2==6||gen_match_1==6) isZJ = true;
	else isZL = true;
	isZLL = isZL || isZJ;
      }
      tree->Fill();
      selEvents++;
    }//end of event
    nFiles++;
    delete _tree;
    file_->Close();
    delete file_;
  }//close file
  
  std::cout << std::endl;
  int allEvents = int(inputEventsH->GetEntries());
  std::cout << "Total number of input events    = " << allEvents << std::endl;
  std::cout << "Total number of events in Tree  = " << nEvents << std::endl;
  std::cout << "Total number of selected events = " << selEvents << std::endl;
  std::cout << std::endl;

  file->cd("");
  file->Write();
  file->Close();
  delete file;
}//close of main
