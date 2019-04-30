#ifndef AC1B_h
#define AC1B_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "map"

class AC1B{
 public:
  TTree *fChain;
  Int_t *fCurrent;

  //List of leaves
  UInt_t  primvertex_count;
  Float_t primvertex_x;
  Float_t primvertex_y;
  Float_t primvertex_z;
  Float_t primvertex_cov[6];
  Int_t   muon_count;
  Float_t muon_pt[100];
  Float_t muon_eta[100];
  Float_t muon_phi[100];
  Int_t   electron_count;
  Float_t electron_pt[100];
  Float_t electron_eta[100];
  Float_t electron_phi[100];
  Int_t   tau_count;
  Float_t tau_charge[100];
  Float_t tau_pt[100];
  Float_t tau_eta[100];
  Float_t tau_phi[100];
  UInt_t  tau_decayMode[100];
  Float_t tau_leadchargedhadrcand_px[100];
  Float_t tau_leadchargedhadrcand_py[100];
  Float_t tau_leadchargedhadrcand_pz[100];
  Float_t tau_leadchargedhadrcand_mass[100];
  Float_t tau_pca3D_x[1000];
  Float_t tau_pca3D_y[1000];
  Float_t tau_pca3D_z[1000];
  Float_t tau_pca2D_x[1000];
  Float_t tau_pca2D_y[1000];
  Float_t tau_pca2D_z[1000];
  UInt_t  genparticles_count;
  Float_t genparticles_vx[100];
  Float_t genparticles_vy[100];
  Float_t genparticles_vz[100];
  Float_t genparticles_px[100];
  Float_t genparticles_py[100];
  Float_t genparticles_pz[100];
  Float_t genparticles_e[100];
  Int_t   genparticles_pdgid[100];
  Int_t   genparticles_isLastCopy[100];
  UInt_t  genparticles_isPrompt[100];
  Int_t   genparticles_isDirectPromptTauDecayProduct[100];
  Int_t   genparticles_status[100];
  UChar_t genparticles_mother[100];
  Int_t   gentau_count;
  Float_t gentau_charge[100];
  Float_t gentau_mother[100];
  Float_t gentau_visible_px[100];
  Float_t gentau_visible_py[100];
  Float_t gentau_visible_pz[100];
  Float_t gentau_visible_e[100];
  Float_t gentau_px[100];
  Float_t gentau_py[100];
  Float_t gentau_pz[100];
  Float_t gentau_e[100];
  UInt_t  gentau_decayMode[100];
  Float_t gentau_visibleNoLep_pt[100];
  Int_t   refitvertex_count;
  Float_t refitvertex_x[100];
  Float_t refitvertex_y[100];
  Float_t refitvertex_z[100];
  Int_t   refitvertex_eleIndex[100][100];
  Int_t   refitvertex_muIndex[100][100];
  Int_t   refitvertex_tauIndex[100][100];
  Float_t refitvertex_cov[100][6];
  Float_t beamspot_x;
  Float_t beamspot_y;
  Float_t beamspot_z;
  Float_t beamspot_cov[6];
  
  Float_t tau_byMediumIsolationMVArun2v1DBoldDMwLT[100];
  Float_t tau_againstElectronTightMVA6[100];
  //vector<float> *tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT = new vector<float>();
  Float_t tau_againstMuonTight3[100];
  //List of branches
  TBranch *b_primvertex_count;
  TBranch *b_primvertex_x;   
  TBranch *b_primvertex_y;   
  TBranch *b_primvertex_z;
  TBranch *b_primvertex_cov;
  TBranch *b_muon_count;
  TBranch *b_muon_pt;
  TBranch *b_muon_eta;
  TBranch *b_muon_phi;
  TBranch *b_electron_count;
  TBranch *b_electron_pt;
  TBranch *b_electron_eta;
  TBranch *b_electron_phi;
  TBranch *b_tau_decayMode;
  TBranch *b_tau_count;
  TBranch *b_tau_charge;
  TBranch *b_tau_pt;
  TBranch *b_tau_eta;
  TBranch *b_tau_phi;
  TBranch *b_tau_leadchargedhadrcand_px;
  TBranch *b_tau_leadchargedhadrcand_py;
  TBranch *b_tau_leadchargedhadrcand_pz;
  TBranch *b_tau_leadchargedhadrcand_mass;
  TBranch *b_genparticles_count;
  TBranch *b_genparticles_vx;
  TBranch *b_genparticles_vy;
  TBranch *b_genparticles_vz;
  TBranch *b_genparticles_px;
  TBranch *b_genparticles_py;
  TBranch *b_genparticles_pz;
  TBranch *b_genparticles_e;
  TBranch *b_genparticles_mother;
  TBranch *b_genparticles_pdgid;
  TBranch *b_genparticles_isLastCopy;
  TBranch *b_genparticles_isPrompt;
  TBranch *b_genparticles_isDirectPromptTauDecayProduct;
  TBranch *b_genparticles_status;
  TBranch *b_gentau_count;
  TBranch *b_gentau_mother;
  TBranch *b_gentau_charge;
  TBranch *b_gentau_visible_px;
  TBranch *b_gentau_visible_py;
  TBranch *b_gentau_visible_pz;
  TBranch *b_gentau_visible_e;
  TBranch *b_gentau_px;
  TBranch *b_gentau_py;
  TBranch *b_gentau_pz;
  TBranch *b_gentau_e;
  TBranch *b_gentau_decayMode;
  TBranch *b_gentau_visibleNoLep_pt;
  TBranch *b_refitvertex_count;
  TBranch *b_refitvertex_x;   
  TBranch *b_refitvertex_y;   
  TBranch *b_refitvertex_z;
  TBranch *b_refitvertex_eleIndex;
  TBranch *b_refitvertex_muIndex;
  TBranch *b_refitvertex_tauIndex;
  TBranch *b_refitvertex_cov;
  TBranch *b_tau_pca2D_x;
  TBranch *b_tau_pca2D_y;
  TBranch *b_tau_pca2D_z;
  TBranch *b_tau_pca3D_x;
  TBranch *b_tau_pca3D_y;
  TBranch *b_tau_pca3D_z;
  TBranch *b_beamspot_x;
  TBranch *b_beamspot_y;
  TBranch *b_beamspot_z;
  TBranch *b_beamspot_cov;
  TBranch *b_tau_byMediumIsolationMVArun2v1DBoldDMwLT;
  TBranch *b_tau_againstMuonTight3;
  TBranch *b_tau_againstElectronTightMVA6;
  Int_t GetEntry(Long64_t entry)
  {
    // Read contents of entry.
    if (!fChain) return 0;
    return fChain->GetEntry(entry);
  }
  Int_t GetEntries()
  {
    if (!fChain) return 0;
    return fChain->GetEntries();
  }
  void Init(TTree *tree, bool isData) 
  {
    tree->SetMaxVirtualSize(3000000);
    fChain = tree;
    fChain->SetMakeClass(1);
    
    fChain->SetBranchAddress("primvertex_count", &primvertex_count, &b_primvertex_count);
    fChain->SetBranchAddress("primvertex_x", &primvertex_x, &b_primvertex_x);
    fChain->SetBranchAddress("primvertex_y", &primvertex_y, &b_primvertex_y);
    fChain->SetBranchAddress("primvertex_z", &primvertex_z, &b_primvertex_z);
    fChain->SetBranchAddress("primvertex_cov", &primvertex_cov, &b_primvertex_cov);
    fChain->SetBranchAddress("muon_count", &muon_count,&b_muon_count);
    fChain->SetBranchAddress("muon_pt", muon_pt,&b_muon_pt);
    fChain->SetBranchAddress("muon_eta", muon_eta,&b_muon_eta);
    fChain->SetBranchAddress("muon_phi", muon_phi,&b_muon_phi);
    fChain->SetBranchAddress("electron_count", &electron_count,&b_electron_count);
    fChain->SetBranchAddress("electron_pt", electron_pt,&b_electron_pt);
    fChain->SetBranchAddress("electron_eta", electron_eta,&b_electron_eta);
    fChain->SetBranchAddress("electron_phi", electron_phi,&b_electron_phi);
    fChain->SetBranchAddress("tau_decayMode", tau_decayMode, &b_tau_decayMode);
    fChain->SetBranchAddress("tau_count", &tau_count,&b_tau_count);
    fChain->SetBranchAddress("tau_charge", tau_charge, &b_tau_charge);
    fChain->SetBranchAddress("tau_pt", tau_pt,&b_tau_pt);
    fChain->SetBranchAddress("tau_eta", tau_eta,&b_tau_eta);
    fChain->SetBranchAddress("tau_phi", tau_phi,&b_tau_phi);
    fChain->SetBranchAddress("tau_pca3D_x", tau_pca3D_x,&b_tau_pca3D_x);
    fChain->SetBranchAddress("tau_pca3D_y", tau_pca3D_y,&b_tau_pca3D_y);
    fChain->SetBranchAddress("tau_pca3D_z", tau_pca3D_z,&b_tau_pca3D_z);
    fChain->SetBranchAddress("tau_pca2D_x", tau_pca2D_x,&b_tau_pca2D_x);
    fChain->SetBranchAddress("tau_pca2D_y", tau_pca2D_y,&b_tau_pca2D_y);
    fChain->SetBranchAddress("tau_pca2D_z", tau_pca2D_z,&b_tau_pca2D_z);
    fChain->SetBranchAddress("tau_leadchargedhadrcand_px", tau_leadchargedhadrcand_px,&b_tau_leadchargedhadrcand_px);
    fChain->SetBranchAddress("tau_leadchargedhadrcand_py", tau_leadchargedhadrcand_py,&b_tau_leadchargedhadrcand_py);
    fChain->SetBranchAddress("tau_leadchargedhadrcand_pz", tau_leadchargedhadrcand_pz,&b_tau_leadchargedhadrcand_pz);
    fChain->SetBranchAddress("tau_leadchargedhadrcand_mass", tau_leadchargedhadrcand_mass,&b_tau_leadchargedhadrcand_mass);
    fChain->SetBranchAddress("genparticles_count",&genparticles_count, &b_genparticles_count);
    fChain->SetBranchAddress("genparticles_vx", genparticles_vx, &b_genparticles_vx);
    fChain->SetBranchAddress("genparticles_vy", genparticles_vy, &b_genparticles_vy);
    fChain->SetBranchAddress("genparticles_vz", genparticles_vz, &b_genparticles_vz);
    fChain->SetBranchAddress("genparticles_px", genparticles_px, &b_genparticles_px);
    fChain->SetBranchAddress("genparticles_py", genparticles_py, &b_genparticles_py);
    fChain->SetBranchAddress("genparticles_pz", genparticles_pz, &b_genparticles_pz);
    fChain->SetBranchAddress("genparticles_e", genparticles_e, &b_genparticles_e);
    fChain->SetBranchAddress("genparticles_status", genparticles_status, &b_genparticles_status);
    fChain->SetBranchAddress("genparticles_pdgid", genparticles_pdgid, &b_genparticles_pdgid);
    fChain->SetBranchAddress("genparticles_isPrompt", genparticles_isPrompt, &b_genparticles_isPrompt);
    fChain->SetBranchAddress("genparticles_isLastCopy", genparticles_isLastCopy, &b_genparticles_isLastCopy);
    fChain->SetBranchAddress("genparticles_isDirectPromptTauDecayProduct", genparticles_isDirectPromptTauDecayProduct, &b_genparticles_isDirectPromptTauDecayProduct);
    fChain->SetBranchAddress("genparticles_mother", genparticles_mother, &b_genparticles_mother);
    fChain->SetBranchAddress("gentau_count",&gentau_count, &b_gentau_count);
    fChain->SetBranchAddress("gentau_mother", gentau_mother, &b_gentau_mother);
    fChain->SetBranchAddress("gentau_charge", gentau_charge, &b_gentau_charge);
    fChain->SetBranchAddress("gentau_px", gentau_px, &b_gentau_px);
    fChain->SetBranchAddress("gentau_py", gentau_py, &b_gentau_py);
    fChain->SetBranchAddress("gentau_pz", gentau_pz, &b_gentau_pz);
    fChain->SetBranchAddress("gentau_e", gentau_e, &b_gentau_e);
    fChain->SetBranchAddress("gentau_decayMode", gentau_decayMode, &b_gentau_decayMode);
    fChain->SetBranchAddress("gentau_visible_px", gentau_visible_px, &b_gentau_visible_px);
    fChain->SetBranchAddress("gentau_visible_py", gentau_visible_py, &b_gentau_visible_py);
    fChain->SetBranchAddress("gentau_visible_pz", gentau_visible_pz, &b_gentau_visible_pz);
    fChain->SetBranchAddress("gentau_visible_e", gentau_visible_e, &b_gentau_visible_e);
    fChain->SetBranchAddress("gentau_visibleNoLep_pt", gentau_visibleNoLep_pt, &b_gentau_visibleNoLep_pt);
    fChain->SetBranchAddress("refitvertex_count", &refitvertex_count, &b_refitvertex_count);
    fChain->SetBranchAddress("refitvertex_x", refitvertex_x, &b_refitvertex_x);
    fChain->SetBranchAddress("refitvertex_y", refitvertex_y, &b_refitvertex_y);
    fChain->SetBranchAddress("refitvertex_z", refitvertex_z, &b_refitvertex_z);
    fChain->SetBranchAddress("refitvertex_cov", refitvertex_cov, &b_refitvertex_cov);
    fChain->SetBranchAddress("refitvertex_eleIndex", refitvertex_eleIndex, &b_refitvertex_eleIndex);
    fChain->SetBranchAddress("refitvertex_muIndex", refitvertex_muIndex, &b_refitvertex_muIndex);
    fChain->SetBranchAddress("refitvertex_tauIndex", refitvertex_tauIndex, &b_refitvertex_tauIndex);
    fChain->SetBranchAddress("beamspot_x",&beamspot_x,&b_beamspot_x);
    fChain->SetBranchAddress("beamspot_y",&beamspot_y,&b_beamspot_y);
    fChain->SetBranchAddress("beamspot_z",&beamspot_z,&b_beamspot_z);
    fChain->SetBranchAddress("beamspot_cov",beamspot_cov,&b_beamspot_cov);
    fChain->SetBranchAddress("tau_byMediumIsolationMVArun2v1DBoldDMwLT", tau_byMediumIsolationMVArun2v1DBoldDMwLT, &b_tau_byMediumIsolationMVArun2v1DBoldDMwLT);
    fChain->SetBranchAddress("tau_againstMuonTight3", tau_againstMuonTight3, &b_tau_againstMuonTight3);
    fChain->SetBranchAddress("tau_againstElectronTightMVA6", tau_againstElectronTightMVA6, &b_tau_againstElectronTightMVA6);
    
  }
  
  AC1B(TTree *tree, bool isData) : fChain(0)
    {
      
      if (tree == 0) {
	TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("output_MC.root");
	if (!f || !f->IsOpen()) {
	  f = new TFile("output_MC.root");
	}
	TDirectory * dir = (TDirectory*)f->Get("output_MC.root:/makeroottree");
	dir->GetObject("AC1B",tree);
	
      }
      Init(tree, isData);
    }
  
  ~AC1B()
    {
      if (!fChain) return;
      delete fChain->GetCurrentFile();
    }
  
};
#endif


