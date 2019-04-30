#include "AC1B.h"
#include "Math/GenVector/Boost.h"
#include <TMath.h>
#include <TROOT.h>
#include "Math/Point3D.h"
#include "Math/LorentzVector.h"

typedef ROOT::Math::XYZPointD Point3D;
typedef ROOT::Math::XYZVectorD PV;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LV;
#define  PI 3.14159265358979312e+00

/*------------------------GenMatch function-----------------*/


float deltaR(float eta_1, float phi_1, float eta_2, float phi_2){
  float deta = fabs(eta_1-eta_2);
  float dphi = fabs(phi_1-phi_2);
  if(dphi > PI)
    dphi = 2*PI -dphi;
  return abs(sqrt(deta*deta+dphi*dphi));
}

//Matching function for gentaus with reco taus
float GenMatch(AC1B *analysisTree, TString particleType, int Index,TString istau,int igen,float *minDR){
  
  float gen_match = 6;
  
  // isZTT = false;
  // isZL  = false;
  // isZJ  = false;
  
 
  float pt,eta,phi;
  if(particleType=="m"){
    pt =      analysisTree->muon_pt[Index]; 
    eta =     analysisTree->muon_eta[Index]; 
    phi =     analysisTree->muon_phi[Index];}
  if(particleType=="e"){         
    pt =      analysisTree->electron_pt[Index];
    eta =     analysisTree->electron_eta[Index]; 
    phi =     analysisTree->electron_phi[Index];}
  if(particleType=="t"){         
    pt =      analysisTree->tau_pt[Index];
    eta =     analysisTree->tau_eta[Index]; 
    phi =     analysisTree->tau_phi[Index];
  }
  
  if (istau=="no"){
    //for (unsigned int igen=0; igen < analysisTree.genparticles_count; ++igen) {
    
    TLorentzVector genLV; genLV.SetXYZT(analysisTree->genparticles_px[igen],
					analysisTree->genparticles_py[igen],
					analysisTree->genparticles_pz[igen],
					analysisTree->genparticles_e[igen]);
    float ptGen = genLV.Pt();
    bool type1 = abs(analysisTree->genparticles_pdgid[igen])==11 && analysisTree->genparticles_isPrompt[igen] && ptGen>8;
    bool type2 = abs(analysisTree->genparticles_pdgid[igen])==13 && analysisTree->genparticles_isPrompt[igen] && ptGen>8;
    bool type3 = abs(analysisTree->genparticles_pdgid[igen])==11 && analysisTree->genparticles_isDirectPromptTauDecayProduct[igen] && ptGen>8;
    bool type4 = abs(analysisTree->genparticles_pdgid[igen])==13 && analysisTree->genparticles_isDirectPromptTauDecayProduct[igen] && ptGen>8;
    
    bool isAnyType = type1 || type2 || type3 || type4;
    if (isAnyType && analysisTree->genparticles_status[igen]==1) {
      float etaGen = genLV.Eta();
      float phigen = genLV.Phi();
      float dR = deltaR(eta,phi,etaGen,phigen);
      if (dR<(*minDR)) {
	minDR = &dR;
	if (type1) gen_match = 1;
	else if (type2) gen_match = 2;
	else if (type3) gen_match = 3;
	else if (type4) gen_match = 4;
      }	
    }
  }
  if(istau=="yes"){
    if (gen_match<1||gen_match>4) {
      //for (unsigned int igen=0; igen < analysisTree.gentau_count; ++igen) {
      
      if (analysisTree->gentau_visibleNoLep_pt[igen]>15.) {
	TLorentzVector genTauLV; genTauLV.SetXYZT(analysisTree->gentau_visible_px[igen],
						  analysisTree->gentau_visible_py[igen],
						  analysisTree->gentau_visible_pz[igen],
						  analysisTree->gentau_visible_e[igen]);
	float dR = deltaR(eta,phi,
			  genTauLV.Eta(),genTauLV.Phi());
	if (dR<(*minDR)) {
	  minDR = &dR;
	  gen_match = 5;
	}
      }
    }
  }
  return gen_match;  
}

void phi_rec_distr(){
  gROOT->LoadMacro("AC1B.h+");
  ifstream inFile;
 
  inFile.open("ggHToTauTau_with.txt");  // The file list calling Ntuple files from nfs directory
  TString infilename_;
  int nf=0;
  while(!inFile.eof()){
 
    cout<<"File No: "<<nf<<" Taking"<<endl;
    inFile >> infilename_;
    TString stri = to_string(nf);
    TFile* file = new TFile(infilename_);
    TTree* tree=(TTree*)file->Get("/makeroottree/AC1B");
    AC1B* t1 = new AC1B(tree,false);
    TFile* f = new TFile("ofile_"+stri+".root","recreate");
    TH1F* h_phistar_0 =  new TH1F("h_phistar_0","",20,0,6);
   
    Int_t TotEntries=t1->GetEntries(),iEv=0;
    while(iEv<TotEntries){
      float kk=0.2;
      float* minDR;minDR=&kk;
      Float_t tauP_pt=0,tauN_pt=0;
      Float_t pionPVtx_x=0,pionPVtx_y=0,pionPVtx_z=0;
      Float_t pionNVtx_x=0,pionNVtx_y=0,pionNVtx_z=0;
      Float_t pionP_px=0,pionP_py=0,pionP_pz=0,pionP_e=0;
      Float_t pionN_px=0,pionN_py=0,pionN_pz=0,pionN_e=0;
      t1->GetEntry(iEv);
    
      for(int gtau=0;gtau<t1->gentau_count;gtau++){
	bool found_piP=false; bool found_piN=false;
	for(int rtau=0;rtau<t1->tau_count;rtau++){
	  
	  if(t1->tau_byMediumIsolationMVArun2v1DBoldDMwLT[rtau]>0.5&&t1->tau_againstMuonTight3[rtau]>0.5 && t1->tau_againstElectronTightMVA6[rtau]>0.5 && t1->tau_pt[rtau]>30. && fabs(t1->tau_eta[rtau])<2.3){
	  
	    if(t1->tau_charge[rtau]>0&&t1->tau_decayMode[rtau]==0 && !found_piP){
	     
	      if(GenMatch(t1,"t",rtau,"yes",gtau,minDR)==5){	
		tauP_pt=sqrt(t1->gentau_visible_px[gtau]*t1->gentau_visible_px[gtau]+t1->gentau_visible_py[gtau]*t1->gentau_visible_py[gtau]);
		pionPVtx_x=t1->tau_pca3D_x[rtau];
		pionPVtx_y=t1->tau_pca3D_y[rtau];
		pionPVtx_z=t1->tau_pca3D_z[rtau];
		pionP_px=t1->tau_leadchargedhadrcand_px[rtau];
		pionP_py=t1->tau_leadchargedhadrcand_py[rtau];
		pionP_pz=t1->tau_leadchargedhadrcand_pz[rtau];
		pionP_e=sqrt(pionP_px*pionP_px+pionP_py*pionP_py+pionP_pz*pionP_pz+(t1->tau_leadchargedhadrcand_mass[rtau]*t1->tau_leadchargedhadrcand_mass[rtau]));
		found_piP=true;
	      }
	    }//tau+
	    
	
	
	    if(t1->tau_charge[rtau]<0 &&t1->tau_decayMode[rtau]==0 && !found_piN){
	      
	      if(GenMatch(t1,"t",rtau,"yes",gtau,minDR)==5){	
		tauN_pt=sqrt(t1->gentau_visible_px[gtau]*t1->gentau_visible_px[gtau]+t1->gentau_visible_py[gtau]*t1->gentau_visible_py[gtau]);
		pionNVtx_x=t1->tau_pca3D_x[rtau];
		pionNVtx_y=t1->tau_pca3D_y[rtau];
		pionNVtx_z=t1->tau_pca3D_z[rtau];
		pionN_px=t1->tau_leadchargedhadrcand_px[rtau];
		pionN_py=t1->tau_leadchargedhadrcand_py[rtau];
		pionN_pz=t1->tau_leadchargedhadrcand_pz[rtau];
		pionN_e=sqrt(pionN_px*pionN_px+pionN_py*pionN_py+pionN_pz*pionN_pz+(t1->tau_leadchargedhadrcand_mass[rtau]*t1->tau_leadchargedhadrcand_mass[rtau]));
		found_piN=true;
	      }
	    }
	  }//tau-
       }
	  
	  //}
      }
     
      LV pionP_lab(pionP_px,pionP_py,pionP_pz,pionP_e);
      LV pionN_lab(pionN_px,pionN_py,pionN_pz,pionN_e);
      if(pionP_lab.pt()>0 && pionN_lab.pt()>0){
        float Higgs_pt=tauP_pt+tauN_pt;

	//compute phi_star in pi-pi+ rf.
	LV pion_pair_lab=pionP_lab+pionN_lab;
	ROOT::Math::Boost boost_to_rf_rec(pion_pair_lab.BoostToCM());
	//boost to rf of two pion system
	LV pionP_rf_2p=boost_to_rf_rec(pionP_lab);
	LV pionN_rf_2p=boost_to_rf_rec(pionN_lab);
	//get unit momentum vectors
	PV pionP_lab_3d=pionP_lab.Vect();
	PV pionN_lab_3d=pionN_lab.Vect();
	PV pionP_lab_3d_u=pionP_lab_3d/sqrt(pionP_lab_3d.mag2());
	PV pionN_lab_3d_u=pionN_lab_3d/sqrt(pionN_lab_3d.mag2());
	//Get PCAs, by extrapolating the line represented by pion vertex and pion momentum
	double tP=pionP_lab_3d_u.x()*(t1->primvertex_x-pionPVtx_x)+pionP_lab_3d_u.y()*(t1->primvertex_y-pionPVtx_y)+pionP_lab_3d_u.z()*(t1->primvertex_z-pionPVtx_z);
	double tN=pionN_lab_3d_u.x()*(t1->primvertex_x-pionNVtx_x)+pionN_lab_3d_u.y()*(t1->primvertex_y-pionNVtx_y)+pionN_lab_3d_u.z()*(t1->primvertex_z-pionNVtx_z);

	Point3D pcaP(pionPVtx_x+pionP_lab_3d_u.x()*tP,pionPVtx_y+pionP_lab_3d_u.y()*tP,pionPVtx_z+pionP_lab_3d_u.z()*tP);
	Point3D pcaN(pionNVtx_x+pionN_lab_3d_u.x()*tN,pionNVtx_y+pionN_lab_3d_u.y()*tN,pionNVtx_z+pionN_lab_3d_u.z()*tN);
	//Get the normalized IP vectors
	PV ipvP_lab_rec_3d(pcaP.x()-t1->primvertex_x,pcaP.y()-t1->primvertex_y,pcaP.z()-t1->primvertex_z);
	PV ipvN_lab_rec_3d(pcaN.x()-t1->primvertex_x,pcaN.y()-t1->primvertex_y,pcaN.z()-t1->primvertex_z);

	LV ipvP_lab_rec(ipvP_lab_rec_3d.x()/sqrt(ipvP_lab_rec_3d.mag2()),ipvP_lab_rec_3d.y()/sqrt(ipvP_lab_rec_3d.mag2()),ipvP_lab_rec_3d.z()/sqrt(ipvP_lab_rec_3d.mag2()),0);
	LV ipvN_lab_rec(ipvN_lab_rec_3d.x()/sqrt(ipvN_lab_rec_3d.mag2()),ipvN_lab_rec_3d.y()/sqrt(ipvN_lab_rec_3d.mag2()),ipvN_lab_rec_3d.z()/sqrt(ipvN_lab_rec_3d.mag2()),0);

	//boost to ZMF
	LV ipvP_rf_rec = boost_to_rf_rec(ipvP_lab_rec);
	LV ipvN_rf_rec = boost_to_rf_rec(ipvN_lab_rec);

	//Get only the position component of the vector
	PV ipvP_rf_rec_3d = ipvP_rf_rec.Vect();
	PV ipvN_rf_rec_3d = ipvN_rf_rec.Vect();
	PV pionP_rf_2p_3d = pionP_rf_2p.Vect();
	PV pionN_rf_2p_3d = pionN_rf_2p.Vect();

	//Get the unit (normalized) vector along pion momentum
	PV pionP_rf_2p_3d_u = pionP_rf_2p_3d/TMath::Sqrt(pionP_rf_2p_3d.mag2());
	PV pionN_rf_2p_3d_u = pionN_rf_2p_3d/TMath::Sqrt(pionN_rf_2p_3d.mag2());
	//Get the longitudinal component of IP vector parallel to pion momenta
	PV ipvP_rf_rec_3d_l = ipvP_rf_rec_3d.Dot(pionP_rf_2p_3d_u)*pionP_rf_2p_3d_u;
	PV ipvN_rf_rec_3d_l = ipvN_rf_rec_3d.Dot(pionN_rf_2p_3d_u)*pionN_rf_2p_3d_u;
	//Get IP vector normal to pion momenta
	PV ipvP_rf_rec_3d_t = ipvP_rf_rec_3d - ipvP_rf_rec_3d_l;
	PV ipvN_rf_rec_3d_t = ipvN_rf_rec_3d - ipvN_rf_rec_3d_l;
	//Get normalized normal IP vector
	PV ipvP_rf_rec_3d_t_u = ipvP_rf_rec_3d_t/TMath::Sqrt(ipvP_rf_rec_3d_t.mag2());
	PV ipvN_rf_rec_3d_t_u = ipvN_rf_rec_3d_t/TMath::Sqrt(ipvN_rf_rec_3d_t.mag2());

	//Get CP angle in ZMF
	double phi_star_rec = TMath::ACos(ipvP_rf_rec_3d_t_u.Dot(ipvN_rf_rec_3d_t_u));
	double OstarCP_rec = pionN_rf_2p_3d_u.Dot(ipvP_rf_rec_3d_t_u.Cross(ipvN_rf_rec_3d_t_u));
	double phi_star_rec_cp = (OstarCP_rec >= 0) ? phi_star_rec : (TMath::TwoPi() - phi_star_rec);
	double phi_star_rec_deg=180*phi_star_rec/PI;
	h_phistar_0->Fill(phi_star_rec_cp);
	
      }
      ++iEv;
    }
    h_phistar_0->Draw();
    ++nf;
  }
  inFile.close();
    
}
