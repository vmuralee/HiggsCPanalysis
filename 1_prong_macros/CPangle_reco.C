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

/*--------------------------Computing Refitted Vertex-------------------------------*/


Float_t refittedVertex(AC1B *t1,TString Cord,Int_t tau_index1,Int_t tau_index2){
 
  //variables
  // Float_t genPV_x=0;Float_t genPV_y=0;Float_t genPV_z=0;
  Float_t CrefitPV;
  Float_t Higgs_e=0;Float_t Higgs_px=0;Float_t Higgs_py=0;Float_t Higgs_pz=0;
  Int_t ele_index=0,mu_index=0;float kk=0.2;
  float* minDR;minDR=&kk;
  
  for(unsigned int igen =0;igen < t1->genparticles_count;igen++){            /*genparticle loop*/
    // if((t1->genparticles_pdgid[igen]==25||t1->genparticles_pdgid[igen]==35||t1->genparticles_pdgid[igen]==36)&&t1->genparticles_isLastCopy[igen]>0.5&&t1->genparticles_isPrompt[igen]>0.5){
    //   genPV_x = t1->genparticles_vx[igen];
    //   genPV_y = t1->genparticles_vy[igen];
    //   genPV_z = t1->genparticles_vz[igen];
    // }
    if(t1->genparticles_pdgid[igen]==25||t1->genparticles_pdgid[igen]==35||t1->genparticles_pdgid[igen]==36){
      Higgs_px= t1->genparticles_px[igen];
      Higgs_py= t1->genparticles_py[igen];
      Higgs_pz= t1->genparticles_pz[igen];
      Higgs_e = t1->genparticles_e[igen];
    }
    
  }       /*end of genparticle loop*/
  
  for(int itau=0;itau<t1->gentau_count;itau++){                               /*itau*/
    for(int itau1=itau+1;itau1<t1->gentau_count;itau1++){                     /*itau1*/
      if((t1->gentau_e[itau]+t1->gentau_e[itau1]==Higgs_e)&&(t1->gentau_px[itau]+t1->gentau_px[itau1]==Higgs_px)&&(t1->gentau_py[itau]+t1->gentau_py[itau1]==Higgs_py)&&(t1->gentau_pz[itau]+t1->gentau_pz[itau1]==Higgs_pz))
	{
	  for(int igen=0;igen<t1->genparticles_count;igen++){
	    for(int ie=0;ie<t1->electron_count;ie++){
	      if(GenMatch(t1,"e",ie,"no",igen,minDR)==3){
		ele_index=igen;
	      }
	      
	    }
	    for(int im=0;im<t1->muon_count;im++){
	      if(GenMatch(t1,"m",im,"no",igen,minDR)==4){
		mu_index=igen;
	      }
	    }
	  }
	  // for(int it=0;it<t1->tau_count;it++){
	  //   if(GenMatch(t1,"t",it,"yes",itau,minDR)==5){
	      
	  //     tau_index1=itau;
	  //   }
	  //   if(GenMatch(t1,"t",it,"yes",itau1,minDR)==5){
	  //     tau_index2=itau1;
	  //   }
	  // }
	  
	}
    }//end of itau
  }//end of itau1
  for(unsigned int i=0; i<t1->refitvertex_count; i++){
      
    if(((ele_index==t1->refitvertex_eleIndex[i][0]||ele_index==t1->refitvertex_eleIndex[i][1])||(mu_index==t1->refitvertex_muIndex[i][0]||mu_index==t1->refitvertex_muIndex[i][1]))||(tau_index1==t1->refitvertex_tauIndex[i][0]||tau_index1==t1->refitvertex_tauIndex[i][1])){
      if(((ele_index==t1->refitvertex_eleIndex[i][0]||ele_index==t1->refitvertex_eleIndex[i][1])||(mu_index==t1->refitvertex_muIndex[i][0]||mu_index==t1->refitvertex_muIndex[i][1]))||(tau_index2==t1->refitvertex_tauIndex[i][0]||tau_index2==t1->refitvertex_tauIndex[i][1])){
	if(Cord=="x"){    
	  CrefitPV=t1->refitvertex_x[i];}
	if(Cord=="y"){
	  CrefitPV=t1->refitvertex_y[i];}
	if(Cord=="z"){
	  CrefitPV=t1->refitvertex_z[i];}	   
      }
    }
  }
  
  return CrefitPV;
}
Float_t RecoVertex_with_bs(AC1B* t1,TString Cord){
  Float_t RecoVertex;
  TMatrixD cov_pv_Matrix(2,2);
  TMatrixD cov_bs_Matrix(2,2);
  TMatrixD cov_pv_InverseMatrix(2,2);
  TMatrixD cov_bs_InverseMatrix(2,2);
  TMatrixD x_pv(1,2);
  TMatrixD x_bs(1,2);
  TMatrixD x_pv_recal(1,2);

  cov_pv_Matrix(0,0)=t1->primvertex_cov[0];
  cov_pv_Matrix(0,1)=t1->primvertex_cov[1];
  cov_pv_Matrix(1,0)=t1->primvertex_cov[1];
  cov_pv_Matrix(1,1)=t1->primvertex_cov[3];

  cov_bs_Matrix(0,0)=t1->beamspot_cov[0];
  cov_bs_Matrix(0,1)=t1->beamspot_cov[1];
  cov_bs_Matrix(1,0)=t1->beamspot_cov[1];
  cov_bs_Matrix(1,1)=t1->beamspot_cov[3];

  x_pv(0,0)=t1->primvertex_x;
  x_pv(0,1)=t1->primvertex_y;

  x_bs(0,0)=t1->beamspot_x;
  x_bs(0,1)=t1->beamspot_y;

  cov_pv_InverseMatrix=cov_pv_Matrix.Invert();
  cov_bs_InverseMatrix=cov_bs_Matrix.Invert();

  TMatrixD cov_combined=cov_pv_InverseMatrix+cov_bs_InverseMatrix;
  cov_combined.Invert();
  x_pv_recal=(x_bs*cov_bs_InverseMatrix+x_pv*cov_pv_InverseMatrix)*cov_combined;

  if(Cord=="x"){
    RecoVertex=x_pv_recal(0,0);
  }
  if(Cord=="y"){
    RecoVertex=x_pv_recal(0,1);
  }
  if(Cord=="z"){
    RecoVertex=t1->primvertex_z;
  }
  return RecoVertex;
}
Float_t RefitVertex_with_bs(AC1B* t1,TString Cord,Int_t tau_index1,Int_t tau_index2){
  Float_t RefitVertex;
  TMatrixD cov_rv_Matrix(2,2);
  TMatrixD cov_bs_Matrix(2,2);
  TMatrixD cov_rv_InverseMatrix(2,2);
  TMatrixD cov_bs_InverseMatrix(2,2);
  TMatrixD x_rv(1,2);
  TMatrixD x_bs(1,2);
  TMatrixD x_rv_recal(1,2);
  for(unsigned int ir=0;ir<t1->refitvertex_count;ir++){
    //if((tau_index1==t1->refitvertex_tauIndex[ir][0]||tau_index1==t1->refitvertex_tauIndex[ir][1])){
    //if((tau_index2==t1->refitvertex_tauIndex[ir][0]||tau_index2==t1->refitvertex_tauIndex[ir][1])){
    cov_rv_Matrix(0,0)=t1->refitvertex_cov[ir][0];
    cov_rv_Matrix(0,1)=t1->refitvertex_cov[ir][1];
    cov_rv_Matrix(1,0)=t1->refitvertex_cov[ir][1];
    cov_rv_Matrix(1,1)=t1->refitvertex_cov[ir][3];
    if(cov_rv_Matrix(0,0)==0 && cov_rv_Matrix(0,1)==0 && cov_rv_Matrix(1,1)==0){
      cov_rv_Matrix(0,0)=t1->primvertex_cov[0];
      cov_rv_Matrix(0,1)=t1->primvertex_cov[1];
      cov_rv_Matrix(1,0)=t1->primvertex_cov[1];
      cov_rv_Matrix(1,1)=t1->primvertex_cov[3];    
    }
   
  }
  
  cov_bs_Matrix(0,0)=t1->beamspot_cov[0];
  cov_bs_Matrix(0,1)=t1->beamspot_cov[1];
  cov_bs_Matrix(1,0)=t1->beamspot_cov[1];
  cov_bs_Matrix(1,1)=t1->beamspot_cov[3];

  x_rv(0,0)=refittedVertex(t1,"x",tau_index1,tau_index2);
  x_rv(0,1)=refittedVertex(t1,"y",tau_index1,tau_index2);

  x_bs(0,0)=t1->beamspot_x;
  x_bs(0,1)=t1->beamspot_y;

  cov_rv_InverseMatrix=cov_rv_Matrix.Invert();
  cov_bs_InverseMatrix=cov_bs_Matrix.Invert();

  TMatrixD cov_combined=(cov_rv_InverseMatrix+cov_bs_InverseMatrix).Invert();
  x_rv_recal=(x_bs*cov_bs_InverseMatrix+x_rv*cov_rv_InverseMatrix)*cov_combined;

  if(Cord=="x"){
    RefitVertex=cov_rv_Matrix(0,0);//x_rv_recal(0,0);
  }
  if(Cord=="y"){
    RefitVertex=cov_rv_Matrix(0,1);//x_rv_recal(0,1);
  }
  if(Cord=="z"){
    RefitVertex=cov_rv_Matrix(1,1);//refittedVertex(t1,"z",tau_index1,tau_index2);
  }
  return RefitVertex;
}
void phi_rec_distr(){
  gROOT->LoadMacro("AC1B.h+");
  ifstream inFile;
 
  inFile.open("ggHToTauTau.txt");
  TString infilename_;
  int nf=0;
  while(!inFile.eof()){
    //for(int ifi=0;ifi<15;ifi++){
    cout<<"File No: "<<nf<<" Taking"<<endl;
    inFile >> infilename_;
    TString stri = to_string(nf);
    TFile* file = new TFile(infilename_);
    TTree* tree=(TTree*)file->Get("/makeroottree/AC1B");
    AC1B* t1 = new AC1B(tree,false);
    TFile* f = new TFile("ofile_"+stri+".root","recreate");
    TH1F* h_phistar = new TH1F("h_phistar","",20,0,6);
    Int_t TotEntries=t1->GetEntries(),iEv=0;
    while(iEv<TotEntries){
      float kk=0.2;
      float* minDR;minDR=&kk;
      Float_t pionPVtx_x=0,pionPVtx_y=0,pionPVtx_z=0;
      Float_t pionNVtx_x=0,pionNVtx_y=0,pionNVtx_z=0;
      Float_t pionP_px=0,pionP_py=0,pionP_pz=0,pionP_e=0;
      Float_t pionN_px=0,pionN_py=0,pionN_pz=0,pionN_e=0;
      t1->GetEntry(iEv);
    
      //for(int gtau=0;gtau<t1->gentau_count;gtau++){
      for(int rtau=0;rtau<t1->tau_count;rtau++){
	if(t1->tau_byMediumIsolationMVArun2v1DBoldDMwLT[rtau]>0.5&&t1->tau_againstMuonTight3[rtau]>0.5 && t1->tau_againstElectronTightMVA6[rtau]>0.5 && t1->tau_pt[rtau]>30.){
	  if(t1->tau_charge[rtau]>0&&t1->tau_decayMode[rtau]==0){
	    pionPVtx_x=t1->tau_pca3D_x[rtau];
	    pionPVtx_y=t1->tau_pca3D_y[rtau];
	    pionPVtx_z=t1->tau_pca3D_z[rtau];
	    pionP_px=t1->tau_leadchargedhadrcand_px[rtau];
	    pionP_py=t1->tau_leadchargedhadrcand_py[rtau];
	    pionP_pz=t1->tau_leadchargedhadrcand_pz[rtau];
	    pionP_e=sqrt(pionP_px*pionP_px+pionP_py*pionP_py+pionP_pz*pionP_pz+(t1->tau_leadchargedhadrcand_mass[rtau]*t1->tau_leadchargedhadrcand_mass[rtau]));
	  }//tau+
	    
	
	
	  if(t1->tau_charge[rtau]<0 &&t1->tau_decayMode[rtau]==0){
	    pionNVtx_x=t1->tau_pca3D_x[rtau];
	    pionNVtx_y=t1->tau_pca3D_y[rtau];
	    pionNVtx_z=t1->tau_pca3D_z[rtau];
	    pionN_px=t1->tau_leadchargedhadrcand_px[rtau];
	    pionN_py=t1->tau_leadchargedhadrcand_py[rtau];
	    pionN_pz=t1->tau_leadchargedhadrcand_pz[rtau];
	    pionN_e=sqrt(pionN_px*pionN_px+pionN_py*pionN_py+pionN_pz*pionN_pz+(t1->tau_leadchargedhadrcand_mass[rtau]*t1->tau_leadchargedhadrcand_mass[rtau]));
	  }
	}//tau-
	//}
      }
      //cout<<"      PV  "<<t1->primvertex_x<<endl;
      //cout<<"Recal PV  "<<RecoVertex_with_bs(t1,"x")<<endl;
      LV pionP_lab(pionP_px,pionP_py,pionP_pz,pionP_e);
      LV pionN_lab(pionN_px,pionN_py,pionN_pz,pionN_e);
      if(pionP_lab.pt()>0 && pionN_lab.pt()>0){
	//cout<<"pi+_e "<<pionP_e<<"  pi-_e "<<pionN_e<<endl;
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
	double OstarCP_gen = pionN_rf_2p_3d_u.Dot(ipvP_rf_rec_3d_t_u.Cross(ipvN_rf_rec_3d_t_u));
	double phi_star_gen_cp = (OstarCP_gen >= 0) ? phi_star_rec : (TMath::TwoPi() - phi_star_rec);
	double phi_star_rec_deg=180*phi_star_rec/PI;
	h_phistar->Fill(phi_star_gen_cp);
      }
      ++iEv;
    }
    double norm=1/h_phistar->Integral();
    h_phistar->Scale(norm);
    h_phistar->Draw();
    f->Write();
    ++nf;
  }
  inFile.close();
    
}
void phi_refit_distr(){
  gROOT->LoadMacro("AC1B.h+");
  ifstream inFile;
  inFile.open("ggHToTauTau.txt");
  TString infilename_;
  int nf=0;
  while(!inFile.eof()){
    cout<<"File No: "<<nf<<" Taking"<<endl;
    inFile >> infilename_;
    TString stri = to_string(nf);
    TFile* file = new TFile(infilename_);
    TTree* tree=(TTree*)file->Get("/makeroottree/AC1B");
    AC1B *t1 = new AC1B(tree,false);
    TFile* f = new TFile("ofile_"+stri+".root","recreate");
    TH1F* h_phistar = new TH1F("h_phistar","",20,0,6);
    Int_t TotEntries=t1->GetEntries(),iEv=0;
    while(iEv<TotEntries){
      float kk=0.2;
      float* minDR;minDR=&kk;
      Int_t tau_indexP=0,tau_indexN=0;
      Float_t refitVtx_x=0,refitVtx_y=0,refitVtx_z=0;
      Float_t pionPVtx_x=0,pionPVtx_y=0,pionPVtx_z=0;
      Float_t pionNVtx_x=0,pionNVtx_y=0,pionNVtx_z=0;
      Float_t pionP_px=0,pionP_py=0,pionP_pz=0,pionP_e=0;
      Float_t pionN_px=0,pionN_py=0,pionN_pz=0,pionN_e=0;
      t1->GetEntry(iEv);
   
      for(int gtau=0;gtau<t1->gentau_count;gtau++){
	for(int rtau=0;rtau<t1->tau_count;rtau++){
	  if(t1->tau_byMediumIsolationMVArun2v1DBoldDMwLT[rtau]>0.5&&t1->tau_againstMuonTight3[rtau]>0.5 && t1->tau_againstElectronTightMVA6[rtau]>0.5){
	    if(t1->tau_charge[rtau]>0 && t1->tau_decayMode[rtau]==0.0){
	      if(GenMatch(t1,"t",rtau,"yes",gtau,minDR)==5){
		pionPVtx_x=t1->tau_pca3D_x[rtau];
		pionPVtx_y=t1->tau_pca3D_y[rtau];
		pionPVtx_z=t1->tau_pca3D_z[rtau];
		pionP_px=t1->tau_leadchargedhadrcand_px[rtau];
		pionP_py=t1->tau_leadchargedhadrcand_py[rtau];
		pionP_pz=t1->tau_leadchargedhadrcand_pz[rtau];
		pionP_e=sqrt(pionP_px*pionP_px+pionP_py*pionP_py+pionP_pz*pionP_pz+(t1->tau_leadchargedhadrcand_mass[rtau]*t1->tau_leadchargedhadrcand_mass[rtau]));
		tau_indexP=rtau;
		}
	    }//tau+
	
	    if(t1->tau_charge[rtau]<0 && t1->tau_decayMode[rtau]==0){
	      if(GenMatch(t1,"t",rtau,"yes",gtau,minDR)==5){
		pionNVtx_x=t1->tau_pca3D_x[rtau];
		pionNVtx_y=t1->tau_pca3D_y[rtau];
		pionNVtx_z=t1->tau_pca3D_z[rtau];
		pionN_px=t1->tau_leadchargedhadrcand_px[rtau];
		pionN_py=t1->tau_leadchargedhadrcand_py[rtau];
		pionN_pz=t1->tau_leadchargedhadrcand_pz[rtau];
		pionN_e=sqrt(pionN_px*pionN_px+pionN_py*pionN_py+pionN_pz*pionN_pz+(t1->tau_leadchargedhadrcand_mass[rtau]*t1->tau_leadchargedhadrcand_mass[rtau]));

		tau_indexN=rtau;
		}
	    }//tau-
		
	  }
	  // double mini=0.00001;
	  //if(abs(RefitVertex_with_bs(t1,"x",tau_indexP,tau_indexN))>mini&&abs(RefitVertex_with_bs(t1,"y",tau_indexP,tau_indexN))>mini&&abs(RefitVertex_with_bs(t1,"z",tau_indexP,tau_indexN))){
	    
	    //}
	  // double mini=0.00001;
	  
	}// tau loop
	if(abs(refittedVertex(t1,"x",tau_indexP,tau_indexN))!=0&&abs(refittedVertex(t1,"y",tau_indexP,tau_indexN))!=0&&abs(refittedVertex(t1,"z",tau_indexP,tau_indexN))!=0){
	  
	  refitVtx_x=/*RefitVertex_with_bs(t1,"x",tau_indexP,tau_indexN);*/refittedVertex(t1,"x",tau_indexP,tau_indexN);
	  refitVtx_y=/*RefitVertex_with_bs(t1,"y",tau_indexP,tau_indexN);*/refittedVertex(t1,"y",tau_indexP,tau_indexN);
	  refitVtx_z=/*RefitVertex_with_bs(t1,"z",tau_indexP,tau_indexN);*/refittedVertex(t1,"z",tau_indexP,tau_indexN);
        }  
	else{
	  refitVtx_x=t1->primvertex_x;
	  refitVtx_y=t1->primvertex_y;
	  refitVtx_z=t1->primvertex_z;
	}
	  
      }
     
      //cout<<"cov_xx  "<<refitVtx_x<<endl;
      //cout<<"cov_xy  "<<refitVtx_y<<endl;
      //cout<<"cov_yy  "<<refitVtx_z<<endl;
      LV pionP_lab(pionP_px,pionP_py,pionP_pz,pionP_e);
      LV pionN_lab(pionN_px,pionN_py,pionN_pz,pionN_e);
      if(pionP_lab.pt()>0 && pionN_lab.pt()>0){
	//compute phi_star in pi-pi+ rf.
	LV pion_pair_lab=pionP_lab+pionN_lab;
	ROOT::Math::Boost boost_to_rf_refit(pion_pair_lab.BoostToCM());
	//boost to rf of two pion system
	LV pionP_rf_2p=boost_to_rf_refit(pionP_lab);
	LV pionN_rf_2p=boost_to_rf_refit(pionN_lab);
	//get unit momentum vectors
	PV pionP_lab_3d=pionP_lab.Vect();
	PV pionN_lab_3d=pionN_lab.Vect();
	PV pionP_lab_3d_u=pionP_lab_3d/sqrt(pionP_lab_3d.mag2());
	PV pionN_lab_3d_u=pionN_lab_3d/sqrt(pionN_lab_3d.mag2());
	//Get PCAs, by extrapolating the line represented by pion vertex and pion momentum
	double tP=pionP_lab_3d_u.x()*(refitVtx_x-pionPVtx_x)+pionP_lab_3d_u.y()*(refitVtx_y-pionPVtx_y)+pionP_lab_3d_u.z()*(refitVtx_z-pionPVtx_z);
	double tN=pionN_lab_3d_u.x()*(refitVtx_x-pionNVtx_x)+pionN_lab_3d_u.y()*(refitVtx_y-pionNVtx_y)+pionN_lab_3d_u.z()*(refitVtx_z-pionNVtx_z);
      
	Point3D pcaP(pionPVtx_x+pionP_lab_3d_u.x()*tP,pionPVtx_y+pionP_lab_3d_u.y()*tP,pionPVtx_z+pionP_lab_3d_u.z()*tP);
	Point3D pcaN(pionNVtx_x+pionN_lab_3d_u.x()*tN,pionNVtx_y+pionN_lab_3d_u.y()*tN,pionNVtx_z+pionN_lab_3d_u.z()*tN);
	//Get the normalized IP vectors
	PV ipvP_lab_refit_3d(pcaP.x()-refitVtx_x,pcaP.y()-refitVtx_y,pcaP.z()-refitVtx_z);
	PV ipvN_lab_refit_3d(pcaN.x()-refitVtx_x,pcaN.y()-refitVtx_y,pcaN.z()-refitVtx_z);
      
	LV ipvP_lab_refit(ipvP_lab_refit_3d.x()/sqrt(ipvP_lab_refit_3d.mag2()),ipvP_lab_refit_3d.y()/sqrt(ipvP_lab_refit_3d.mag2()),ipvP_lab_refit_3d.z()/sqrt(ipvP_lab_refit_3d.mag2()),0);
	LV ipvN_lab_refit(ipvN_lab_refit_3d.x()/sqrt(ipvN_lab_refit_3d.mag2()),ipvN_lab_refit_3d.y()/sqrt(ipvN_lab_refit_3d.mag2()),ipvN_lab_refit_3d.z()/sqrt(ipvN_lab_refit_3d.mag2()),0);
      
	//boost to ZMF
	LV ipvP_rf_refit = boost_to_rf_refit(ipvP_lab_refit);
	LV ipvN_rf_refit = boost_to_rf_refit(ipvN_lab_refit);
      
	//Get only the position component of the vector
	PV ipvP_rf_refit_3d = ipvP_rf_refit.Vect();
	PV ipvN_rf_refit_3d = ipvN_rf_refit.Vect();
	PV pionP_rf_2p_3d = pionP_rf_2p.Vect();
	PV pionN_rf_2p_3d = pionN_rf_2p.Vect();
      
	//Get the unit (normalized) vector along pion momentum
	PV pionP_rf_2p_3d_u = pionP_rf_2p_3d/TMath::Sqrt(pionP_rf_2p_3d.mag2());
	PV pionN_rf_2p_3d_u = pionN_rf_2p_3d/TMath::Sqrt(pionN_rf_2p_3d.mag2());
	//Get the longitudinal component of IP vector parallel to pion momenta
	PV ipvP_rf_refit_3d_l = ipvP_rf_refit_3d.Dot(pionP_rf_2p_3d_u)*pionP_rf_2p_3d_u;
	PV ipvN_rf_refit_3d_l = ipvN_rf_refit_3d.Dot(pionN_rf_2p_3d_u)*pionN_rf_2p_3d_u;
	//Get IP vector normal to pion momenta
	PV ipvP_rf_refit_3d_t = ipvP_rf_refit_3d - ipvP_rf_refit_3d_l;
	PV ipvN_rf_refit_3d_t = ipvN_rf_refit_3d - ipvN_rf_refit_3d_l;
	//Get normalized normal IP vector
	PV ipvP_rf_refit_3d_t_u = ipvP_rf_refit_3d_t/TMath::Sqrt(ipvP_rf_refit_3d_t.mag2());
	PV ipvN_rf_refit_3d_t_u = ipvN_rf_refit_3d_t/TMath::Sqrt(ipvN_rf_refit_3d_t.mag2());
      
	//Get CP angle in ZMF
	double phi_star_refit = TMath::ACos(ipvP_rf_refit_3d_t_u.Dot(ipvN_rf_refit_3d_t_u));
	double OstarCP_gen = pionN_rf_2p_3d_u.Dot(ipvP_rf_refit_3d_t_u.Cross(ipvN_rf_refit_3d_t_u));
	double phi_star_gen_cp = (OstarCP_gen >= 0) ? phi_star_refit : (TMath::TwoPi() - phi_star_refit);
	double phi_star_refit_deg=180*phi_star_refit/PI;
	h_phistar->Fill(phi_star_gen_cp);
      }
      ++iEv;
    }
    double norm=1/h_phistar->Integral();
    h_phistar->Scale(norm);
    h_phistar->Draw();
    f->Write();
    ++nf;
    
  }
  inFile.close();
}

