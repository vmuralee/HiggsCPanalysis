#include<iostream> 
#include "AC1B.h"
#include "Math/GenVector/Boost.h"
#include <TMath.h>
#include <TROOT.h>
#include "Math/Point3D.h"
#include "Math/LorentzVector.h"
#include<string>
#include <fstream>
typedef ROOT::Math::XYZPointD Point3D;
typedef ROOT::Math::XYZVectorD PV;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LV;
#define  PI 3.14159265358979312e+00

float deltaR(float eta_1, float phi_1, float eta_2, float phi_2){
  float deta = fabs(eta_1-eta_2);
  float dphi = fabs(phi_1-phi_2);
  if(dphi > PI)
    dphi = 2*PI -dphi;
  return abs(sqrt(deta*deta+dphi*dphi));
}

float GenPion(AC1B* t1,int ig1,string Cord){
  float genpion=-6;
  float pion_px=0,pion_py=0,pion_pz=0,pion_e=0;
  if(t1->genparticles_pdgid[ig1]==211 || t1->genparticles_pdgid[ig1]==-211){
    float _px = t1->genparticles_px[ig1];
    float _py = t1->genparticles_py[ig1];
    float _pz = t1->genparticles_pz[ig1];
    float _e  = t1->genparticles_e[ig1];
    
    LV pion(_px,_py,_pz,_e);float mindR=0.03;
    float pion_eta=pion.Eta();float pion_phi = pion.Phi();
    for(unsigned it=0;it<t1->gentau_count;it++){
      if(t1->gentau_decayMode[it]!=1)continue;
      LV gentau(t1->gentau_px[it],t1->gentau_py[it],t1->gentau_pz[it],t1->gentau_e[it]);
      float gentau_eta = gentau.Eta();
      float gentau_phi = gentau.Phi();
      float dR = deltaR(gentau_eta,gentau_phi,pion_eta,pion_phi);
      if(dR<mindR){
	mindR=dR;
	pion_px=t1->gentau_visible_px[it]-_px;
	pion_py=t1->gentau_visible_py[it]-_py;
	pion_pz=t1->gentau_visible_pz[it]-_pz;
	pion_e =t1->gentau_visible_e[it]-_e;
      }
    }
  }
  if(Cord=="px")genpion=pion_px;
  if(Cord=="py")genpion=pion_py;
  if(Cord=="pz")genpion=pion_pz;
  if(Cord=="e") genpion=pion_e;
  return genpion;
}
void phi_gen_distr(){
  gROOT->LoadMacro("AC1B.h+");
  ifstream inFile;
  inFile.open("file_list.txt");
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
    TH1F* h_phistar_plus = new TH1F("h_phistar_plus","",20,0,6);
    TH1F* h_phistar_neg = new TH1F("h_phistar_neg","",20,0,6);
    Int_t TotEntries=t1->GetEntries(),iEv=0;
    while(iEv<TotEntries){
      Float_t genPV_x=0;Float_t genPV_y=0;Float_t genPV_z=0;
      Float_t rhoP_px=0,rhoP_py=0,rhoP_pz=0,rhoP_e=0,rhoP_m=0;
      Float_t rhoN_px=0,rhoN_py=0,rhoN_pz=0,rhoN_e=0,rhoN_m=0;
      Float_t genpion0P_px_lab=0,genpion0P_py_lab=0,genpion0P_pz_lab=0,genpion0P_e_lab=0;
      Float_t genpion0N_px_lab=0,genpion0N_py_lab=0,genpion0N_pz_lab=0,genpion0N_e_lab=0;
      Float_t genPiN_px_lab=0,genPiN_py_lab=0,genPiN_pz_lab=0,genPiN_e_lab=0;
      Float_t genpionNVtx_x=0,genpionNVtx_y=0,genpionNVtx_z=0;
      Float_t genPiP_px_lab=0,genPiP_py_lab=0,genPiP_pz_lab=0,genPiP_e_lab=0;
      Float_t genpionPVtx_x=0,genpionPVtx_y=0,genpionPVtx_z=0;
      t1->GetEntry(iEv);
      for(unsigned int igen2=0;igen2<t1->genparticles_count;igen2++){
	if(t1->genparticles_pdgid[igen2]==25&&t1->genparticles_pdgid[igen2]==36&&t1->genparticles_isLastCopy[igen2]>0.5&&t1->genparticles_isPrompt[igen2]>0.5){
	  genPV_x = t1->genparticles_vx[igen2];
	  genPV_y = t1->genparticles_vy[igen2];
	  genPV_z = t1->genparticles_vz[igen2];
	
	}
	LV genpion(GenPion(t1,igen2,"px"),GenPion(t1,igen2,"py"),GenPion(t1,igen2,"pz"),GenPion(t1,igen2,"e"));
	for(unsigned int igen1=0;igen1<t1->genparticles_count;igen1++){
	  if(t1->genparticles_pdgid[igen1]!=211)continue;
	  for(unsigned int itauP=0;itauP<t1->gentau_count;itauP++){
	    if(t1->gentau_charge[itauP]==-1)continue;
	    if(t1->gentau_visible_px[itauP]==(t1->genparticles_px[igen1]+genpion.x()) && t1->gentau_visible_py[itauP]==t1->genparticles_py[igen1]+genpion.y()&&t1->gentau_visible_pz[itauP]==t1->genparticles_pz[igen1]+genpion.z()&&t1->gentau_visible_e[itauP]==t1->genparticles_e[igen1]+genpion.t()){
	      genPiP_px_lab=t1->genparticles_px[igen1];
	      genPiP_py_lab=t1->genparticles_py[igen1];
	      genPiP_pz_lab=t1->genparticles_pz[igen1];
	      genPiP_e_lab=t1->genparticles_e[igen1];
	      genpionPVtx_x=t1->genparticles_vx[igen1];
	      genpionPVtx_y=t1->genparticles_vy[igen1];
	      genpionPVtx_z=t1->genparticles_vz[igen1];
	      genpion0P_px_lab=genpion.x();
	      genpion0P_py_lab=genpion.y();
	      genpion0P_pz_lab=genpion.z();
	      genpion0P_e_lab= genpion.t();
	    }
	  }
	
	}//tau+
	for(unsigned int igen1=0;igen1<t1->genparticles_count;igen1++){
	  if(t1->genparticles_pdgid[igen1]!=-211)continue;
	  for(unsigned int itauN=0;itauN<t1->gentau_count;itauN++){
	    if(t1->gentau_charge[itauN]==1)continue;
	    if(t1->gentau_visible_px[itauN]==(t1->genparticles_px[igen1]+genpion.x())&& t1->gentau_visible_py[itauN]==t1->genparticles_py[igen1]+genpion.y()&&t1->gentau_visible_pz[itauN]==t1->genparticles_pz[igen1]+genpion.z()&&t1->gentau_visible_e[itauN]==t1->genparticles_e[igen1]+genpion.t()){
	      genPiN_px_lab=t1->genparticles_px[igen1];
	      genPiN_py_lab=t1->genparticles_py[igen1];
	      genPiN_pz_lab=t1->genparticles_pz[igen1];
	      genPiN_e_lab=t1->genparticles_e[igen1];
	      genpionNVtx_x=t1->genparticles_vx[igen1];
	      genpionNVtx_y=t1->genparticles_vy[igen1];
	      genpionNVtx_z=t1->genparticles_vz[igen1];
	      genpion0N_px_lab=genpion.x();
	      genpion0N_py_lab=genpion.y();
	      genpion0N_pz_lab=genpion.z();
	      genpion0N_e_lab=genpion.t();
	    }
	  
	  }
	}//tau-
      }//genpion loop
      LV genpionP_lab(genPiP_px_lab,genPiP_py_lab,genPiP_pz_lab,genPiP_e_lab);
      LV genpionN_lab(genPiN_px_lab,genPiN_py_lab,genPiN_pz_lab,genPiN_e_lab);
      LV genpion0P_lab(genpion0P_px_lab,genpion0P_py_lab,genpion0P_pz_lab,genpion0P_e_lab);
      LV genpion0N_lab(genpion0N_px_lab,genpion0N_py_lab,genpion0N_pz_lab,genpion0N_e_lab);
      
      double y_plus=(genpionP_lab.e()-genpion0P_lab.e())/(genpionP_lab.e()+genpion0P_lab.e());
      double y_neg =(genpionN_lab.e()-genpion0N_lab.e())/(genpionN_lab.e()+genpion0N_lab.e());
      if(genpionP_lab.pt() > 0 && genpionN_lab.pt() > 0 && genpion0P_lab.pt() > 0 && genpion0N_lab.pt() > 0 ){
	//cout<<genpion0P_lab.pt()<<endl;
      	if(y_plus*y_neg >0){
	  LV gen_pion_pair_lab = genpionP_lab + genpionN_lab;
	  ROOT::Math::Boost boost_to_rf_gen(gen_pion_pair_lab.BoostToCM());
      
	  //boost to rf of two pion system
	  LV genpionP_rf_2p = boost_to_rf_gen(genpionP_lab);
	  LV genpionN_rf_2p = boost_to_rf_gen(genpionN_lab);
	  cout<<genpionP_rf_2p.pt()<<endl;
	  //get unit momentum vectors
	  PV genpionP_lab_3d = genpionP_lab.Vect();
	  PV genpionN_lab_3d = genpionN_lab.Vect();
	  PV genpionP_lab_3d_u = genpionP_lab_3d/sqrt(genpionP_lab_3d.mag2());
	  PV genpionN_lab_3d_u = genpionN_lab_3d/sqrt(genpionN_lab_3d.mag2());
      
	  //Get PCAs, by extrapolating the line represented by pion vertex and pion momentum
	  double tP = genpionP_lab_3d_u.x()*(genPV_x-genpionPVtx_x)+genpionP_lab_3d_u.y()*(genPV_y-genpionPVtx_y)+genpionP_lab_3d_u.z()*(genPV_z-genpionPVtx_z);
	  double tN = genpionN_lab_3d_u.x()*(genPV_x-genpionNVtx_x)+genpionN_lab_3d_u.y()*(genPV_y-genpionNVtx_y)+genpionN_lab_3d_u.z()*(genPV_z-genpionNVtx_z);
	  Point3D pcaP(genpionPVtx_x + genpionP_lab_3d_u.x()*tP, genpionPVtx_y + genpionP_lab_3d_u.y()*tP, genpionPVtx_z + genpionP_lab_3d_u.z()*tP);
	  Point3D pcaN(genpionNVtx_x + genpionN_lab_3d_u.x()*tN, genpionNVtx_y + genpionN_lab_3d_u.y()*tN, genpionNVtx_z + genpionN_lab_3d_u.z()*tN);
      
	  //Get the normalized IP vectors
	  PV ipvP_lab_gen_3d(pcaP.x() - genPV_x, pcaP.y() - genPV_y, pcaP.z() - genPV_z);
	  PV ipvN_lab_gen_3d(pcaN.x() - genPV_x, pcaN.y() - genPV_y, pcaN.z() - genPV_z);
	  LV ipvP_lab_gen(ipvP_lab_gen_3d.x()/sqrt(ipvP_lab_gen_3d.mag2()), ipvP_lab_gen_3d.y()/sqrt(ipvP_lab_gen_3d.mag2()), ipvP_lab_gen_3d.z()/sqrt(ipvP_lab_gen_3d.mag2()), 0);
	  LV ipvN_lab_gen(ipvN_lab_gen_3d.x()/sqrt(ipvN_lab_gen_3d.mag2()), ipvN_lab_gen_3d.y()/sqrt(ipvN_lab_gen_3d.mag2()), ipvN_lab_gen_3d.z()/sqrt(ipvN_lab_gen_3d.mag2()), 0);
	  //cout<<ipvP_lab_gen_3d.r()<<endl;
	  //boost to ZMF
	  LV ipvP_rf_gen = boost_to_rf_gen(ipvP_lab_gen);
	  LV ipvN_rf_gen = boost_to_rf_gen(ipvN_lab_gen);
      
	  //Get only the position component of the vector
	  PV ipvP_rf_gen_3d = ipvP_rf_gen.Vect();
	  PV ipvN_rf_gen_3d = ipvN_rf_gen.Vect();
      
	  PV genpionP_rf_2p_3d = genpionP_rf_2p.Vect();
	  PV genpionN_rf_2p_3d = genpionN_rf_2p.Vect();
      
	  //Get the unit (normalized) vector along pion momentum
	  PV genpionP_rf_2p_3d_u = genpionP_rf_2p_3d/TMath::Sqrt(genpionP_rf_2p_3d.mag2());
	  PV genpionN_rf_2p_3d_u = genpionN_rf_2p_3d/TMath::Sqrt(genpionN_rf_2p_3d.mag2());
      
	  //Get the longitudinal component of IP vector parallel to pion momenta
	  PV ipvP_rf_gen_3d_l = ipvP_rf_gen_3d.Dot(genpionP_rf_2p_3d_u)*genpionP_rf_2p_3d_u;
	  PV ipvN_rf_gen_3d_l = ipvN_rf_gen_3d.Dot(genpionN_rf_2p_3d_u)*genpionN_rf_2p_3d_u;
     
	  //Get IP vector normal to pion momenta
	  PV ipvP_rf_gen_3d_t = ipvP_rf_gen_3d - ipvP_rf_gen_3d_l;
	  PV ipvN_rf_gen_3d_t = ipvN_rf_gen_3d - ipvN_rf_gen_3d_l;
      
	  //Get normalized normal IP vector
	  PV ipvP_rf_gen_3d_t_u = ipvP_rf_gen_3d_t/TMath::Sqrt(ipvP_rf_gen_3d_t.mag2());
	  PV ipvN_rf_gen_3d_t_u = ipvN_rf_gen_3d_t/TMath::Sqrt(ipvN_rf_gen_3d_t.mag2());
      
	  //Get CP angle in ZMF
	  double phi_star_gen = acos(ipvP_rf_gen_3d_t_u.Dot(ipvN_rf_gen_3d_t_u));
	  double phi_star_deg=180*phi_star_gen/PI;
	  double OstarCP_gen = genpionN_rf_2p_3d_u.Dot(ipvP_rf_gen_3d_t_u.Cross(ipvN_rf_gen_3d_t_u));
	  double phi_star_gen_cp = (OstarCP_gen >= 0) ? phi_star_gen : (TMath::TwoPi() - phi_star_gen);
	  h_phistar_plus->Fill(phi_star_gen_cp);

      	  }
      	if(y_plus*y_neg <0){
	  LV gen_pion_pair_lab = genpionP_lab + genpionN_lab;
	  ROOT::Math::Boost boost_to_rf_gen(gen_pion_pair_lab.BoostToCM());
      
	  //boost to rf of two pion system
	  LV genpionP_rf_2p = boost_to_rf_gen(genpionP_lab);
	  LV genpionN_rf_2p = boost_to_rf_gen(genpionN_lab);
	  //cout<<genpionP_rf_2p.pt()<<endl;
	  //get unit momentum vectors
	  PV genpionP_lab_3d = genpionP_lab.Vect();
	  PV genpionN_lab_3d = genpionN_lab.Vect();
	  PV genpionP_lab_3d_u = genpionP_lab_3d/sqrt(genpionP_lab_3d.mag2());
	  PV genpionN_lab_3d_u = genpionN_lab_3d/sqrt(genpionN_lab_3d.mag2());
      
	  //Get PCAs, by extrapolating the line represented by pion vertex and pion momentum
	  double tP = genpionP_lab_3d_u.x()*(genPV_x-genpionPVtx_x)+genpionP_lab_3d_u.y()*(genPV_y-genpionPVtx_y)+genpionP_lab_3d_u.z()*(genPV_z-genpionPVtx_z);
	  double tN = genpionN_lab_3d_u.x()*(genPV_x-genpionNVtx_x)+genpionN_lab_3d_u.y()*(genPV_y-genpionNVtx_y)+genpionN_lab_3d_u.z()*(genPV_z-genpionNVtx_z);
	  Point3D pcaP(genpionPVtx_x + genpionP_lab_3d_u.x()*tP, genpionPVtx_y + genpionP_lab_3d_u.y()*tP, genpionPVtx_z + genpionP_lab_3d_u.z()*tP);
	  Point3D pcaN(genpionNVtx_x + genpionN_lab_3d_u.x()*tN, genpionNVtx_y + genpionN_lab_3d_u.y()*tN, genpionNVtx_z + genpionN_lab_3d_u.z()*tN);
      
	  //Get the normalized IP vectors
	  PV ipvP_lab_gen_3d(pcaP.x() - genPV_x, pcaP.y() - genPV_y, pcaP.z() - genPV_z);
	  PV ipvN_lab_gen_3d(pcaN.x() - genPV_x, pcaN.y() - genPV_y, pcaN.z() - genPV_z);
	  LV ipvP_lab_gen(ipvP_lab_gen_3d.x()/sqrt(ipvP_lab_gen_3d.mag2()), ipvP_lab_gen_3d.y()/sqrt(ipvP_lab_gen_3d.mag2()), ipvP_lab_gen_3d.z()/sqrt(ipvP_lab_gen_3d.mag2()), 0);
	  LV ipvN_lab_gen(ipvN_lab_gen_3d.x()/sqrt(ipvN_lab_gen_3d.mag2()), ipvN_lab_gen_3d.y()/sqrt(ipvN_lab_gen_3d.mag2()), ipvN_lab_gen_3d.z()/sqrt(ipvN_lab_gen_3d.mag2()), 0);
	  //cout<<ipvP_lab_gen_3d.r()<<endl;
	  //boost to ZMF
	  LV ipvP_rf_gen = boost_to_rf_gen(ipvP_lab_gen);
	  LV ipvN_rf_gen = boost_to_rf_gen(ipvN_lab_gen);
      
	  //Get only the position component of the vector
	  PV ipvP_rf_gen_3d = ipvP_rf_gen.Vect();
	  PV ipvN_rf_gen_3d = ipvN_rf_gen.Vect();
      
	  PV genpionP_rf_2p_3d = genpionP_rf_2p.Vect();
	  PV genpionN_rf_2p_3d = genpionN_rf_2p.Vect();
      
	  //Get the unit (normalized) vector along pion momentum
	  PV genpionP_rf_2p_3d_u = genpionP_rf_2p_3d/TMath::Sqrt(genpionP_rf_2p_3d.mag2());
	  PV genpionN_rf_2p_3d_u = genpionN_rf_2p_3d/TMath::Sqrt(genpionN_rf_2p_3d.mag2());
      
	  //Get the longitudinal component of IP vector parallel to pion momenta
	  PV ipvP_rf_gen_3d_l = ipvP_rf_gen_3d.Dot(genpionP_rf_2p_3d_u)*genpionP_rf_2p_3d_u;
	  PV ipvN_rf_gen_3d_l = ipvN_rf_gen_3d.Dot(genpionN_rf_2p_3d_u)*genpionN_rf_2p_3d_u;
     
	  //Get IP vector normal to pion momenta
	  PV ipvP_rf_gen_3d_t = ipvP_rf_gen_3d - ipvP_rf_gen_3d_l;
	  PV ipvN_rf_gen_3d_t = ipvN_rf_gen_3d - ipvN_rf_gen_3d_l;
      
	  //Get normalized normal IP vector
	  PV ipvP_rf_gen_3d_t_u = ipvP_rf_gen_3d_t/TMath::Sqrt(ipvP_rf_gen_3d_t.mag2());
	  PV ipvN_rf_gen_3d_t_u = ipvN_rf_gen_3d_t/TMath::Sqrt(ipvN_rf_gen_3d_t.mag2());
      
	  //Get CP angle in ZMF
	  double phi_star_gen = acos(ipvP_rf_gen_3d_t_u.Dot(ipvN_rf_gen_3d_t_u));
	  double phi_star_deg=180*phi_star_gen/PI;
	  double OstarCP_gen = genpionN_rf_2p_3d_u.Dot(ipvP_rf_gen_3d_t_u.Cross(ipvN_rf_gen_3d_t_u));
	  double phi_star_gen_cp = (OstarCP_gen >= 0) ? phi_star_gen : (TMath::TwoPi() - phi_star_gen);
	  h_phistar_neg->Fill(phi_star_gen_cp);
	  
      	  }
      }
      ++iEv;
    }//Event loop
    Double_t norm = 1;
    Double_t scale = norm/(h_phistar_plus->Integral());
    Double_t scale1 = norm/(h_phistar_neg->Integral());
    h_phistar_plus->Scale(scale);
    h_phistar_neg->Scale(scale1);
    h_phistar_plus->Draw();
    h_phistar_neg->Draw();
    f->Write();
    ++nf;
  }// file close
  inFile.close();
}
