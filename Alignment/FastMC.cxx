/////////////////////////////////////////////////
// Fast MC for MVTX geometry and hits generation
//    Multiple Scattering + Single Hit Resolution
////////////////////////////////////////////////
#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class PlotFile;
#endif

//#ifndef __CINT__
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>

#include "math.h"
#include "string.h"

#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TSystem.h"
#include "TUnixSystem.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TGeoMatrix.h"
//#endif

#include "mTree.h"
#include "StPhysicalHelix.h"
#include "SystemOfUnits.h"
#include "PhysicalConstants.h"

// Stave Geometry
const Double_t ZMAX = 27.1/2.;   // Total Length 2*ZMax
const Double_t XMAX = 1.5/2.;    // Total Width 2*XMAX
const Int_t NL = 3;
const Int_t NSMAX = 20;
const Int_t NS[NL] = {12, 16, 20};
const Double_t R0[NL] = {2.46, 3.23, 4.00};  // radial middle point
const Double_t PitchSize = 0.0028;   //

bool onStave(TVector3 pos)
{
  return fabs(pos.x())<XMAX && fabs(pos.y())<1.e-4 && fabs(pos.z())<ZMAX;
}

double theta0(double p, double m, double dx)
{
  double beta = p/TMath::Sqrt(p*p+m*m);
  return 0.0136/beta/p*TMath::Sqrt(dx)*(1.+0.038*TMath::Log(dx/beta/beta));
}

int main(int argc, char **argv)
{
  const Int_t Nevt = 100000;

  std::cout << " Reading in Stave Geometry Tables ... " << std::endl;
  TGeoHMatrix *stave[NL][NSMAX];
  TFile *fin = new TFile("StaveGeoMatrix_Ideal.root");
  for(int i=0;i<NL;i++) {
    for(int j=0;j<NS[i];j++) {
      std::cout << Form("geoM_%d_%d",i,j) << std::endl;
      stave[i][j] = (TGeoHMatrix *)fin->Get(Form("geoM_%d_%d",i,j));
      stave[i][j]->Print();
    }
  }
  fin->Close();

  TFile f("test.root","RECREATE");
  MVTXTREE mT;
  TTree *mTree = new TTree("mTree","MVTX Hit Tree");
  mTree->Branch("ox",&mT.ox,"ox/F");
  mTree->Branch("oy",&mT.oy,"oy/F");
  mTree->Branch("oz",&mT.oz,"oz/F");
  mTree->Branch("px",&mT.px,"px/F");
  mTree->Branch("py",&mT.py,"py/F");
  mTree->Branch("pz",&mT.pz,"pz/F");
  mTree->Branch("nh",&mT.nh,"nh/I");
  mTree->Branch("id",mT.id,"id[nh]/I");
  mTree->Branch("xL_mc",mT.xL_mc,"xL_mc[nh]/F");
  mTree->Branch("yL_mc",mT.yL_mc,"yL_mc[nh]/F");
  mTree->Branch("zL_mc",mT.zL_mc,"zL_mc[nh]/F");
  mTree->Branch("xG_mc",mT.xG_mc,"xG_mc[nh]/F");
  mTree->Branch("yG_mc",mT.yG_mc,"yG_mc[nh]/F");
  mTree->Branch("zG_mc",mT.zG_mc,"zG_mc[nh]/F");
  mTree->Branch("xL_rc",mT.xL_rc,"xL_rc[nh]/F");
  mTree->Branch("yL_rc",mT.yL_rc,"yL_rc[nh]/F");
  mTree->Branch("zL_rc",mT.zL_rc,"zL_rc[nh]/F");
  mTree->Branch("xG_rc",mT.xG_rc,"xG_rc[nh]/F");
  mTree->Branch("yG_rc",mT.yG_rc,"yG_rc[nh]/F");
  mTree->Branch("zG_rc",mT.zG_rc,"zG_rc[nh]/F");

  std::cout << " Generate Initial Cosmic tracks .... " << std::endl;  
  //////////////////////////////////////////////////  
  // Generating the cosmic ray track
  //////////////////////////////////////////////////  
  const TVector3 origin(0., 10., 0.);
  const Double_t thetaP = 0.;  // angle w.r.t Y-axis going down wards - cosmic ray
  const Double_t momP = 3.0; // 3.0 GeV
  Double_t pyP = - momP * TMath::Cos(thetaP);

  TRandom3 *gRandom = new TRandom3();
  Double_t phi_xz = gRandom->Rndm()*TMath::Pi()*2.;  // phi in x-z plane
  Double_t pxP = momP * TMath::Sin(thetaP) * TMath::Cos(phi_xz);
  Double_t pzP = momP * TMath::Sin(thetaP) * TMath::Sin(phi_xz);
  TVector3 pmom(pxP, pyP, pzP);

  StPhysicalHelix cosmicRay(pmom, origin, 0, 0);
  //////////////////////////////////////////////////  
  std::cout << " Calculating projections to stave planes " << std::endl;

  mT.ox = origin.x();
  mT.oy = origin.y();
  mT.oz = origin.z();
  mT.px = pmom.x();
  mT.py = pmom.y();
  mT.pz = pmom.z();
  
  int nh = 0;
  for(int i=NL-1;i>=0;i--) {  // going from most outer layer
    for(int j=0;j<NS[i];j++) {
      double *tra = stave[i][j]->GetTranslation();
      double *rot = stave[i][j]->GetRotationMatrix();

      TVector3 ori(tra[0], tra[1], tra[3]);
      TVector3 norm(rot[1], rot[4], rot[7]); // normal direction of the plane
      double s = cosmicRay.pathLength(ori, norm);
      if(s<0) continue;  // cannot go backwards
      double hitG_mc[3] = {cosmicRay.x(s), cosmicRay.y(s), cosmicRay.z(s)};
      double hitL_mc[3] = {-999., -999., -999.};
      stave[i][j]->MasterToLocal(hitG_mc, hitL_mc);
      TVector3 hit_mc(hitL_mc[0], hitL_mc[1], hitL_mc[2]);
      
      if(!onStave(hit_mc)) continue;

      std::cout << "MC Hit on stave:   layer# " << i << " sector# " << j << std::endl;
      std::cout << "     xyz = " << hitL_mc[0] << " " << hitL_mc[1] << " " << hitL_mc[2] << std::endl;

      double hitL_rc[3];;
      hitL_rc[0] = gRandom->Gaus(hitL_mc[0], PitchSize/sqrt(12.));
      hitL_rc[1] = 0.;
      hitL_rc[2] = gRandom->Gaus(hitL_mc[2], PitchSize/sqrt(12.));
      double hitG_rc[3] = {-999., -999., -999.};
      stave[i][j]->LocalToMaster(hitL_rc, hitG_rc);

      mT.id[nh] = i * 100 + j;
      mT.xL_mc[nh] = hitL_mc[0];
      mT.yL_mc[nh] = hitL_mc[1];
      mT.zL_mc[nh] = hitL_mc[2];
      mT.xG_mc[nh] = hitG_mc[0];
      mT.yG_mc[nh] = hitG_mc[1];
      mT.zG_mc[nh] = hitG_mc[2];
      mT.xL_rc[nh] = hitL_rc[0];
      mT.yL_rc[nh] = hitL_rc[1];
      mT.zL_rc[nh] = hitL_rc[2];
      mT.xG_rc[nh] = hitG_rc[0];
      mT.yG_rc[nh] = hitG_rc[1];
      mT.zG_rc[nh] = hitG_rc[2];

      nh++;
    }
  }
  std::cout << " # of Hits in this event " << nh << std::endl;
  mT.nh = nh;
  mTree->Fill();

  
  mTree->Write();
  f.Close();

  return 0;
}
