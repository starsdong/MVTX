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

#ifndef __CINT__
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
#include "TGeoHMatrix.h"
#endif

#include "StPhysicalHelix.h"
#include "SystemOfUnits.h"
#include "PhysicalConstants.h"

/*
void FreeStreaming(TVector3 pos_0, TVector3 mom_0, float s_z, TVector3 pos_1, TVector3 mom_1)
{
}
*/

double theta0(double p, double m, double dx)
{
  double beta = p/TMath::Sqrt(p*p+m*m);
  return 0.0136/beta/p*TMath::Sqrt(dx)*(1.+0.038*TMath::Log(dx/beta/beta));
}

int main(int argc, char **argv)
{
  const Int_t Nevt = 100000.;

  const Int_t NL = 3;
  const Int_t NSMAX = 20;
  const Int_t NS[NL] = {12, 16, 20};
  const Double_t ZMAX = 27.1/2.;   // Total Length 2*ZMax
  const Double_t XMAX = 1.5/2.;    // Total Width 2*XMAX
  const Double_t R0[NL] = {2.46, 3.23, 4.00};  // radial middle point
  const Double_t PitchSize = 0.0028;   //

  TGeoHMatrix *stave[NL][NSMAX];
  TFile *fin = new TFile("StaveGeoMatrix_Ideal.root");
  for(int i=0;i<NL;i++) {
    for(int j=0;j<NS[i];j++) {
      stave[i][j] = (TGeoHMatrix *)fin->Get(Form("geoM_%d_%d",i,j));
    }
  }
  
  //////////////////////////////////////////////////  
  // Generating the cosmic ray track
  //////////////////////////////////////////////////  
  const TVector3 origin(0., 10., 0.);
  const Double_t thetaP = 0.;  // angle w.r.t Y-axis going down wards - cosmic ray
  const Double_t momP = 3.0; // 3.0 GeV
  doublex pyP = - momP * TMath::Cos(thetaP);

  TRandom3 *gRandom = new TRandom3();
  Double_t phi_xz = gRandom->Rndm()*TMath::Pi()*2.;  // phi in x-z plane
  double pxP = momP * TMath::Sin(thetaP) * TMath::Cos(phi_xz);
  double pzP = momP * TMath::Sin(thetaP) * TMath::Sin(phi_xz);
  TVector3 pmom(pxP, pyP, pzP);

  StPhysicalHelix cosmicRay(pmom, origin, 0, 0);
  //////////////////////////////////////////////////  
  

}
