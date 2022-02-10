/////////////////////////////////////////////////
// analyze MC trees for alingment calibration
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
#include "MvtxConstants.h"
#include "StPhysicalHelix.h"
#include "SystemOfUnits.h"
#include "PhysicalConstants.h"

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
  std::cout << " Reading in Stave Geometry Tables ... " << std::endl;
  TGeoHMatrix *stave[NL][NSMAX];
  TGeoHMatrix *staveP[NL][NSMAX];
  TGeoHMatrix corr[NL][NSMAX];
  TGeoHMatrix corrL[NL][NSMAX];
  TFile *fin = new TFile("StaveGeoMatrix.root");
  for(int i=0;i<NL;i++) {
    for(int j=0;j<NS[i];j++) {
      std::cout << Form("geoM_%d_%d",i,j) << std::endl;
      stave[i][j] = (TGeoHMatrix *)fin->Get(Form("geoM_%d_%d",i,j));  // default ideal geometry
      stave[i][j]->Print();
      staveP[i][j] = (TGeoHMatrix *)fin->Get(Form("geoM_mis_%d_%d",i,j));  // mis-aligned geometry
      staveP[i][j]->Print();

      corr[i][j] = (*staveP[i][j]) * (*stave[i][j]).Inverse();
      corr[i][j].SetName(Form("geoM_corr_%d_%d",i,j));
      corr[i][j].Print();

      corrL[i][j] = (*stave[i][j]).Inverse() * (*staveP[i][j]);
      corrL[i][j].SetName(Form("geoM_corrL_%d_%d",i,j));
      corrL[i][j].Print();
    }
  }
  fin->Close();

  TH2F *hRes[3][3][NL][NSMAX];  // Global residuals
  for(int ixyz=0;ixyz<9;ixyz++) {
    for(int i=0;i<NL;i++) {
      for(int j=0;j<NS[i];j++) {
	if(ixyz%3<2) {
	  hRes[ixyz/3][ixyz%3][i][j] = new TH2F(Form("Res_%d_%d_%d_%d",ixyz/3,ixyz%3,i,j),"",200,-5.,5.,1000,-0.25,0.25);
	} else{
	  hRes[ixyz/3][ixyz%3][i][j] = new TH2F(Form("Res_%d_%d_%d_%d",ixyz/3,ixyz%3,i,j),"",200,-ZMAX,ZMAX,1000,-0.25,0.25);
	}
      }
    }
  }
  TH2F *hResL[2][2][NL][NSMAX]; // Local residuals X-Z only
  for(int ixz=0;ixz<4;ixz++) {
    for(int i=0;i<NL;i++) {
      for(int j=0;j<NS[i];j++) {
	if(ixz%2<1) {
	  hResL[ixz/2][ixz%2][i][j] = new TH2F(Form("ResL_%d_%d_%d_%d",ixz/2,ixz%2,i,j),"",200,-XMAX,XMAX,1000,-0.25,0.25);
	} else {
	  hResL[ixz/2][ixz%2][i][j] = new TH2F(Form("ResL_%d_%d_%d_%d",ixz/2,ixz%2,i,j),"",200,-ZMAX,ZMAX,1000,-0.25,0.25);
	}
      }
    }
  }

  MVTXTREE mT;
  TFile *f1 = new TFile("test.root");
  TTree *mTree = (TTree *)f1->Get("mTree");
  mTree->SetBranchAddress("ox",&mT.ox);
  mTree->SetBranchAddress("oy",&mT.oy);
  mTree->SetBranchAddress("oz",&mT.oz);
  mTree->SetBranchAddress("px",&mT.px);
  mTree->SetBranchAddress("py",&mT.py);
  mTree->SetBranchAddress("pz",&mT.pz);
  mTree->SetBranchAddress("nh",&mT.nh);
  mTree->SetBranchAddress("id",mT.id);
  mTree->SetBranchAddress("xL_mc",mT.xL_mc);
  mTree->SetBranchAddress("yL_mc",mT.yL_mc);
  mTree->SetBranchAddress("zL_mc",mT.zL_mc);
  mTree->SetBranchAddress("xG_mc",mT.xG_mc);
  mTree->SetBranchAddress("yG_mc",mT.yG_mc);
  mTree->SetBranchAddress("zG_mc",mT.zG_mc);
  mTree->SetBranchAddress("xL_rc",mT.xL_rc);
  mTree->SetBranchAddress("yL_rc",mT.yL_rc);
  mTree->SetBranchAddress("zL_rc",mT.zL_rc);
  mTree->SetBranchAddress("xG_rc",mT.xG_rc);
  mTree->SetBranchAddress("yG_rc",mT.yG_rc);
  mTree->SetBranchAddress("zG_rc",mT.zG_rc);

  //////////////////////////////////////////////////  
  std::cout << " Number of entries in the tree " << mTree->GetEntries() << std::endl;

  for(int i=0;i<mTree->GetEntries();i++) {
    mTree->GetEntry(i);
    if(i%10000==0) std::cout << " nh = " << mT.nh << std::endl;

    if(mT.nh<3) continue;   // at least 3 hits

    /*
    int indexSort[20];  // store indices of sorted hits (according descending yG)
    bool sorted[20] = {0};
    for(int j=0;j<mT.nh;j++) {
      double y_tmp = -999.;
      for(int k=0;k<mT.nh;k++) {
	if(!sorted[k] && mT.yG_rc[k]>y_tmp) {
	  y_tmp = mT.yG_rc[k];
	  indexSort[j] = k;
	}
      }
      sorted[indexSort[j]] = true;      
    }
    std::cout << "Sorted indecies = ";
    for(int j=0;j<mT.nh;j++) std::cout << indexSort[j] << " ";
    std::cout << std::endl;
    */

    for(int j=0;j<mT.nh-1;j++) {
      int i_p1 = j; // indexSort[j];
      for(int k=j+1;k<mT.nh;k++) {
	int i_p2 = k; // indexSort[k];

	TVector3 p1(mT.xG_rc[i_p1], mT.yG_rc[i_p1], mT.zG_rc[i_p1]);
	TVector3 p2(mT.xG_rc[i_p2], mT.yG_rc[i_p2], mT.zG_rc[i_p2]);
	if((p2-p1).Mag()<0.5) continue;  // don't use hits from the same layer

	StPhysicalHelix track(p2-p1, p2, 0, 0);

	for(int l=0;l<mT.nh;l++) {
	  if(l==j || l==k) continue;  // project to any stave
	  int i_p3 = l; // indexSort[l];
	  TVector3 p3(mT.xG_rc[i_p3], mT.yG_rc[i_p3], mT.zG_rc[i_p3]);
	  if((p3-p2).Mag()<0.5) continue;  // don't use hits from the same layer
	  if((p3-p1).Mag()<0.5) continue;  // don't use hits from the same layer

	  int il = mT.id[i_p3]/100;  // decode the detector ID
	  int is = mT.id[i_p3]%100;
	  
	  //	  if(il==2 && is==5) std::cout << " test " << std::endl;
	  double *tra = stave[il][is]->GetTranslation();
	  double *rot = stave[il][is]->GetRotationMatrix();
	  TVector3 ori(tra[0], tra[1], tra[3]);
	  TVector3 norm(rot[1], rot[4], rot[7]); // normal direction of the plane

	  double s = track.pathLength(ori, norm);
	  double hitG_proj[3] = {track.x(s), track.y(s), track.z(s)};
	  double hitL_proj[3] = {-999., -999., -999.};
	  stave[il][is]->MasterToLocal(hitG_proj, hitL_proj);
	  //	  if(il==2 && is==5)
	  //  std::cout << " Projection " << hitL_proj[0] << " " << hitL_proj[1] << " " << hitL_proj[2] << std::endl;
	  TVector3 hit_proj(hitL_proj[0], hitL_proj[1], hitL_proj[2]);
	  if(!onStave(hit_proj)) continue;
	  

	  hRes[0][0][il][is]->Fill(mT.xG_rc[i_p3], hitG_proj[0] - mT.xG_rc[i_p3]);
	  hRes[0][1][il][is]->Fill(mT.yG_rc[i_p3], hitG_proj[0] - mT.xG_rc[i_p3]);
	  hRes[0][2][il][is]->Fill(mT.zG_rc[i_p3], hitG_proj[0] - mT.xG_rc[i_p3]);
	  hRes[1][0][il][is]->Fill(mT.xG_rc[i_p3], hitG_proj[1] - mT.yG_rc[i_p3]);
	  hRes[1][1][il][is]->Fill(mT.yG_rc[i_p3], hitG_proj[1] - mT.yG_rc[i_p3]);
	  hRes[1][2][il][is]->Fill(mT.zG_rc[i_p3], hitG_proj[1] - mT.yG_rc[i_p3]);
	  hRes[2][0][il][is]->Fill(mT.xG_rc[i_p3], hitG_proj[2] - mT.zG_rc[i_p3]);
	  hRes[2][1][il][is]->Fill(mT.yG_rc[i_p3], hitG_proj[2] - mT.zG_rc[i_p3]);
	  hRes[2][2][il][is]->Fill(mT.zG_rc[i_p3], hitG_proj[2] - mT.zG_rc[i_p3]);

	  hResL[0][0][il][is]->Fill(mT.xL_rc[i_p3], hitL_proj[0] - mT.xL_rc[i_p3]);
	  hResL[0][1][il][is]->Fill(mT.zL_rc[i_p3], hitL_proj[0] - mT.xL_rc[i_p3]);
	  hResL[1][0][il][is]->Fill(mT.xL_rc[i_p3], hitL_proj[2] - mT.zL_rc[i_p3]);
	  hResL[1][1][il][is]->Fill(mT.zL_rc[i_p3], hitL_proj[2] - mT.zL_rc[i_p3]);
	}
      }
    }
  }

  TFile *fout = new TFile("residualHist.root","recreate");
  for(int ixyz=0;ixyz<9;ixyz++) {
    for(int i=0;i<NL;i++) {
      for(int j=0;j<NS[i];j++) {
	hRes[ixyz/3][ixyz%3][i][j]->Write();
      }
    }
  }
  for(int ixz=0;ixz<4;ixz++) {
    for(int i=0;i<NL;i++) {
      for(int j=0;j<NS[i];j++) {
	hResL[ixz/2][ixz%2][i][j]->Write();
      }
    }
  }
  fout->Close();


  return 0;
}
