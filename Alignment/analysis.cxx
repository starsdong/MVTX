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
#endif

//#include "MuDst.h"
#include "StPhysicalHelix.h"
#include "SystemOfUnits.h"
#include "PhysicalConstants.h"

bool passEvent();
void writeHistograms(char *outFile);
void deleteHistograms();

double twoPi = 2.*3.1415927;
double eMass = 0.000511;
int iran = 0;


int main(int argc, char **argv)
{
  TChain *chain = new TChain("MuDst");

  if(argc!=3 && argc!=1) return 0;

  char *FileInput;
  char *FileOutput;

  if(argc==1){
    FileInput  = "test.list";
    FileOutput = "test.root";
  }

  if(argc==3){
    FileInput = argv[1];
    FileOutput = argv[2];
  }  

  int fileNumber = 0;

  char FileList[512];

  ifstream* inputStream = new ifstream;
  inputStream->open(FileInput);
  if (!(inputStream)) {
    printf("can not open list file\n");
    return 0;
  }
  for (;inputStream->good();) {
    inputStream->getline(FileList,512);
    if  ( inputStream->good() ) {
      printf(" read in file %s\n",FileList);
      chain->Add(FileList);
      fileNumber++;
    }
  }

  printf("  %d files read in\n",fileNumber);
/*
  MuDst* muDst = new MuDst();
  muDst->SetReadPrimary(kTRUE);
  muDst->SetReadGlobal(kTRUE);
  muDst->InitData(chain);

  int mNEvents = (int)chain->GetEntries();
  int nPTracks = 0;
  int nGTracks = 0;

  printf("total events: %d\n", mNEvents);

  printf("reading phi wgt files ...\n");
  readphiwgt();

  for(int mNEventCount = 0; mNEventCount < mNEvents; mNEventCount++){
    chain->GetEvent(mNEventCount);

    nPTracks = muDst->nPrimaryTracks();
    nGTracks = muDst->nGlobalTracks();

    CurrentEvent_nTags = 0;
    CurrentEvent_nPartners = 0;

    if(passEvent(muDst)){
      // determine the phi wgt file
      
      for(int i=0; i<nPTracks; i++){
        passTagTrack(muDst,i);
      }

      makeRealPairs();

      makeMixedPairs(BufferPointer);

      copyCurrentToBuffer(BufferPointer);      

    }
    
    if(mNEventCount%1000==0) {
      printf("Working on event #%d\n",mNEventCount);
      printf(" there are %d primary tracks and %d global tracks in current event.\n",nPTracks,nGTracks);
      printf(" there are %d tag tracks and %d partner tracks in current event.\n",CurrentEvent_nTags,CurrentEvent_nPartners);
    }
    
    if(mNEventCount%5000==0) writeHistograms(FileOutput);
  }
*/  
  writeHistograms(FileOutput);

  printf("EXIT!\n");

  deleteHistograms();  
  delete chain;

  printf("RETURN!\n");
  return 1;
}

//________________________________________________________________________
bool passEvent()
{
}
//_____________________________________________________________________________
void writeHistograms(char *outFile)
{
  TFile f(outFile,"RECREATE");
  f.cd();
  f.Close();
}
//___________________________________________________________________________
void deleteHistograms()
{
}

