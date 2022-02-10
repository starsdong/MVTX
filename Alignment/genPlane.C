/////////////////////////////////////////////////
// Fast MC for MVTX geometry and hits generation
//    Multiple Scattering + Single Hit Resolution
////////////////////////////////////////////////
void genPlane()
{
  const Int_t NL = 3;
  const Int_t NSMAX = 20;
  const Int_t NS[NL] = {12, 16, 20};
  const Double_t ZMAX = 27.1/2.;   // Total Length 2*ZMax
  const Double_t R0[NL] = {2.46, 3.23, 4.00};  // radial middle point

  TGeoHMatrix stave[NL][NSMAX];  // ideal geometry
  TGeoHMatrix stave_p[NL][NSMAX];  // mis-aligned geometry

  TRandom3 *gRandom = new TRandom3();
  
  ofstream outData;
  outData.open("StavePlane_Ideal.txt");
  for(int i=0;i<NL;i++) {
    for(int j=0;j<NS[i];j++) {
      double phi = TMath::Pi()*2./NS[i] * (j+0.5);
      double nx = TMath::Cos(phi);
      double ny = TMath::Sin(phi);
      double nz = 0.;
      double x0 = R0[i] * nx;
      double y0 = R0[i] * ny;
      double z0 = 0.;
      outData << setw(6) << i << setw(6) << j << setw(15) << x0 << setw(15) << y0 << setw(15) << z0 << setw(15) << nx << setw(15) << ny << setw(15) << nz << endl;

      double tra[3] = {x0, y0, z0};
      double rot[9] = { -ny, -nx, 0, nx, -ny, 0, 0, 0, 1 };
      stave[i][j].SetName(Form("geoM_%d_%d",i,j));
      stave[i][j].SetTranslation(tra);
      stave[i][j].SetRotation(rot);

      stave_p[i][j].SetName(Form("geoM_mis_%d_%d",i,j));
      double tra_p[3] = {x0, y0, z0};
      if(i==2 && j==5) {
	tra_p[0] += gRandom->Gaus(0,0.1);
	tra_p[1] += gRandom->Gaus(0,0.1);
	tra_p[2] += gRandom->Gaus(0,0.1);
      }
      stave_p[i][j].SetTranslation(tra_p);
      stave_p[i][j].SetRotation(rot);
      if(i==2 && j==5) {
	stave_p[i][j].RotateX(gRandom->Gaus(0,0.001*180./TMath::Pi()));
	stave_p[i][j].RotateY(gRandom->Gaus(0,0.001*180./TMath::Pi()));
	stave_p[i][j].RotateZ(gRandom->Gaus(0,0.001*180./TMath::Pi()));
      }
      
      stave[i][j].Print();
      stave_p[i][j].Print();
    }
  }
  outData.close();

  TFile *fout = new TFile("StaveGeoMatrix.root","recreate");
  for(int i=0;i<NL;i++) {
    for(int j=0;j<NS[i];j++) {
      stave[i][j].Write();
    }
  }
  for(int i=0;i<NL;i++) {
    for(int j=0;j<NS[i];j++) {
      stave_p[i][j].Write();
    }
  }
  fout->Close();
  
  
}
