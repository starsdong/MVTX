#include "draw.C+"
#include "style.C+"
#include "/Users/starsdong/work/work/fitfun/GaussFunctions.C"

void plotRes()
{
  style();

  TFile *fin = new TFile("outputFile.root");
  TH2D *hDcaXYP = (TH2D *)fin->Get("h3");
  TH2D *hDcaZP = (TH2D *)fin->Get("h4");

  const Double_t XMAX = 0.05;
  const Int_t NP = 60;
  const Double_t PMax = 3.0;

  TF1 *fitfun = new TF1("fitfun",StudentT,-XMAX,XMAX,4);
  fitfun->SetParameters(1000,0,0.005,3.0);

  //  hDcaXYPt->FitSlicesY(fitfun);


  TH1D *hDcaXY[NP];
  TH1D *hDcaZ[NP];
  for(int i=0;i<NP;i++) {
    hDcaXY[i] = (TH1D *)hDcaXYP->ProjectionY(Form("DcaXY_%d",i),i+1,i+1);
    hDcaZ[i] = (TH1D *)hDcaZP->ProjectionY(Form("DcaZ_%d",i),i+1,i+1);
  }

  
  // // Example: Student's t PDF with mean and sigma as parameters
  // TF1* myStudentT = new TF1("myStudentT", "[0]*ROOT::Math::tdistribution_pdf(x,[1], [2])", -XMAX, XMAX);
  // myStudentT->SetParameters(1000, 0,
  //  TF1 *fitfun = new TF1("fitfun","gaus",-XMAX,XMAX);
  double p[NP], sig1[NP], sig1_e[NP], sig2[NP], sig2_e[NP];
  
  TCanvas *c1 = new TCanvas("c1","",1000,800);
  c1->Draw();
  c1->Divide(5,4);

  for(int i=0;i<NP;i++) {
    p[i] = (i+0.5)*PMax/NP;
    c1->cd(i+1);
    hDcaXY[i]->Draw();


    fitfun->SetParameters(hDcaXY[i]->GetMaximum(), 0, 0.005, 3.0);
    fitfun->FixParameter(1,0.);
    fitfun->SetRange(-XMAX,XMAX);
    hDcaXY[i]->Fit("fitfun","R");

    sig1[i] = fitfun->GetParameter(2)*1e4;
    sig1_e[i] = fitfun->GetParError(2)*1e4;
    c1->cd();
  }

  for(int i=0;i<NP;i++) {
    c1->cd(i+1);
    hDcaZ[i]->Draw();


    fitfun->SetParameters(hDcaXY[i]->GetMaximum(), 0, 0.005, 3.0);
    fitfun->FixParameter(1,0.);
    fitfun->SetRange(-XMAX,XMAX);
    hDcaZ[i]->Fit("fitfun","R");

    sig2[i] = fitfun->GetParameter(2)*1e4;
    sig2_e[i] = fitfun->GetParError(2)*1e4;
    c1->cd();
  }
  
  c1->Update();

  const Int_t N0 = 3;
  TGraphErrors *gr_xy = new TGraphErrors(NP-N0, p+N0, sig1+N0, 0, sig1_e+N0);
  TGraphErrors *gr_z = new TGraphErrors(NP-N0, p+N0, sig2+N0, 0, sig2_e+N0);

  const Double_t x1 = 0.0;
  const Double_t x2 = 3.0;
  const Double_t y1 = 0.0;
  const Double_t y2 = 150.;
  TCanvas *c2 = new TCanvas("c2","",1200,800);
  c2->Draw();

  TPad *p1 = new TPad("p1","",0.0,0.13,0.5,1.0);
  p1->SetLeftMargin(0.2);
  p1->SetRightMargin(0.03);
  p1->SetTopMargin(0.03);
  p1->SetBottomMargin(0.05);
  p1->SetGridx();
  p1->SetGridy();
  p1->Draw();
  p1->cd();
  
  TH1D *h1 = new TH1D("h1","",1,x1,x2);
  h1->SetMinimum(y1);
  h1->SetMaximum(y2);
  //  h1->SetXTitle("Total Momentum p (GeV/c)");
  h1->SetYTitle("#sigma_{XY} (#mum)");
  h1->Draw();

  gr_xy->Draw("p");

  drawHistBox(x1,x2,y1,y2);
  p1->Modified();
  c2->Update();
  c2->cd();

  TPad *p2 = new TPad("p1","",0.5,0.13,1.0,1.0);
  p2->SetLeftMargin(0.2);
  p2->SetRightMargin(0.03);
  p2->SetTopMargin(0.03);
  p2->SetBottomMargin(0.05);
  p2->SetGridx();
  p2->SetGridy();
  p2->Draw();
  p2->cd();
  
  TH1D *h2 = new TH1D("h2","",1,x1,x2);
  h2->SetMinimum(y1);
  h2->SetMaximum(y2);
  //  h2->SetXTitle("Total Momentum p (GeV/c)");
  h2->SetYTitle("#sigma_{Z} (#mum)");
  h2->Draw();

  gr_z->Draw("p");

  drawHistBox(x1,x2,y1,y2);
  p2->Modified();
  c2->Update();
  c2->cd();

  TPad *p3 = new TPad("p1","",0.0,0.0,1.0,0.13);
  p3->Draw();
  p3->cd();

  drawText(0.3,0.3, "Total Momentum p (GeV/c)", 42, 0.5);

  p3->Modified();
  c2->Update();
  c2->SaveAs("fig/PointingResolution.pdf");
  c2->SaveAs("fig/PointingResolution.png");
  
}
