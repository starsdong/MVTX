#include "style.C+"
#include "draw.C+"

void drawBLvsRun()
{
  style();
  const Int_t NMAX = 10000;
  const Int_t NL = 3;
  Int_t Run[NMAX];
  Double_t index[NMAX], runId[NMAX];
  Double_t x0[NL][NMAX], x0_err[NL][NMAX];
  Double_t y0[NL][NMAX], y0_err[NL][NMAX];
  Double_t d0[NL][NMAX], d0_err[NL][NMAX];
  Double_t phi0[NL][NMAX], phi0_err[NL][NMAX];
  Int_t N = 0;

  ifstream inData;
  inData.open("run.list");
  while (!inData.eof()) {
    inData >> Run[N];
    //    inData >> x0[N] >> x0_err[N] >> y0[N] >> y0_err[N] >> d0[N] >> d0_err[N] >> phi0[N] >> phi0_err[N];
    if(!inData.eof())
      N++;
  }
  inData.close();
  cout << " # of runs read in " << N << endl;
  cout << " Run IDs: ";
  for(int i=0;i<N;i++) {
    cout << " " << Run[i];
  }
  cout << endl;

  for(int i=0;i<N;i++) {
    index[i] = i+1;
    runId[i] = Run[i];
    inData.open(Form("dat/MVTX_BL_%d.txt",Run[i]));
    for(int j=0;j<NL;j++) {
      inData >> x0[j][i] >> x0_err[j][i] >> y0[j][i] >> y0_err[j][i] >> d0[j][i] >> d0_err[j][i] >> phi0[j][i] >> phi0_err[j][i];
    }
    inData.close();
  }

  TGraph *gr_run = new TGraph(N, index, runId);
  gr_run->SetName("gr_RunId");

  TGraphErrors *gr_x0[NL], *gr_y0[NL], *gr_d0[NL], *gr_phi0[NL];
  for(int i=0;i<NL;i++) {
    gr_x0[i] = new TGraphErrors(N, index, x0[i], 0, x0_err[i]);
    gr_x0[i]->SetName("gr_X0");
    gr_y0[i] = new TGraphErrors(N, index, y0[i], 0, y0_err[i]);
    gr_y0[i]->SetName("gr_Y0");
    gr_d0[i] = new TGraphErrors(N, index, d0[i], 0, d0_err[i]);
    gr_d0[i]->SetName("gr_d0");
    gr_phi0[i] = new TGraphErrors(N, index, phi0[i], 0, phi0_err[i]);
    gr_phi0[i]->SetName("gr_phi0");
  }
  
  TCanvas *c1 = new TCanvas("c1","",1200,600);
  c1->SetLeftMargin(0.1);
  c1->Draw();
  
  double x1 = -N*0.1;
  double x2 = N*1.25;
  double y1 = -6.;
  double y2 = 3.;
  TH1D *h0 = new TH1D("h0","",1,x1,x2);
  h0->SetMaximum(y2);
  h0->SetMinimum(y1);
  h0->GetXaxis()->SetTickLength(0.001);
  h0->GetXaxis()->SetLabelOffset(999.);
  h0->GetXaxis()->SetTitleOffset(0.5);
  h0->GetXaxis()->SetTitle("Run ID");
  h0->GetYaxis()->SetTitleOffset(0.6);
  h0->GetYaxis()->SetTitle("Beam Offset (mm)");
  h0->Draw("c");

  drawLine(x1, 0, x2, 0, 2, 9);

  const Int_t markerStyle[NL] = {20, 21, 22};
  const Int_t markerColor[NL] = {1, 2, 4};
  for(int i=0;i<NL;i++) {
    gr_x0[i]->SetMarkerStyle(markerStyle[i]);
    gr_x0[i]->SetMarkerColor(markerColor[i]);
    gr_x0[i]->SetLineColor(markerColor[i]);
    gr_x0[i]->SetLineWidth(2);
    gr_x0[i]->Draw("p");
  
    gr_y0[i]->SetMarkerStyle(markerStyle[i]+4);
    gr_y0[i]->SetMarkerColor(markerColor[i]);
    gr_y0[i]->SetLineColor(markerColor[i]);
    gr_y0[i]->SetLineWidth(2);
    gr_y0[i]->Draw("p");
  }

  for(int i=0;i<N;i++) {
    if(i%10==0) drawText(index[i]-0.1, y1-1.6, Form("%d",Run[i]), 42, 0.04, 75);
  }

  drawText(-N*0.07, x0[0][0], "X_{0}", 52, 0.055); 
  drawText(-N*0.07, y0[0][0], "Y_{0}", 52, 0.055);

  TLegend *leg1 = new TLegend(0.82, 0.45, 0.86, 0.65);
  leg1->SetLineColor(10);
  for(int i=NL-1;i>=0;i--) {
    leg1->AddEntry(gr_x0[i], " ", "p");
  }
  leg1->Draw();
  TLegend *leg2 = new TLegend(0.84, 0.45, 0.96, 0.65);
  leg2->SetLineColor(10);
  leg2->SetTextSize(0.04);
  for(int i=NL-1;i>=0;i--) {
    leg2->AddEntry(gr_y0[i], Form(" Layer %d", i), "p");
  }
  leg2->Draw();
  
  c1->Update();
  c1->SaveAs("fig/MVTX_BL_run.pdf");
  c1->SaveAs("fig/MVTX_BL_run.png");

}
