#include <stdio>
#include "iomanip.h"

Double_t CorrBkgd2SignalRatio(Double_t pT)  // obtained from data 
{
  if(pT<1.5) return 0.;
  else if(pT<6.5) return (pT-1.5)/(6.5-1.5)*0.4;
  else return 0.4;
}

void makePlot_preproposal()
{
  gROOT->Reset();

  const Int_t N = 10; // number of pT bins
  const Int_t NC = 3; // number of configurations
  //  const Char_t *Name[NC] = {"noPid","hybrid_ideal","hybrid","clean_ideal","clean"};
  const Char_t *Name[NC] = {"noPid","clean_ideal","clean"};

  const Double_t Nbin_cen = 941.; // Nbin 0-10%
  const Double_t Nbin_mb = 292.;  // Nbin 0-80%
  const Double_t Nbin_peri = 21.6;  // Nbin 60-80%
  const Double_t Npart_cen = 325.; // Npart 0-10%
  const Double_t Npart_mb = 127.;  // Npart 0-80%
  const Double_t Npart_peri = 21.0; // Npart 60-80%
  
  // 240B mb event total
  const Double_t Nevt_mb = 240;  //
  const Double_t Nevt_cen = Nevt_mb*0.1;  // 0-10%
  const Double_t Nevt_peri = Nevt_mb*0.2; // 60-80%
  
  // input files are for 100M 0-10% central
  const Double_t Nevt_cen_0 = 0.1;  // input 100M
  
  const Double_t pt[N] = {0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5};
  Double_t d0[NC][N], d0_B[NC][N], bkgd[NC][N], bkgd_B[NC][N];
  Double_t bkgd_B_corr[NC][N]; // estimated correlated background
  Double_t sig_d0[NC][N], sig_d0_B[NC][N];

  Double_t sig_d0_mb[NC][N], sig_d0_B_mb[NC][N];
  Double_t sig_d0_peri[NC][N], sig_d0_B_peri[NC][N];
  
  TGraph *gr_d0[NC];
  TGraph *gr_d0_B[NC];

  TGraph *gr_d0_mb[NC];
  TGraph *gr_d0_B_mb[NC];
  
  TGraph *gr_d0_peri[NC];
  TGraph *gr_d0_B_peri[NC];

  ifstream inData;
  for(int i=NC-1;i>=0;i--) {
    inData.open(Form("dat/count_%s.txt",Name[i]));
    for(int j=0;j<N;j++) {
      inData >> d0[i][j] >> bkgd[i][j] >> d0_B[i][j] >> bkgd_B[i][j];

      double scale = bkgd_B[i][j]/bkgd_B[2][j]; // scale from the clean PID case
      double bkgd_corr = d0_B[i][j]*CorrBkgd2SignalRatio(pt[j])*scale;
	
      sig_d0[i][j] = d0[i][j]/sqrt(d0[i][j]+bkgd[i][j])*sqrt(Nevt_cen/Nevt_cen_0);
      sig_d0_B[i][j] = d0_B[i][j]/sqrt(d0_B[i][j]+bkgd_B[i][j]+bkgd_corr)*sqrt(Nevt_cen/Nevt_cen_0);

      sig_d0_mb[i][j] = d0[i][j]*(Nbin_mb/Nbin_cen)/sqrt(d0[i][j]*(Nbin_mb/Nbin_cen)+bkgd[i][j]*(Npart_mb*Npart_mb/Npart_cen/Npart_cen))*sqrt(Nevt_mb/Nevt_cen_0);
      sig_d0_B_mb[i][j] = d0_B[i][j]*(Nbin_mb/Nbin_cen)/sqrt(d0_B[i][j]*(Nbin_mb/Nbin_cen)+(bkgd_B[i][j]+bkgd_corr)*(Npart_mb*Npart_mb/Npart_cen/Npart_cen))*sqrt(Nevt_mb/Nevt_cen_0);
      
      sig_d0_peri[i][j] = d0[i][j]*(Nbin_peri/Nbin_cen)/sqrt(d0[i][j]*(Nbin_peri/Nbin_cen)+bkgd[i][j]*(Npart_peri*Npart_peri/Npart_cen/Npart_cen))*sqrt(Nevt_peri/Nevt_cen_0);
      sig_d0_B_peri[i][j] = d0_B[i][j]*(Nbin_peri/Nbin_cen)/sqrt(d0_B[i][j]*(Nbin_peri/Nbin_cen)+(bkgd_B[i][j]+bkgd_corr)*(Npart_peri*Npart_peri/Npart_cen/Npart_cen))*sqrt(Nevt_peri/Nevt_cen_0);
    }
    inData.close();
    gr_d0[i] = new TGraph(N, pt, sig_d0[i]);
    gr_d0[i]->SetName(Form("sig_d0_Cen_%s",Name[i]));
    gr_d0[i]->SetMarkerSize(1.8);
    gr_d0[i]->SetMarkerColor(i+1);
    gr_d0[i]->SetMarkerStyle(24+i);
    gr_d0[i]->Print();
      
    gr_d0_B[i] = new TGraph(N, pt, sig_d0_B[i]);
    gr_d0_B[i]->SetName(Form("sig_d0_B_Cen_%s",Name[i]));
    gr_d0_B[i]->SetMarkerSize(1.8);
    gr_d0_B[i]->SetMarkerColor(i+1);
    gr_d0_B[i]->SetMarkerStyle(20+i);
    gr_d0_B[i]->Print();

    gr_d0_mb[i] = new TGraph(N, pt, sig_d0_mb[i]);
    gr_d0_mb[i]->SetName(Form("sig_d0_MB_%s",Name[i]));
    gr_d0_mb[i]->SetMarkerSize(1.8);
    gr_d0_mb[i]->SetMarkerColor(i+1);
    gr_d0_mb[i]->SetMarkerStyle(24+i);
    gr_d0_mb[i]->Print();
      
    gr_d0_B_mb[i] = new TGraph(N, pt, sig_d0_B_mb[i]);
    gr_d0_B_mb[i]->SetName(Form("sig_d0_B_MB_%s",Name[i]));
    gr_d0_B_mb[i]->SetMarkerSize(1.8);
    gr_d0_B_mb[i]->SetMarkerColor(i+1);
    gr_d0_B_mb[i]->SetMarkerStyle(20+i);
    gr_d0_B_mb[i]->Print();

    gr_d0_peri[i] = new TGraph(N, pt, sig_d0_peri[i]);
    gr_d0_peri[i]->SetName(Form("sig_d0_Peri_%s",Name[i]));
    gr_d0_peri[i]->SetMarkerSize(1.8);
    gr_d0_peri[i]->SetMarkerColor(i+1);
    gr_d0_peri[i]->SetMarkerStyle(24+i);
    gr_d0_peri[i]->Print();
      
    gr_d0_B_peri[i] = new TGraph(N, pt, sig_d0_B_peri[i]);
    gr_d0_B_peri[i]->SetName(Form("sig_d0_B_Peri_%s",Name[i]));
    gr_d0_B_peri[i]->SetMarkerSize(1.8);
    gr_d0_B_peri[i]->SetMarkerColor(i+1);
    gr_d0_B_peri[i]->SetMarkerStyle(20+i);
    gr_d0_B_peri[i]->Print();

  }


  TCanvas *c1 = new TCanvas("c1", "c1",0,0,800,600);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0.01);
  c1->SetFillColor(10);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetFrameFillColor(0);
  c1->SetFrameBorderMode(0);

  c1->SetLogy();
  c1->SetGridx();
  c1->SetGridy();
  
  c1->SetLeftMargin(0.16);
  c1->SetBottomMargin(0.14);
  c1->SetTopMargin(0.02);
  c1->SetRightMargin(0.02);
  
  double x1 = 0;
  double x2 = 10.;
  double y1 = 0.5;
  double y2 = 9000;
  
  TH1D *h0 = new TH1D("h0","",1,x1,x2);
  h0->SetMinimum(y1);
  h0->SetMaximum(y2);
  h0->GetXaxis()->SetNdivisions(205);
  h0->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h0->GetXaxis()->SetTitleOffset(1.0);
  h0->GetXaxis()->SetTitleSize(0.06);
  h0->GetXaxis()->SetLabelOffset(0.005);
  h0->GetXaxis()->SetLabelSize(0.045);
  h0->GetXaxis()->SetLabelFont(42);
  h0->GetXaxis()->SetTitleFont(42);
  h0->GetYaxis()->SetNdivisions(202);
  h0->GetYaxis()->SetTitle("S/#sqrt{S+B}");
  h0->GetYaxis()->SetTitleOffset(1.2);
  h0->GetYaxis()->SetTitleSize(0.06);
  h0->GetYaxis()->SetLabelOffset(0.005);
  h0->GetYaxis()->SetLabelSize(0.045);
  h0->GetYaxis()->SetLabelFont(42);
  h0->GetYaxis()->SetTitleFont(42);
  h0->Draw();
  
  TLine *l1 = new TLine(x1,y1,x2,y1);
  l1->SetLineWidth(3);
  l1->Draw("same");
  TLine *l2 = new TLine(x1,y2,x2,y2);
  l2->SetLineWidth(3);
  l2->Draw("same");
  TLine *l3 = new TLine(x1,y1,x1,y2);
  l3->SetLineWidth(3);
  l3->Draw("same");
  TLine *l4 = new TLine(x2,y1,x2,y2);
  l4->SetLineWidth(3);
  l4->Draw("same");
  

  for(int i=0;i<2;i++) {
    if(i==0) gr_d0[i]->Draw("p");
    gr_d0_B[i]->Draw("p");    
  }

  TLegend *leg = new TLegend(0.45, 0.16, 0.96, 0.4);
  leg->SetFillColor(10);
  leg->SetFillStyle(10);
  leg->SetLineStyle(4000);
  leg->SetLineColor(10);
  leg->SetLineWidth(0.);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->AddEntry(gr_d0[0],"prompt D^{0}","p");
  leg->AddEntry(gr_d0_B[1], "non-prompt D^{0} - w/ TOF", "p");
  leg->AddEntry(gr_d0_B[0], "non-prompt D^{0}", "p");
  leg->Draw();

  TLatex *tex = new TLatex(6.0, 3000, "Au + Au @ 200 GeV");
  tex->SetTextFont(42);
  tex->SetTextSize(0.05);
  tex->Draw("same");

  TLatex *tex = new TLatex(7.5, 1500, "10B 0-10%");
  tex->SetTextFont(42);
  tex->SetTextSize(0.05);
  tex->Draw("same");

  c1->Update();
  c1->SaveAs("fig/D0_10B_wide.eps");
  c1->SaveAs("fig/D0_10B_wide.png");

  TFile *fin = new TFile("root/significance.root","recreate");
  for(int i=0;i<NC;i++) {
    gr_d0[i]->Write();
    gr_d0_B[i]->Write();
    gr_d0_mb[i]->Write();
    gr_d0_B_mb[i]->Write();
    gr_d0_peri[i]->Write();
    gr_d0_B_peri[i]->Write();
  }
  fin->Close();
  
}
