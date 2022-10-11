void makePlot(const char* cent = "0_10")
{
  gROOT->Reset();

  TFile *f1 = new TFile("B_and_D0_ptSpectra.root");
  TH1F *hNonD0Pt_In = (TH1F *)f1->Get("hNonD0Pt");
  hNonD0Pt_In->Scale(1e3);

  TFile *f2 = new TFile("significance.root");
  TGraphErrors *gr_noPID = (TGraphErrors *)f2->Get(Form("nonProD0_noPID_%s",cent));
  TGraphErrors *gr_cleanIdeal = (TGraphErrors *)f2->Get(Form("nonProD0_cleanIdeal_%s",cent));

  TFile *f3 = new TFile("RAA_Ave.root");
  TGraph *gr_RAA_D0_B = (TGraph *)f3->Get("RAA_D0_B");

  TH1F *hNonD0Pt = hNonD0Pt_In->Clone();
  for(int i=0;i<hNonD0Pt->GetNbinsX();i++) {
    double pt = hNonD0Pt->GetBinCenter(i+1);
    double yield = hNonD0Pt->GetBinContent(i+1)*gr_RAA_D0_B->Eval(pt);
    hNonD0Pt->SetBinContent(i+1, yield);
  }
  
  double x[20], x_noPID[20], x_cleanIdeal[20], y[20],ye_noPID[20], ye_cleanIdeal[20];
  int n = 0;
  double sum = 0.;
  double sum_e_noPID = 0.;
  double sum_e_cleanIdeal = 0.;
  for(int i=0;i<gr_noPID->GetN();i++) {
    x[i] = gr_noPID->GetX()[i];
    x_noPID[i] = x[i] - 0.05;
    x_cleanIdeal[i] = x[i] + 0.05;
    int bin = hNonD0Pt->FindBin(x[i]);
    y[i] = hNonD0Pt->GetBinContent(bin);
    ye_noPID[i] = y[i]/gr_noPID->GetY()[i];
    ye_cleanIdeal[i] = y[i]/gr_cleanIdeal->GetY()[i];

    sum += y[i];
    sum_e_noPID += ye_noPID[i]*ye_noPID[i];
    sum_e_cleanIdeal += ye_cleanIdeal[i]*ye_cleanIdeal[i];
    n++;
  }

  sum_e_noPID = sqrt(sum_e_noPID);
  sum_e_cleanIdeal = sqrt(sum_e_cleanIdeal);

  cout << " total yield error = " << sum_e_noPID/sum << "\t" << sum_e_cleanIdeal/sum << endl;
  
   TCanvas *c1 = new TCanvas("c1", "c1",0,0,800,600);
   gStyle->SetOptFit(0);
   gStyle->SetOptStat(0);
   gStyle->SetEndErrorSize(0.01);
   c1->SetFillColor(10);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetFrameFillColor(0);
   c1->SetFrameBorderMode(0);

   c1->SetLeftMargin(0.13);
   c1->SetBottomMargin(0.15);
   c1->SetTopMargin(0.02);
   c1->SetRightMargin(0.02);

   double x1 = 0;
   double x2 = 8;
   double y1 = 0;
   double y2 = hNonD0Pt->GetMaximum()/0.8;
   
   TH1D *d0 = new TH1D("d0","",1,x1,x2);
   d0->SetMinimum(y1);
   d0->SetMaximum(y2);
   d0->GetXaxis()->SetNdivisions(208);
   d0->GetXaxis()->SetTitle("Non-prompt D^{0} p_{T} (GeV/c)");
   d0->GetXaxis()->SetTitleOffset(1.0);
   d0->GetXaxis()->SetTitleSize(0.06);
   d0->GetXaxis()->SetLabelOffset(0.01);
   d0->GetXaxis()->SetLabelSize(0.045);
   d0->GetXaxis()->SetLabelFont(42);
   d0->GetXaxis()->SetTitleFont(42);
   d0->GetYaxis()->SetNdivisions(204);
   d0->GetYaxis()->SetTitle("d#sigma/dp_{T}dy/N_{bin} [#mub/(GeV/c)]");
   d0->GetYaxis()->SetTitleOffset(1.0);
   d0->GetYaxis()->SetTitleSize(0.06);
   d0->GetYaxis()->SetLabelOffset(0.005);
   d0->GetYaxis()->SetLabelSize(0.045);
   d0->GetYaxis()->SetLabelFont(42);
   d0->GetYaxis()->SetTitleFont(42);
   d0->Draw();

   //   hNonD0Pt_In->Draw("same");
   hNonD0Pt->Draw("same");

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
  
   TGraphErrors *gr_noPID = new TGraphErrors(n, x_noPID, y, 0, ye_noPID);
   gr_noPID->SetMarkerStyle(20);
   gr_noPID->SetMarkerSize(1.5);
   gr_noPID->SetLineWidth(2);
   gr_noPID->Draw("p");

   TGraphErrors *gr_cleanIdeal = new TGraphErrors(n, x_cleanIdeal, y, 0, ye_cleanIdeal);
   gr_cleanIdeal->SetMarkerStyle(24);
   gr_cleanIdeal->SetMarkerSize(1.5);
   gr_cleanIdeal->SetLineWidth(2);
   gr_cleanIdeal->Draw("p");


   TLegend *leg = new TLegend(0.60, 0.66, 0.96, 0.8);
   leg->SetFillColor(10);
   leg->SetFillStyle(10);
   leg->SetLineStyle(4000);
   leg->SetLineColor(10);
   leg->SetLineWidth(0.);
   leg->SetTextFont(42);
   leg->SetTextSize(0.045);
   leg->AddEntry(gr_noPID, "No PID,  (7.7%)","p");
   leg->AddEntry(gr_cleanIdeal, "Clean PID,  (2.6%)","p");   
   leg->Draw();

   TLatex *tex = new TLatex(3.0, y2*0.88, "Au+Au 200 GeV, 10B 0-10%");
   tex->SetTextFont(42);
   tex->SetTextSize(0.055);
   tex->Draw("same");

   c1->SaveAs("fig/Xsec.eps");
   c1->SaveAs("fig/Xsec.png");
   
}



