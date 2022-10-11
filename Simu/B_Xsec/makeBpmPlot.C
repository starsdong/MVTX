void makeBpmPlot(const char* cent = "0_10")
{
  gROOT->LoadMacro("sPhenixStyle.C");
  SetsPhenixStyle();

  TFile *f1 = new TFile("B_and_D0_ptSpectra.root");
  TH1F *hBPt_In = (TH1F *)f1->Get("hBPt_fonll");
  hBPt_In->Scale(1e3);

  TFile *f2 = new TFile("Bpm_significance.root");
  TGraph *gr_noPID = (TGraph *)(f2->Get(Form("Bpm_significance_%s",cent)));
  //  TGraphErrors *gr_noPID = new TGraphErrors((TH1F *)f2->Get(Form("Bpm_significance_noPid_%s",cent)));
  gr_noPID->Print();
  TGraph *gr_cleanIdeal = (TGraph *)(f2->Get(Form("Bpm_significance_%s",cent)));
  // TGraphErrors *gr_cleanIdeal = new TGraphErrors((TH1F *)f2->Get(Form("Bpm_significance_clean_%s",cent)));

  TFile *f3 = new TFile("RAA_Ave.root");
  TGraph *gr_RAA_B = (TGraph *)f3->Get("RAA_B");

  TH1F *hBPt = hBPt_In->Clone();
  for(int i=0;i<hBPt->GetNbinsX();i++) {
    double pt = hBPt->GetBinCenter(i+1);
    double yield = hBPt->GetBinContent(i+1)*gr_RAA_B->Eval(pt);
    hBPt->SetBinContent(i+1, yield);
  }

  cout << " AAAA " << endl;
  
  double x[20], x_noPID[20], x_cleanIdeal[20], y[20],ye_noPID[20], ye_cleanIdeal[20];
  int n = 0;
  double sum = 0.;
  double sum_e_noPID = 0.;
  double sum_e_cleanIdeal = 0.;
  for(int i=0;i<gr_noPID->GetN();i++) {
    x[i] = gr_noPID->GetX()[i];
    x_noPID[i] = x[i] - 0.05;
    x_cleanIdeal[i] = x[i] + 0.05;
    int bin = hBPt->FindBin(x[i]);
    y[i] = hBPt->GetBinContent(bin);
    ye_noPID[i] = y[i]/gr_noPID->GetY()[i];
    ye_cleanIdeal[i] = y[i]/gr_cleanIdeal->GetY()[i];

    sum += y[i];
    sum_e_noPID += ye_noPID[i]*ye_noPID[i];
    sum_e_cleanIdeal += ye_cleanIdeal[i]*ye_cleanIdeal[i];
    n++;
  }

  cout << " BBBB " << endl;
  
  sum_e_noPID = sqrt(sum_e_noPID);
  sum_e_cleanIdeal = sqrt(sum_e_cleanIdeal);

  cout << " total yield error = " << sum_e_noPID/sum << "\t" << sum_e_cleanIdeal/sum << endl;
  
   TCanvas *c1 = new TCanvas("c1", "c1",0,0,800,600);
   // gStyle->SetOptFit(0);
   // gStyle->SetOptStat(0);
   // gStyle->SetEndErrorSize(0.01);
   // c1->SetFillColor(10);
   // c1->SetBorderMode(0);
   // c1->SetBorderSize(2);
   // c1->SetFrameFillColor(0);
   // c1->SetFrameBorderMode(0);

   // c1->SetLeftMargin(0.13);
   // c1->SetBottomMargin(0.15);
   // c1->SetTopMargin(0.02);
   // c1->SetRightMargin(0.02);

   double x1 = 0;
   double x2 = 8;
   double y1 = 0;
   double y2 = hBPt->GetMaximum()/0.8;
   
   TH1D *d0 = new TH1D("d0","",1,x1,x2);
   d0->SetMinimum(y1);
   d0->SetMaximum(y2);
   // d0->GetXaxis()->SetNdivisions(208);
   d0->GetXaxis()->SetTitle("B^{#pm} p_{T} (GeV/c)");
   // d0->GetXaxis()->SetTitleOffset(1.0);
   // d0->GetXaxis()->SetTitleSize(0.06);
   // d0->GetXaxis()->SetLabelOffset(0.01);
   // d0->GetXaxis()->SetLabelSize(0.045);
   // d0->GetXaxis()->SetLabelFont(42);
   // d0->GetXaxis()->SetTitleFont(42);
   // d0->GetYaxis()->SetNdivisions(204);
   d0->GetYaxis()->SetTitle("d#sigma/dp_{T}dy/N_{bin} [#mub/(GeV/c)]");
   // d0->GetYaxis()->SetTitleOffset(1.0);
   // d0->GetYaxis()->SetTitleSize(0.06);
   // d0->GetYaxis()->SetLabelOffset(0.005);
   // d0->GetYaxis()->SetLabelSize(0.045);
   // d0->GetYaxis()->SetLabelFont(42);
   // d0->GetYaxis()->SetTitleFont(42);
   d0->Draw();

   //   hBPt_In->Draw("same");
   hBPt->Draw("same");

   // TLine *l1 = new TLine(x1,y1,x2,y1);
   // l1->SetLineWidth(3);
   // l1->Draw("same");
   // TLine *l2 = new TLine(x1,y2,x2,y2);
   // l2->SetLineWidth(3);
   // l2->Draw("same");
   // TLine *l3 = new TLine(x1,y1,x1,y2);
   // l3->SetLineWidth(3);
   // l3->Draw("same");
   // TLine *l4 = new TLine(x2,y1,x2,y2);
   // l4->SetLineWidth(3);
   // l4->Draw("same");
  
   //   TGraphErrors *gr_noPID = new TGraphErrors(n, x_noPID, y, 0, ye_noPID);
   TGraphErrors *gr_noPID_proj = new TGraphErrors(n, x, y, 0, ye_noPID);
   gr_noPID_proj->SetMarkerStyle(20);
   gr_noPID_proj->SetMarkerSize(1.5);
   gr_noPID_proj->SetLineWidth(2);
   gr_noPID_proj->Draw("p");

   TGraphErrors *gr_cleanIdeal_proj = new TGraphErrors(n, x_cleanIdeal, y, 0, ye_cleanIdeal);
   gr_cleanIdeal_proj->SetMarkerStyle(24);
   gr_cleanIdeal_proj->SetMarkerSize(1.5);
   gr_cleanIdeal_proj->SetLineWidth(2);
   //   gr_cleanIdeal_proj->Draw("p");


   // leg->SetFillColor(10);
   // leg->SetFillStyle(10);
   // leg->SetLineStyle(4000);
   // leg->SetLineColor(10);
   // leg->SetLineWidth(0.);
   // leg->SetTextFont(42);
   TLegend *leg = new TLegend(0.5, 0.68, 0.92, 0.9);
   leg->SetTextSize(0.045);
   leg->AddEntry("","#it{#bf{sPHENIX}} Simulation","");
   leg->AddEntry("","Au+Au #sqrt{s_{NN}}=200 GeV","");
   leg->AddEntry("","24B 0-10%","");
   leg->Draw();

   // TLatex *tex = new TLatex(3.0, y2*0.88, "Au+Au 200 GeV, 24B 0-10%");
   // tex->SetTextFont(42);
   // tex->SetTextSize(0.055);
   // tex->Draw("same");

   c1->SaveAs(Form("fig/Xsec_Bpm_%s_sPHENIX.eps",cent));
   c1->SaveAs(Form("fig/Xsec_Bpm_%s_sPHENIX.png",cent));
   
}



