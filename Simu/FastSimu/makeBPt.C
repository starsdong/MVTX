void makeBPt()
{
  gROOT->Reset();

  TFile *fin = new TFile("B0_clean.root");
  TH1D *hPt = (TH1D *)fin->Get("hptWg");
  double N = hPt->Integral(1,hPt->GetNbinsX());
  double dpt = hPt->GetXaxis()->GetBinWidth(1);
  hPt->Scale(1./dpt/N);
  double tot = hPt->Integral(1,200);
  double int_1 = hPt->Integral(hPt->GetXaxis()->FindBin(1.001),200)/tot;
  double int_2 = hPt->Integral(hPt->GetXaxis()->FindBin(2.001),200)/tot;
  cout << int_1 << " " << int_2 << endl;
  
  TGraph *gr = new TGraph(hPt);
  
  TCanvas *c1 = new TCanvas("c1", "c1",0,0,700,600);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0.01);
  c1->SetFillColor(10);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetFrameFillColor(0);
  c1->SetFrameBorderMode(0);

  // c1->SetGridx();
  // c1->SetGridy();
  
  c1->SetLeftMargin(0.16);
  c1->SetBottomMargin(0.14);
  c1->SetTopMargin(0.02);
  c1->SetRightMargin(0.02);
  
  double x1 = 0.01;
  double x2 = 5.9;
  double y1 = 0.0;
  double y2 = 0.599;
  
  TH1D *h0 = new TH1D("h0","",1,x1,x2);
  h0->SetMinimum(y1);
  h0->SetMaximum(y2);
  h0->GetXaxis()->SetNdivisions(208);
  h0->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h0->GetXaxis()->SetTitleOffset(1.0);
  h0->GetXaxis()->SetTitleSize(0.06);
  h0->GetXaxis()->SetLabelOffset(0.005);
  h0->GetXaxis()->SetLabelSize(0.045);
  h0->GetXaxis()->SetLabelFont(42);
  h0->GetXaxis()->SetTitleFont(42);
  h0->GetYaxis()->SetNdivisions(204);
  h0->GetYaxis()->SetTitle("1/N dN/dp_{T}");
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

  gr->SetLineWidth(2);
  gr->Draw("c");

  TLatex * tex = new TLatex(2.8, 0.53, "FONLL #sqrt{s} = 200 GeV");
  tex->SetTextFont(42);
  tex->SetTextSize(0.045);
  tex->Draw("same");

  TLatex * tex = new TLatex(3.8, 0.47, "D^{0} from B");
  tex->SetTextFont(42);
  tex->SetTextSize(0.045);
  tex->Draw("same");
  
  TLine *la = new TLine(1, y1, 1, y2*0.8);
  la->SetLineStyle(2);
  la->SetLineWidth(2);
  la->Draw("same");
  TArrow *aa = new TArrow(1, y2*0.75, 1.6, y2*0.75, 0.05, "|>");
  aa->SetLineStyle(1);
  aa->SetAngle(30);
  aa->SetLineWidth(2);
  aa->Draw();
  TLatex * tex = new TLatex(1.4, y2*0.77, "75%");
  tex->SetTextFont(22);
  tex->SetTextSize(0.04);
  tex->Draw("same");
  
  TLine *lb = new TLine(2, y1, 2, y2*0.8);
  lb->SetLineStyle(2);
  lb->SetLineWidth(2);
  lb->Draw("same");
  TArrow *ab = new TArrow(2, y2*0.75, 2.6, y2*0.75, 0.05, "|>");
  ab->SetLineStyle(1);
  ab->SetLineWidth(2);
  ab->SetAngle(30);
  ab->Draw();
  TLatex * tex = new TLatex(2.4, y2*0.77, "35%");
  tex->SetTextFont(22);
  tex->SetTextSize(0.04);
  tex->Draw("same");

  c1->Update();
  c1->SaveAs("fig/D0_B_Pt.eps");
  c1->SaveAs("fig/D0_B_Pt.png");
  
}
