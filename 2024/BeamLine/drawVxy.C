#include "draw.C+"
#include "style.C+"

void drawVxy()
{
  style();

  TFile *fin = new TFile("HIST_PHYSICS_run2pp_new_2024p001-00046667-0000.root");

  TH2F *hVxy = (TH2F *)fin->Get("h_SiliconSeedsQA_vx_vy");
  TH1F *hVx = (TH1F *)fin->Get("h_SiliconSeedsQA_vx");
  TH1F *hVy = (TH1F *)fin->Get("h_SiliconSeedsQA_vy");

  TCanvas *c1 = new TCanvas("c1","c1",1000,600);
  c1->Draw();

  TF1 *fg[2];
  for(int i=0;i<2;i++) {
    fg[i] = new TF1(Form("gauss_%d",i),"gaus",-2.,2.);
    fg[i]->SetParameters(1., 0., 0.1);
  }
  TPad *p1 = new TPad("p1","",0,0,0.6,1.0);
  p1->SetLeftMargin(0.15);
  p1->SetBottomMargin(0.15);
  p1->SetTopMargin(0.04);
  p1->SetRightMargin(0.04);
  p1->Draw();
  p1->cd();
  hVxy->SetTitle("");
  hVxy->GetXaxis()->SetNdivisions(105);
  hVxy->GetYaxis()->SetNdivisions(105);
  hVxy->Draw();
  
  p1->Modified();
  c1->Update();
  c1->cd();

  TPad *p2 = new TPad("p2","",0.6,0.5,1.0,1.0);
  p2->SetLeftMargin(0.2);
  p2->SetBottomMargin(0.2);
  p2->SetTopMargin(0.05);
  p2->SetRightMargin(0.05);
  p2->Draw();
  p2->cd();

  hVx->SetTitle("");
  hVx->GetXaxis()->SetNdivisions(105);
  hVx->GetXaxis()->SetTitleOffset(0.9);
  hVx->GetXaxis()->SetTitleSize(0.07);
  hVx->GetYaxis()->SetNdivisions(105);
  hVx->GetYaxis()->SetTitleOffset(0.9);
  hVx->GetYaxis()->SetTitleSize(0.07);
  hVx->Draw();
  hVx->Fit(fg[0],"R");

  drawText(0.5, hVx->GetMaximum()*0.8, Form("Vx_{0} = %4.2f mm", fg[0]->GetParameter(1)*10.), 42, 0.06);
  p2->Modified();
  c1->Update();
  c1->cd();
  
  TPad *p3 = new TPad("p3","",0.6,0.0,1.0,0.5);
  p3->SetLeftMargin(0.2);
  p3->SetBottomMargin(0.2);
  p3->SetTopMargin(0.05);
  p3->SetRightMargin(0.05);
  p3->Draw();
  p3->cd();
  
  hVy->SetTitle("");
  hVy->GetXaxis()->SetNdivisions(105);
  hVy->GetXaxis()->SetTitleOffset(0.9);
  hVy->GetXaxis()->SetTitleSize(0.07);
  hVy->GetYaxis()->SetNdivisions(105);
  hVy->GetYaxis()->SetTitleOffset(0.9);
  hVy->GetYaxis()->SetTitleSize(0.07);
  hVy->Draw();
  hVy->Fit(fg[1],"R");

  drawText(0.5, hVy->GetMaximum()*0.8, Form("Vy_{0} = %4.2f mm", fg[1]->GetParameter(1)*10.), 42, 0.06);
  p3->Modified();
  c1->Update();
  c1->cd();

  c1->SaveAs("fig/Vxy.pdf");
  c1->SaveAs("fig/Vxy.png");
}
