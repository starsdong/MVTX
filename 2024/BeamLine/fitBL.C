#include "draw.C+"
#include "style.C+"


void fitBL()
{
  style();

  const Int_t NL = 3;
  const Double_t R[NL] = {25.233, 33.355, 41.478};

  TH1F *hOccu[NL];
  TFile *fin = new TFile("HIST_PHYSICS_run2pp_new_2024p001-00046667-0000.root");
  for(int i=0;i<NL;i++) {
    hOccu[i] = (TH1F*)fin->Get(Form("h_MvtxClusterQA_clusterPhi_l%d",i));
  }

  TF1 *func[NL];
  for(int i=0;i<NL;i++) {
    func[i] = new TF1(Form("func_%d",i),"[0]*(1.+2*[1]/[3]*cos(x-[2]))",-TMath::Pi(),TMath::Pi());
    //    func[i] = new TF1(Form("func_%d",i),"[0]/(1.-2*[1]/[3]*cos(x-[2])+[1]*[1]/[3]/[3]*cos(x-[2])*cos(x-[2]))",-TMath::Pi(),TMath::Pi());   // second order
  }

  TCanvas *c1 = new TCanvas("c1","c1",600,1000);
  c1->Draw();
  TPad *pad[NL];

  double d0[NL], d0_err[NL];
  double phi0[NL], phi0_err[NL];
  for(int i=0;i<NL;i++) {
    c1->cd();
    pad[i] = new TPad(Form("pad_%d",i),"",0.,(2-i)*1./NL,1.0,(3-i)*1./NL);
    pad[i]->SetLeftMargin(0.12);
    pad[i]->Draw();
    pad[i]->cd();
    
    hOccu[i]->SetMaximum(2000);
    hOccu[i]->SetMarkerSize(0.8);
    hOccu[i]->GetXaxis()->CenterTitle();
    hOccu[i]->GetXaxis()->SetTitleOffset(0.9);
    hOccu[i]->GetXaxis()->SetTitleSize(0.08);
    hOccu[i]->GetYaxis()->SetNdivisions(105);
    hOccu[i]->GetYaxis()->SetTitleOffset(0.7);
    hOccu[i]->GetYaxis()->SetTitleSize(0.08);
    hOccu[i]->SetTitle("");
    hOccu[i]->Draw();
    func[i]->SetParameters(hOccu[i]->GetBinContent(1), 0, TMath::Pi(), R[i],0);
    func[i]->FixParameter(3, R[i]);
    hOccu[i]->Fit(Form("func_%d", i),"R");

    d0[i] = func[i]->GetParameter(1);
    d0_err[i] = func[i]->GetParError(1);
    phi0[i] = func[i]->GetParameter(2);
    phi0_err[i] = func[i]->GetParError(2);

    drawText(-2, 400, Form("d = %4.2f #pm %4.2f mm", d0[i], d0_err[i]), 42, 0.06);
    drawText(-2, 200, Form("#phi_{0} = %4.2f #pm %4.2f rad", phi0[i], phi0_err[i]), 42, 0.06);
    drawText(0.5, 200, Form("(X_{0}, Y_{0}) = (%4.2f, %4.2f) mm", d0[i]*cos(phi0[i]), d0[i]*sin(phi0[i])), 42, 0.06);

    pad[i]->Modified();
    c1->Update();
  }

  c1->Update();
  c1->SaveAs("fig/MVTX_occupancy_BL.pdf");
  c1->SaveAs("fig/MVTX_occupancy_BL.png");
}
