#include "draw.C+"
#include "style.C+"


void fitBL(const Int_t run = 46667)
{
  style();

  const Int_t NL = 3;
  const Double_t R[NL] = {25.233, 33.355, 41.478};

  TH1F *hOccu[NL];
  //  TFile *fin = new TFile(Form("root/HIST_PHYSICS_run2pp_new_2024p001-000%d-0000.root",run));
  TFile *fin = new TFile(Form("root/HIST_DST_TRKR_CLUSTER_run2pp_new_2024p004-000%d-9000.root",run));
  for(int i=0;i<NL;i++) {
    hOccu[i] = (TH1F*)fin->Get(Form("h_MvtxClusterQA_clusterPhi_l%d",i));
  }

  TF1 *func[NL];
  TF1 *funcm[NL];
  for(int i=0;i<NL;i++) {
    func[i] = new TF1(Form("func_%d",i),"[0]*(1.+2*[1]/[3]*cos(x-[2]))",-TMath::Pi(),TMath::Pi());
    //    func[i] = new TF1(Form("func_%d",i),"[0]/(1.-2*[1]/[3]*cos(x-[2])+[1]*[1]/[3]/[3]*cos(x-[2])*cos(x-[2]))",-TMath::Pi(),TMath::Pi());   // second order
    funcm[i] = new TF1(Form("funcm_%d",i),"[0]*(1.+2*[1]/[3]*cos(x)+2*[2]/[3]*sin(x))",-TMath::Pi(),TMath::Pi());
  }

  TCanvas *c1 = new TCanvas("c1","c1",600,1000);
  c1->Draw();
  TPad *pad[NL];

  double d0[NL], d0_err[NL];
  double phi0[NL], phi0_err[NL];
  double x0[NL], x0_err[NL];
  double y0[NL], y0_err[NL];
  for(int i=0;i<NL;i++) {
    c1->cd();
    pad[i] = new TPad(Form("pad_%d",i),"",0.,(2-i)*1./NL,1.0,(3-i)*1./NL);
    pad[i]->SetLeftMargin(0.12);
    pad[i]->Draw();
    pad[i]->cd();

    if(!hOccu[i]) continue;
    double ymax = (floor(hOccu[0]->GetEntries()/hOccu[0]->GetNbinsX()*2./1000)+0.4)*1000;
    hOccu[i]->SetMaximum(ymax);
    hOccu[i]->SetMarkerSize(0.8);
    hOccu[i]->GetXaxis()->CenterTitle();
    hOccu[i]->GetXaxis()->SetTitleOffset(0.9);
    hOccu[i]->GetXaxis()->SetTitleSize(0.08);
    hOccu[i]->GetYaxis()->SetNdivisions(105);
    hOccu[i]->GetYaxis()->SetTitleOffset(0.7);
    hOccu[i]->GetYaxis()->SetTitleSize(0.08);
    hOccu[i]->SetTitle("");
    hOccu[i]->Draw();

    // pre-fit, remove ineffecient entries
    func[i]->SetParameters(hOccu[i]->GetBinContent(1), 4.0, TMath::Pi(), R[i],0);
    func[i]->FixParameter(1, 4.0);
    func[i]->FixParameter(2, 3.0);
    func[i]->FixParameter(3, R[i]);
    hOccu[i]->Fit(Form("func_%d", i), "R", "", -3.1, 3.1);

    // first iteration filtering - window (0.5, 2.0)
    for(int j=0;j<hOccu[i]->GetNbinsX();j++) {
      if(hOccu[i]->GetBinContent(j+1) < 0.5*func[i]->Eval(hOccu[i]->GetBinCenter(j+1)) ||
	 hOccu[i]->GetBinContent(j+1) > 2*func[i]->Eval(hOccu[i]->GetBinCenter(j+1)) ) {
	hOccu[i]->SetBinContent(j+1, 0);
	hOccu[i]->SetBinError(j+1, 0);
      }
    }
    //    continue;

    func[i]->SetParameters(hOccu[i]->GetBinContent(1), 4.0, TMath::Pi(), R[i],0);
    func[i]->ReleaseParameter(1);
    func[i]->ReleaseParameter(2);
    func[i]->FixParameter(3, R[i]);
    hOccu[i]->Fit(Form("func_%d", i), "R", "", -3.1, 3.1);

    // second iteration filtering - window (0.8, 1.2)
    for(int j=0;j<hOccu[i]->GetNbinsX();j++) {
      if(hOccu[i]->GetBinContent(j+1) < 0.8*func[i]->Eval(hOccu[i]->GetBinCenter(j+1)) ||
	 hOccu[i]->GetBinContent(j+1) > 1.2*func[i]->Eval(hOccu[i]->GetBinCenter(j+1)) ) {
	hOccu[i]->SetBinContent(j+1, 0);
	hOccu[i]->SetBinError(j+1, 0);
      }
    }
    
    func[i]->SetParameters(hOccu[i]->GetBinContent(1), 4.0, TMath::Pi(), R[i],0);
    func[i]->ReleaseParameter(1);
    func[i]->ReleaseParameter(2);
    func[i]->FixParameter(3, R[i]);
    hOccu[i]->Fit(Form("func_%d", i), "R", "", -3.1, 3.1);

    funcm[i]->SetParameters(hOccu[i]->GetBinContent(1), -4.0, 0, R[i],0);
    funcm[i]->SetParLimits(1, -5.0, -3.0);
    funcm[i]->SetParLimits(2, -0.5, 1.5);
    funcm[i]->FixParameter(3, R[i]);
    hOccu[i]->Fit(Form("funcm_%d", i), "R", "", -3.1, 3.1);

    d0[i] = func[i]->GetParameter(1);
    d0_err[i] = func[i]->GetParError(1);
    phi0[i] = func[i]->GetParameter(2);
    phi0_err[i] = func[i]->GetParError(2);

    x0[i] = funcm[i]->GetParameter(1);
    x0_err[i] = funcm[i]->GetParError(1);
    y0[i] = funcm[i]->GetParameter(2);
    y0_err[i] = funcm[i]->GetParError(2);
    

    drawText(-0.7, ymax*0.88, Form("Run # %d", run), 42, 0.07);
    drawText(-2, ymax*0.16, Form("d = %4.2f #pm %4.2f mm", d0[i], d0_err[i]), 42, 0.06);
    drawText(-2, ymax*0.08, Form("#phi_{0} = %4.2f #pm %4.2f rad", phi0[i], phi0_err[i]), 42, 0.06);
    drawText(0.5, ymax*0.16, Form("X_{0} = %4.2f #pm %4.2f mm", x0[i], x0_err[i]), 42, 0.06);
    drawText(0.5, ymax*0.08, Form("Y_{0} = %4.2f #pm %4.2f mm", y0[i], y0_err[i]), 42, 0.06);
    //    drawText(0.5, 200, Form("(X_{0}, Y_{0}) = (%4.2f, %4.2f) mm", d0[i]*cos(phi0[i]), d0[i]*sin(phi0[i])), 42, 0.06);

    pad[i]->Modified();
    c1->Update();
  }

  c1->Update();
  c1->SaveAs(Form("fig/MVTX_BL_%d.pdf",run));
  c1->SaveAs(Form("fig/MVTX_BL_%d.png",run));

  ofstream outData;
  outData.open(Form("dat/MVTX_BL_%d.txt",run));
  for(int i=0;i<NL;i++) {
    outData << setw(12) << x0[i] << setw(12) << x0_err[i]
	    << setw(12) << y0[i] << setw(12) << y0_err[i]
	    << setw(12) << d0[i] << setw(12) << d0_err[i]
	    << setw(12) << phi0[i] << setw(12) << phi0_err[i] << endl;
  }
  outData.close();
}
