#include "draw.C+"
#include "style.C+"
#include "/Users/starsdong/work/work/fitfun/GaussFunctions.C"

void plotTwoDca()
{
  style();

  const Double_t DcaMax = 200;  // +/- 200 micron fit window
  const Double_t um2cm = 1e-4;  // data in cm, simu in um
  const Int_t cm2um = 1e4;
  const Double_t BinWidth = 0.1;  // pT bin desired width
  const Double_t PtMax = 2.0;
  const Int_t NPt = floor(PtMax/BinWidth+0.1);

  TH1D *h1d_data_tmp[NPt], *h1d_data[NPt];
  TFile *f_data = new TFile("pp200_twoTrackDca.root");
  TH2D *h_data = (TH2D *)f_data->Get("h_two_track_pt_dcaxy");
  int NRebin = floor(BinWidth/h_data->GetXaxis()->GetBinWidth(1)+0.1);
  cout << "NRebin = " << NRebin << endl;
  h_data->RebinX(NRebin);
  int i1 = h_data->GetYaxis()->FindBin((-DcaMax+0.1)*um2cm);
  int i2 = h_data->GetYaxis()->FindBin((DcaMax-0.1)*um2cm);
  cout << " i1 = " << i1 << " i2 = " << i2 << endl;
  h_data->GetYaxis()->SetRange(i1,i2);
  for(int i=0;i<NPt;i++) {
    h1d_data_tmp[i] = h_data->ProjectionY(Form("data_tmp_%d",i+1),i+1,i+1);

    h1d_data[i] = new TH1D(Form("data_%d",i+1),"",h1d_data_tmp[i]->GetNbinsX(),h1d_data_tmp[i]->GetXaxis()->GetXmin()*cm2um,h1d_data_tmp[i]->GetXaxis()->GetXmax()*cm2um);
    for(int j=0;j<h1d_data_tmp[i]->GetNbinsX();j++) {
      h1d_data[i]->SetBinContent(j+1, h1d_data_tmp[i]->GetBinContent(j+1));
    }
  }

  const Int_t NMC = 2; // 2 simu
  const Char_t *MCName[20] = {"Gauss","StudentT_nu_1.0"};
  const Char_t *MCLabel[20] = {"Gauss","stu.-t #nu=1"};
  TFile *f_mc[NMC];
  TH2D *h_mc[NMC];
  TH1D *h1d_mc[NMC][NPt];
  for(int i=0;i<NMC;i++) {
    f_mc[i] = new TFile(Form("output_%s.root",MCName[i]));
    h_mc[i] = (TH2D *)f_mc[i]->Get("DCA2");
    int NRebinMC = floor(BinWidth/h_mc[i]->GetXaxis()->GetBinWidth(1)+0.1);
    cout << "NRebinMC = " << NRebinMC << endl;
    h_mc[i]->RebinX(NRebinMC);
    int i1mc = h_mc[i]->GetYaxis()->FindBin(-DcaMax+0.1);
    int i2mc = h_mc[i]->GetYaxis()->FindBin(DcaMax-0.1);
    cout << " i1mc = " << i1mc << " i2mc = " << i2mc << endl;
    h_mc[i]->GetYaxis()->SetRange(i1mc, i2mc);
    for(int j=0;j<NPt;j++) {
      h1d_mc[i][j] = (TH1D *)h_mc[i]->ProjectionY(Form("mc_%d_%d",i,j+1),j+1,j+1);
    }
  }


  /// Now start plotting and fitting
  TCanvas *c1 = new TCanvas("c1","",1000,800);
  c1->Divide(5,4);

  const Int_t kColor[NMC] = {kRed, kBlue};
  TH1D *h1d_mc_p[NMC][NPt]; // scaled histogram for plotting
  for(int i=0;i<NPt;i++) {
    c1->cd(i+1);
    h1d_data[i]->SetTitle("");
    h1d_data[i]->SetMarkerSize(0.8);
    h1d_data[i]->Draw("pe");

    for(int j=0;j<NMC;j++) {
      
      h1d_mc_p[j][i] = (TH1D *)h1d_mc[j][i]->Clone(Form("%s_new",h1d_mc[j][i]->GetName()));
      h1d_mc_p[j][i]->Scale(h1d_data_tmp[i]->GetEntries()/h1d_mc[j][i]->GetEntries());
      h1d_mc_p[j][i]->SetLineColor(kColor[j]);
      h1d_mc_p[j][i]->Draw("histsame");
    }
    c1->Update();
  }
  c1->SaveAs("fig/TwoTrackDcaComp.pdf");
  c1->SaveAs("fig/TwoTrackDcaComp.png");

  // start fitting
  double pT[NPt], sig_data[NPt], sige_data[NPt];
  double sig_mc[NMC][NPt], sige_mc[NMC][NPt];

  TCanvas *c2 = new TCanvas("c2","",1000,800);
  c2->Divide(5,4);
  c2->Draw();
  TF1 *fitfun = new TF1("fitfun",StudentT,-DcaMax,DcaMax,4);
  fitfun->SetParameters(100, 0, 50, 5.);
  for(int i=0;i<NPt;i++) {
    pT[i] = (i+0.5)*PtMax/NPt;
    c2->cd(i+1);
    h1d_data[i]->Draw("pe");
    h1d_data[i]->Fit("fitfun","R");
    sig_data[i] = fitfun->GetParameter(2);
    sige_data[i] = fitfun->GetParError(2);
  }
  c2->Update();
  c2->SaveAs("fig/TwoTrackDca_DataFit.pdf");
  c2->SaveAs("fig/TwoTrackDca_DataFit.png");

  cout << "==================" << endl;
  TCanvas *c3 = new TCanvas("c3","",1000,800);
  c3->Divide(5,4);
  c3->Draw();
  for(int i=0;i<NPt;i++) {
      c3->cd(i+1);
      for(int j=0;j<NMC;j++) {
	fitfun->SetParameters(1000,0.,50.,5.);
	h1d_mc[j][i]->Draw("histsame");
	h1d_mc[j][i]->Fit("fitfun","R");
	fitfun->SetLineColor(kColor[j]);
	fitfun->Clone()->Draw("same");
	sig_mc[j][i] = fitfun->GetParameter(2);
	sige_mc[j][i] = fitfun->GetParError(2);
      }
  }

  c3->Update();
  c3->SaveAs("fig/TwoTrackDca_McFit.pdf");
  c3->SaveAs("fig/TwoTrackDca_McFit.png");
  
  TCanvas *c4 = new TCanvas("c4","",800,600);
  //  c4->SetLogy();
  c4->SetGridx();
  c4->SetGridy();
  c4->Draw();

  TH1D *d0 = new TH1D("d0","",1,0,PtMax);
  d0->SetMaximum(159);
  d0->SetMinimum(0.);
  d0->SetXTitle("Transerverse Momentum p_{T} (GeV/c)");
  d0->SetYTitle("Two-Track #sigma_{DCA} (#mum)");  
  d0->Draw();

  TGraphErrors *gr_mc[NMC];
  for(int i=0;i<NMC;i++) {
    gr_mc[i] = new TGraphErrors(NPt-3, pT+3, sig_mc[i]+3, 0, sige_mc[i]+3);
    gr_mc[i]->SetLineWidth(2);
    gr_mc[i]->SetLineColor(kColor[i]);
    gr_mc[i]->SetFillColor(kColor[i]);
    gr_mc[i]->Draw("e3");
  }
  
  TGraphErrors *gr_data = new TGraphErrors(NPt-2, pT+2, sig_data+2, 0, sige_data+2);
  gr_data->SetMarkerSize(1.5);
  gr_data->SetLineWidth(2);
  gr_data->Draw("p");
  
  TLegend *leg = new TLegend(0.65, 0.7, 0.95, 0.92);
  leg->SetTextSize(0.045);
  leg->AddEntry(gr_data," data @pp200 ", "p");
  for(int i=0;i<NMC;i++)
    leg->AddEntry(gr_mc[i],Form(" MC %s",MCLabel[i]), "l");
  leg->Draw();

  drawText(0.1,10,Form("Track: #sigma_{p_{T}}/p_{T} = 10 %%"));

  c4->Update();
  c4->SaveAs("fig/ResTwoTrackDca_Comp.pdf");
  c4->SaveAs("fig/ResTwoTrackDca_Comp.png");
    
}
