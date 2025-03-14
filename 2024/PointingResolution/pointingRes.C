#include "draw.C+"
#include "style.C+"
#include "/Users/starsdong/work/work/fitfun/Levy.C"
#include "TRandom3.h"

void pointingRes(const Int_t flag = 1, const Int_t nu = 1)
// flag: 1 - student-T sampling,
// flag: 0 - Gaussian sampling
{
  style();
  gSystem->Load("libMathCore");
  const Int_t NConf = 2;
  const Char_t *ConfigName[NConf] = {"Gauss","StudentT"};
  TString OutString;
  if(flag) {
    OutString = Form("%s_nu_%d",ConfigName[flag],nu);
  } else {
    OutString = Form("%s",ConfigName[flag]);
  }

  const Int_t NP = 2; // pi and proton spectra in pp
  const Char_t *Name[NP] = {"pion","proton"};
  const Int_t index[NP] = {2, 12};
  const Double_t Mass[NP] = {0.13957, 0.93827};
  const Int_t kColor[NP] = {kRed, kBlue};

  const Int_t Nevt = 1e5;
  const Int_t Ntrk = 4;  // mean number of tracks per event
  const Int_t NTMax = 50;  // 100 tracks maximum per event
  const Double_t pTres = 0.1;  // 1/pT or pT resolution
  const Double_t pT_th = 0.5;  // pT threshold for vtx
  const Int_t nTrkVtx_th = 3;  // number of tracks for vtx
  gRandom = new TRandom3;
  
  ////////////////////////////////
  // Read in spectra data
  ////////////////////////////////
  TGraphAsymmErrors *gr_data[NP];
  for(int i=0;i<NP;i++) {
    TFile *fin = new TFile(Form("data/HEPData-ins709170-v1-Table_%d.root",index[i]));
    gr_data[i] = (TGraphAsymmErrors *)fin->Get(Form("Table %d/Graph1D_y1",index[i]));
    fin->Close();
  }

  TF1 *func = new TF1("func",Levy,0,5.0,4);
  TF1 *spec[NP];

  ////////////////////////////////
  // Two options for pT and DCA smearing
  ////////////////////////////////
  TF1 *fStudentT = new TF1("studentT",Form("ROOT::Math::tdistribution_pdf(x,%d)",nu),-10.,10.);
  TF1 *fGauss = new TF1("gauss","gaus",-10.,10.);
  fGauss->SetParameters(1,0,1);

  ////////////////////////////////
  // Levy fit to spectra
  ////////////////////////////////
  TCanvas *c1 = new TCanvas("c1","",600,600);
  c1->SetLogy();
  c1->Draw();
  c1->cd();

  TH1D *d1 = new TH1D("d1","",1,0,5);
  d1->SetMinimum(1e-8);
  d1->SetMaximum(10);
  d1->SetXTitle("p_{T} (GeV/c)");
  d1->SetYTitle("1/2#pid^{2}N/p_{T}dp_{T}dy (c^{2}/GeV^{2})");
  d1->Draw();

  for(int i=0;i<NP;i++) {
    gr_data[i]->SetMarkerColor(kColor[i]);
    gr_data[i]->SetLineColor(kColor[i]);
    gr_data[i]->Draw("p");

    func->SetParameters(1, 0.2, 6, Mass[i]);
    func->FixParameter(3, Mass[i]);
    func->SetLineColor(kColor[i]);
    gr_data[i]->Fit("func","R");

    spec[i] = new TF1(Form("spec_%s",Name[i]),Levy_dNdpT,0,5.0,4);
    spec[i]->SetParameters(func->GetParameters());
  }

  TLegend *leg1 = new TLegend(0.7, 0.7, 0.95, 0.9);
  for(int i=0;i<NP;i++) {
    leg1->AddEntry(gr_data[i], Name[i], "p");
  }
  leg1->Draw();
  
  c1->Update();
  c1->SaveAs(Form("fig/invyield_pip_%s.pdf",OutString.Data()));
  c1->SaveAs(Form("fig/invyield_pip_%s.png",OutString.Data()));
  
  ////////////////////////////////
  // Data_rphi resolution from MVTX proposal
  ////////////////////////////////
  const Int_t NPt_DCA = 12;
  const Double_t PT[NPt_DCA] = {0.1, 0.3, 0.5, 0.7, 0.9,
				1.1, 1.3, 1.5, 1.9, 2.9,
				3.9, 4.9};
  const Double_t DCA[NPt_DCA] = {1021, 143, 50.1, 35.8, 29.4,
				 25.4, 23.9, 21.7, 19.1, 15.5,
				 13.5, 12.2};
  TGraph *gr_DCA = new TGraph(NPt_DCA, PT, DCA);
  TCanvas *c1a = new TCanvas("c1a","",200,0,600,600);
  c1a->SetLogy();
  c1a->Draw();
  c1a->cd();

  TH1D *d1a = new TH1D("d1a","",1,0,5);
  d1a->SetMinimum(1);
  d1a->SetMaximum(5000);
  d1a->SetXTitle("p_{T} (GeV/c)");
  d1a->SetYTitle("Resolution r-#phi (#mu m)");
  d1a->Draw();

  gr_DCA->Draw("cp");

  c1a->Update();
  c1a->SaveAs(Form("fig/resolution_simu_%s.pdf",OutString.Data()));
  c1a->SaveAs(Form("fig/resolution_simu_%s.png",OutString.Data()));
  
  ////////////////////////////////
  // plot dN/dpT
  ////////////////////////////////
  TCanvas *c2 = new TCanvas("c2","",400,0,600,600);
  c2->SetLogy();
  c2->Draw();
  c2->cd();

  TH1D *d2 = new TH1D("d2","",1,0,5);
  d2->SetMinimum(1e-8);
  d2->SetMaximum(10);
  d2->SetXTitle("p_{T} (GeV/c)");
  d2->SetYTitle("dN/dp_{T} (a.u.)");
  d2->Draw();
  
  for(int i=0;i<NP;i++) {
    spec[i]->SetLineColor(kColor[i]);
    spec[i]->Draw("csame");
  }

  TLegend *leg2 = new TLegend(0.7, 0.7, 0.95, 0.9);
	     for(int i=0;i<NP;i++) {
    leg2->AddEntry(spec[i], Name[i], "l");
  }
  leg2->Draw();
	     
  c2->Update();
  c2->SaveAs(Form("fig/dNdpT_pip_%s.pdf",OutString.Data()));
  c2->SaveAs(Form("fig/dNdpT_pip_%s.png",OutString.Data()));
  
  
  ////////////////////////////////
  // pT and DCA smearing
  ////////////////////////////////
  
  TH1D *hPtMc[NP], *hPtRc[NP];
  TH2D *hDCAMc[NP], *hDCARc[NP];
  TH2D *hDCAVtxMc[NP], *hDCAVtxRc[NP];
  TH2D *hDCAVtxWtMc[NP], *hDCAVtxWtRc[NP];
  for(int i=0;i<NP;i++) {
    hPtMc[i] = new TH1D(Form("PtMc_%s",Name[i]),"",500,0.,5.0);
    hPtRc[i] = new TH1D(Form("PtRc_%s",Name[i]),"",500,0.,5.0);
    hDCAMc[i] = new TH2D(Form("DCAMc_%s",Name[i]),"",500,0.,5.0,1000,-2500,2500);
    hDCARc[i] = new TH2D(Form("DCARc_%s",Name[i]),"",500,0.,5.0,1000,-2500,2500);

    hDCAVtxMc[i] = new TH2D(Form("DCAVtxMc_%s",Name[i]),"",500,0.,5.0,1000,-2500,2500);
    hDCAVtxRc[i] = new TH2D(Form("DCAVtxRc_%s",Name[i]),"",500,0.,5.0,1000,-2500,2500);

    hDCAVtxWtMc[i] = new TH2D(Form("DCAVtxWtMc_%s",Name[i]),"",500,0.,5.0,1000,-2500,2500);
    hDCAVtxWtRc[i] = new TH2D(Form("DCAVtxWtRc_%s",Name[i]),"",500,0.,5.0,1000,-2500,2500);
  }

  TH2D *hDCAvsNTrk = new TH2D("DCAvsNTrk","",50,0,50,1000,-2500,2500);
  TH2D *hDCAWtvsNTrk = new TH2D("DCAWtvsNTrk","",50,0,50,1000,-2500,2500);
  
  for(int i=0;i<Nevt;i++) {
    Int_t nTrk = 0;
    do {
      nTrk = gRandom->Poisson(Ntrk);
    } while (nTrk>=NTMax);

    if(i%(Nevt/10)==0)
      cout << " Processing " << i << "-th event:  # of tracks = " << nTrk << endl;

    double pT[NP][NTMax], pT_rec[NP][NTMax];
    double dca[NP][NTMax];
    for(int it=0;it<nTrk;it++) {  
      for(int j=0;j<NP;j++) {
	pT[j][it] = spec[j]->GetRandom();
	do {
	    //	pT_rec = 1./gRandom->Gaus(1/pT, 1/pT*pTres);
	    //	pT_rec = gRandom->Gaus(pT, pT*pTres);
	    
	    if(flag) {
	      pT_rec[j][it] = 1./(fStudentT->GetRandom()*1/pT[j][it]*pTres + 1/pT[j][it]);
	    } else {
	      pT_rec[j][it] = 1./(fGauss->GetRandom()*1/pT[j][it]*pTres + 1/pT[j][it]);
	    }
	} while (pT_rec[j][it]<0);
	hPtMc[j]->Fill(pT[j][it]);
	hPtRc[j]->Fill(pT_rec[j][it]);
	    
	double dcaRes = gr_DCA->Eval(pT[j][it]);
	if(flag) {
	  dca[j][it] = fStudentT->GetRandom()*dcaRes;
	} else {
	  dca[j][it] = gRandom->Gaus(0, dcaRes);
	}
	  hDCAMc[j]->Fill(pT[j][it], dca[j][it]);
	  hDCARc[j]->Fill(pT_rec[j][it], dca[j][it]);      
      } // end j
    } // end it
	  
    double vtx = 0.;
    double vtxwt = 0., wt = 0.;  // 1/sigma^2 weighted 
    int ntrk_used = 0;
    for(int it=0;it<nTrk;it++) {
      if(pT_rec[0][it]>pT_th) {
	vtx += dca[0][it];
	double dcaRes = gr_DCA->Eval(pT[0][it]);
	vtxwt += dca[0][it]/dcaRes/dcaRes;
	wt += 1./dcaRes/dcaRes;
	ntrk_used++;
      }
    }
    if(ntrk_used>=nTrkVtx_th) {
      vtx /= ntrk_used;
      hDCAvsNTrk->Fill(ntrk_used, vtx);

      vtxwt /= wt;
      hDCAWtvsNTrk->Fill(ntrk_used, vtxwt);
    } else {
      vtx = -999.;
      vtxwt = -999.;
    }


    if(ntrk_used>=nTrkVtx_th) {
      for(int it=0;it<nTrk;it++) {
	for(int j=0;j<NP;j++) {
	  hDCAVtxMc[j]->Fill(pT[j][it], dca[j][it] - vtx);
	  hDCAVtxRc[j]->Fill(pT_rec[j][it], dca[j][it] - vtx);      

	  hDCAVtxWtMc[j]->Fill(pT[j][it], dca[j][it] - vtxwt);
	  hDCAVtxWtRc[j]->Fill(pT_rec[j][it], dca[j][it] - vtxwt);      
	}
      }
    }
  }

  ////////////////////////////////
  // Check MC and RC pT spectra
  ////////////////////////////////
  TCanvas *c3 = new TCanvas("c3","",800,0,600,600);
  c3->SetLogy();
  c3->Draw();
  c3->cd();

  TH1D *d3 = new TH1D("d3","",1,0,5);
  d3->SetMinimum(1e-1);
  d3->SetMaximum(hPtMc[0]->GetMaximum()*2.0);
  d3->SetXTitle("p_{T} (GeV/c)");
  d3->SetYTitle("Counts");
  d3->Draw();
  
  for(int i=0;i<NP;i++) {
    hPtMc[i]->SetLineStyle(2);
    hPtMc[i]->SetLineColor(kColor[i]);
    hPtMc[i]->Draw("csame");

    hPtRc[i]->SetLineColor(kColor[i]);
    hPtRc[i]->Draw("csame");
  }

  drawText(3.0,hPtMc[0]->GetMaximum()*0.7,Form("#sigma_{p_{T}}/p_{T} = %d %%", (int)(pTres*100)));
  c3->Update();
  c3->SaveAs(Form("fig/pTdist_McRc_%s.pdf",OutString.Data()));
  c3->SaveAs(Form("fig/pTdist_McRc_%s.png",OutString.Data()));

  ////////////////////////////////
  // Check MC and RC DCA distributions, w.r.t. 0
  ////////////////////////////////
  TCanvas *c4 = new TCanvas("c4","",1000,0,600,800);
  c4->Draw();
  c4->cd();
  
  c4->Divide(1,2);
  
  c4->cd(1);
  TH1D *d41 = new TH1D("d41","",1,0,5);
  d41->SetMinimum(-1000);
  d41->SetMaximum(1000);
  d41->Draw();
  
  hDCAMc[0]->SetXTitle("MC p_{T} (GeV/c)");
  hDCAMc[0]->SetYTitle("DCA w.r.t Mc Vtx (#mu m)");
  hDCAMc[0]->Draw("colz");
  hDCAMc[0]->RebinX(10);
  hDCAMc[0]->FitSlicesY();
  
  c4->Update();
  
  c4->cd(2);
  TH1D *d42 = new TH1D("d42","",1,0,5);
  d42->SetMinimum(-1000);
  d42->SetMaximum(1000);
  d42->Draw();
  
  hDCARc[0]->SetXTitle("Rec p_{T} (GeV/c)");
  hDCARc[0]->SetYTitle("DCA w.r.t Mc Vtx (#mu m)");
  hDCARc[0]->Draw("colz");
  hDCARc[0]->RebinX(10);
  hDCARc[0]->FitSlicesY();
  
  drawText(3.0,-2000,Form("#sigma_{p_{T}}/p_{T} = %d %%", (int)(pTres*100)),42,0.065);
  c4->Update();
  c4->SaveAs(Form("fig/DCA_McVtx_%s.pdf",OutString.Data()));
  c4->SaveAs(Form("fig/DCA_McVtx_%s.png",OutString.Data()));
  

  
  TGraphErrors *gr_dca_mc = new TGraphErrors((TH1D *)gDirectory->Get(Form("DCAMc_pion_2")));
  TGraphErrors *gr_dca_rc = new TGraphErrors((TH1D *)gDirectory->Get(Form("DCARc_pion_2")));
  

  ////////////////////////////////
  // DCA resolution comparison
  ////////////////////////////////
  TCanvas *c5 = new TCanvas("c5","",1200,0,600,600);
  c5->SetLogy();
  c5->Draw();
  c5->cd();

  TH1D *d5 = new TH1D("d5","",1,0,5);
  d5->SetMinimum(1);
  d5->SetMaximum(1000);
  d5->SetXTitle("p_{T} (GeV/c)");
  d5->SetYTitle("#sigma_{DCA} w.r.t. MC Vtx (#mu m)");
  d5->Draw();

  gr_dca_mc->SetMarkerStyle(24);
  gr_dca_mc->Draw("p");
  gr_dca_rc->Draw("p");

  TLegend *leg5 = new TLegend(0.7, 0.7, 0.95, 0.9);
  leg5->AddEntry(gr_dca_rc," Rec p_{T}", "p");
  leg5->AddEntry(gr_dca_mc," MC p_{T}", "p");
  leg5->Draw();

  drawText(0.5,2.0,Form("#sigma_{p_{T}}/p_{T} = %d %%", (int)(pTres*100)));

  c5->Update();
  c5->SaveAs(Form("fig/DCA_McVtx_res_%s.pdf",OutString.Data()));
  c5->SaveAs(Form("fig/DCA_McVtx_res_%s.png",OutString.Data()));
  
  ////////////////////////////////
  // Check MC and RC DCA distributions, w.r.t. the Rc vtx
  ////////////////////////////////
  TCanvas *c6 = new TCanvas("c6","",1400,0,600,600);
  c6->Draw();
  c6->cd();

  c6->Divide(1,2);

  c6->cd(1);
  TH1D *d61 = new TH1D("d61","",1,0,5);
  d61->SetMinimum(-1000);
  d61->SetMaximum(1000);
  d61->Draw();

  hDCAVtxMc[0]->SetXTitle("MC p_{T} (GeV/c)");
  hDCAVtxMc[0]->SetYTitle("DCA w.r.t Rec Vtx (#mu m)");
  hDCAVtxMc[0]->Draw("colz");
  hDCAVtxMc[0]->RebinX(10);
  hDCAVtxMc[0]->FitSlicesY();

  TGraphErrors *gr_dcavtx_mc = new TGraphErrors((TH1D *)gDirectory->Get(Form("DCAVtxMc_pion_2")));
  c6->Update();

  c6->cd(2);
  TH1D *d62 = new TH1D("d62","",1,0,5);
  d62->SetMinimum(-1000);
  d62->SetMaximum(1000);
  d62->Draw();

  hDCAVtxRc[0]->SetXTitle("Rec p_{T} (GeV/c)");
  hDCAVtxRc[0]->SetYTitle("DCA w.r.t Rec Vtx (#mu m)");
  hDCAVtxRc[0]->Draw("colz");
  hDCAVtxRc[0]->RebinX(10);
  hDCAVtxRc[0]->FitSlicesY();

  TGraphErrors *gr_dcavtx_rc = new TGraphErrors((TH1D *)gDirectory->Get(Form("DCAVtxRc_pion_2")));


  drawText(3.0,-1500,Form("Vtx: NTrk>= %d, p_{T}>%4.2f", nTrkVtx_th, pT_th),42,0.065);
  drawText(3.0,-2100,Form("Track: #sigma_{p_{T}}/p_{T} = %d %%", (int)(pTres*100)),42,0.065);
  
  c6->Update();
  c6->SaveAs(Form("fig/DCA_RcVtx_%s.pdf",OutString.Data()));
  c6->SaveAs(Form("fig/DCA_RcVtx_%s.png",OutString.Data()));

  ////////////////////////////////
  // DCA resolution comparison
  ////////////////////////////////
  TCanvas *c7 = new TCanvas("c7","",1600,0,600,600);
  c7->SetLogy();
  c7->Draw();
  c7->cd();

  TH1D *d7 = new TH1D("d7","",1,0,5);
  d7->SetMinimum(1);
  d7->SetMaximum(1000);
  d7->SetXTitle("p_{T} (GeV/c)");
  d7->SetYTitle("#sigma_{DCA} w.r.t. Rec Vtx (#mu m)");
  d7->Draw();

  gr_dcavtx_mc->SetMarkerStyle(24);
  gr_dcavtx_mc->Draw("p");
  gr_dcavtx_rc->Draw("p");

  TLegend *leg7 = new TLegend(0.7, 0.7, 0.95, 0.9);
  leg7->AddEntry(gr_dca_rc," Rec p_{T}", "p");
  leg7->AddEntry(gr_dca_mc," MC p_{T}", "p");
  leg7->Draw();

  drawText(0.5,5.0,Form("Vtx: NTrk>= %d, p_{T}>%4.2f", nTrkVtx_th, pT_th));
  drawText(0.5,2.0,Form("Track: #sigma_{p_{T}}/p_{T} = %d %%", (int)(pTres*100)));

  c7->Update();
  c7->SaveAs(Form("fig/DCA_RcVtx_res_%s.pdf",OutString.Data()));
  c7->SaveAs(Form("fig/DCA_RcVtx_res_%s.png",OutString.Data()));

  ////////////////////////////////
  // Vtx resolution
  ////////////////////////////////
  TCanvas *c8 = new TCanvas("c8","",1800,0,600,600);
  //  c8->SetLogy();
  c8->Draw();
  c8->cd();
  
  TH1D *d8 = new TH1D("d8","",1,0,20);
  d8->SetMinimum(-1000);
  d8->SetMaximum(1000);
  d8->SetXTitle("Number of Tracks");
  d8->SetYTitle("#Delta Vtx (#mu m)");
  d8->Draw();

  hDCAvsNTrk->Draw("same");

  drawText(5,-900,Form("Vtx: NTrk>= %d, p_{T}>%4.2f", nTrkVtx_th, pT_th));
  
  c8->Update();
  c8->SaveAs(Form("fig/VtxDiff_NTrk_%s.pdf",OutString.Data()));
  c8->SaveAs(Form("fig/VtxDiff_NTrk_%s.png",OutString.Data()));

  hDCAvsNTrk->FitSlicesY();
  TGraphErrors *gr_vtx = new TGraphErrors((TH1D *)gDirectory->Get(Form("DCAvsNTrk_2")));
  
  TCanvas *c9 = new TCanvas("c9","",2000,0,600,600);
  c9->SetLogy();
  c9->Draw();
  c9->cd();
  
  TH1D *d9 = new TH1D("d9","",1,0,20);
  d9->SetMinimum(1);
  d9->SetMaximum(1000);
  d9->SetXTitle("Number of Tracks");
  d9->SetYTitle("#sigma_{vtx} (#mu m)");
  d9->Draw();

  gr_vtx->Draw("p");

  drawText(5,3,Form("Vtx: NTrk>= %d, p_{T}>%4.2f", nTrkVtx_th, pT_th));

  c9->Update();
  c9->SaveAs(Form("fig/VtxRes_NTrk_%s.pdf",OutString.Data()));
  c9->SaveAs(Form("fig/VtxRes_NTrk_%s.png",OutString.Data()));
  

  ////////////////////////////////
  // Check MC and RC DCA distributions, w.r.t. the weighted Rc vtx
  ////////////////////////////////
  TCanvas *c10 = new TCanvas("c10","",2200,0,600,600);
  c10->Draw();
  c10->cd();

  c10->Divide(1,2);

  c10->cd(1);
  TH1D *d101 = new TH1D("d101","",1,0,5);
  d101->SetMinimum(-1000);
  d101->SetMaximum(1000);
  d101->Draw();

  hDCAVtxWtMc[0]->SetXTitle("MC p_{T} (GeV/c)");
  hDCAVtxWtMc[0]->SetYTitle("DCA to Weighted Rec Vtx (#mu m)");
  hDCAVtxWtMc[0]->Draw("colz");
  hDCAVtxWtMc[0]->RebinX(10);
  hDCAVtxWtMc[0]->FitSlicesY();

  TGraphErrors *gr_dcavtxwt_mc = new TGraphErrors((TH1D *)gDirectory->Get(Form("DCAVtxWtMc_pion_2")));
  c10->Update();

  c10->cd(2);
  TH1D *d102 = new TH1D("d102","",1,0,5);
  d102->SetMinimum(-1000);
  d102->SetMaximum(1000);
  d102->Draw();

  hDCAVtxWtRc[0]->SetXTitle("Rec p_{T} (GeV/c)");
  hDCAVtxWtRc[0]->SetYTitle("DCA to Weighted Rec Vtx (#mu m)");
  hDCAVtxWtRc[0]->Draw("colz");
  hDCAVtxWtRc[0]->RebinX(10);
  hDCAVtxWtRc[0]->FitSlicesY();

  TGraphErrors *gr_dcavtxwt_rc = new TGraphErrors((TH1D *)gDirectory->Get(Form("DCAVtxWtRc_pion_2")));


  drawText(3.0,-1500,Form("Wt. Vtx: NTrk>= %d, p_{T}>%4.2f", nTrkVtx_th, pT_th),42,0.065);
  drawText(3.0,-2100,Form("Track: #sigma_{p_{T}}/p_{T} = %d %%", (int)(pTres*100)),42,0.065);
  
  c10->Update();
  c10->SaveAs(Form("fig/DCA_RcVtxWt_%s.pdf",OutString.Data()));
  c10->SaveAs(Form("fig/DCA_RcVtxWt_%s.png",OutString.Data()));

  ////////////////////////////////
  // DCA resolution comparison w.r.t. weighted vtx
  ////////////////////////////////
  TCanvas *c11 = new TCanvas("c11","",2400,0,600,600);
  c11->SetLogy();
  c11->Draw();
  c11->cd();

  TH1D *d11 = new TH1D("d11","",1,0,5);
  d11->SetMinimum(1);
  d11->SetMaximum(1000);
  d11->SetXTitle("p_{T} (GeV/c)");
  d11->SetYTitle("#sigma_{DCA} to Wt. Rec Vtx (#mu m)");
  d11->Draw();

  gr_dcavtxwt_mc->SetMarkerStyle(24);
  gr_dcavtxwt_mc->Draw("p");
  gr_dcavtxwt_rc->Draw("p");

  TLegend *leg11 = new TLegend(0.7, 0.7, 0.95, 0.9);
  leg11->AddEntry(gr_dcavtxwt_rc," Rec p_{T}", "p");
  leg11->AddEntry(gr_dcavtxwt_mc," MC p_{T}", "p");
  leg11->Draw();

  drawText(0.5,5.0,Form("Wt. Vtx: NTrk>= %d, p_{T}>%4.2f", nTrkVtx_th, pT_th));
  drawText(0.5,2.0,Form("Track: #sigma_{p_{T}}/p_{T} = %d %%", (int)(pTres*100)));

  c11->Update();
  c11->SaveAs(Form("fig/DCA_RcVtxWt_res_%s.pdf",OutString.Data()));
  c11->SaveAs(Form("fig/DCA_RcVtxWt_res_%s.png",OutString.Data()));

  ////////////////////////////////
  // Weighted Vtx resolution
  ////////////////////////////////
  TCanvas *c12 = new TCanvas("c12","",2600,0,600,600);
  //  c12->SetLogy();
  c12->Draw();
  c12->cd();
  
  TH1D *d12 = new TH1D("d12","",1,0,20);
  d12->SetMinimum(-1000);
  d12->SetMaximum(1000);
  d12->SetXTitle("Number of Tracks");
  d12->SetYTitle("#Delta Vtx (#mu m)");
  d12->Draw();

  hDCAWtvsNTrk->Draw("same");

  drawText(5,-900,Form("Wt. Vtx: NTrk>= %d, p_{T}>%4.2f", nTrkVtx_th, pT_th));
  
  c12->Update();
  c12->SaveAs(Form("fig/VtxWtDiff_NTrk_%s.pdf",OutString.Data()));
  c12->SaveAs(Form("fig/VtxWtDiff_NTrk_%s.png",OutString.Data()));

  hDCAWtvsNTrk->FitSlicesY();
  TGraphErrors *gr_vtxwt = new TGraphErrors((TH1D *)gDirectory->Get(Form("DCAWtvsNTrk_2")));
  
  TCanvas *c13 = new TCanvas("c13","",2800,0,600,600);
  c13->SetLogy();
  c13->Draw();
  c13->cd();
  
  TH1D *d13 = new TH1D("d13","",1,0,20);
  d13->SetMinimum(1);
  d13->SetMaximum(1000);
  d13->SetXTitle("Number of Tracks");
  d13->SetYTitle("#sigma_{vtx} (#mu m)");
  d13->Draw();

  gr_vtxwt->Draw("p");

  drawText(5,3,Form("Wt. Vtx: NTrk>= %d, p_{T}>%4.2f", nTrkVtx_th, pT_th));

  c13->Update();
  c13->SaveAs(Form("fig/VtxWtRes_NTrk_%s.pdf",OutString.Data()));
  c13->SaveAs(Form("fig/VtxWtRes_NTrk_%s.png",OutString.Data()));

  ////////////////////////////////
  // Output histograms
  ////////////////////////////////
  TFile *fout = new TFile("output.root","recreate");
  for(int i=0;i<NP;i++) {
    hPtMc[i]->Write();
    hPtRc[i]->Write();
    hDCAMc[i]->Write();
    hDCARc[i]->Write();
    hDCAVtxMc[i]->Write();
    hDCAVtxRc[i]->Write();
    hDCAVtxWtMc[i]->Write();
    hDCAVtxWtRc[i]->Write();
  }
  hDCAvsNTrk->Write();
  hDCAWtvsNTrk->Write();
  gr_dca_mc->Write();
  gr_dca_rc->Write();
  gr_dcavtx_mc->Write();
  gr_dcavtx_rc->Write();
  
  fout->Close();
  
  
}
