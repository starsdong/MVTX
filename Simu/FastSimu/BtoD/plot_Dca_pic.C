#include "myFunction.h"
#include "myConst.h"

void plot_Dca_pic()
{
    globalSetting();
    char dir[250];
    char name[250];
    char CMD[250];
    char lname[16][250];
    TPaveStats* ptstates;
    TLegend* legend;
    TH1F* h0;
    
    int   const npt = 8;
    double const nptbin[npt+1] = {0, 1.0, 2.0, 3.0, 4.0, 5.0, 8.0, 9.0, 10.0};
    int ptbin_lw, ptbin_up;

    //read
    //D0
    TFile* fD0 = new TFile("data/D0_clean.root");
    TH1D* hD0 = (TH1D*)fD0->Get("hD0");
    hD0->SetDirectory(0);
    TH2F* h2D0DcaSig = (TH2F*)fD0->Get("h2D0DcaSig");
    h2D0DcaSig->SetDirectory(0);
    TH2F* h2DcaBg = (TH2F*)fD0->Get("h2D0DcaBg");
    h2DcaBg->SetDirectory(0);
    //B0
    TFile* fB0 = new TFile("data/B0_clean.root");
    TH1D* hB0 = (TH1D*)fB0->Get("hD0");
    hB0->SetDirectory(0);
    TH2F* h2B0DcaSig = (TH2F*)fB0->Get("h2D0DcaSig");
    h2B0DcaSig->SetDirectory(0);
    fB0->Close();
    //Bpm
    TFile* fBpm = new TFile("data/Bpm_clean.root");
    TH1D* hBpm = (TH1D*)fBpm->Get("hD0");
    hBpm->SetDirectory(0);
    TH2F* h2BpmDcaSig = (TH2F*)fBpm->Get("h2D0DcaSig");
    h2BpmDcaSig->SetDirectory(0);
    fBpm->Close();
    
    //scale background (++,--,+-,-+)
    h2DcaBg->Scale(0.5);
    
    //combine B and scale it properly
    //no any cuts
    double neventD0 = hD0->GetEntries();
    double neventB0 = hB0->GetEntries();
    double neventBpm = hBpm->GetEntries();
    cout << "Number of D0:\t" << neventD0 << endl;
    cout << "Number of B+/-:\t" << neventB0 << endl;
    cout << "Number of B0:\t" << neventBpm << endl;
    double FR1 = 0.4; //b->B+
    double BR1 = 0.086+0.79;//0.086: B+ -> D0 X; 0.79: B+ -> D0bar X
    double FR2 = 0.4; //b->B0
    double BR2 = 0.081+0.474; // 0.081: B0->D0 X, 0.474: B0->D0bar X
    double scale_Bpm = FR1*BR1 * neventD0/neventBpm; //Bpm
    double scale_B0  = FR2*BR2 * neventD0/neventB0;  // B0
    double scale_D0  = 1.0;  //D0
    cout << "scale for D0:\t" << scale_D0 << endl;
    cout << "scale for B+/-:\t" << scale_Bpm << endl;
    cout << "scale for B0:\t" << scale_B0 << endl;
    h2D0DcaSig->Scale(scale_D0);
    h2B0DcaSig->Scale(scale_B0);
    h2BpmDcaSig->Scale(scale_Bpm);
    TH2F* h2BDcaSig = new TH2F(*h2B0DcaSig);
    h2BDcaSig->Add(h2B0DcaSig,h2BpmDcaSig,1,1);
    
    //plot
    sprintf(dir,"pic/Dca_clean");
    sprintf(CMD,"[ -d %s ] || mkdir -p %s", dir, dir);
    gSystem->Exec(CMD);
    
    TH1F* hD0Dca;
    TH1F* hBDca;
    TH1F* hDcaBg;
    
    TCanvas* c1 = new TCanvas("c1", "A Canvas",10,10,800,800);
    setPad(c1);
    
    for(int ipt=0; ipt<npt; ipt++) {
        ptbin_lw = h2D0DcaSig->GetXaxis()->FindBin(nptbin[ipt]+1e-6);
        ptbin_up = h2D0DcaSig->GetXaxis()->FindBin(nptbin[ipt+1]-1e-6);
        hD0Dca = (TH1F*)h2D0DcaSig->ProjectionY(Form("hDcaSigD0%i",ipt),ptbin_lw,ptbin_up);
        hBDca = (TH1F*)h2BDcaSig->ProjectionY(Form("hDcaSigB%i",ipt),ptbin_lw,ptbin_up);
        hDcaBg = (TH1F*)h2DcaBg->ProjectionY(Form("hDcaBg%i",ipt),ptbin_lw,ptbin_up);
        
        //Set
        setHisto(hD0Dca, "", "Dca_{D^{0}} (cm)", "Counts");
        hD0Dca->SetLineColor(COLOR[0]);
        hD0Dca->SetLineWidth(2);
        hD0Dca->SetMaximum(10*hD0Dca->GetMaximum());
        
        hBDca->SetLineColor(COLOR[1]);
        hBDca->SetLineWidth(2);
        
        hDcaBg->SetLineColor(COLOR[2]);
        hDcaBg->SetLineWidth(2);
        
        //Draw
        gPad->SetLogy();
        hD0Dca->Draw("HIST");
        hBDca->Draw("HISTSAME");
        hDcaBg->Draw("HISTSAME");
        sprintf(name,"%.0f #leq p_{T} < %.0f",nptbin[ipt],nptbin[ipt+1]);
        drawLatex(0.2,0.9,name,22,0.05,4);
        sprintf(name,"0-10%s","%");
        drawLatex(0.2,0.83,name,22,0.05,4);
        legend = new TLegend(0.45,0.8,0.9,0.95);
        legend->SetFillStyle(0);
        legend->SetFillColor(10);
        legend->SetBorderSize(0);
        legend->SetTextSize(0.04);
        //legend->SetHeader("0-10%");
        legend->AddEntry(hD0Dca,"prompt D^{0}","l");
        legend->AddEntry(hBDca,"non-prompt D^{0}","l");
        legend->AddEntry(hDcaBg,"background","l");
        legend->Draw();
        
        sprintf(name,"%s/Dca_pT_%.0f_%.0f.gif",dir,nptbin[ipt],nptbin[ipt+1]);
        c1->SaveAs(name);
    }
}
