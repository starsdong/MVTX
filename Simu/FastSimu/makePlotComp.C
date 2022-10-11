void makePlotComp()
{
  TFile *fin = new TFile("root/significance.root");
  TGraph *gr_old_D_MB = (TGraph *)fin->Get("sig_d0_MB_noPid");
  TGraph *gr_old_D_B_MB = (TGraph *)fin->Get("sig_d0_B_MB_noPid");
  TGraph *gr_old_D_Cen = (TGraph *)fin->Get("sig_d0_Cen_noPid");
  TGraph *gr_old_D_B_Cen = (TGraph *)fin->Get("sig_d0_B_Cen_noPid");
  TGraph *gr_old_D_Peri = (TGraph *)fin->Get("sig_d0_Peri_noPid");
  TGraph *gr_old_D_B_Peri = (TGraph *)fin->Get("sig_d0_B_Peri_noPid");
  fin->Close();

  
  TGraph *gr_new_D_MB = new TGraph("dat/SigPromptD0_0_80.txt","%lg %lg");
  TGraph *gr_new_D_B_MB = new TGraph("dat/SigNonPromptD0_0_80.txt","%lg %lg");
  TGraph *gr_new_D_Cen = new TGraph("dat/SigPromptD0_0_10.txt","%lg %lg");
  TGraph *gr_new_D_B_Cen = new TGraph("dat/SigNonPromptD0_0_10.txt","%lg %lg");
  TGraph *gr_new_D_Peri = new TGraph("dat/SigPromptD0_60_80.txt","%lg %lg");
  TGraph *gr_new_D_B_Peri = new TGraph("dat/SigNonPromptD0_60_80.txt","%lg %lg");

}
