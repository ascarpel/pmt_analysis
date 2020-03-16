#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"

void compare_rms()
{

  //gStyle->SetPalette(kRainBow);

  TFile *tfile = TFile::Open("noise_run1381_001.root");
  TCanvas *c = new TCanvas( "c", "", 500, 400 );

  TH1D *hist_any = (TH1D*)tfile->Get("hist_board0_channel1_noise");
  hist_any->SetName("Channel 1");
  hist_any->SetTitle("Channel 1");
  hist_any->SetLineColor(kRed+1);
  hist_any->Draw("hist");

  TH1D *hist_last = (TH1D*)tfile->Get("hist_board0_channel15_noise");
  hist_last->SetName("Channel 15");
  hist_last->SetTitle("Channel 15");
  hist_last->SetLineColor(kBlue+1);
  hist_last->Draw("hist SAME");

  TLegend *l = new TLegend(0.162651, 0.744, 0.461847, 0.896);
  //l->SetNColumns(2);
  l->AddEntry(hist_any, "Channel 1", "l");
  l->AddEntry(hist_last, "Channel 15", "l");
  l->Draw();

  //c.Update();
}

void look_last_channel()
{
  compare_rms();


}
