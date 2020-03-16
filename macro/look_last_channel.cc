#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"

void look_last_channel()
{

  int nboards=12;
  int nchannels=16;

  vector<double> rms;
  TFile *tfile = TFile::Open("noise_run1286_001.root");

  TCanvas *c = new TCanvas( "c", "", 500, 400 );

  for(int board=0; board<nboards; board++)
  {
    char hname[100]; sprintf(hname, "hist_board%d_channel15_noise", board);
    printf("%s",hname);
    TH1D *hist = (TH1D*)tfile->Get(hname);
    cout << board << " " << hist->GetRMS() << endl;
    hist->Draw("hist SAME PLC");
  }
}
