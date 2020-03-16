#include "TFile.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TMultiGraph.h"

#include "Utils.h"

void compare_noise()
{

  gStyle->SetPalette(kCMYK);

  TMultiGraph *mgr = new TMultiGraph();


   vector<string> filenames = {
    "noise_run1375_001.root",
    "noise_run1377_001.root",
    "noise_run1379_001.root",
    "noise_run1381_001.root",
  };

  for( string filename : filenames )
  {

    cout << "Processing file: " + filename << endl;

    int run=0.0; int subrun=0.0;
    utils::get_rundata( filename.substr(6, filename.size()), run, subrun );

    TFile *_file = new TFile(filename.c_str(), "READ");

    char gname[100]; sprintf(gname, "g_rms_run%d_subrun%d", run, subrun );
    TGraph *gRMSBoard = (TGraph*)_file->Get(gname);

    mgr->Add(gRMSBoard);

  }


  TCanvas *c = new TCanvas("c", "c", 700, 400);

  mgr->GetXaxis()->SetTitle("Channel number");
  mgr->GetXaxis()->CenterTitle();
  mgr->GetXaxis()->SetRangeUser(0, 15);
  mgr->GetYaxis()->SetTitle("Channel RMS");
  mgr->GetYaxis()->CenterTitle();

  mgr->Draw("APL");

  gPad->BuildLegend();

}
