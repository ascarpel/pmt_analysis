////////////////////////////////////////////////////////////////////////////////
// Macro that interfaces the output of the CAEN decoder script and produces
// a file format compatible with this small analysis framework
//
// mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////////////

#include "Waveform.h"
#include "Pmt.h"
#include "Utils.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"

#include <Eigen/Core>
#include <unsupported/Eigen/FFT>

using namespace std;



void plot_waveform(
  string filename="../data/noise/run1398_001.root",
  int event=0,
  int board = 0,
  int channel = 0)
{

  //****************************************************************************
  // General definitions ( ..import from DB )

  RUN my_run(filename);
  int run=my_run.getRun();
  int subrun=my_run.getSubrun();
  
  int nboards=12;
  int nchannels=16;

  //****************************************************************************
  // Input

  // Open TFile
  TFile* ifile = new TFile(filename.c_str(), "READ");
  cout << "Open TFile"+filename << endl;

  // Get the TTres
  TTree* tevents = (TTree*)ifile->Get("caenv1730dump/events");
  int nevents = tevents->GetEntries();

  // Set Branch address
  std::vector<std::vector<uint16_t> > *data=0; //unsigned short
  tevents->SetBranchAddress("fWvfmsVec", &data);

  //****************************************************************************
  // Event loop

  cout << "Get event: " << event << endl;
  tevents->GetEvent(event);

  const int n_channels = 16;
  const int n_samples = (*data)[0].size();
  const int n_boards = (*data).size()/n_channels;

  // Create the PMT object
  Waveform *waveform = new Waveform(run, subrun ,event, board, channel);
  waveform->loadData((*data).at(channel+n_channels*board));

  char cname[100];

  // Waveform - baseline
  sprintf(cname, "Run%d-Subrun%d-Event%d-Board%d-Channel%d_time",
                                            run, subrun ,event, board, channel);
  printf("%s \n", cname);
  TH1D* hist = waveform->getWaveformHist();
  TCanvas *c = new TCanvas(cname, cname, 1200, 400); c->cd();
  hist->GetXaxis()->SetRangeUser(4000, 5000);
  hist->GetXaxis()->CenterTitle();
  hist->GetYaxis()->SetTitleOffset(0.33);
  hist->GetYaxis()->CenterTitle();
  hist->GetYaxis()->SetRangeUser(-20, 20);
  hist->Draw("hist");

  double baselineMean  = waveform->getBaselineMean();
  double baselineWidth = waveform->getBaselineWidth();
  cout << baselineMean << " " << baselineWidth << endl;

  // Raw waveform
  sprintf(cname, "Run%d-Subrun%d-Event%d-Board%d-Channel%d_raw",
                                            run, subrun ,event, board, channel);
  printf("%s \n", cname);
  TH1D* hraw = waveform->getRawWaveformHist();
  TCanvas *c2 = new TCanvas(cname, cname, 1200, 400); c2->cd();
  hraw->GetXaxis()->SetRangeUser(4000, 5000);
  hraw->GetXaxis()->CenterTitle();
  hraw->GetYaxis()->SetTitleOffset(0.33);
  hraw->GetYaxis()->CenterTitle();
  hraw->GetYaxis()->SetRangeUser(baselineMean-3*baselineWidth,
                                                  baselineMean+3*baselineWidth);
  hraw->Draw("hist");

  // Power spectrum
  TH1D* hfreq = waveform->getPowerSpectrum();
  sprintf(cname, "Run%d-Subrun%d-Event%d-Board%d-Channel%d_freq",
                                            run, subrun ,event, board, channel);
  TCanvas *c1 = new TCanvas(cname, cname, 1200, 400); c1->cd();
  hfreq->GetXaxis()->SetRangeUser(0, 80);
  hfreq->GetXaxis()->CenterTitle();
  hfreq->GetYaxis()->SetTitleOffset(0.52);
  hfreq->GetYaxis()->CenterTitle();
  hfreq->Draw("hist");

} //end macro
