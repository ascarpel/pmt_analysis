////////////////////////////////////////////////////////////////////////////////
// Macro to analyze the PMT noise
// Process the noise analyis for the filename given as argument. Edit the input
// filename with the file you prefer. Note: the output will always be in the format
// "noise_run%d_00%d.root". you must change manually the variables run and subrun for now.
//
// Usage: root -l loadLib.cc noise_ana.cc
// Out of the box, it will produce a root file in the macro folder called "noise_run1264_1.root"
// with the ADC distribution of all the pmts in data_dl1_run1264_1_20200227T235326-decoded.root
// after baseline subtraction.
//
// mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////////////

#include "Run.h"
#include "Waveform.h"
#include "Pmt.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"

using namespace std;

void noise_ana( string filename="../data/noise/run1377_001.root" )
{

  // CHANGE THE FILENAME BELOW IF WANT TO USE ONE OTHER FILE FOR YOUR ANALYSIS
  RUN my_run(filename);
  int run=my_run.getRun();
  int subrun=my_run.getSubrun();

  const int nboards=1;
  const int nchannels=16;

  //****************************************************************************
  // Input
  //

  TChain *tchain =  new TChain("caenv1730dump/events");

  // HERE you filename is loaded inside the macro to be processed!
  tchain->Add(filename.c_str());

  std::vector<std::vector<uint16_t> > *data=0; //unsigned short
  tchain->SetBranchAddress("fWvfmsVec", &data);


  TH1D *h_rms_before[nboards][nchannels];
  TH1D *h_rms_after[nboards][nchannels];

  for(int board=0; board<nboards; board++)
  {
    for(int channel=0; channel<nchannels; channel++)
    {
      char hname[100];
      sprintf(hname, "board%d_channel%d_rms_before", board, channel);
      h_rms_before[board][channel]= new TH1D(hname, hname, 40, -20, 20);

      sprintf(hname, "board%d_channel%d_rms_after", board, channel);
      h_rms_after[board][channel]= new TH1D(hname, hname, 40, -20, 20);
    }
  }


  for(int e=0; e<tchain->GetEntries(); e++)
  {
    cout << "Processing event: " << e << endl;

    // WE TAKE THE EVENT
    tchain->GetEvent(e);

    for(int board=0; board<nboards; board++)
    {
      for(int channel=0; channel<nchannels; channel++)
      {

        Waveform *waveform = new Waveform(run, subrun ,e, board, channel);

        waveform->loadData((*data).at(channel+nchannels*board));

        if(!waveform->isValidWaveform()){ continue; }
        for( float entry : waveform->getWaveform() )
        {
          h_rms_before[board][channel]->Fill( entry );
        }

        waveform->filterNoise(200, false, 100);
        for( float entry : waveform->getWaveform() )
        {
          h_rms_after[board][channel]->Fill( entry );
        }

      } // channel
    } // boards
  } // event



  //****************************************************************************
  // Output: create an output TFile and write the histogram in it
  //

  char ofilename[100]; sprintf(ofilename, "noise_run%d_%d.root", run, subrun);
  TFile ofile(ofilename, "RECREATE"); ofile.cd();

  TGraph *g_rms_before = new TGraph( nboards*nchannels );
  g_rms_before->SetTitle("RMS Before filter");
  g_rms_before->SetName("RMS_Before_filter");

  TGraph *g_rms_after = new TGraph( nboards*nchannels );
  g_rms_after->SetTitle("RMS Before after");
  g_rms_after->SetName("RMS_Before_after");

  for(int board=0; board<nboards; board++)
  {
    for(int channel=0; channel<nchannels; channel++)
    {
      h_rms_before[board][channel]->Write();
      h_rms_after[board][channel]->Write();

      double rms_before = h_rms_before[board][channel]->GetRMS();
      double rms_after = h_rms_after[board][channel]->GetRMS();

      int number=channel+board*nchannels;
      g_rms_before->SetPoint( number, (double)number, rms_before);
      g_rms_after->SetPoint( number, (double)number, rms_after);
    }
  }

  ofile.Close();

  TCanvas *g = new TCanvas("g", "g", 500, 400);
  TMultiGraph *mgr = new TMultiGraph();
  mgr->Add(g_rms_before); mgr->Add(g_rms_after);

  mgr->Draw("APL");

  cout << "All done" << endl;

} //end macro
