////////////////////////////////////////////////////////////////////////////////
// This macro will fit waveforms on a single pmt ( given by a board + channel )
// configuration over a number of events specified. the input file can be changed
// in the argument of the macro as well as the board number, channel number, The
// start event and the number of events to process.
//
// It works best with a file from november data. Note that pulses are not present
// on all boards and all channels, some were in dark.
// As exercise: you can add a loop over many channels and boards, and skip all
// waveforms without a pulse.
//
// The class expGauss contains the fit function
// You can add more functions to test in a similar fashon...
//
// The fit is performed in the function fitWaveform
//
// the output histograms are saved in the output rootfilename whose name can be
// changed in the arguments of the macro
//
// usage:
// $ cd macro/
// $ root -l loadLib.cc fit_waveform.cc

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

// Fit function for the pulse
class expGaus {
  public:
    // use constructor to customize the function object
    double operator() (double* x, double * par){
      double t0 = par[0];
      //double mu = par[1];
      double w = par[1];
      double c = par[2];
      double a = par[3];

      double t = x[0];

      return a*c/2.0*TMath::Exp(c*c*w*w/2.0)*TMath::Exp(-1.0*c*(t-t0))
                                   * TMath::Erfc( 1.0/1.414* (c*w-(t-t0)/w) );
    }
};




void fitWaveform( TH1D *m_wave, double t_min, double t_max  ) {

  // NB: t_min and t_max are given in time domain, not in bin

  // Get The position of the peak and the amplitude;
  double ampl=0; double t_peak = t_min;
  for( int bin=t_min/2.0; bin<t_max/2.0; bin++  ) {
    if( ampl < m_wave->GetBinContent(bin) ){
      ampl = m_wave->GetBinContent(bin);
      t_peak = m_wave->GetBinCenter(bin);
    }
  }

  double area = 2.0*m_wave->Integral( int((t_min)/2), int((t_max)/2) );

  char funcname[100]; sprintf(funcname, "func_pulse_direct");


  expGaus function_obj;
  TF1* func = new TF1(funcname, function_obj, t_min, t_max, 4);
  func->SetParNames("t0","#mu","#sigma","a");
  func->SetParameters(t_peak - 5, 2, 0.1, area );
  func->SetLineColor(kRed);
  func->SetLineStyle(2);
  func->SetLineWidth(2);

  // DO the fit: also do refit two more times if minimization didn't succeed.
  int status = m_wave->Fit(funcname,"NR+","",t_min, t_max);
  double par[5];
  for(int k=0; k<2; k++){
    for(int j=0; j<4; j++){
      par[j] = func->GetParameter(j);
      func->SetParameter(j,par[j]);
    }
    status = m_wave->Fit(funcname,"NR+","",t_min, t_max);
  }

}


//------------------------------------------------------------------------------


void fit_waveform(
  string filename="",
  string ofilename="",
  int startevent=0,
  int eventmax=10,
  int board = 0,
  int channel = 2)
{



  int nboards=12;
  int nchannels=16;

  //****************************************************************************
  // Input

  // Open TFile
  TFile* ifile = new TFile(filename.c_str(), "READ");
  cout << "Open TFile"+filename << endl;

  // Get the TTres and set the number of entries to process
  TTree* tevents = (TTree*)ifile->Get("caenv1730dump/events");
  std::vector<std::vector<uint16_t> > *data=0; //unsigned short
  tevents->SetBranchAddress("fWvfmsVec", &data);
  int nevents = tevents->GetEntries();
  if( eventmax<0 || eventmax>nevents ){ eventmax=nevents; }



  //Here we create an output where to store the result of the fit
  TFile* ofile = new TFile(ofilename.c_str(), "RECREATE");


  Waveform *myWaveAna = new Waveform();

  for( int event=startevent; event < eventmax; event++ ) {

    cout << "Get event: " << event << endl;
    tevents->GetEvent(event);

    const int n_channels = 16;
    const int n_samples = (*data)[0].size();
    const int n_boards = (*data).size()/n_channels;

    // Load data into the waveform analyzer ( remove automatically baseline )
    myWaveAna->loadData((*data).at(channel+n_channels*board));

    // Optional noise removal ( Uncomment if you want to use it )
    // myWaveAna->filterNoise( 200, false, 100 );


    // For the fit we need the binned distribution,
    // otherwise we cannot use the function TH1::Fit();
    TH1D * m_wave = myWaveAna->getWaveformHist();
    m_wave->SetName( ("Wave_event_"+to_string(event)).c_str() );


    fitWaveform( m_wave, 1900, 2000);



    m_wave->Write();


    // We prepare for the next waveform
    myWaveAna->clean();

  } //end event loop



} //end macro
