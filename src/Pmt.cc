////////////////////////////////////////////////////////////////////////////////
// PMT Methods implementation
//
// mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////////////

#include "Pmt.h"

PMT::PMT( )
{ };

//------------------------------------------------------------------------------

PMT::PMT( int run, int board, int channel )
  : m_run(run)
  , m_board(board)
  , m_channel(channel)
{
  this->initHist();
};

//------------------------------------------------------------------------------

PMT::~PMT()
{};

//------------------------------------------------------------------------------

void PMT::initHist()
{
  char hname[100];

  // Noise
  sprintf(hname, "hist_board%d_channel%d_noise", m_board,  m_channel);
  h_noise = new TH1D(hname, hname, 40, -20, 20);

  // Amplitude
  sprintf(hname, "hist_board%d_channel%d_amplitude", m_board,  m_channel);
  h_amplitude = new TH1D(hname, hname, 100, 0, 2000);

  sprintf(hname, "hist_board%d_channel%d_amplitude_low", m_board,  m_channel);
  h_amplitude_low = new TH1D(hname, hname, 20, 0, 150);

};

//------------------------------------------------------------------------------

void PMT::loadWaveform( Waveform *waveform )
{

  // Fill nosie (just consider the last 200 bins)
  Waveform::Waveform_t wave = waveform->getWaveform();
  for(size_t t=0; t<wave.size(); t++)
  {
    h_noise->Fill( wave.at(t) );
  }

  // Fill the amplitude histograms
  double amplitude = abs( waveform->getAmplitude() );
  h_amplitude->Fill(amplitude);
  h_amplitude_low->Fill(amplitude);

  // Fill the arrays
  // m_amplitude_array.push_back(amplitude);

};

//------------------------------------------------------------------------------

void PMT::writeHist()
{

  // Write the histograms to file
  h_noise->Write();

  h_amplitude->Write();
  h_amplitude_low->Write();

};

//------------------------------------------------------------------------------

void PMT::clean()
{
  // m_amplitude_array.clear();
};

//------------------------------------------------------------------------------
