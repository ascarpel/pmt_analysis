////////////////////////////////////////////////////////////////////////////////
// Waveform methods implementation
//
// mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////////////

#include "Waveform.h"



Waveform::Waveform(){};


//------------------------------------------------------------------------------


Waveform::Waveform(int run, int subrun, int event, int board, int channel)
  : m_run(run)
  , m_subrun(subrun)
  , m_event(event)
  , m_board(board)
  , m_channel(channel)
{ };


//------------------------------------------------------------------------------


Waveform::Waveform( int run, int subrun, int event ,int board, int channel,
                                                      Rawdigits_t raw_waveform )
  : Waveform(run, subrun, event, board, channel)
{
  this->loadData(raw_waveform);
};


//------------------------------------------------------------------------------


Waveform::~Waveform(){};


//------------------------------------------------------------------------------


void Waveform::loadData( Rawdigits_t raw_waveform )
{
  m_nsamples=int(raw_waveform.size());
  m_raw_waveform = raw_waveform;
  for( auto w : raw_waveform ){ m_waveform.push_back( double(w) ); }
  removeBaseline();
};


//------------------------------------------------------------------------------


void Waveform::removeBaseline()
{
  // Calculate the baseline as the mean values on the first part of the spectrum
  n_sample_baseline = m_nsamples;

  for(int t=0; t<n_sample_baseline; t++)
  {
    baseline_mean += m_waveform.at(t);
  }
  baseline_mean /= n_sample_baseline;

  // Calculate the stdev of the baseline
  for(int t=0; t<n_sample_baseline; t++)
  {
    baseline_width += pow(m_waveform.at(t)-baseline_mean, 2);
  }
  baseline_width = sqrt( baseline_width / (n_sample_baseline-1) );

  // Subtract the baseline from the signal waveform
  std::transform(m_waveform.begin(), m_waveform.end(), m_waveform.begin(),
                                [ & ] (double x) { return x - baseline_mean; });

};


//------------------------------------------------------------------------------


bool Waveform::find(int run, int subrun, int event, int board, int channel )
{
  return (run==m_run && subrun==m_subrun && event==m_event
                                       && board==m_board && channel==m_channel);
};


//------------------------------------------------------------------------------


bool Waveform::hasSignal(double n_sigma)
{
  // Find if the waveform has n consecutive counts above a threshold expressed
  // in number of sigmas

  bool has_signal=false;
  int n_counts=5.0;
  int counts=0;

  for( double value : m_waveform )
  {
    if( abs(value) > n_sigma*baseline_width ){ counts++; }
    else{ counts=0; } // Reset the counts

    if( counts > n_counts ){ has_signal=true; break; }
  }

  return has_signal;
};


//------------------------------------------------------------------------------


bool Waveform::hasPulse( double n_sigma )
{
  // Define the Pulse region of the signal as the region with >n consecutive
  // counts above a threshold expressed in number of sigmas.

  bool has_pulse=false;
  int counts=0;
  int n_counts=5.0;

  int t_start=0, t_end=0;

  for(int t=0; t<m_nsamples; t++)
  {
    double value = m_waveform.at(t);

    if( abs(value) > n_sigma*baseline_width )
    {
      counts++;
      if( counts==1 ){ t_start=t-2; }
    }
    else
    {
      // if it is a pulse region, then this is the ending clause
      if( counts>n_counts && t_start>0 )
      {
        t_end=t+3;
        has_pulse=true;
        break; // Force to have only one pulse
      }
      else
      {
        // reset everything
        counts=0; t_start=0; t_end=0;
      }
    }
  }

  // If no pulse is found, we end the games here..
  if( !has_pulse ) { return has_pulse; }

  // Here we define the characteristics of the pulse
  double amp = 0., charge=0;
  for( int t=t_start; t<t_end; t++ )
  {
    if( m_waveform.at(t) < amp ){ amp = m_waveform.at(t); }
    charge += abs(m_waveform.at(t));
  }

  m_start_time = t_start;
  m_width = (t_end-t_start);
  m_amplitude = amp;
  m_integral = charge;

  return has_pulse;
};


//------------------------------------------------------------------------------


bool Waveform::isValidWaveform()
{
  bool isValid;

  if(m_raw_waveform.size()==0){
    isValid=false;
  }
  else {
    isValid=true;
  }

  return isValid;
}


//------------------------------------------------------------------------------


Waveform::Complex_t Waveform::doFFT(Waveform::Waveform_t m_time_domain)
{
    Eigen::FFT<double> fft;
    Waveform::Complex_t  m_frequency_domain;
    fft.fwd(m_frequency_domain, m_time_domain);
    return m_frequency_domain;
}


Waveform::Waveform_t Waveform::doIFFT(Waveform::Complex_t m_frequency_domain)
{
    Eigen::FFT<double> fft;
    Waveform::Waveform_t  m_time_domain;
    fft.inv(m_time_domain, m_frequency_domain);
    return m_time_domain;
}

//------------------------------------------------------------------------------


void Waveform::filterNoise(int window_size=200, bool reverse=false,
                                                            float threshold=100)
{
  //Noise filter algorithm. It will produce a new waveform object after noise
  //mitigation. Noise patterns to be mitigated are selected on a given window at
  //the beginning or at end of the waveform.
  //Arguments:
  //  window_size: set the window to produce the nosie model
  //  reverse: if true, the noise window is calculated at the end of the waveform

  // Select a noise window sufficiently large to be representative of the noise
  // patters of the original waveform
  Waveform::Waveform_t tmp_noise(window_size);
  copy(m_waveform.begin(), m_waveform.begin()+window_size, tmp_noise.begin());

  if( reverse ){
    copy(m_waveform.end()-window_size, m_waveform.end(), tmp_noise.begin());
  }

  // Get the noise spectra and get rid of some of the most nasty frequencies
  // which survive above a threshold
  vector<complex<double>> spec = doFFT(tmp_noise);
  vector<complex<double>> tmp_spec(spec.size());
  for(int i=0; i<int(spec.size()); i++ )
  {
    double pwr = sqrt( pow(spec[i].real(), 2) + pow(spec[i].imag(), 2) );
    if( pwr  > threshold ){
      tmp_spec[i] = spec[i];
    }
  }

  // Produce the filtered waveform in time domain
  Waveform::Waveform_t tmp_filter = doIFFT(tmp_spec);

  // Now we mirror the tmp_filter waveform to match the same length of the
  // original waveform. The mirroring is reasonable in this case since the
  // base noise we would like to remove is periodic.
  tmp_filter.resize(m_nsamples);

  int groups = ceil(float(m_nsamples)/float(window_size));
  for(int group=1; group<groups; group++ )
  {
    for( int i=0; i<window_size; i++ )
    {
      if( group*window_size+i == m_nsamples ){ break; }
      tmp_filter[group*window_size+i] = tmp_filter[i];
    }
  }

  // Finally subtract the noise filter from the original waveform
  for(int t=0; t<m_nsamples; t++)
  {
    m_waveform[t] -= tmp_filter[t];
  }
}


//------------------------------------------------------------------------------


TH1D* Waveform::getPowerSpectrum()
{
  Waveform::Complex_t m_spectrum = this->Waveform::doFFT(m_waveform);
  int m_fft_size = int(m_spectrum.size())/2;
  double max_sampling = m_sampling_freq/2;
  double freq_res = m_sampling_freq/m_spectrum.size();

  TH1D *h_power = new TH1D("", ";Frequency [MHz];Power",
                             m_spectrum.size(), 0, m_spectrum.size()*freq_res );

  for(int f=0; f<int(m_spectrum.size()); f++)
  {
    double ampls = pow(m_spectrum.at(f).real(), 2)
                                              + pow(m_spectrum.at(f).imag(), 2);
    h_power->Fill(f*freq_res, ampls );
  }

  return h_power;
}

//------------------------------------------------------------------------------

TH1D* Waveform::getWaveformHist()
{
  char hname[100];
  sprintf(hname, "Run%d-Subrun%d-Event%d-Board%d-Channel%d_hist", m_run,
                                         m_subrun, m_event, m_board, m_channel);

  TH1D *hist = new TH1D(hname, ";Time (ns);ADC", m_nsamples,
                                               0, m_nsamples*m_sampling_period);

  for(int t=0; t<m_nsamples; t++){ hist->Fill( t*m_sampling_period,
                                                           m_waveform.at(t) ); }

  return hist;
};

//------------------------------------------------------------------------------

TH1D* Waveform::getRawWaveformHist()
{
  char hname[100];
  sprintf(hname, "Run%d-Subrun%d-Event%d-Board%d-Channel%d_raw_hist", m_run,
                                         m_subrun, m_event, m_board, m_channel);

  TH1D *hist = new TH1D(hname, ";Time [ns];ADC", m_nsamples,
                                               0, m_nsamples*m_sampling_period);

  for(int t=0; t<m_nsamples; t++){ hist->Fill( t*m_sampling_period,
                                                       m_raw_waveform.at(t) ); }

  return hist;
};
