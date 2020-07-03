////////////////////////////////////////////////////////////////////////////////
// Waveform methods implementation
//
// mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////////////

#include "Waveform.h"



Waveform::Waveform(){};

//------------------------------------------------------------------------------


Waveform::~Waveform() { this->clean(); };


//------------------------------------------------------------------------------


void Waveform::loadData( Rawdigits_t raw_waveform )
{

  m_nsamples= raw_waveform.size();
  m_raw_waveform = raw_waveform;

  for( auto w : raw_waveform ) {

    double value = -1*double(w); // Reverse polarity
    m_waveform.push_back( value );

  }

  removeBaseline();

}


//------------------------------------------------------------------------------


void Waveform::calculateWaveformMeanAndRMS( double &mean, double &width  )
{

  // Calculate the baseline as the mean values

  for(int t=0; t<n_sample_baseline; t++){
    mean += m_waveform.at(t);
  }
  mean /= n_sample_baseline;

  // Calculate the stdev of the baseline
  for(int t=0; t<n_sample_baseline; t++){
    width += pow(m_waveform.at(t)-mean, 2);
  }
  width = sqrt( width / (n_sample_baseline-1) );

  return;

}


//------------------------------------------------------------------------------


void Waveform::removeBaseline()
{

  m_baseline_mean=0; m_baseline_width=0;
  this->calculateWaveformMeanAndRMS(m_baseline_mean, m_baseline_width);

  // Subtract the baseline from the signal waveform
  std::transform(m_waveform.begin(), m_waveform.end(), m_waveform.begin(),
                            [ & ] (double x) { return x - m_baseline_mean; });

  }


//------------------------------------------------------------------------------




//------------------------------------------------------------------------


void Waveform::resetPulse(Waveform::Pulse &pulse)
{
  pulse.start_time = 0;
  pulse.end_time=0;
  pulse.time_peak = 0;
  pulse.width = 0;
  pulse.amplitude = 0;
  pulse.integral = 0;
  pulse.fit_start_time = 0;
  pulse.error_start_time = 0;
  pulse.fit_sigma = 0;
  pulse.error_sigma = 0;
  pulse.fit_mu = 0;
  pulse.error_mu = 0;
  pulse.fit_amplitude = 0;
  pulse.error_amplitude = 0;
  pulse.chi2 = 0;
  pulse.ndf = 0;
  pulse.fitstatus = 999;
}



//------------------------------------------------------------------------



std::vector<Waveform::Pulse> Waveform::findPulses()
{
  // Simple thresold pulse finder
  std::vector<Waveform::Pulse> pulse_v;

  auto start_threshold = m_start_adc_thres;
  auto end_threshold   = m_end_adc_thres;

  bool fire = false;
  double counter = 0;

  Waveform::Pulse pulse;

  for( auto const &value : m_waveform ){

    // Start logic
    if( !fire && value >= start_threshold ){
      // Found a new pulse
      fire = true;
      pulse.start_time = counter - 1 > 0 ? counter - 1 : counter;
    }

    // End logic
    if( fire && value < end_threshold ){
        fire = false;
        pulse.end_time = counter < m_waveform.size() ? counter : counter - 1;
        pulse.width = pulse.end_time-pulse.start_time;
        pulse_v.push_back(pulse);
        this->resetPulse(pulse);
    }

    // Find max and integral
    if(fire){
      pulse.integral += value;
      if( pulse.amplitude < value ){
        pulse.amplitude = value;
        pulse.time_peak = counter;
      }
    }

    counter++;

  } //end waveform

  if(fire){
    // Take care of a pulse that did not finish within the readout window.
    fire = false;
    pulse.end_time = counter - 1;
    pulse.width = pulse.end_time-pulse.start_time;
    pulse_v.push_back(pulse);
    this->resetPulse(pulse);
  }

  return pulse_v;

} // end findPulse()



//--------------------------------------------------------------------------


Waveform::Pulse Waveform::getLaserPulse() {


//Signal integral over a fixed time window: used for direct light
//calibration if do fit is true, a fit is performed to find the t0

Waveform::Pulse temp_pulse;

double charge=0.0;
for( int i = m_startbin; i<m_startbin+m_nbins; i++ ) {
  charge += m_waveform.at(i);
}

std::vector<Waveform::Pulse> pulses;
for( auto & pulse : this->findPulses()) {
  if( (pulse.time_peak > m_trigger_time[0]) && (pulse.time_peak < m_trigger_time[1]) ) {
    pulses.push_back( pulse );
  }
}


// If pulses are not found, just return the integration window
if( pulses.size() == 0 ) {

  temp_pulse.start_time = 0;
  temp_pulse.time_peak = 0;
  temp_pulse.width = 0;
  temp_pulse.amplitude = 0;
  temp_pulse.integral = charge * adc_to_pC;

  return temp_pulse;
}

// Sort pulses cronologically and choose the first
Waveform::Pulse laserPulse;

if( pulses.size() > 1 ){
  std::sort( pulses.begin(), pulses.end(),
        []( Waveform::Pulse & a, Waveform::Pulse & b ) -> bool {
                    return a.start_time > b.start_time;
        } );
}

temp_pulse = pulses.at(0);

// Now we add the fit information
if( m_dofit ) {

  auto m_wave = this->getWaveformHist();

  double t_min = temp_pulse.time_peak - m_fitrange[0];
  double t_max = temp_pulse.time_peak + m_fitrange[1];

  char funcname[100]; sprintf(funcname, "func_pulse_direct");

  PulseShapeFunction_ExpGaus function_obj;
  TF1* func = new TF1(funcname, function_obj, t_min, t_max, 4);
  func->SetParNames("t0","#mu","#sigma","a");
  func->SetParameters(temp_pulse.time_peak - 5, 2, 0.1,
                2.0*m_wave->Integral( int((t_min)/2), int((t_max)/2) ));
  int status = m_wave->Fit(funcname,"RQN","",t_min, t_max);

  // Refit two more times with the latest version of the parameters
  double par[5];
  for(int k=0; k<2; k++){
    for(int j=0; j<4; j++){
      par[j] = func->GetParameter(j);
      func->SetParameter(j,par[j]);
    }
    status = m_wave->Fit(funcname,"RQN","",t_min, t_max);
  }

  //If the fit status is ok we calculated the rising time as the time
  // when the fitted function has value 10% of its max

  double first_spe_time = -999;
  double dt=-999;

  if( status ==0 ) {

    //No point in doing that if fit is garbage

    int npoints = 5000; // should grant decent resolution
    dt = (m_trigger_time[1]-m_trigger_time[0])/npoints;
    double max=0.0;
    for(int i=0; i<npoints; i++){
      double t=m_trigger_time[0] + dt*i;
      if(max < func->Eval(t) ){ max = func->Eval(t); };
    }

    double startval = 0.1 * max;
    for(int i=0; i<npoints; i++){
      double t=m_trigger_time[0] + dt*i;
      if(startval < func->Eval(t) ){ first_spe_time = t; break; };
    }
  }

  // Save the fit paramteters to the pulse object
  temp_pulse.fit_start_time = first_spe_time;
  temp_pulse.error_start_time = dt; // TODO: if needed should use the par errors
  temp_pulse.fit_mu = func->GetParameter(1);
  temp_pulse.error_mu = func->GetParError(1);
  temp_pulse.fit_sigma = func->GetParameter(2);
  temp_pulse.error_sigma = func->GetParError(2);
  temp_pulse.fit_amplitude = func->GetParameter(3);
  temp_pulse.error_amplitude = func->GetParError(3);
  temp_pulse.chi2 = func->GetChisquare();
  temp_pulse.ndf = func->GetNDF();
  temp_pulse.fitstatus = status;

}

  return temp_pulse;


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


void Waveform::filterNoise(size_t window_size=200, bool reverse=false, float threshold=100)
{
  //Noise filter algorithm. It will produce a new waveform object after noise
  //mitigation. Noise patterns to be mitigated are selected on a given window at
  //the beginning or at end of the waveform.
  //Arguments:
  //  window_size: set the window to produce the nosie model
  //  reverse: if true, the noise window is calculated at the end of the waveform


  Waveform::Waveform_t tmp_noise(window_size);
  copy(m_waveform.begin(), m_waveform.begin()+window_size, tmp_noise.begin());

  if( reverse ){
    copy(m_waveform.end()-window_size, m_waveform.end(), tmp_noise.begin());
  }

  // Get the noise spectra and get rid of some of the most nasty frequencies
  // which survive above a threshold
  vector<complex<double>> spec = doFFT(tmp_noise);
  vector<complex<double>> tmp_spec(spec.size());
  for(size_t i=0; i<spec.size(); i++ )
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

  size_t groups = ceil(float(m_nsamples)/float(window_size));
  for(size_t group=0; group<groups+1; group++ )
  {
    for( size_t i=0; i<window_size; i++ ){
      if( int(group*window_size+i) < int(m_nsamples) ){

        m_waveform[group*window_size+i] -= tmp_filter[i];

      }
      else{ break; }
    }
  }
};


//------------------------------------------------------------------------------


TH1D* Waveform::getPowerSpectrum()
{
  Waveform::Complex_t m_spectrum = this->Waveform::doFFT(m_waveform);
  int m_fft_size = int(m_spectrum.size())/2;
  double max_sampling = m_sampling_freq/2;
  double freq_res = m_sampling_freq/m_spectrum.size();

  TH1D *h_power = new TH1D("freq_domain", ";Frequency [MHz];Power",
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
  sprintf(hname, "binned_wave");

  TH1D *hist = new TH1D(hname, ";Time (ns);ADC", m_nsamples,
                                               0, m_nsamples*m_sampling_period);

  for(int t=0; t<m_nsamples; t++) {
    hist->Fill( t*m_sampling_period, m_waveform.at(t) );
    hist->SetBinError(t, 1.0);
  }

  return hist;
};


//------------------------------------------------------------------------------


TH1D* Waveform::getRawWaveformHist()
{
  char hname[100];
  sprintf(hname, "binned_wave_raw");

  TH1D *hist = new TH1D(hname, ";Time [ns];ADC", m_nsamples,
                                               0, m_nsamples*m_sampling_period);

  for(int t=0; t<m_nsamples; t++) {
    hist->Fill( t*m_sampling_period, m_raw_waveform.at(t) );
    hist->SetBinError(t, 1.0);
  }

  return hist;
};


//------------------------------------------------------------------------------


void Waveform::clean() {

  m_waveform.clear();
  m_raw_waveform.clear();

  m_baseline_mean=0;
  m_baseline_width=0;

};
