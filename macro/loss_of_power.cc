////////////////////////////////////////////////////////////////////////////////
// Use Milind's slides to understand the Eigen FFT
// mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////////////

#include "Waveform.h"
#include "Pmt.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"

#include <Eigen/Core>
#include <unsupported/Eigen/FFT>

//------------------------------------------------------------------------------

void makehist(vector<float> x, vector<float> y, string cname, string histname)
{
  TCanvas *c = new TCanvas(cname.c_str(), cname.c_str(), 500, 400);

  TH1F *hist = new TH1F(histname.c_str(), histname.c_str(), x.size(), 0, x.at(y.size()-1));
  for(int i=0; i<x.size(); i++){ hist->Fill(x.at(i), y.at(i) ); }

  hist->Draw("hist");
}

//------------------------------------------------------------------------------

vector<complex<float>> dft(vector<float> time)
{
    Eigen::FFT<float> gEigenFFT;
    vector<complex<float>>  spec;
    gEigenFFT.fwd(spec, time);
    return spec;
}

vector<float> idft(vector<complex<float>> spec)
{
    Eigen::FFT<float> gEigenFFT;
    auto v = Eigen::Map<Eigen::VectorXcf>(spec.data(), spec.size());
    Eigen::VectorXf ret;
    gEigenFFT.inv(ret, v);
    return vector<float>(ret.data(), ret.data()+ret.size());
}

//------------------------------------------------------------------------------

float filter( float x, float param=4 )
{
  double factor=1;
  for(float i=param;i>0; i-- ){ factor *=i; }
  float y = exp(-x)*pow( x, param )/factor;

  return y;
}

//------------------------------------------------------------------------------

void loss_of_power()
{
  // init constants
  const int m_samples = 257;
  const float m_sampling_time = 0.1; // in sec
  const float m_sampling_freq = 1./m_sampling_time; //in Hertz
  const float m_freq_res = m_sampling_freq/m_samples; // in Hertz
  const float m_norm=sqrt(2*TMath::Pi()/m_samples)*(1./m_sampling_time);

  // Calculate the shaping response in time domain
  vector<float> time; time.resize(m_samples);
  vector<float> resp; resp.resize(m_samples);
  for( int i=0; i<m_samples; i++ )
  {
    time[i] = i*m_sampling_time;
    resp[i] = filter(time[i]);
  }
  makehist(time, resp, "c", "base_shaping");

  // Now perform the FFT using eigen
  vector<complex<float>> spec = dft(resp);
  vector<complex<float>> norm_spec; norm_spec.resize(m_samples);
  vector<float> freq; freq.resize(m_samples);
  vector<float> pwrs; pwrs.resize(m_samples);
  vector<float> norm_pwrs; norm_pwrs.resize(m_samples);
  vector<float> phse; phse.resize(m_samples);

  for( int i=0; i<m_samples; i++ )
  {
    freq[i] = i*m_freq_res;
    norm_spec[i] = m_norm*spec[i];
    pwrs[i] = log10( sqrt(pow(spec[i].real(),2) + pow(spec[i].imag(), 2)) );
    norm_pwrs[i] = log10( sqrt(pow(norm_spec[i].real(),2) + pow(norm_spec[i].imag(), 2)) );
    phse[i] = spec[i].imag();
  }
  makehist(freq, pwrs, "c1", "power_shaping");
  makehist(freq, norm_pwrs, "c2", "norm_power_shaping");

  //re-trasform this in time domain
  vector<float> no_norm_resp = idft(spec);
  vector<float> norm_resp = idft(norm_spec);
  vector<float> diff; diff.resize(m_samples);
  for( int i=0; i<m_samples; i++ )
  {
    diff[i] = norm_resp[i]-no_norm_resp[i];
  }
  makehist(time, no_norm_resp, "c3", "base_shaping_return");
  makehist(time, norm_resp, "c4", "norm_base_shaping_return");
  makehist(time, diff, "c5", "time_shaping_diff");

  // Post-note: eigen claims scailing is performed: https://eigen.tuxfamily.org/index.php?title=EigenFFT

}
