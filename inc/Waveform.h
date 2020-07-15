////////////////////////////////////////////////////////////////////////////////
// Waveform class definition
//
// mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////////////

#ifndef  __WAVEFORM_H
#define __WAVEFORM_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <Eigen/Core>
#include <unsupported/Eigen/FFT>
#include <numeric>
#include <complex>

#include "TH1D.h"
#include "TF1.h"
#include "TMath.h"

using namespace std;

// Fit function for the pulse
class PulseShapeFunction_ExpGaus {
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


class Waveform
{

    public:

      struct Pulse
      {
        // Pulse characteristics
        double start_time = 0;
        double end_time=0;
        double time_peak = 0;
        double width = 0;
        double amplitude = 0;
        double integral = 0;

        // Fit paramteters
        double fit_start_time = 0;
        double error_start_time = 0;
        double fit_sigma = 0;
        double error_sigma = 0;
        double fit_mu = 0;
        double error_mu = 0;
        double fit_amplitude = 0;
        double error_amplitude = 0;
        double chi2 = 0;
        double ndf = 0;
        double fitstatus = 999;
      };

      typedef vector<unsigned short> Rawdigits_t;
      typedef vector<double> Waveform_t;
      typedef vector<complex<double>> Complex_t;

      Waveform();
      ~Waveform();

      // Import data
      void loadData( Rawdigits_t raw_waveform );


      Rawdigits_t getRawWaveform(){ return m_raw_waveform; };
      Waveform_t getWaveform(){ return m_waveform; };


      // Baseline removal
      void calculateWaveformMeanAndRMS( double &mean, double &width  );
      void removeBaseline();
      double getBaselineMean(){return m_baseline_mean;};
      double getBaselineWidth(){return m_baseline_width;}

      // Noise filter
      Complex_t doFFT(Waveform_t m_time_domain);
      Waveform_t doIFFT(Waveform::Complex_t m_frequency_domain);
      void filterNoise(size_t window_size, bool reverse,float threshold);

      // Pulse
      Pulse getLaserPulse();
      Pulse getIntegral();
      std::vector<Pulse> findPulses( double startTh, double endTh );
      void resetPulse(Pulse &pulse);


      // Histograms
      TH1D* getRawWaveformHist();
      TH1D* getWaveformHist();
      TH1D* getPowerSpectrum();

      void clean();

    private:


      // NB! PROVIDE EXTERNAL CONFIGURATION FOR SOME PARAMTERS << TO BE DONE

      int m_nsamples;
      double m_sampling_period=2.; //in ns
      double m_sampling_freq = (1./m_sampling_period)*1e3; //in MHz

      int n_sample_baseline=200;
      double m_baseline_mean=0.0;
      double m_baseline_width=0.0;

      Rawdigits_t m_raw_waveform;
      Waveform_t m_waveform;

      double m_start_time=0.0;
      double m_width=0.0;
      double m_amplitude=0.0;
      double m_integral=0.0;

      bool m_dofit;
      std::vector<double> m_fitrange;

      int m_startbin; // start of the integration window for the charge
      int m_nbins; // integration window for the charge
      std::vector<double> m_trigger_time; //trigger window
      double adc_to_pC = 0.122*2.0*0.02; // TODO: make it configurable

      double m_start_adc_thres; //start Threshold for a pulse
      double m_end_adc_thres; // end threshold for a pulse


};


#endif //__WAVEFORM_H
