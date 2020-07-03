////////////////////////////////////////////////////////////////////////////////
// PMT Methods implementation
//
// mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////////////

#include "Pmt.h"

PMT::PMT( )
{ };

//------------------------------------------------------------------------------

PMT::PMT( int channel )
  : m_channel(channel)
  , m_hv(0)
{ };

//------------------------------------------------------------------------------

PMT::PMT( int channel, double hv )
  : m_channel(channel)
  , m_hv(hv)
{ };

//------------------------------------------------------------------------------

PMT::PMT( int channel, string name )
  : m_channel(channel)
  , pmtname(name)
{ };

//------------------------------------------------------------------------------

PMT::~PMT()
{
  clean();
};

//------------------------------------------------------------------------------


void PMT::findExtremes( std::vector<double> vec, double &min, double &max )
{

  max = 0;
  min =9999;

  for( double v : vec )
  {
    if( v > max  ){ max = v;  }
    if( v < min ){ min = v; }
  }
};


void PMT::initChargeHist( int nbins = -999, double min=-999, double max=-999 )
{

  // Charge hist: dynamical binning if I don't put anything
  findExtremes(m_charge_array, m_charge_min, m_charge_max);
  if (min != -999 ){ m_charge_min = min; }
  if( max !=-999 ){ m_charge_max = max; }
  if( nbins == -999 ){ m_nbins = sqrt( m_charge_array.size() ); }
  else{ m_nbins = nbins; }

  h_charge = new TH1D(pmtname.c_str(), pmtname.c_str(), m_nbins, m_charge_min, m_charge_max);

  for(auto entry : m_charge_array){ h_charge->Fill(entry); }

};


//------------------------------------------------------------------------------

void PMT::clean()
{
  m_charge_array.clear();
  m_amplitude_array.clear();
};

//------------------------------------------------------------------------------
