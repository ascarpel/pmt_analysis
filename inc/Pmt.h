////////////////////////////////////////////////////////////////////////////////
// PMT Class definition
//
// mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////////////

#ifndef  __PMT_H
#define __PMT_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

#include "TH1D.h"

using namespace std;

class PMT
{

  public:

    PMT();
    PMT( int channel );
    PMT( int channel, double hv);
    PMT( int channel, std::string name );
    ~PMT();

    void loadCharge( double charge ){ m_charge_array.push_back(charge); };
    void loadAmplitude( double amp ){ m_amplitude_array.push_back(amp); }

    //getters
    double getHV(){ return m_hv; }\
    vector<double> getArrayCharge(){ return m_charge_array; }
    TH1D* getHistCharge(){ return h_charge; }

    // Helpers
    void findExtremes( std::vector<double> vec, double &min, double &max );
    void initChargeHist( int nbins, double min, double max );
    void clean();


  private:

    std::string pmtname="hist";

    int m_channel;
    double m_hv;

    TH1D *h_charge;

    vector<double> m_charge_array;
    vector<double> m_amplitude_array;

    double m_charge_max;
    double m_charge_min;
    int m_nbins;

};

#endif //__PMT_H
