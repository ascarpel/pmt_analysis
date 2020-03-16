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

#include "Waveform.h"
#include "TH1D.h"

using namespace std;

class PMT
{

  public:

    PMT();
    PMT(int run, int board, int channel );
    ~PMT();

    // Loaders
    void loadWaveform( Waveform *waveform );

    // Getters
    int getRun(){ return m_run; };
    int getBoard(){ return m_board; };
    int getChannel(){ return m_channel; }

    TH1D* getNoiseHist(){ return h_noise; };

    // Helpers
    void clean();
    void initHist();
    void writeHist();


  private:

    int m_run;
    int m_board;
    int m_channel;

    TH1D *h_amplitude;
    TH1D *h_amplitude_low;
    TH1D *h_noise;


    // vector<double> m_amplitude_array;

};

#endif //__PMT_H
