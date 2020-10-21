////////////////////////////////////////////////////////////////////////////////
// ChargeFit definition
//
// mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////////////

#ifndef  __CHARGEFIT_H
#define __CHARGEFIT_H

#include "TCanvas.h"
#include "TH1D.h"
#include "TF1.h"
#include "TMath.h"
#include "TMinuit.h"

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

class ChargeFit
{

  public:

    ChargeFit();
    ChargeFit( vector<TH1D*> hist_array, bool norm );
    ~ChargeFit();

    void loadHist( TH1D *hist ){ m_hist_array.push_back(hist); }

    void loadInitialConditions();
    void multiHistFit();

    double jointFit(double npe,
                    double Q1, double w1, double a1,
                    double Q2, double w2, double a2,
                    double Q3, double w3, double a3,
                    double Q4, double w4, double a4,
                    double pedm1, double pedw1, double pedamp1, double pedtau1,
                    double pedm2, double pedw2, double pedamp2, double pedtau2,
                    double pedm3, double pedw3, double pedamp3, double pedtau3,
                    double pedm4, double pedw4, double pedamp4, double pedtau4 );


    static void jointFitFunc(int &npar, double * deriv, double &f, double * par, int iflag);
    static double singleFitFunc(double* x, double * par);

    void getParameters( int i, double &par, double &err )
    {
      par = m_result[i];
      err = m_errors[i];
    };

    void GetChisquare( double &chi2, int &ndf )
    {
      chi2 = m_chi2;
      ndf = m_ndf;
    }

    int GetFitstatus(){ return m_fitstatus; }

    void getCanvas();

    void clean(){ m_hist_array.clear(); }

  private:

    vector<TH1D*> m_hist_array;
    bool m_normalize=false;

    static const int m_parameters = 29;

    double m_chi2=0;
    int m_ndf=0;
    int m_fitstatus;

    double m_vstart[m_parameters];
    double m_step[m_parameters];
    double m_minval[m_parameters];
    double m_maxval[m_parameters];
    string m_parname[m_parameters];

    double m_result[m_parameters];
    double m_errors[m_parameters];

};

#endif //__CHARGEFIT_H
