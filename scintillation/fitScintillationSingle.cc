////////////////////////////////////////////////////////////////////////////////
// Macro that interfaces the output of the CAEN decoder script and produces
// a file format compatible with this small analysis framework
//
// mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////////////


#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TMinuit.h"
#include "TMath.h"
#include "TString.h"
#include "TF1.h"
#include "TCanvas.h"


class Scintillation {

  public:
    Scintillation();
    ~Scintillation();

    void fitResponse( TProfile *prof );

  private:

    TProfile *m_profile;
    static const int m_pars = 8;
    int m_ndf;
    double m_chi2;

    double singleComponentResponse(double x, double tm, double s, double tau, double a );

    static double singleComponentFcn( double* x, double * par );
    static double allComponentFcn( double* x, double * par );

    double larResponse(  double tm, double s,
                         double tauf, double taui, double taus,
                         double af,   double ai,  double as     );

    static void larResponseFcn(int &npar, double * deriv, double &f, double * par, int iflag);

};

static Scintillation* fitobj;


Scintillation::Scintillation() {
  fitobj=this;
};

Scintillation::~Scintillation(){ };


double Scintillation::singleComponentResponse(double x, double tm, double s, double tau, double a ) {

  double expPart = exp( (1./2.) * pow( (s / tau), 2 ) - ( x - tm ) / tau );

  double erfPart = 1-TMath::Erf( sqrt(1./2.) * ( (s / tau) - ( x - tm ) / s ) );

  double expNorm = ( 1./ 2. ) * ( a / tau ) * expPart * erfPart;

  return expNorm;

};


double Scintillation::singleComponentFcn( double* x, double *par ) {

  return fitobj->singleComponentResponse(x[0], par[0], par[1], par[2], par[3]);

};

double Scintillation::allComponentFcn( double* x, double * par ) {

  double slowc = fitobj->singleComponentResponse(x[0], par[0], par[1], par[2], par[5]);
  double intmc = fitobj->singleComponentResponse(x[0], par[0], par[1], par[3], par[6]);
  double fastc = fitobj->singleComponentResponse(x[0], par[0], par[1], par[4], par[7]);

  return slowc + intmc + fastc;

};


double Scintillation::larResponse(  double tm, double s,
                                  double tauf, double taui, double taus,
                                  double af,   double ai,  double as     ) {

  // Minimization of the fit response using the three components
  double tau[3]; tau[0]=tauf; tau[1]=taui; tau[2]=taus;
  double a[3]; a[0]=af; a[1]=ai; tau[2]=as;

  double chi2=0;
  int ndf = 0;

  //double baseline = this->m_profile->Integral(0, 50) / 50.;

  int nbins = this->m_profile->GetNbinsX();
  for( int bin=0; bin<nbins; bin++ ) {

      double x = this->m_profile->GetBinCenter(bin);
      double yo = this->m_profile->GetBinContent(bin);
      double ey = this->m_profile->GetBinError(bin);

      if(yo < 1e-3){ continue; }

      if(x < -20 || x > 1500){ continue; }

      //Loop over the components
      double ye=0.0;
      for( int i=0; i<3; i++ ){

        double expNorm = this->singleComponentResponse(x, tm, s, tau[i], a[i] );
        ye += expNorm;
      }

      if( ey > 0 ) { chi2 += pow( (yo-ye), 2 ) / ( ey*ey );   ndf++; }

  } // end loop over bins

  fitobj->m_ndf = ndf;
  fitobj->m_chi2 = chi2;
  return chi2;

};


void Scintillation::larResponseFcn(int &npar, double * deriv, double &f, double * par, int iflag){
  f = fitobj->larResponse(par[0], par[1], par[2], par[3], par[4], par[5], par[6], par[7] );
};





void Scintillation::fitResponse( TProfile *prof ) {


  // Pass the argument hist to the class member;
  m_profile = prof;

  /*

  double vstart[m_pars]  = { 0.0, 5.0, 6.0, 50, 1400, 3.34, 15, 120  };
  double vminval[m_pars] = { -5, 0.0, 5.9, 0.0, 1200, 0.0, 0.0, 0.0 };
  double vmaxval[m_pars] = {  5, 15.0, 6.1, 1000, 2000, 100, 100, 350 };
  double vsteps[m_pars] = { 0.1, 0.1, 0.1, 5, 10, 30, 10, 10 };

  string vnames[m_pars]={"#t_m","#sigma",
                          "#tau_{f}", "#tau_{i}", "#tau_{s}",
                          "#a_{f}", "#a_{i}", "#a_{s}" };

 double vresult[m_pars];
 double verrors[m_pars];

 TMinuit* minimizer = new TMinuit(m_pars);
 minimizer->SetFCN(Scintillation::larResponseFcn);

 for(int i=0; i<m_pars; i++){
   minimizer->DefineParameter(i, vnames[i].c_str(), vstart[i], vsteps[i],
                                                        vminval[i], vmaxval[i]);
 }

 // do the fit
 minimizer->Command("SET PRINT 2");//
 minimizer->mnscan();

 Int_t ierrs=0;
 minimizer->mnexcm("SET NOWarnings", 0, 0, ierrs);
 int fitstatus= minimizer->Migrad();

 // get parameters
 for(int i=0; i<m_pars; i++){
   minimizer->GetParameter(i, vresult[i], verrors[i]);

   std::cout << vnames[i] << " " << vresult[i] << " " << verrors[i] << std::endl;

 }

 std::cout << fitstatus << " " << m_chi2 << " " <<  m_ndf-m_pars  << std::endl;
 */


 int nbins = m_profile->GetNbinsX();
 double tlow = m_profile->GetBinLowEdge(0);
 double thigh = m_profile->GetBinLowEdge(nbins) + m_profile->GetBinWidth(nbins);


 double vthispar[m_pars]  = { 0.0, 5.0, 6.0, 50, 1400, 3.34, 3.4, 23 };
 TF1 *fitf = new TF1("allComponents", Scintillation::allComponentFcn, -50, 2400, 8);
 fitf->SetParameters( vthispar );
 fitf->SetLineColor(kOrange);
 fitf->SetLineStyle(2);


 fitf->SetParLimits( 0, -5, 5 );
 fitf->SetParLimits( 1, 3.0, 20 );
 //fitf->SetParLimits( 2, 4.5, 6.5 );
 //fitf->SetParLimits( 3, 550, 700 );
 fitf->FixParameter( 2, 6.0 );
 //fitf->FixParameter( 3, 50 );
 //fitf->SetParLimits( 4, 1350, 1600 );
 //fitf->FixParameter( 5, 4.0 );
 //fitf->FixParameter( 6, 4.0 );
 //fitf->FixParameter( 6, 28 );

 TCanvas *c = new TCanvas("c", "c", 600, 500);
 m_profile->Draw("hist E1 same");

 int fitstatus = m_profile->Fit("allComponents", "RN", "", -25, 1200);

 int counts = 0;

 /*
 while( fitstatus < 4 && counts < 3  ){

   for( int i=0; i<m_pars; i++ ){
     fitf->SetParameters( i, fitf->GetParameter(i) );
   }

   fitf->SetParLimits( 0, -5, 5 );
   fitf->SetParLimits( 1, 3.0, 20 );
   fitf->FixParameter( 2, 6.0 );

   m_profile->Fit("allComponents", "RN", "", -25, 1200);


   counts++;
 }
 */

 std::cout << counts << " " << fitstatus << " " << fitf->GetChisquare() / fitf->GetNDF() << std::endl;


 TF1 *fcomp[3];

 for(int i = 0; i<3; i++) {

   char fname[100]; sprintf(fname, "component_%d", i );

   fcomp[i] = new TF1(fname, singleComponentFcn, -50, 2400, 4);

   fcomp[i]->SetParameter(0, fitf->GetParameter(0) );
   fcomp[i]->SetParameter(1, fitf->GetParameter(1) );
   fcomp[i]->SetParameter(2, fitf->GetParameter(i+2) );
   fcomp[i]->SetParameter(3, fitf->GetParameter(i+5) );

   fcomp[i]->SetLineColor(i+1);
   fcomp[i]->SetLineStyle(2);

   fcomp[i]->Draw("same");
   //m_profile->GetListOfFunctions()->Add(fcomp[i]);

 }

 //m_profile->GetListOfFunctions()->Add(fitf);

 fitf->Draw("same");

 c->SetLogy();

 // Define the cumulative repsonse function and the three independed components
 // for drawing the fit



};


//------------------------------------------------------------------------------



void fitScintillationSingle() {


  std::string ifilename = "./sample/scintillation_fit_run1717.root";
  std::string histname = "prof_daq160";

  TFile *tfile = new TFile(ifilename.c_str(), "READ");

  TProfile *prof = (TProfile*)tfile->Get( histname.c_str() );



  // feed the profile to the fitter
  Scintillation s;
  s.fitResponse( prof );

}
