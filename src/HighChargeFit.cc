////////////////////////////////////////////////////////////////////////////////
// HighChargeFit Methods implementation
//
// mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////////////

#include "TMinuit.h"
#include "TFitter.h"
#include "HighChargeFit.h"
#include "TPaveText.h"

static HighChargeFit* fitobj_hc;


HighChargeFit::HighChargeFit()
{
  fitobj_hc=this;
};


//------------------------------------------------------------------------------


HighChargeFit::HighChargeFit( vector<TH1D*> hist_array, bool norm )
  : m_hist_array(hist_array)
  , m_normalize(norm)
{
  fitobj_hc=this;

  if(norm){
    for( TH1D* hist : m_hist_array ){
      hist->Scale(1./hist->GetEntries());
    }
  }



};



//------------------------------------------------------------------------------


HighChargeFit::~HighChargeFit()
{
  clean();
};


//------------------------------------------------------------------------------


void HighChargeFit::loadInitialConditions()
{

  double vstart[m_parameters]={ 18,
                            0.5, 0.5, m_hist_array[0]->Integral()*0.5,
                            0.6, 0.5, m_hist_array[1]->Integral()*0.5,
                            0.7, 0.5, m_hist_array[2]->Integral()*0.5,
                            0.8, 0.5, m_hist_array[3]->Integral()*0.5 };
  double step[m_parameters]={ 0.1,
                              0.01, 0.01, 0.1,
                              0.01, 0.01, 0.1,
                              0.01, 0.01, 0.1,
                              0.01, 0.01, 0.1 };
  double minVal[m_parameters]={ 8,
                                0.1,  0.1,  0.0,
                                0.2,  0.1,  0.0,
                                0.3,  0.1,  0.0,
                                0.4,  0.1,  0.0};
  double maxVal[m_parameters]={ 28,
                                2.0, 3, m_hist_array[0]->Integral()*2.0,
                                2.2, 3, m_hist_array[1]->Integral()*2.0,
                                2.4, 3, m_hist_array[2]->Integral()*2.0,
                                2.6, 3, m_hist_array[3]->Integral()*2.0 };

  string parName[m_parameters]={ "npe",
                             "q1", "w1", "a1",
                             "q2", "w2", "a2",
                             "q3", "w3", "a3",
                             "q4", "w4", "a4" };

  for(int i=0; i<m_parameters; i++)
  {
    m_vstart[i] = vstart[i];
    m_step[i] = step[i];
    m_minval[i] = minVal[i];
    m_maxval[i] = maxVal[i];
    m_parname[i] = parName[i];
  }

};



//------------------------------------------------------------------------------


double HighChargeFit::jointFit(double npe,
                           double Q1, double w1, double a1,
                           double Q2, double w2, double a2,
                           double Q3, double w3, double a3,
                           double Q4, double w4, double a4 )
{

  const int npoints = 4;

  // I have these parameters:
  //double npe = mu; // mean npe
  double q[npoints]; q[0]=Q1; q[1]=Q2; q[2]=Q3; q[3]=Q4;
  double sigma[npoints]; sigma[0]=w1; sigma[1]=w2; sigma[2]=w3; sigma[3]=w4;
  double amp[npoints]; amp[0]=a1; amp[1]=a2; amp[2]=a3; amp[3]=a4;

  // compute the chi2
  int ndf=0;
  double chi2[4];

  for(int i=0; i<npoints; i++)
  {
    chi2[i]=0.0;
    int nbins = (int)this->m_hist_array[i]->GetNbinsX();//200;

    for(int j=0; j<nbins-1; j++)
    {
      double x = this->m_hist_array[i]->GetBinCenter(j+1);
      double y = this->m_hist_array[i]->GetBinContent(j+1);
      double ey = this->m_hist_array[i]->GetBinError(j+1);
      double prediction = 0;
      double p2 = 0;

      for(int k=1;k<50;k++){
        p2 += (TMath::Power(npe,k)*TMath::Exp(-1.0*npe)/TMath::Factorial(k)
          *TMath::Exp(-1.0*(x-q[i]*k)*(x-q[i]*k)
            /(2.0*k*sigma[i]*sigma[i]))/(sigma[i]*TMath::Sqrt(2.0*TMath::Pi()*k)) ) ;
      }

      prediction = p2*amp[i];

      if(ey!=0){
        chi2[i] += TMath::Power( y-prediction , 2)/(ey*ey); // the chi2 function
        ndf++;
      }
    }
  }

  fitobj_hc->m_ndf=ndf;
  return chi2[0]+chi2[1]+chi2[2]+chi2[3];

};


//------------------------------------------------------------------------------


void HighChargeFit::jointFitFunc(int &npar, double * deriv, double &f, double * par, int iflag)
{
  f = fitobj_hc->jointFit(par[0],
               par[1], par[2], par[3],
               par[4], par[5], par[6],
               par[7], par[8], par[9],
               par[10], par[11], par[12] );
};


//------------------------------------------------------------------------------


double HighChargeFit::singleFitFunc(double* x, double * par)
{
  double npe = par[0];
  double q = par[1];
  double sigma = par[2];
  double amplitude = par[3];

  double prediction = 0;
  for(int k=1;k<100;k++){
    prediction += (TMath::Power(npe,k)*TMath::Exp(-1.0*npe)/TMath::Factorial(k)
      *TMath::Exp(-1.0*(x[0] - q*k)*(x[0] - q*k)
        /(2.0*k*sigma*sigma))/(sigma*TMath::Sqrt(2.0*TMath::Pi()*k)));
  }

  return prediction*amplitude ;
};


//------------------------------------------------------------------------------


void HighChargeFit::multiHistFit( )
{

  loadInitialConditions();

  // Use TMinuit for the minimization
  TMinuit* minimizer = new TMinuit(m_parameters);

  minimizer->SetFCN(HighChargeFit::jointFitFunc);

  for(int i=0; i<m_parameters; i++){
    minimizer->DefineParameter(i, m_parname[i].c_str(), m_vstart[i],
                                          m_step[i],  m_minval[i], m_maxval[i]);
  }

  // do the Fit
  minimizer->Command("SET PRINT 1");//
  minimizer->mnscan();
  minimizer->SetPrintLevel(0);

  Int_t ierrs=0;
  minimizer->mnexcm("SET NOWarnings", 0, 0, ierrs);
  int fitstatus= minimizer->Migrad();


  for(int i=0; i<m_parameters; i++){
    minimizer->GetParameter(i, m_result[i], m_errors[i]);
  }

  // Continue minimization, but slower, if if didn't converge
  if(fitstatus==4)
  {
    for(int i=0; i<m_parameters; i++){
      minimizer->DefineParameter(i, m_parname[i].c_str(), m_result[i],
                                      m_step[i]/2.0,  m_minval[i], m_maxval[i]);
    }
    fitstatus = minimizer->Migrad();
    cout << "status of the second fit " << fitstatus << endl;
    for(int i=0; i<m_parameters; i++){
      minimizer->GetParameter(i, m_result[i], m_errors[i]);
    }
  }

  // We might need to minimize again, now move close to the found npe value
  if(fitstatus==4 || ((int)m_result[0] == (int)m_minval[0]) ){

      double newpars[m_parameters];
      double newsteps[m_parameters];
      double newmax[m_parameters];
      double newmin[m_parameters];

      for(int i=0; i<m_parameters; i++){
        newpars[i] = m_result[i];
        newsteps[i] = m_step[i];
        newmax[i] = m_minval[i];
        newmin[i] = m_maxval[i];
      }

      newmax[0] = newpars[0]+4;
      newmin[0] = newpars[0]-4;
      newsteps[0] = 0.01;

      for(int i=0; i<m_parameters; i++){
        minimizer->DefineParameter(i, m_parname[i].c_str(), newpars[i], newsteps[i],
                               newmin[i], newmax[i]);
      }
      fitstatus = minimizer->Migrad();
      cout<<"status of the third fit "<<fitstatus<<endl;
      for(int i=0; i<m_parameters; i++){
        minimizer->GetParameter(i, m_result[i], m_errors[i]);
      }
  }

  // If still gives troubles
  if(fitstatus==4)
  {
    for(int i=0; i<m_parameters; i++){
      minimizer->DefineParameter(i, m_parname[i].c_str(), m_result[i],
                                      m_step[i]/2.0,  m_minval[i], m_maxval[i]);
    }
    fitstatus = minimizer->Migrad();
    cout << "status of the last fit " << fitstatus << endl;
    for(int i=0; i<m_parameters; i++){
      minimizer->GetParameter(i, m_result[i], m_errors[i]);
    }
  }


  // Here we get the result from the fit
  double amin,edm,errdef;
  int nvpar,nparx,icstat;
  minimizer->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  cout << amin << "\t" << m_ndf-m_parameters << "\t"
       << edm << "\t" << errdef << "\t" << nvpar << "\t"
       << nparx << "\t" << icstat << endl;

  cout<<"---fit results "<<endl;

  for(int i=0; i<m_parameters; i++){
    cout<<m_parname[i]<<"\t"<<m_result[i]<<"\t"<<m_errors[i]<<endl;
  }

  m_chi2=amin;
  m_ndf -= m_parameters;
  m_fitstatus = fitstatus;
};


//------------------------------------------------------------------------------


void HighChargeFit::getCanvas()
{

  string cname(m_hist_array[0]->GetName());

  //Dereference the histogram to a static copy
  TH1D hist_array[4];
  hist_array[0] = *m_hist_array[0];
  hist_array[1] = *m_hist_array[1];
  hist_array[2] = *m_hist_array[2];
  hist_array[3] = *m_hist_array[3];

  TCanvas *c = new TCanvas(cname.c_str(), cname.c_str(), 600, 400);
  hist_array[0].SetLineColor(1);
  hist_array[0].Draw("hist");
  hist_array[1].SetLineColor(2);
  hist_array[1].Draw("hist same");
  hist_array[2].SetLineColor(3);
  hist_array[2].Draw("hist same");
  hist_array[3].SetLineColor(4);
  hist_array[3].Draw("hist same");

  TF1* func[4];
  char text[50];
  for(int i=0; i<4; i++)
  {
    sprintf(text, "func%d", i);
    func[i] = new TF1(text, HighChargeFit::singleFitFunc, 0, 150, 4);
    func[i]->SetNpx(2000);
    func[i]->SetParameters(m_result[0], m_result[i*3+1], m_result[i*3+2], m_result[i*3+3] );
    func[i]->SetLineColor(i+1);
    func[i]->SetLineStyle(2);
    func[i]->Draw("same");
  }

  char gname[100];
  TPaveText* ptGainPar = new TPaveText(0.6, 0.55, 0.9, 0.9,"brNDC");
  sprintf(gname,"#chi^2 / ndf = %.0f/%d", m_chi2, m_ndf-13);
  ptGainPar->AddText(gname);
  sprintf(gname,"#mu = %.2f #pm %.2e", m_result[0], m_errors[0]);
  ptGainPar->AddText(gname);
  sprintf(gname,".................................");
  ptGainPar->AddText(gname);
  sprintf(gname,"#color[1]{q_{1} = %.2f #pm %.2e}", m_result[1], m_errors[1]);
  ptGainPar->AddText(gname);
  sprintf(gname,"#color[2]{q_{2} = %.2f #pm %.2e}", m_result[4], m_errors[4]);
  ptGainPar->AddText(gname);
  sprintf(gname,"#color[3]{q_{3} = %.2f #pm %.2e}", m_result[7], m_errors[7]);
  ptGainPar->AddText(gname);
  sprintf(gname,"#color[4]{q_{4} = %.2f #pm %.2e}", m_result[10], m_errors[10]);
  ptGainPar->AddText(gname);

  ptGainPar->SetBorderSize(0);
  ptGainPar->SetFillColor(0);
  ptGainPar->SetTextFont(42);
  //ptGainPar->SetTextSize(14);
  ptGainPar->SetTextAlign(12);
  ptGainPar->Draw();

  c->Write();

}
