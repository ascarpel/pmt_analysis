////////////////////////////////////////////////////////////////////////////////
// ChargeFit Methods implementation
//
// mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////////////

#include "ChargeFit.h"

static ChargeFit* fitobj;


ChargeFit::ChargeFit()
{
  fitobj=this;
};


//------------------------------------------------------------------------------


ChargeFit::ChargeFit( vector<TH1D*> hist_array, bool norm )
  : m_hist_array(hist_array)
  , m_normalize(norm)
{
  fitobj=this;

  if(norm){
    for( TH1D* hist : m_hist_array ){
      hist->Scale(1./hist->GetEntries());
    }
  }

  loadInitialConditions();
};



//------------------------------------------------------------------------------


ChargeFit::~ChargeFit()
{
  clean();
};


//------------------------------------------------------------------------------


void ChargeFit::loadInitialConditions()
{
  // Basic inital condition loader for the joint SPE fit

  // Select the pedestal boundaries
  int bin1[3]; int bin2[3];
  for(int i=0;i<3;i++){
    bin1[i] = m_hist_array[i]->GetXaxis()->FindBin(-0.2);
    bin2[i] = m_hist_array[i]->GetXaxis()->FindBin(0.2);
  }

  double vstart[m_parameters]={m_hist_array[0]->GetMean()/1.6,
                     1.0, 0.5, m_hist_array[0]->Integral()/400.0*40,
                     1.0, 0.5, m_hist_array[1]->Integral()/400.0*40,
                     1.0, 0.5, m_hist_array[2]->Integral()/400.0*40,
                     0.0, 0.1, m_hist_array[0]->Integral(bin1[0],bin2[0])/10.0*4, 0.05,
                     0.0, 0.1, m_hist_array[1]->Integral(bin1[1],bin2[1])/10.0*4, 0.05,
                     0.0, 0.1, m_hist_array[2]->Integral(bin1[2],bin2[2])/10.0*4, 0.05};
  double step[m_parameters]={0.02,
                     0.05, 0.1, 0.005,
                     0.05, 0.1, 0.005,
                     0.05, 0.1, 0.005,
                     0.005, 0.01, 0.0010, 0.001,
                     0.005, 0.01, 0.0010, 0.001,
                     0.005, 0.01, 0.0010, 0.001};
  double minVal[m_parameters]={0,
                     0.01, 0.001, 0,
                     0.01, 0.001, 0,
                     0.01, 0.001, 0,
                     -0.5, 0.001, 0, 0.04,
                     -0.5, 0.001, 0, 0.04,
                     -0.5, 0.001, 0, 0.04};
  double maxVal[m_parameters]={10,
                      10, 10, m_hist_array[0]->Integral()/4.0,
                      10, 10, m_hist_array[1]->Integral()/4.0,
                      10, 10, m_hist_array[2]->Integral()/4.0,
                      0.5, 0.5, m_hist_array[0]->Integral(bin1[0],bin2[0])*10.0, 0.2,
                      0.5, 0.5, m_hist_array[1]->Integral(bin1[1],bin2[1])*10.0, 0.2,
                      0.5, 0.5, m_hist_array[2]->Integral(bin1[2],bin2[2])*10.0, 0.2 };
  string parName[m_parameters]={"npe",
                      "q1", "w1", "a1",
                      "q2", "w2", "a2",
                      "q3", "w3", "a3",
                      "pedm1", "pedw1", "peda1", "pedt1",
                      "pedm2", "pedw2", "peda2", "pedt2",
                      "pedm3", "pedw3", "peda3", "pedt3" };

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


double ChargeFit::jointFit(double npe,
                           double Q1, double w1, double a1,
                           double Q2, double w2, double a2,
                           double Q3, double w3, double a3,
                           double pedm1, double pedw1, double pedamp1, double pedtau1,
                           double pedm2, double pedw2, double pedamp2, double pedtau2,
                           double pedm3, double pedw3, double pedamp3, double pedtau3)
{

  const int npoints = 3;

  // I have these parameters:
  //double npe = mu; // mean npe
  double q[npoints]; q[0]=Q1; q[1]=Q2; q[2]=Q3;
  double sigma[npoints]; sigma[0]=w1; sigma[1]=w2; sigma[2]=w3;
  double amp[npoints]; amp[0]=a1; amp[1]=a2; amp[2]=a3;
  double q0[npoints]; q0[0]=pedm1; q0[1]=pedm2; q0[2]=pedm3;
  double w0[npoints]; w0[0]=pedw1; w0[1]=pedw2; w0[2]=pedw3;
  double a0[npoints]; a0[0]=pedamp1; a0[1]=pedamp2; a0[2]=pedamp3;
  double tau0[npoints]; tau0[0]=pedtau1; tau0[1]=pedtau2; tau0[2]=pedtau3;

  // compute the chi2
  int ndf=0;
  double chi2[3];
  for(int i=0; i<npoints; i++){
    chi2[i]=0.0;
    int nbins = (int)this->m_hist_array[i]->GetNbinsX();//200;

    for(int j=0; j<nbins-1; j++){
      double x = this->m_hist_array[i]->GetBinCenter(j+1);
      double y = this->m_hist_array[i]->GetBinContent(j+1);
      double ey = this->m_hist_array[i]->GetBinError(j+1);
      double prediction = 0;

      // ideal response term
      double p2 = 0;
      for(int k=1;k<100;k++){
          p2 +=    (TMath::Power(npe,k)*TMath::Exp(-1.0*npe)/TMath::Factorial(k)
            *TMath::Exp(-1.0*(x-q[i]*k)*(x-q[i]*k)
              /(2.0*k*sigma[i]*sigma[i]))/(sigma[i]*TMath::Sqrt(2.0*TMath::Pi()*k)) ) ;
      }

      prediction = p2*amp[i]
        + a0[i]*w0[i]/tau0[i]*TMath::Sqrt(TMath::Pi()/2.0)
           *TMath::Exp(0.5*(w0[i]*w0[i])/(tau0[i]*tau0[i])-(x-q0[i])/tau0[i])
              *TMath::Erfc(TMath::Sqrt(1./2.)*(w0[i]/tau0[i]-(x-q0[i])/w0[i]));

      if(ey!=0){
        chi2[i] += TMath::Power( y-prediction , 2)/(ey*ey); // the chi2 function
        ndf++;
      }
    }
  }

  fitobj->m_ndf=ndf;
  return chi2[0]+chi2[1]+chi2[2];

};


//------------------------------------------------------------------------------


void ChargeFit::jointFitFunc(int &npar, double * deriv, double &f, double * par, int iflag)
{
  f = fitobj->jointFit(par[0],
               par[1], par[2], par[3],
               par[4], par[5], par[6],
               par[7], par[8], par[9],
               par[10], par[11], par[12], par[13],
               par[14], par[15], par[16], par[17],
               par[18], par[19], par[20], par[21] );
};


//------------------------------------------------------------------------------


double ChargeFit::singleFitFunc(double* x, double * par)
{
  double mu = par[0];
  double q = par[1];
  double sigma = par[2];
  double amplitude = par[3];
  double sum=0;
  for(int n=1; n<100; n++){
    sum += amplitude* (TMath::Power(mu,n)*TMath::Exp(-1.0*mu)/TMath::Factorial(n)
            *TMath::Exp(-1.0*(x[0]-q*n)*(x[0]-q*n)
              /(2.0*n*sigma*sigma))/(sigma*TMath::Sqrt(2.0*TMath::Pi()*n))) ;

  }

  sum = sum + par[6]*par[5]/par[7]*TMath::Sqrt(TMath::Pi()/2.0)*TMath::Exp(0.5*(par[5]*par[5])/(par[7]*par[7])-(x[0]-par[4])/par[7])*TMath::Erfc(TMath::Sqrt(1./2.)*(par[5]/par[7]-(x[0]-par[4])/par[5]));

  return sum ;
};


//------------------------------------------------------------------------------


void ChargeFit::multiHistFit( )
{

  loadInitialConditions();

  // Use TMinuit for the minimization
  TMinuit* minimizer = new TMinuit(m_parameters);
  minimizer->SetFCN(ChargeFit::jointFitFunc);

  for(int i=0; i<m_parameters; i++){
    minimizer->DefineParameter(i, m_parname[i].c_str(), m_vstart[i],
                                          m_step[i],  m_minval[i], m_maxval[i]);
  }


  // do the Fit
  minimizer->Command("SET PRINT -3");//
  Int_t e=0;
  minimizer->mnexcm("SET NOWarnings", 0, 0, e);

  int fitstatus= minimizer->Migrad();


  // Refit if doesn't converge
  if(fitstatus==4)
  {
    for(int i=0; i<m_parameters; i++)
    {
      minimizer->GetParameter(i, m_result[i], m_errors[i]);
      minimizer->DefineParameter(i, m_parname[i].c_str(), m_result[i],
                                      m_step[i]/2.0,  m_minval[i], m_maxval[i]);
    }
    fitstatus = minimizer->Migrad();
    cout << "status of the second fit " << fitstatus << endl;
  }

  // Now we get the paramters
  for(int i=0; i<m_parameters; i++){
    minimizer->GetParameter(i, m_result[i], m_errors[i]);
  }

  // Here we get the result from the fit
  double amin,edm,errdef;
  int nvpar,nparx,icstat;
  minimizer->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  cout<<amin<<"\t"<<m_ndf-m_parameters<<"\t"
              <<edm<<"\t"<<errdef<<"\t"<<nvpar<<"\t"<<nparx<<"\t"<<icstat<<endl;
  cout<<"---fit results "<<endl;
  for(int i=0; i<m_parameters; i++){
    cout<<m_parname[i]<<"\t"<<m_result[i]<<"\t"<<m_errors[i]<<endl;
  }

  m_chi2=amin;
  m_ndf -= m_parameters;
  m_fitstatus = fitstatus;

};


//------------------------------------------------------------------------------


void ChargeFit::getCanvas()
{

  string cname(m_hist_array[0]->GetName());

  TCanvas* c = new TCanvas(cname.c_str(), cname.c_str(), 600, 400);
  m_hist_array[0]->SetLineColor(1);
  m_hist_array[0]->Draw("hist");
  m_hist_array[0]->GetXaxis()->SetRangeUser(-2, 30);
  m_hist_array[1]->SetLineColor(2);
  m_hist_array[1]->Draw("hist same");
  m_hist_array[2]->SetLineColor(3);
  m_hist_array[2]->Draw("hist same");

  TF1* func[3];
  char text[50];
  for(int i=0; i<3; i++){
     sprintf(text, "func%d", i);
     func[i] = new TF1(text, ChargeFit::singleFitFunc, -2, 30, 8);
     func[i]->SetNpx(2000);
     func[i]->SetParameters(m_result[0],
                            m_result[i*3+1], m_result[i*3+2], m_result[i*3+3],
                            m_result[i*4+10], m_result[i*4+11], m_result[i*4+12], m_result[i*4+13]);

    func[i]->SetLineColor(i+1);
    func[i]->SetLineStyle(2);
    func[i]->Draw("same");

  }

  c->Write();

  return;

}
