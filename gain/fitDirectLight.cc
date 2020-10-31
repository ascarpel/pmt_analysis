#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <iterator>
#include <string>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TGraphErrors.h"
#include "TPaveText.h"
#include "Pmt.h"
#include "ChargeFit.h"
#include "HighChargeFit.h"


//------------------------------------------------------------------------------

void loadPMTData(std::map< int, std::vector<double> > & hvmap)
{

  string filename="../dbase/hvpoint/hvdbase.csv";
  string line="";

  ifstream file(filename);

  while(getline(file, line))
  {
    istringstream s_stream(line);
    string token;
    vector<string> vec;

    while( getline(s_stream, token, ',') ){
      vec.push_back( token );
    }

    int pmt = int(stod(vec[0]));
    for(int i=1; i<7; i++){
      hvmap[pmt].push_back( stod(vec[i]) );
    }
  }

  return;
}

//------------------------------------------------------------------------------

void fitGainCurve( int ch,  double hv[4], double q[4], double eq[4],
                        std::map< int, std::vector<double> > & m_results_map ) {


    double hvmax=0.0; double hvmin=9999.0;
    for( int i=0; i<4; i++ ){
      if( hv[i] > hvmax ){
        hvmax = hv[i];
      }
      if(hv[i] < hvmin){
        hvmin = hv[i];
      }
    }

    char gname[100];
    sprintf(gname, "gain_ch%d", ch);
    TCanvas* cGain = new TCanvas(gname, gname, 500, 400);
    sprintf(gname, "c_gain_ch%d",ch);
    cGain->SetName(gname);
    TGraphErrors* gGain = new TGraphErrors();
    gGain->SetName(gname); gGain->SetTitle("");
    gGain->GetXaxis()->SetTitle("HV (V)");
    gGain->GetYaxis()->SetTitle("Gain (10^{7})");
    for(int i=0; i < 4; i++){
        gGain->SetPoint(i, hv[i], q[i]/1.6);
        gGain->SetPointError(i, 0, eq[i]/1.6);
    }
    gGain->Draw("ALP");
    gGain->SetMarkerStyle(8);
    // fit the gain curve
    TF1* fGain = new TF1("fGain","[0]*TMath::Power(x,[1])",1,2000);
    fGain->SetParNames("a","b");
    fGain->SetParameters(1e-24,7.5);
    double a=0, b=0;
    for(int i=0; i<4;i++){
      a=fGain->GetParameter(0);
      b=fGain->GetParameter(1);
      fGain->SetParameters(a,b);
      gGain->Fit("fGain","QR","",hvmin-10, hvmax+10);
    }
    a=fGain->GetParameter(0);
    b=fGain->GetParameter(1);
    cGain->Update();
    TPaveText* ptGainPar = new TPaveText(0.2,0.7,0.6,0.9,"brNDC");
    sprintf(gname,"PMT: %d", 360-ch);
    ptGainPar->AddText(gname);
    sprintf(gname,"chi2/ndf = %.2f/%d",fGain->GetChisquare(), fGain->GetNDF());
    ptGainPar->AddText(gname);
    sprintf(gname,"a = %e #pm %e", a, fGain->GetParError(0));
    ptGainPar->AddText(gname);
    sprintf(gname,"k = %.3f #pm %.3f", b, fGain->GetParError(1));
    ptGainPar->AddText(gname);
    double V0 = TMath::Power(1.0/a, 1.0/b);
    double dvda = TMath::Power(1.0/a, -1.0+1.0/b)/(a*a*b);
    double dvdb = TMath::Power(1.0/a, 1.0/b)*TMath::Log(1.0/a)/b/b;
    double V0err = TMath::Sqrt( dvda*dvda*fGain->GetParError(0)*fGain->GetParError(0) + dvdb*dvdb*fGain->GetParError(1)*fGain->GetParError(1) );
    sprintf(gname,"V(10^{7})=%.1f #pm %.1f", V0, V0err);
    ptGainPar->AddText(gname);
    ptGainPar->SetBorderSize(0);
    ptGainPar->SetFillColor(0);
    ptGainPar->SetTextFont(42);
    ptGainPar->SetTextAlign(12);

    ptGainPar->Draw();
    cGain->Update();
    cGain->Write();

    m_results_map[ch].push_back( a*1e24 );
    m_results_map[ch].push_back( b );
    m_results_map[ch].push_back( fGain->GetParError(0)*1e24 );
    m_results_map[ch].push_back( fGain->GetParError(1) );
    m_results_map[ch].push_back( V0 );
    m_results_map[ch].push_back( V0err );
    m_results_map[ch].push_back( fGain->GetChisquare()/fGain->GetNDF() );

}



//------------------------------------------------------------------------------


bool isHighCharge( vector<TH1D*> hist_array )
{

  double xlow = hist_array[0]->GetXaxis()->GetBinCenter(0);

  if(xlow > 0){
    return true;
  }

  int bin1 = hist_array[0]->GetXaxis()->FindBin(-0.1);
  int bin2 = hist_array[0]->GetXaxis()->FindBin(0.1);

  if( bin1 == 0 || bin2 ==0  ){
    return true;
  }

  double ratio = hist_array[3]->Integral(bin1, bin2)/hist_array[3]->Integral();

  return (ratio <= 0.01) ? true : false;

}


//------------------------------------------------------------------------------

void saveToFile(string m_filename, std::map< int, std::vector<double> > m_results_map)
{
  std::ofstream myfile;
  myfile.open(m_filename.c_str());

  //Write the file header
  myfile << "pmt,n1,n2,n3,n4,type,npe,q1,q2,q3,q4,eq1,eq2,eq3,eq4,cgchi2,ndf,fitstatus,warmV,a,k,ea,ek,V,eV,chi2\n";

  // Write the lines
  for( auto item : m_results_map ){
    int pmtid = item.first;
    std::string line = to_string(pmtid);
    for( auto v : item.second ){
      line += ","+to_string(v);
    }
    myfile << line+"\n";
  }

  myfile.close();

}


////////////////////////////////////////////////////////////////////////////////


void parseKeyName( std::string name, int & ch, int & filenum ) {

  std::stringstream ss(name);
  std::string token;
  while (std::getline(ss, token, '_')) {
    std::cout << token << std::endl;
  }


}


int main( int argc, char* argv[] ) {

  // As usual define some handy variables in here
  std::map< int, std::vector<double> > m_results_map;
  const int nhvpoints = 4;
  int startch = 180; // <<< Edit to select a start channel
  int endch = 359; // << Edit to selec an end channel

  std::string histfilename, outfilename, databasename;
  for ( int i=1; i<argc; i=i+2 )
  {
    if      ( std::string(argv[i]) == "-i" ) histfilename = argv[i+1];
    else if ( std::string(argv[i]) == "-o" ) outfilename = argv[i+1];
    else if ( std::string(argv[i]) == "-d" ) databasename = argv[i+1];
    else if ( std::string(argv[i]) == "-s" ) startch = atoi(argv[i+1]);
    else if ( std::string(argv[i]) == "-e" ) endch = atoi(argv[i+1]);
    else {
      std::cout << "Unknown option " << argv[i+1] << std::endl;
      return 1;
    }
  }


  //----------------------------------------------------------------------------


  std::cout << "\nLOAD METADATA: " << std::endl;

  std::map< int, std::vector<double> > hvmap;
  loadPMTData(hvmap);


  //----------------------------------------------------------------------------


  std::cout << "\nLOAD HISTOGRAMS " << std::endl;

  std::map< int, std::vector<TH1D*> > m_hists;

  TFile *tfile = new TFile(histfilename.c_str(), "READ" );
  TList *list = tfile->GetListOfKeys(); TIter iter(list->MakeIterator());

  while(TObject* obj = iter() ) {

    TKey* theKey = (TKey*)obj;

    std::string keyname(theKey->GetName());

    // Parse the histname: assumes format hist_<filenum>_pmtid_<ch>
    std::vector<std::string> stringbuff;
    std::stringstream ss(keyname); std::string token;
    while (std::getline(ss, token, '_')) {
      stringbuff.push_back( token );
    }

    int ch = stoi( stringbuff[1] );
    m_hists[ch].push_back( (TH1D*)tfile->Get(theKey->GetName()) );

  }


  //----------------------------------------------------------------------------

  std::cout << "\nFIT DATA in range [ " << startch << ", " << endch << "]:" << std::endl;

  TFile *outfile = new TFile(outfilename.c_str(), "RECREATE");

  for(int ch=startch; ch<endch+1; ch++) {

    m_results_map[ch].clear();

    cout << " >>> Fit channel_id: " << ch << " pmt_id: " << 360-ch << endl;

    double nomhv = hvmap[360-ch][0];

    double hv[nhvpoints];
    double q[4] = {0.0, 0.0, 0.0, 0.0};
    double eq[4] = {0.0, 0.0, 0.0, 0.0};

    // This part is crucual to pass the histogram from the file to the fitter in good order
    std::vector<int> points = {2,3,4,5}; // this is the fixed order which appears in the file
    std::vector<int> index = {0,1,2,3}; // this is the actual index to sort the histograms in ascending order
    std::vector<TH1D*> m_fit_hists; // << Allocates only the actual histograms for the fit

    cout << "Nominal (warm): " << nomhv << " " << "Voltage points: ";
    for(size_t i=0; i<points.size(); i++ ){

      cout << hvmap[360-ch][index[i]+1] << " ";

      hv[i] = hvmap[360-ch][index[i]+1];
      m_fit_hists.push_back(m_hists[ch][index[i]]);
    }
    cout << "\n";


    //check if it is low or high charge
    if ( isHighCharge( m_fit_hists ) && ch != 89 )
    {

      std::cout << "IS HIGH CHARGE" << std::endl;

      for(int i=0; i<nhvpoints; i++) {
        m_results_map[ch].push_back( m_fit_hists[i]->GetEntries() );
      }

      HighChargeFit myfitter( m_fit_hists, true );
      myfitter.multiHistFit();

      //outfile->cd();
      myfitter.getCanvas(ch);

      double npe, enpe;
      myfitter.getParameters(0, npe, enpe);
      myfitter.getParameters(1, q[0], eq[0]);
      myfitter.getParameters(4, q[1], eq[1]);
      myfitter.getParameters(7, q[2], eq[2]);
      myfitter.getParameters(10, q[3], eq[3]);

      double chi2=0;
      int ndf=0;
      myfitter.GetChisquare( chi2, ndf );


      m_results_map[ch].push_back( 0 );
      m_results_map[ch].push_back(npe);
      m_results_map[ch].push_back( q[0] );
      m_results_map[ch].push_back( q[1] );
      m_results_map[ch].push_back( q[2] );
      m_results_map[ch].push_back( q[3] );
      m_results_map[ch].push_back( eq[0] );
      m_results_map[ch].push_back( eq[1] );
      m_results_map[ch].push_back( eq[2] );
      m_results_map[ch].push_back( eq[3] );
      m_results_map[ch].push_back( chi2 );
      m_results_map[ch].push_back( ndf );
      m_results_map[ch].push_back( myfitter.GetFitstatus() );
      m_results_map[ch].push_back( nomhv );

    }
    else if( !isHighCharge( m_fit_hists ) )
    {

      std::cout << "IS LOW CHARGE" << std::endl;

      for(int i=0; i<nhvpoints; i++) {
        m_results_map[ch].push_back( m_fit_hists[i]->GetEntries() );
      }


      ChargeFit myfitter( m_fit_hists, true );

      //outfile->cd();
      myfitter.multiHistFit();

      myfitter.getCanvas(ch);

      double npe, enpe;
      myfitter.getParameters(0, npe, enpe);
      myfitter.getParameters(1, q[0], eq[0]);
      myfitter.getParameters(4, q[1], eq[1]);
      myfitter.getParameters(7, q[2], eq[2]);
      myfitter.getParameters(10, q[3], eq[3]);

      double chi2=0;
      int ndf=0;
      myfitter.GetChisquare( chi2, ndf );

      m_results_map[ch].push_back( 1 );
      m_results_map[ch].push_back(npe);
      m_results_map[ch].push_back( q[0] );
      m_results_map[ch].push_back( q[1] );
      m_results_map[ch].push_back( q[2] );
      m_results_map[ch].push_back( q[3] );
      m_results_map[ch].push_back( eq[0] );
      m_results_map[ch].push_back( eq[1] );
      m_results_map[ch].push_back( eq[2] );
      m_results_map[ch].push_back( eq[3] );
      m_results_map[ch].push_back( chi2 );
      m_results_map[ch].push_back( ndf );
      m_results_map[ch].push_back( myfitter.GetFitstatus() );
      m_results_map[ch].push_back( nomhv );

    }

    std::cout << "Fit gain curve" << std::endl;
    fitGainCurve( ch,  hv, q, eq, m_results_map);

  }

  tfile->Close();
  outfile->Close();


  std::cout << " Save results to file " << std::endl;

  saveToFile(databasename, m_results_map);


  std::cout << "All done!" << std::endl;

  return 0;

}
