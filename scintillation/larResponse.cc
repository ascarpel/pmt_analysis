////////////////////////////////////////////////////////////////////////////////
// Macro that interfaces the output of the CAEN decoder script and produces
// a file format compatible with this small analysis framework
//
// mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////////////


#include "TFile.h"
#include "TChain.h"
#include "TProfile.h"
#include "TProfile.h"
#include "TMinuit.h"
#include "TMath.h"


class Scintillation {

  public:
    Scintillation();
    ~Scintillation();

    int fitResponse( TProfile *prof );

  private:

    TProfile *m_profile;
    static const int m_pars = 8;
    int m_ndf;
    double m_chi2;

    double singleComponentResponse(double x, double tm, double s, double tau, double a );

    static double singleComponentFcn( double* x, double * par );
    static double allComponentFcn( double* x, double * par );

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

  double fastc = fitobj->singleComponentResponse(x[0], par[0], par[1], par[2], par[5]);
  double intmc = fitobj->singleComponentResponse(x[0], par[0], par[1], par[3], par[6]);
  double slowc = fitobj->singleComponentResponse(x[0], par[0], par[1], par[4], par[7]);



  return par[8]+fastc + intmc + slowc ;

};



int Scintillation::fitResponse( TProfile *prof ) {


  // Pass the argument hist to the class member;
  m_profile = prof;

 int nbins = m_profile->GetNbinsX();
 double tlow = m_profile->GetBinLowEdge(0);
 double thigh = m_profile->GetBinLowEdge(nbins) + m_profile->GetBinWidth(nbins);


 double vthispar[9]  = { 0.0, 5.0, 6.0, 50, 1400, 3.34, 3.4, 23, 1e-4 };
 TF1 *fitf = new TF1("allComponents", Scintillation::allComponentFcn, tlow, thigh, 9);
 fitf->SetParameters( vthispar );
 fitf->SetLineColor(kOrange);
 fitf->SetLineStyle(2);


 fitf->SetParLimits( 0, -10, 5 );
 fitf->SetParLimits( 1, 1.0, 20 );
 fitf->FixParameter( 2, 6.0  );
 //fitf->SetParLimits( 2, 4.0, 30 );
 fitf->SetParLimits( 5, 0, 500 );
 fitf->SetParLimits( 6, 0, 500 );
 fitf->SetParLimits( 7, 0, 500 );

 fitf->SetParLimits( 8, 0.0 , 5e-1 );
 //fitf->FixParameter( 8, 0.0  );


 int fitstatus = m_profile->Fit("allComponents", "RN", "", -100, 3000);

 std::cout << fitstatus << " " << fitf->GetChisquare() / fitf->GetNDF() << std::endl;

 TF1 *fcomp[3];

 for(int i = 0; i<3; i++) {

   char fname[100]; sprintf(fname, "component_%d", i );

   fcomp[i] = new TF1(fname, singleComponentFcn, tlow, thigh, 4);

   fcomp[i]->SetParameter(0, fitf->GetParameter(0) );
   fcomp[i]->SetParameter(1, fitf->GetParameter(1) );
   fcomp[i]->SetParameter(2, fitf->GetParameter(i+2) );
   fcomp[i]->SetParameter(3, fitf->GetParameter(i+5) );

   fcomp[i]->SetLineColor(i+1);
   fcomp[i]->SetLineStyle(2);

   m_profile->GetListOfFunctions()->Add( fcomp[i] );

 }

 m_profile->GetListOfFunctions()->Add( fitf );


 //m_res->


 TCanvas *c = new TCanvas("c", "c", 600, 500);
 m_profile->Draw("E1 same");
 c->SetLogy();

 //TCanvas *c1 = new TCanvas("c1", "c1", 600, 200);
 //m_res->Draw("SAME");
 //c1->SetLogy();

 // Define the cumulative repsonse function and the three independed components
 // for drawing the fit

 return fitstatus;

};




//------------------------------------------------------------------------------




void getMetadata( std::string name, int &event, int &daqid, int &idx ) {

  std::string each;
  std::istringstream splitbuff(name);
  std::vector<std::string> tokens;
  while ( std::getline(splitbuff, each, '_') ) { tokens.push_back(each); }

  event = stoi( tokens[1].substr( tokens[1].find("e")+1) );
  daqid = stoi( tokens[2].substr( tokens[2].find("daq")+3) );
  idx = stoi( tokens[3].substr( tokens[3].find("idx")+3) );

}


//------------------------------------------------------------------------------


void saveToFile(string m_filename, std::map< int, std::vector<double> > m_results_map)
{
  std::ofstream myfile;
  myfile.open(m_filename.c_str());

  //Write the file header
  myfile << "daqid,tm,etm,sigma,esigma,tauf,etauf,taui,etaui,taus,etaus,af,eaf,ai,eai,as,eas,p,ps,fitstatus,chi2,ndof,nevents\n";

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



//------------------------------------------------------------------------------


void binning( double tlow, double thigh, double tsep, double narrowW, double coarseW, int nbins, double *binning ) {

  double lowEdge = tlow;
  for( int bin=0; bin < nbins+1; bin++){

    binning[bin] = lowEdge;

    if( lowEdge<tsep ){
      lowEdge += narrowW;
    }
    else{
      lowEdge += coarseW;
    }

  }

  return;
}


int findBin( double time, double *binning, int nbins ){

  int ret=0;
  for( int bin=0; bin<nbins-1; bin++ ){
    if( time >  binning[bin] && time < binning[bin+1]) { ret = bin; break; }
    else{ continue; }
  }

  return ret;

}


void mergeHist( TFile* tfile, std::map<int , TProfile*> & m_profile_map, std::map<int , TH2D*> & m_color_map ){

  std::cout << "MERGE HISTOGRAMS" << std::endl;

  TList *list = tfile->GetListOfKeys();

  // Check the first item and extract the boundaries
  TKey *firstKey = (TKey*)list->First();
  TProfile *firstHist = (TProfile*)tfile->Get( firstKey->GetName() );

  int nbins = firstHist->GetNbinsX();
  double width = firstHist->GetBinWidth(nbins);
  double tlow = firstHist->GetBinLowEdge(0)+width;
  double thigh = firstHist->GetBinLowEdge(nbins)+width;

  std::vector<int> activeBoards = { 1,2,3,4,5,6,7,8,9,10,11 };;
  std::vector<int> activeChannels = { 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 };


  //Get binning vector
  // degine the number of bins I want
  double tsep = -100; // moment when we start to apply the rebinning
  double narroWidth = 2.0;
  double coarseWidth = 4.0; // the new bin width

  int nbinsNarrow = (tsep-tlow)/narroWidth ;
  int nbinsCoarse = (thigh-tsep)/coarseWidth ;
  int nbinsNew = nbinsNarrow+nbinsCoarse;

  //std::cout << nbinsNarrow << " " << nbinsCoarse << " " << nbinsNew << std::endl;

  double newBinning[nbinsNew+1];

  binning( tlow, thigh, tsep, narroWidth, coarseWidth, nbinsNew, newBinning );

  for( int board : activeBoards ){
    for( int channel : activeChannels ){

      int daqid = board*16+channel;

      char pname[20]; sprintf(pname, "prof_daq%d", daqid);
      TProfile *profile = new TProfile(pname, pname, nbinsNew, newBinning);
      m_profile_map[daqid] = profile;

      char cname[20]; sprintf(cname, "cmap_daq%d", daqid);
      TH2D *cmap = new TH2D(cname, cname, 10, 0, 10, 100, 0.001, 1);
      m_color_map[daqid] = cmap;

    }
  }


  // Now we start iterating on the entries
  TIter iter(list->MakeIterator());

  while( TObject* obj = iter() ) {


    TKey* theKey = (TKey*)obj;

    std::string histname( theKey->GetName() );

    int event, daqid, idx;
    getMetadata( histname, event, daqid, idx );


    // We get the raw waveform
    TProfile *tmphist = (TProfile*)tfile->Get( theKey->GetName() );
    //tmphist->Scale(1./tmphist->GetMaximum());


    // Fill the correct bin
    for( int bin=0; bin<tmphist->GetNbinsX(); bin++ ) {

        double value = tmphist->GetBinContent(bin);
        double time = tmphist->GetBinCenter(bin);

        //cout << time << " " << m_profile_map[daqid]->FindBin(time) << endl;

        m_profile_map[daqid]->Fill(time, value);
        //m_color_map[daqid]->Fill(thisBin, value );

      }
  }

  return ;

}



void larResponse() {


  int run = 1882;

  std::string fittype = "settau";
  std::string ifilename="./sample/coincidence_waveforms_run"+to_string(run)+".root";
  std::string outfilename="./sample/scintillation_"+fittype+"_run"+to_string(run)+".root";
  std::string dbfilename =  "./sample/db_"+fittype+"_run"+to_string(run)+".csv";

  //Read from file
  TFile *tfile = new TFile(ifilename.c_str());


  // HERE WE MERGE THE HISTOGRAMS ----------------------------------------------
  // A tprofile is creaded for each daqid. If more pulses are found on one event
  // on the same tube, they are merged togheter

  std::map<int , TProfile*> m_profile_map;
  std::map<int , TH2D*> m_color_map;
  std::map<int , std::vector<double> > m_results_map;
  mergeHist( tfile,  m_profile_map, m_color_map );


  // HERE WE DO THE FIT --------------------------------------------------------
  std::cout << "FITTING: " << std::endl;


  TFile *ofile = new TFile(outfilename.c_str(), "RECREATE");
  for( auto & profileIt : m_profile_map ) {

    int daqid = profileIt.first;
    TProfile *prof = profileIt.second;
    TH2D *cmap = m_color_map[daqid];

    // Just fit profile with the entries
    if(prof->GetEntries() > 0 && daqid == 120) {

      std::cout << " >>>> " << daqid << " " << prof->GetName() << " " << prof->GetEntries() / 2000 << std::endl;

      prof->Scale(1./prof->GetMaximum());

      // Here we do the fit
      Scintillation myScintillation;
      int fitstatus = myScintillation.fitResponse( prof );

      // Here we read the parameters <<< this is a good exit point if you want to save the paramters to file
      std::cout << "  Status of the fit: " << fitstatus << std::endl;

      TF1 *fitf = prof->GetFunction("allComponents");


      std::vector<string> vnames={ "#t_m","#sigma", "#tau_{f}", "#tau_{i}", "#tau_{s}", "#a_{f}",   "#a_{i}",   "#a_{s}", "p"  };

      for( int i=0; i<vnames.size(); i++ ) {
        std::cout << "  " << vnames[i] << " " << fitf->GetParameter(i) << " " << fitf->GetParError(i) << endl;
        m_results_map[daqid].push_back( fitf->GetParameter(i) );
        m_results_map[daqid].push_back( fitf->GetParError(i) );

      }

      std::cout << "  " << fitf->GetChisquare() << " " <<  fitf->GetNDF() << std::endl;
      std::cout << "\n";

      m_results_map[daqid].push_back( fitstatus );
      m_results_map[daqid].push_back( fitf->GetChisquare() );
      m_results_map[daqid].push_back( fitf->GetNDF() );
      m_results_map[daqid].push_back( prof->GetEntries() );


      prof->Write();
      cmap->Write();


    }
  }

  std::cout << "Save to file" << std::endl;
  saveToFile(dbfilename, m_results_map);

  std::cout << "ALL DONE" << std::endl;

}
