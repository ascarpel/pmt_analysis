 ////////////////////////////////////////////////////////////////////////////////
//
// Process the ROOT files produced during the calibration and returns a file
// with the appropriate histograms. The input is a list of files sorted in
// increasing order of HV
// Usage: ./makeHist -i filelist.txt -o output/path/file.root -t treename
// the threename allows to choose between direct and indirect light
//
// mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"

#include "Pmt.h"


//------------------------------------------------------------------------------

/*
bool isGoodFile( int m_run, int m_subrun, int m_optical_channel )
{

  bool isGood=false;

  // Want to exclude side files due to unwanted pedestal events
  string label = make_string( m_run, m_optical_channel );
  std::pair<int, int> edgefiles = m_optch_map[label];

  if( (edgefiles.first < m_subrun) & (edgefiles.second > m_subrun) ){
    isGood = true;
  }

  return isGood;
}
*/



//------------------------------------------------------------------------------


int main( int argc, char* argv[] ) {

  // Define here some quantities to have them handy in case they have to change
  // One may also consider to have them configured exterenally
  int npmts = 360;
  double hmin = -2; double hmax = 100;
  int nbins = 300;


  // Read the input
  std::string ifilename, ofilename, ttreename;
  for ( int i=1; i<argc; i=i+2 )
  {
    if      ( std::string(argv[i]) == "-i" ) ifilename = argv[i+1];
    else if ( std::string(argv[i]) == "-o" ) ofilename = argv[i+1];
    else if ( std::string(argv[i]) == "-t" ) ttreename = argv[i+1];
    else {
      std::cout << "Unknown option " << argv[i+1] << std::endl;
      return 1;
    }
  }


  //Read the list of input files and create an array
  std::vector<std::string> files;
  std::ifstream filep(ifilename.c_str()); std::string line;
  while( std::getline( filep, line ) ){
    files.push_back(line);
  }


  size_t nfiles = files.size();
  std::map < int, std::vector<PMT*> > m_pmts;

  std::cout << "\nREAD THE FILES " << std::endl;

  //Process each file and create histograms
  for(size_t i=0; i<files.size(); i++ ){

    std::cout << " >> Processing file: " << files[i] << std::endl;

    for( int ch=1; ch<npmts+1; ch++ )
    {
      std::string name = "hist_"+to_string(i)+"_pmtid_"+to_string(ch);
      PMT *mypmt = new PMT(ch, name);
      m_pmts[ch].push_back( mypmt );
    }


    TFile *tfile = new TFile( files[i].c_str(), "OPEN" );
    TTree *ttree = (TTree*)tfile->Get(ttreename.c_str());

    double m_integral;
    int m_optical_channel;
    int m_pmt_illuminated;
    int m_board_illuminated;
    int m_run, m_subrun;

    ttree->SetBranchAddress("m_integral", & m_integral);
    ttree->SetBranchAddress("m_optical_channel", & m_optical_channel);
    ttree->SetBranchAddress("m_pmt_illuminated", & m_pmt_illuminated);
    ttree->SetBranchAddress("m_run", & m_run);
    ttree->SetBranchAddress("m_subrun", & m_subrun);

    // NOW we loop over the entries. the total number of entries is:
    // - For the direct light case nevents*10 ( laser illuminates groups of 10 )
    // - For the indirect ligth case nevents*(180-10).
    // i don't know nevents with precision, for the direct light case is around
    // 10k events per PMT.
    // Here we process all events. The total number of events can be capped by
    // checking how many entries a pmt has accumulated at every iteration
    // Uncomment the lines with the symbol >> to apply the cap

    int nentries = ttree->GetEntries();
    // >> int entriesPms = 100;

    for(int entry=0; entry<nentries; entry++  ) {

      ttree->GetEntry(entry);

      // Cap the total number of events
      // >>> if( m_pmts[m_pmt_illuminated][i]->getArrayCharge().size() > entriesPms) {
      // >>   continue;
      // >> }

      m_pmts[m_pmt_illuminated][i]->loadCharge( m_integral );

    } // end loop over entries

    tfile->Close();

  } // end loop over files



  std::cout << "\nMAKE HISTOGRAMS" << std::endl;

  TFile *outfile = new TFile(ofilename.c_str(), "RECREATE");

  for(int ch=1; ch<npmts+1; ch++)
  {

    cout << " >>> Make histograms for channel: " << ch << endl;

    if( m_pmts[ch][0]->getArrayCharge().size() == 0){ continue; }

    // Find the extremes
    double m_min[3]; double m_max[3];
    for(size_t i=0; i<files.size(); i++) {

      std::vector<double> m_tmp_array = m_pmts[ch][i]->getArrayCharge();
      m_pmts[ch][i]->findExtremes(m_tmp_array, m_min[i], m_max[i]);
      if( m_min[i] < -2 ){ m_min[i] = hmin; }
      if( m_max[i] > 100 ){ m_max[i] = hmax; }

    }

    // Create and save the histogram with a valid range for all the HV points
    for(size_t i=0; i<files.size(); i++) {
      m_pmts[ch][i]->initChargeHist(nbins, m_min[0], m_max[2]);
      TH1D *tmphist = m_pmts[ch][i]->getHistCharge();
      tmphist->Write();
    }
  }

  outfile->Close();

  std::cout << "\nALL DONE" << std::endl;

  return 0;
}
