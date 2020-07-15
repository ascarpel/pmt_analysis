////////////////////////////////////////////////////////////////////////////////
//
// usage makePulseHist.cc inputfile.root outputfile.root nevents applycut
//  - inputfile.root is a file with _pulses.root, it must have a ttree called calo/pulsetree
//  - outputfile.root is your own choice. It contains the histograms needed for the
// analysis
//  - nevents is an integer form -1 to inf. if the real number events of the file
// is lower than the number of events set, the total number of events is then
// processed. if nevents < 0, also all events are processed
// - applycut: select some pulses based on some criteria specified in the function
// passCuts.
//
// mailto: ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////////////

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include <stdio.h>
#include <iostream>


//------------------------------------------------------------------------------


struct Pulse{
  double peaktime;
  double amplitude;
  double integral;
  int pmt;
};


//------------------------------------------------------------------------------
// Class PMTDark here is just a placeholder to have all the histograms in the
// same place.

class PMTDark {

  public:

    PMTDark();
    PMTDark( int num );
    ~PMTDark();

    void loadPulse(Pulse pulse);
    void loadEvents(int nevents){ m_nevents = nevents; };
    void writeHist();

  private:
    void initHist();
    void clean();
    int m_num =0;
    int m_nevents =0;

    char m_hname[100];
    TH1D *hamplitude;
    TH1D *hamplitude_norm;
    TH1D *hintegral;
    TH1D *hintegral_norm;
    TH1D *htime;

};

PMTDark::PMTDark(){ this->initHist(); }
PMTDark::PMTDark(int num) : m_num(num){ this->initHist(); }
PMTDark::~PMTDark(){ this->clean(); }

void PMTDark::initHist(){

  sprintf(m_hname, "hamplitude_%d", m_num);
  hamplitude = new TH1D( m_hname, m_hname, 100, 0, 100 );

  sprintf(m_hname, "hamplitude_norm_%d", m_num);
  hamplitude_norm = new TH1D( m_hname, m_hname, 100, 0, 100 );

  sprintf(m_hname, "hintegral_%d", m_num);
  hintegral = new TH1D( m_hname, m_hname, 100, 0, 500 );

  sprintf(m_hname, "hintegral_norm_%d", m_num);
  hintegral_norm = new TH1D( m_hname, m_hname, 100, 0, 500 );

  sprintf(m_hname, "htime_%d", m_num);
  htime = new TH1D( m_hname, m_hname, 100, 0, 5000 );

}


void PMTDark::loadPulse(Pulse pulse)
{
  hamplitude->Fill( pulse.amplitude );
  hintegral->Fill( pulse.integral );
  htime->Fill( pulse.peaktime );
}


void PMTDark::writeHist()
{
  (*hamplitude).Write();
  (*hintegral).Write();
  (*htime).Write();

  if(m_nevents>0)
  {
    hamplitude->Scale( 1./ m_nevents );
    hintegral->Scale( 1./ m_nevents );

    sprintf(m_hname, "hamplitude_norm_%d", m_num);
    hamplitude->SetName(m_hname); hamplitude->SetTitle(m_hname);
    (*hamplitude).Write();

    sprintf(m_hname, "hintegral_norm_%d", m_num);
    hintegral->SetName(m_hname); hintegral->SetTitle(m_hname);
    (*hintegral).Write();
  }

}

void PMTDark::clean(){

  m_num = 0;
  m_nevents = 0;

  delete hamplitude;
  delete hintegral;
  delete htime;
}


//------------------------------------------------------------------------------


bool activePMT(int num, std::vector<int> array ){
  auto findIt = std::find( array.begin(), array.end(), num );
  if( findIt != array.end() ){ return true; }
  else{ return false; }
}


//------------------------------------------------------------------------------


bool passCuts( std::map<int, std::vector<Pulse>> map ){

  bool is_good=true; int ncounts=0;

  for( auto & pmt : map ){
    if( (pmt.second).size() > 2 ){ is_good=false; break;  }
    if( (pmt.second).size() > 1 ){ ncounts++; }
    if( ncounts > 3 ){ is_good=false; break; }
  }

  return is_good;

}

//------------------------------------------------------------------------------


bool hasValue( std::vector<size_t> array, size_t value ) {

  auto checkit = std::find( array.begin(), array.end(), value );

  bool found = checkit != array.end();

  return found;

};

// selection function for coincidence. I would recommend to use window 10 (bins in this case so 20 ns) and thresold 100, but you can select differnet values
bool hasCoincidence(double window, double threshold,  std::vector<Pulse> m_pulses ) {

    // This index hold the entries of m_pulse which are found in coincidence
    std::vector<size_t>  m_indx;

    for( size_t i=0; i<m_pulses.size(); i++ ) {

      for( size_t j=i+1; j<m_pulses.size(); j++ ) {

        //Do not let small noise peak spoil the selection
        if( m_pulses[i].amplitude < threshold || m_pulses[j].amplitude < threshold ){ continue; }


        // Check if there is a match inside the window
        if( abs(m_pulses[i].peaktime - m_pulses[j].peaktime) < window ){

          // Check if the index has been matched already, if yet continue
          if( !hasValue(m_indx, i) && !hasValue(m_indx, j) ) {

            m_indx.push_back(i);
            m_indx.push_back(j);

          }
        }

      } // end j
    } //end i


    // Super dumb selection: waveforms are skipped if there are more than two PMTs in coincidence,
    // you can put in your preferred logic:
    // m_indx holds all the index of the pulse array found in coincidence.
    // the struct Pulse also has pmt information, you can access it in a loop like this and then do further manipulations
    // Remember that even after the thresold cut there might be more than one pulse associated to the same pmt
    // for( int & idx : m_indx ){
    //  int pmtid = m_pulses[idx].pmt;
    // }

    bool isCoincidence = false;
    if( m_indx.size() > 0 ){ isCoincidence = true; }

    return isCoincidence;

};

//------------------------------------------------------------------------------


int main( int argc, char **argv ){

  if(argc != 5) {

      std::cout << "ERROR: " << argc << "/5 arguments" << std::endl;
      return -1;
  }


  // TODO: make it configurable!
  std::vector<int> activepmts = { 156, 157, 158, 159, 160, 161, 163, 164 };

  std::string input = argv[1];
  std::cout << ">>> Process input: " << input << std::endl;

  TChain tchain( "calo/pulsetree" );
  tchain.Add( input.c_str() );

  int nevents = atoi(argv[3]);
  if( nevents<0 || nevents>tchain.GetEntries() ){
    nevents = tchain.GetEntries();
  }

  std::vector<int> *m_pmt_array = 0;
  std::vector<float> *m_pulse_time = 0;
  std::vector<float> *m_pulse_amplitude = 0;
  std::vector<float> *m_pulse_integral = 0;

  tchain.SetBranchAddress("m_pmt_array", & m_pmt_array);
  tchain.SetBranchAddress("m_pulse_time", & m_pulse_time);
  tchain.SetBranchAddress("m_pulse_amplitude", & m_pulse_amplitude);
  tchain.SetBranchAddress("m_pulse_integral", & m_pulse_integral);


  // Event loop ................................................................

  int m_apply_cuts = atoi(argv[4]);
  if(m_apply_cuts != 0){ std::cout << "Using cuts" << std::endl; }

  if(nevents<0 || nevents>tchain.GetEntries()){ nevents = tchain.GetEntries(); }
  std::cout << ">>> Process number of events: " << nevents << std::endl;


  std::map<int, PMTDark*> m_pmts;
  for( int pmt : activepmts ){
    PMTDark *mypmt = new PMTDark(pmt); m_pmts[pmt] = mypmt;
  }

  int ngoodevent=0;
  for( int event=0; event<nevents; event++ ){

    tchain.GetEntry(event);

    // Group the pulses on each PMT
    std::map<int, std::vector<Pulse>> m_pmt_pulse;
    std::vector<Pulse> m_pulses;
    for( size_t pmtindex=0; pmtindex<m_pmt_array->size(); pmtindex++ ){

      int pmtid = (*m_pmt_array).at(pmtindex);

      if(activePMT( pmtid, activepmts ) ){
        Pulse pulse;
        pulse.peaktime = (*m_pulse_time).at(pmtindex);
        pulse.amplitude = (*m_pulse_amplitude).at(pmtindex);
        pulse.integral = (*m_pulse_integral).at(pmtindex);
        pulse.pmt = pmtid;
        m_pmt_pulse[pmtid].push_back( pulse );
        m_pulses.push_back( pulse );
      }

    }

    // If applicable skip bad events -- replace this line using the appropriate selection function if you want to use the coincidence selection
    // m_pmt_pulse is a map associating pulses and pmts; m_pulses is just an array with all the pulses
    if(  !passCuts( m_pmt_pulse ) && m_apply_cuts != 0 ){ continue; }


    // Save to files
    for( auto & pulses : m_pmt_pulse ){
      int pmtid = pulses.first;
      for( auto & pulse : pulses.second  )
        m_pmts[pmtid]->loadPulse( pulse );
    }

    ngoodevent++;

  }

  std::cout<< ">>> Good events: "<< ngoodevent << std::endl;


  // Write the output ..........................................................


  const std::string output = argv[2];
  std::cout << "Output name: " << output << std::endl;

  TFile ofile(output.c_str(), "RECREATE");

  for( auto & m_pmt : m_pmts ){

    PMTDark pmt = *m_pmt.second;

    pmt.loadEvents( ngoodevent );
    pmt.writeHist();

  }

  ofile.Close();

  std::cout << "All DONE!" << std::endl;

  return 0;
}
