////////////////////////////////////////////////////////////////////////////////
// Macro to analyze the PMT signal
// Process the pulse analyis for the filename given as argument. Edit the input
// filename with the file you prefer. Note: the output will always be in the format
// "pmt_run%d_00%d.root". you must change manually the variables run and subrun for now.
//
// Usage: root -l loadLib.cc load_pmt_tchain.cc
// Out of the box, it will produce a root file in the macro folder called "pmt_run1067_1.root"
// with the ADC distribution of all the pmts in
// after baseline subtraction.
//
// mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////////////

#include "Geometry.h"
#include "Run.h"
#include "Waveform.h"
#include "Pmt.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"

#include <fstream>
#include <sstream>

using namespace std;

void load_pmt_tchain()
{

  // CHANGE THE FILENAME BELOW IF WANT TO USE ONE OTHER FILE FOR YOUR ANALYSIS
  string filename="/media/disk_a/ICARUS/PMT_gainData_20200205_decoded/data_dl1_run1067_1_20200204T225822_dl3_decoded.root";

  // THE RUN CLASS SHOULD ASSIGN THE RIGTH RUN AND SUBRUN NUMBERS ACCORDING TO
  // THE FILENAME GIVEN
  RUN my_run(filename);
  int run=my_run.getRun();
  int subrun=my_run.getSubrun();

  // We get the optical channel and find if it is straighforwardly invalid
  int optical_channel = my_run.getOpticalChannel();
  if(optical_channel < 0)
  {
    cout<<"ERROR! Check if the optical channel is correctly identified!"<< endl;
  }

  // Print out what the macro is doing so far
  cout << "Processing run: " << run << " and subrun: " << subrun << endl;
  cout << "Optical channel is: " << optical_channel << endl;

  //****************************************************************************
  // Input
  //

  TChain *tchain =  new TChain("caenv1730dump/events");

  // HERE you filename is loaded inside the macro to be processed!
  tchain->Add(filename.c_str());

  std::vector<std::vector<uint16_t> > *data=0; //unsigned short
  tchain->SetBranchAddress("fWvfmsVec", &data);


  //****************************************************************************
  // Here we create and initialize a vector with the class PMT. If you want you
  // can define also other objects (like TH1D classes) in here. See macro
  // noise_ana.cc with an example
  //

  PMT* pmts[geo::nboards][geo::nchannels];
  for( int board=0; board<geo::nboards; board++ )
  {
    for( int channel=0; channel<geo::nchannels; channel++ )
    {
      pmts[board][channel] = new PMT( run, board, channel );
    }
  }


  //****************************************************************************
  // Now it is time to loop over the events
  //

  int nentries = tchain->GetEntries();
  for(int e=0; e<nentries; e++)
  {
    cout << "Processing event: " << e << endl;

    // WE TAKE THE EVENT
    tchain->GetEntry(e);

    for(int board=0; board<geo::nboards; board++)
    {
      for(int channel=0; channel<geo::nchannels; channel++)
      {

        // THIS CONDITION EXCLUDES PMTS WHICH ARE DIRECTLY ILLUMINATED BY
        // THE LASER
        if( pmts[board][channel]->isIlluminated(optical_channel) ){ continue; }

        // Create the WAVEFORM OBJECT
        Waveform *waveform = new Waveform(run, subrun ,e, board, channel);


        // ONCE THE DATA ARE LOADED IN THE WAVEFORM OBJECT, BASELINE IS
        // AUTOMATICALLY SUBSTRACTED. YOU CAN USE THE FUNCTION "getWaveform()"
        // TO EXTRACT THE WAVEFORM VECTOR AFTER BASELINE SUBTRACTION OR
        // "getRawWaveform()" BEFORE BASELINE SUBTRACTION.
        waveform->loadData((*data).at(channel+geo::nchannels*board));

        // THIS PART IS SOMETHING THAT YOU PROBABLY HAVE TO MODIFY:
        // WHAT THE CODE DOES HERE IS CHECKING IF THE WAFEFORM SATISFY MY CRITERIA
        // TO FIND A PULSE. IF YES, THIS WAVEFORM IS SAVED INSIDE THE PMT OBJECT
        //if(waveform->hasPulse(1.0))
        //{
         pmts[board][channel]->loadWaveform(waveform);
        //}

      } // end channel
    } // end boards

  } // end events


  //****************************************************************************
  // Output: create an output TFile and write the histogram in it
  //

  char ofilename[100]; sprintf(ofilename, "pmt_run%d_%d.root", run, subrun);
  TFile ofile(ofilename, "RECREATE"); ofile.cd();

  for( int board=0; board<geo::nboards; board++ )
  {
    for( int channel=0; channel<geo::nchannels; channel++ )
    {
      pmts[board][channel]->writeHist();
    }
  }

  ofile.Close();

  cout << "All done" << endl;

  //****************************************************************************
  //
  // For now you want to focus mostly on the waveform class and change the logic
  // behind the pulse finder.
  //
  //****************************************************************************

} //end macro
