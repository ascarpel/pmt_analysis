////////////////////////////////////////////////////////////////////////////////
// Macro that interfaces the output of the CAEN decoder script and produces
// a file format compatible with this small analysis framework
//
// mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////////////

#include "Waveform.h"
#include "Pmt.h"
#include "Utils.h"

#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"
#include "TProfile.h"

#include <Eigen/Core>
#include <unsupported/Eigen/FFT>


bool hasValue( std::vector<size_t> array, size_t value ) {

  auto checkit = std::find( array.begin(), array.end(), value );

  bool found = checkit != array.end();

  return found;

};


void findCoincidencePulses(double window, std::vector<Waveform::Pulse> m_pulses, std::vector<size_t>  & m_indx ) {


    for( size_t i=0; i<m_pulses.size(); i++ ) {

      for( size_t j=i+1; j<m_pulses.size(); j++ ) {


        // Check if there is a match
        if( abs(m_pulses[i].time_peak - m_pulses[j].time_peak) < window ){

          // Check if the index has been matched already, if yet continue
          if( !hasValue(m_indx, i) && !hasValue(m_indx, j) ) {

            m_indx.push_back(i);
            m_indx.push_back(j);

          }
        }

      } // end j
    } //end i



};



void findCoincidences(){

  int startevent=0;
  int nevents = -1;

  double window = 3; // Pulse coincidence window

  double startTh = 300; // ADC
  double endTh = 200; //ADC

  // ROI where the response of the LAr Scintillation should be contained
  // startROI is the part before the 0 ( 0 is set on the peak )
  // endROI is the end after the 0 ( large enough to contain the 1.6 us slow component )
  int startROI = 250; // in bins
  int endROI = 1750; // in bin


  std::string filePattern = "../data/coldTest/run1882/cold-commissioning-run1882.root";
  std::string outfilename = "./sample/coincidence_waveforms_run1882.root";
  std::vector<int> activeBoards = { 1,2,3,4,5,6,7,8,9,10,11 };
  std::vector<int> activeChannels = { 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14, 15 };


  //----------------------------------------------------------------------------
  int nchannels=16;

  TChain* tchain = new TChain("caenv1730dump/events");
  std::vector<std::vector<uint16_t> > *data=0; //unsigned short
  tchain->SetBranchAddress("fWvfmsVec", &data);

  tchain->Add( filePattern.c_str() );

  std::cout << "Total number of entries: " << tchain->GetEntries() << std::endl;

  if(nevents < 0 || nevents > tchain->GetEntries()){ nevents = tchain->GetEntries(); }
  std::cout << "Total number of events: " << nevents << std::endl;


  TFile *outfile = new TFile(outfilename.c_str(), "RECREATE");


  //----------------------------------------------------------------------------


  // Book the waveAna object to process the wavefrom Analysis
  Waveform *waveAna = new Waveform();


  for( int event=startevent; event<startevent+nevents; event++  ) {

    tchain->GetEntry(event);

    std::vector<Waveform::Pulse> m_pulses;
    std::vector<int> m_daqid;
    std::vector<Waveform::Waveform_t> m_waves;
    std::vector<TH1D*> m_hists;



    // We do a first loop to find the pulses
    for( int & board : activeBoards ) {
      for( int & channel : activeChannels ){

        int daqid = channel+nchannels*board;

        std::vector<uint16_t> dataArray = (*data).at(daqid);
        if( dataArray.size() == 0 ){ continue; }

        //As always load the waveform and subtract baseline
        waveAna->loadData(dataArray);

        //Check pulses and keep waveforms with only one pulse
        auto pulses = waveAna->findPulses(startTh, endTh);

        // Want to study only events max 1 pulses in the 10 us window
        // this would reduce pileup hopefully
        if( pulses.size() > 0 && pulses.size() < 3) {

          for( auto & pulse : pulses ){

            m_pulses.push_back( pulse );
            m_daqid.push_back( daqid );
            m_waves.push_back( waveAna->getWaveform() );

            TH1D *tmpwfm = waveAna->getWaveformHist();
            char hname[100]; sprintf(hname, "h_e%d_id%d", event, daqid);
            tmpwfm->SetName(hname); tmpwfm->SetTitle(hname);

            m_hists.push_back( tmpwfm );
          }

        }

        // Clean waveAna for the next iteration
        waveAna->clean();

      }
    }


    if( m_pulses.size() ==0 ){ continue; }

    // Book an array to hold the location of the pulses
    std::vector< size_t > m_indx;

    findCoincidencePulses( window, m_pulses, m_indx );

    if( m_indx.size() > 0 ) {

      std::cout << "Event: " << event;
      std::cout << " coincidences on ch: ";

      for( size_t & idx : m_indx ) {

         int daqid = m_daqid[idx];
         std::cout <<  daqid << " ";

         // Here we process the waveform to select the best part of the
         // scintillation signal around the peak found in coincidence

         Waveform::Pulse m_pulse = m_pulses[idx];
         Waveform::Waveform_t m_wave = m_waves[idx];


         char name[100]; sprintf(name, "roi_e%d_daq%d_idx%zu", event, daqid, idx);
         TH1D* m_roi = new TH1D(name, name, endROI+startROI, -startROI*2.0, endROI*2.0);

         int maxBin = m_pulse.time_peak;
         double maxAmpl = m_pulse.amplitude;

         // Select pulses only in a good range
         if(maxAmpl > 14000 || maxAmpl < 100 ){ continue; }

         // Jump if the peak happens too close to the edges of the waveform
         if( ((maxBin-startROI) < 0) || ((maxBin+endROI) > m_wave.size()) ) {
           continue;
         }

         int counts=0;
         for( size_t bin=maxBin-startROI; bin<maxBin+endROI; bin++, counts++) {

           double time = (counts-startROI)*2.0-1.0;
           m_roi->Fill(time, m_wave[bin]);
           m_roi->SetBinError(counts, 1.0);

         }


         m_roi->Write();
         //m_hists[idx]->Write();

      }


      std::cout << std::endl;

    }

  } // end loop over events


  outfile->Close();

  std::cout << "All done" << std::endl;

}
