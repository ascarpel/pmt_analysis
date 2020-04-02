////////////////////////////////////////////////////////////////////////////////
// RUN Methods implementation
//
// mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////////////

#include "Utils.h"
#include "Run.h"

RUN::RUN()
{
  //this->loadDatabase();
};

//------------------------------------------------------------------------------

RUN::RUN(int run, int subrun)
  : m_run(run)
  , m_subrun(subrun)
{
  this->loadDatabase();
};

//------------------------------------------------------------------------------

RUN::RUN(string filename)
  : m_run_file(filename)
{
  m_run=0; m_subrun=0;
  this->getRundata(m_run_file, m_run, m_subrun);

  this->loadDatabase();
};

//------------------------------------------------------------------------------

RUN::~RUN()
{ };

//------------------------------------------------------------------------------


void RUN::loadDatabase()
{
  // Load the database with the association between run-subrun and the
  // optical channel


  utils::CSVReader reader("../dbase/Calibration_20200204.csv");
  m_dataList = reader.getData();


  return;
}

//------------------------------------------------------------------------------

int RUN::getOpticalChannel()
{
  // Produce the association between run-subrun and optical channel

  if( m_dataList.size() == 0)
  {
    cout << "ERROR! load a valid CSV file!" << endl;
  }

  // Explanation of the entries:
  // vec[0] = run number
  // vec[1] = real optical channel (the value set in the switch controller)
  // vec[2] = file start
  // vec[3] = file end
  // vec[4] = optical channel list ( a reconversion of the real optical channel)
  //          value into a list of numbers from 0-17, identical on both east and
  //          west side. This is correct since I associate the optical channel
  //          to the absolute pmt number and we have only 180 pmts per file.
  //          NB: in future, if data are taken from all boards, the list must
  //          be from 0-35 for all the 36 possible switcher configurations.

  // Skip the first line since it has just the titles of the columns
  const int nskip = 1;
  for(size_t i=nskip; i<m_dataList.size(); i++)
  {
    vector<string> vec = m_dataList.at(i);

    int run = stoi(vec[0]);
    int real_opt_ch = stoi(vec[1]);
    int file_start = stoi(vec[2]);
    int file_end = stoi(vec[3]);
    int opt_ch_list =  stoi(vec[4]);

    if(run > m_run){ break; }
    else if( m_run == run )
    {
      if( file_start <= m_subrun & file_end >= m_subrun )
      {
        m_optical_channel = opt_ch_list;
      }
    }
  }

  return m_optical_channel;
};


//------------------------------------------------------------------------------

void RUN::getRundata(string filename, int & run, int & subrun)
{
  // NB: It assumes filename with structure *runRUN_SUBRUN_*.root
  int suffix[3];
  int end = filename.find(".root");

  suffix[0] = filename.find("run")+3;
  suffix[1] = suffix[0]+filename.substr(suffix[0], end-suffix[0]).find("_");
  suffix[2] = suffix[1]+filename.substr(suffix[1]+1, end-suffix[1]).find("_");

  run = std::stoi(filename.substr(suffix[0], suffix[1]-suffix[0]));

  // In some files the subrun might not be in the middle of the word, but at the
  // end, before the extension .root
  string tmp_str = filename.substr(suffix[1]+1, suffix[2]-suffix[1]);
  if( tmp_str.find(".root") < tmp_str.size() )
  {
    int pos = tmp_str.find(".root");
    subrun = std::stoi(filename.substr(suffix[1]+1, tmp_str.size()-pos-2));
  }
  else
  {
    subrun = std::stoi(filename.substr(suffix[1]+1, suffix[2]-suffix[1]));
  }

  return ;
};

//------------------------------------------------------------------------------
