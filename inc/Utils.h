#include <iostream>
#include <stdio.h>


namespace utils
{

////////////////////////////////////////////////////////////////////////////////

  void get_rundata( string filename, int & run, int & subrun )
  {
	  // NB: assumes file with structure run*_***.root
 	  int suffix[3];
    suffix[0] = filename.find("run")+3;
    suffix[1] = filename.find("_");
    suffix[2] = filename.find(".root");

    run = stoi(filename.substr(suffix[0], suffix[1]-suffix[0]));
    subrun = stoi(filename.substr(suffix[1]+1, suffix[2]-suffix[1]-1));

    return ;
  }

////////////////////////////////////////////////////////////////////////////////

} //end namespace
