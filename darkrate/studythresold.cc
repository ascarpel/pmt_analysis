////////////////////////////////////////////////////////////////////////////////
// this macro computes the darkrate for the active pmts as fuction of a
// thresold in amplitiude
//
// usage root -l studythresold.cc
// make sure to replace (if necessary):
// 1) the file with the output of makePulseHist.cc
// 2) the window with the effective waveform window ( if you apply a cut on it )
// 3) the list of active pmts
//
// ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////////////

#include "TH1D.h"
#include "TFile.h"

void studythresold()
{
  TFile *tfile = new TFile("myfile_darkrate.root", "READ");
  double window = 0.1; // 1 us

  std::vector<TGraphErrors*> graphs;
  std::vector<int> pmts = { 156, 157, 158, 159, 160, 161, 163, 164 };

  std::vector<double> rates;

  for( int pmt : pmts )
  {

    std::cout << "Process PMT: " << pmt << std::endl;

    TH1D *hist = (TH1D*)tfile->Get( ("hamplitude_norm_"+to_string(pmt)).c_str() );

    int nbins = hist->GetXaxis()->GetNbins();

    std::cout << "Integrated rate [kHz]: " << hist->Integral(17, nbins+1) / window << std::endl;

    rates.push_back( hist->Integral(17, nbins+1) / window );

    TGraphErrors *graph = new TGraphErrors( nbins );
    graph->SetName( ("graph_"+to_string(pmt)).c_str() );
    graph->SetTitle( ("PMT: "+to_string(pmt)).c_str() );

    for( int bin = 0; bin<nbins; bin++ )
    {
      double ncounts = hist->Integral(bin, nbins+1); // should include the overflow bin too

      double rate = ncounts /  window;
      double erate = TMath::Sqrt((double)ncounts) / window;
      graph->SetPoint(bin, (double)bin*0.122, rate);

    }

    graphs.push_back(graph);
    
  }


  TCanvas *c = new TCanvas("c", "c", 600, 500);


  c->SetGridx(true); c->SetGridy(true);

  graphs[0]->Draw("AL PLC");
  for(size_t i=1; i<graphs.size(); i++){ graphs[i]->Draw("L PLC"); }
  gPad->SetLogy();

  TLegend *legend = new TLegend(0.58, 0.58, 0.92, 0.92);
  legend->SetHeader("#bf{Rate @ thr: > 2 mV}", "C");

  for(size_t i=0; i<graphs.size(); i++){

    graphs[i]->SetLineWidth(2);
    graphs[i]->GetXaxis()->SetRangeUser(1.464, 12);
    graphs[i]->GetYaxis()->SetRangeUser(0.01, 50);
    graphs[i]->GetYaxis()->SetTitle("Rate [kHz]");
    graphs[i]->GetXaxis()->SetTitle("Threshold [mV]");

    char label[100];
    sprintf(label, "PMT %d: %.02f [kHz]", pmts[i], rates[i] );

    legend->AddEntry(graphs[i]->GetName(), label, "L");

  }

  legend->Draw("");
  //c->Update();

}
