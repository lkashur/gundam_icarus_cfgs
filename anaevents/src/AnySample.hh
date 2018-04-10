//////////////////////////////////////////////////////////
//
//  A class for event samples for for any analysis
//
//
//
//  Created: Nov 17 2015
//
//////////////////////////////////////////////////////////
#ifndef __AnySample_hh__
#define __AnySample_hh__

#include <string>
#include <vector>

#include <TH2D.h>
#include <TDirectory.h>
#include <TRandom3.h>
#include <TTree.h>

#include "AnaEvent.hh"
#include "AnaSample.hh"

///////////////////////////////////////
// Class definition
///////////////////////////////////////
class AnySample : public AnaSample
{
public:
  AnySample(int sample_id, std::string name,
       std::vector<std::pair <double,double> > v_d1edges,
       std::vector<std::pair <double,double> > v_d2edges, TTree *data, bool isBuffer, bool isEmpty=false, bool isIngrid=false);
  ~AnySample();

  //binning for various histograms
  void SetD1Binning(int nbins, double *bins);
  void SetD2Binning(int nbins, double *bins);
  void SetEnuBinning(int nbins, double *bins);
  void MakeHistos(); //must be called after binning is changed

  //histogram for event distributions
  void SetData(TObject *hdata);
  void FillEventHisto(int datatype);
  double CalcChi2();

  TH1D* GetPredHisto(){ return m_hpred; }
  TH1D* GetDataHisto(){ return m_hdata; }
  TH1D* GetMCHisto(){ return m_hmc; }
  TH1D* GetMCTruthHisto(){ return m_hmc_true; }

  void GetSampleBreakdown(TDirectory *dirout, std::string tag, bool save);
  TH1D* GetSignalHisto(){ return m_sig; }
  void Write(TDirectory *dirout, const char *bsname, int fititer);
  int GetSampleId(){ return sample_id; }

private:
  int sample_id;
  TH1D *m_hmc_true;
  TH1D *m_hmc;
  TH1D *m_hpred; //n(pRec_mu, thetaRec_mu)
  TH1D *m_hdata;
  TTree * m_data_tree;
  TH1D *m_sig;
  int nbins_D1, nbins_D2, nbins_enu, nAnybins, nbinsD1_toPlot;
  double *bins_D1, *bins_D2, *bins_enu, *bins_Any, *bins_D1toPlot;
  std::vector<std::pair<double, double> > m_D1edges;
  std::vector<std::pair<double, double> > m_D2edges;
  bool m_empty; // If true, we won't include any events in this sample (useful for testing the effect of removing samples)
  bool m_BufferBin; // Should we bother plotting the last bin (dat, dphit), or is it just a buffer (dpt)?
};

#endif
