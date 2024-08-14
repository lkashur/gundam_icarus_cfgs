#ifndef GUNDAM_STAT_COVARIANCE_H
#define GUNDAM_STAT_COVARIANCE_H

#include "JointProbabilityBase.h"
#include "TMatrixTSym.h"


namespace JointProbability{

  class StatCovariance : public JointProbabilityBase {

  public:
    [[nodiscard]] std::string getType() const override { return "StatCovariance"; }

    void initializeCache(const std::vector<Sample>& vec_samples){

    }

    // For StatCovariance
    [[nodiscard]] virtual double eval(const std::vector<Sample>& vec_samples) const override {

      // Binning

      int NTotalBins = 0;
      int NBinsForEachSample[vec_samples.size()];
      for(unsigned int i_sample=0; i_sample<vec_samples.size(); i_sample++){
        NBinsForEachSample[i_sample] = vec_samples[i_sample].getDataContainer().getHistogram().binList.size();
        NTotalBins += NBinsForEachSample[i_sample];
      }

      int SampleIndicesForEachBin[NTotalBins];
      int LocalBinIndicesForEachBin[NTotalBins];
      int binIndex = 0;
      for(unsigned int i_sample=0; i_sample<vec_samples.size(); i_sample++){
        for(int i = 0; i < NBinsForEachSample[i_sample]; ++i) {
          SampleIndicesForEachBin[binIndex] = i_sample;
          LocalBinIndicesForEachBin[binIndex] = i;
          binIndex++;
        }
      }

      TMatrixTSym<double> Cov_Data(NTotalBins);
      Cov_Data.Zero(); // Make sure it is initialized to zero

      TMatrixTSym<double> Cov_MC(NTotalBins);
      Cov_MC.Zero(); // Make sure it is initialized to zero

      double DataArr[NTotalBins];
      double MCArr[NTotalBins];

      for(unsigned int idx_global_i=0; idx_global_i<NTotalBins; idx_global_i++){

        int idx_sample_i = SampleIndicesForEachBin[idx_global_i];
        const auto& sample_i = vec_samples[idx_sample_i];
        int idx_local_i = LocalBinIndicesForEachBin[idx_global_i];

        std::vector<Event*> vec_DataEvtList_i = sample_i.getDataContainer().getHistogram().binList[idx_local_i].eventPtrList;
        std::vector<Event*> vec_MCEvtList_i = sample_i.getMcContainer().getHistogram().binList[idx_local_i].eventPtrList;

        DataArr[idx_global_i] = sample_i.getDataContainer().getHistogram().binList[idx_local_i].content;
        MCArr[idx_global_i] = sample_i.getMcContainer().getHistogram().binList[idx_local_i].content;

        for(unsigned int idx_global_j=idx_global_i; idx_global_j<NTotalBins; idx_global_j++){

          const auto& sample_j = vec_samples[ SampleIndicesForEachBin[idx_global_j] ];
          int idx_local_j = LocalBinIndicesForEachBin[idx_global_j];

          std::vector<Event*> vec_DataEvtList_j = sample_j.getDataContainer().getHistogram().binList[idx_local_j].eventPtrList;
          std::vector<Event*> vec_MCEvtList_j = sample_j.getMcContainer().getHistogram().binList[idx_local_j].eventPtrList;

          // Data
          double DataBinContent = 0.;
          for(Event* EvtList_i: vec_DataEvtList_i){
            EventUtils::Indices& EvtIndices_i = EvtList_i->getIndices();
            for(Event* EvtList_j: vec_DataEvtList_j){
              EventUtils::Indices& EvtIndices_j = EvtList_j->getIndices();
              if(EvtIndices_i.entry==EvtIndices_j.entry){
                DataBinContent += EvtList_i->getEventWeight();
              }
            }
          }
          Cov_Data(idx_global_i, idx_global_j) = DataBinContent;
          Cov_Data(idx_global_j, idx_global_i) = DataBinContent; 

          // MC
          double MCBinContent = 0.;
          for(Event* EvtList_i: vec_MCEvtList_i){
            EventUtils::Indices& EvtIndices_i = EvtList_i->getIndices();
            for(Event* EvtList_j: vec_MCEvtList_j){
              EventUtils::Indices& EvtIndices_j = EvtList_j->getIndices();
              if(EvtIndices_i.entry==EvtIndices_j.entry){
                // We need to take the stat of the raw-MC events, not on the weighted sum;
                // so it is squared
                MCBinContent += EvtList_i->getEventWeight() * EvtList_i->getEventWeight();
              }
            }
          }
          Cov_MC(idx_global_i, idx_global_j) = MCBinContent;
          Cov_MC(idx_global_j, idx_global_i) = MCBinContent;


        }

      } // done constructin covariance

      // Now calculate chi2
      // inverse the chi2
      TMatrixTSym<double> Cov_Sum = Cov_Data+Cov_MC;
      Cov_Sum.Invert(); // inverted

      double chi2 = 0.;
      for(unsigned int i=0; i<NTotalBins; i++){
        for(unsigned int j=0; j<NTotalBins; j++){
          double this_chi2 = ( DataArr[i] - MCArr[i] ) * Cov_Sum(i, j) * ( DataArr[j] - MCArr[j] );
          chi2 += this_chi2;
        }
      }

      return chi2;

    } // END eval(std::vector<Sample> vec_samples)

  };

}

#endif //GUNDAM_STAT_COVARIANCE_H
