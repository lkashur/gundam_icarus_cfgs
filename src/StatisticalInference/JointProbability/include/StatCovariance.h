#ifndef GUNDAM_STAT_COVARIANCE_H
#define GUNDAM_STAT_COVARIANCE_H

#include "JointProbabilityBase.h"
#include "TMatrixTSym.h"
#include <TDecompSVD.h>
#include "TFile.h"
#include "Logger.h"

namespace JointProbability{

  class StatCovariance : public JointProbabilityBase {

  public:
    [[nodiscard]] std::string getType() const override { return "StatCovariance"; }

    bool _isInitialized{false};
    int _nTotalBins;
    std::vector<int> _samplepairIndicesForEachBin;
    std::vector<int> _localBinIndicesForEachBin;

    std::vector< std::vector< std::vector<Event*> > > _arr_DataPtrs;
    std::vector< std::vector< std::vector<Event*> > > _arr_MCPtrs;

    TMatrixTSym<double> Cov_Data_Nominal;
    TMatrixTSym<double> Cov_MC_Nominal;

    void fillEventPtrs(const std::vector<SamplePair>& vec_samplepairs){

      // Binning
      int NBinsForEachSamplePair[vec_samplepairs.size()];
      for(unsigned int i_samplepair=0; i_samplepair<vec_samplepairs.size(); i_samplepair++){
        NBinsForEachSamplePair[i_samplepair] = vec_samplepairs[i_samplepair].data->getHistogram().getBinContentList().size();
        _nTotalBins += NBinsForEachSamplePair[i_samplepair];
      }

      for(unsigned int i_samplepair=0; i_samplepair<vec_samplepairs.size(); i_samplepair++){
        for(int i = 0; i < NBinsForEachSamplePair[i_samplepair]; ++i) {
          _samplepairIndicesForEachBin.push_back( i_samplepair );
          _localBinIndicesForEachBin.push_back( i );
        }
      }

      for (int i = 0; i < _nTotalBins; ++i) {
        _arr_DataPtrs.resize(_nTotalBins);
        _arr_DataPtrs[i].resize(_nTotalBins);
        _arr_MCPtrs.resize(_nTotalBins);
        _arr_MCPtrs[i].resize(_nTotalBins);
      }

      for(unsigned int idx_global_i=0; idx_global_i<_nTotalBins; idx_global_i++){

        int idx_samplepair_i = _samplepairIndicesForEachBin[idx_global_i];
        const auto& samplepair_i = vec_samplepairs[idx_samplepair_i];
        int idx_local_i = _localBinIndicesForEachBin[idx_global_i];

        std::vector<Event*> vec_DataEvtList_i = samplepair_i.data->getHistogram().getBinContextList()[idx_local_i].eventPtrList;
        std::vector<Event*> vec_MCEvtList_i = samplepair_i.model->getHistogram().getBinContextList()[idx_local_i].eventPtrList;

        for(unsigned int idx_global_j=idx_global_i; idx_global_j<_nTotalBins; idx_global_j++){

          const auto& samplepair_j = vec_samplepairs[ _samplepairIndicesForEachBin[idx_global_j] ];
          int idx_local_j = _localBinIndicesForEachBin[idx_global_j];

          std::vector<Event*> vec_DataEvtList_j = samplepair_j.data->getHistogram().getBinContextList()[idx_local_j].eventPtrList;
          std::vector<Event*> vec_MCEvtList_j = samplepair_j.model->getHistogram().getBinContextList()[idx_local_j].eventPtrList;

          // Data
          for(Event* EvtList_i: vec_DataEvtList_i){
            EventUtils::Indices& EvtIndices_i = EvtList_i->getIndices();
            for(Event* EvtList_j: vec_DataEvtList_j){
              EventUtils::Indices& EvtIndices_j = EvtList_j->getIndices();
              if(EvtIndices_i.entry==EvtIndices_j.entry){
                _arr_DataPtrs[idx_global_i][idx_global_j].push_back( EvtList_i );
              }
            }
          }

          // MC
          for(Event* EvtList_i: vec_MCEvtList_i){
            EventUtils::Indices& EvtIndices_i = EvtList_i->getIndices();
            for(Event* EvtList_j: vec_MCEvtList_j){
              EventUtils::Indices& EvtIndices_j = EvtList_j->getIndices();
              if(EvtIndices_i.entry==EvtIndices_j.entry){
                _arr_MCPtrs[idx_global_i][idx_global_j].push_back( EvtList_i );
              }
            }
          }

        }
      }

      printf("[JSKIMDEBUG] Evaluating nominal stat covarinace matrices\n");
      Cov_Data_Nominal.ResizeTo(_nTotalBins, _nTotalBins);
      Cov_Data_Nominal.Zero(); // Make sure it is initialized to zero
      Cov_MC_Nominal.ResizeTo(_nTotalBins, _nTotalBins);
      Cov_MC_Nominal.Zero(); // Make sure it is initialized to zero

      for(unsigned int idx_global_i=0; idx_global_i<_nTotalBins; idx_global_i++){

        int idx_samplepair_i = _samplepairIndicesForEachBin[idx_global_i];
        const auto& samplepair_i = vec_samplepairs[idx_samplepair_i];
        int idx_local_i = _localBinIndicesForEachBin[idx_global_i];

        for(unsigned int idx_global_j=idx_global_i; idx_global_j<_nTotalBins; idx_global_j++){
          // Data
          double DataBinContent = 0.;
          for(Event* EvtList_i: _arr_DataPtrs[idx_global_i][idx_global_j]){
            DataBinContent += EvtList_i->getEventWeight();
          }
          Cov_Data_Nominal(idx_global_i, idx_global_j) = DataBinContent;
          Cov_Data_Nominal(idx_global_j, idx_global_i) = DataBinContent;
          // MC
          double MCBinContent = 0.;
          for(Event* EvtList_i: _arr_MCPtrs[idx_global_i][idx_global_j]){
            // We need to take the stat of the raw-MC events, not on the weighted sum;
            // so it is squared
            MCBinContent += EvtList_i->getEventWeight() * EvtList_i->getEventWeight();
          }
          Cov_MC_Nominal(idx_global_i, idx_global_j) = MCBinContent;
          Cov_MC_Nominal(idx_global_j, idx_global_i) = MCBinContent;
        }
      } // done constructin covariance

      _isInitialized = true;

    }

    // For StatCovariance
    [[nodiscard]] virtual double eval(const std::vector<SamplePair>& vec_samplepairs) const override {

      bool DoDebug = false;

      LogThrowIf(!_isInitialized, "StatCovariance is not initialized");

      double DataArr[_nTotalBins];
      double MCArr[_nTotalBins];

      for(unsigned int idx_global_i=0; idx_global_i<_nTotalBins; idx_global_i++){

        int idx_samplepair_i = _samplepairIndicesForEachBin[idx_global_i];
        const auto& samplepair_i = vec_samplepairs[idx_samplepair_i];
        int idx_local_i = _localBinIndicesForEachBin[idx_global_i];

        DataArr[idx_global_i] = samplepair_i.data->getHistogram().getBinContentList()[idx_local_i].sumWeights;
        MCArr[idx_global_i] = samplepair_i.model->getHistogram().getBinContentList()[idx_local_i].sumWeights;

      }

      TMatrixTSym<double> Cov_Data(Cov_Data_Nominal);
      TMatrixTSym<double> Cov_MC(Cov_MC_Nominal);

      // Now calculate chi2
      // inverse the chi2
      TMatrixT<double> Cov_Sum = Cov_Data+Cov_MC;
      //TMatrixT<double> Cov_Sum = Cov_Data;

      if(DoDebug){
        std::cout << "[JSKIMDEBUG] Cov_Data:" << std::endl;
        Cov_Data.Print();
        std::cout << "[JSKIMDEBUG] Cov_MC:" << std::endl;
        Cov_MC.Print();
        std::cout << "[JSKIMDEBUG] Cov_Sum:" << std::endl;
        Cov_Sum.Print();
      }

      double det = Cov_Sum.Determinant();
      //if( std::abs(det) >= 1e-10 ){
      if( false ){ // force SVD
        Cov_Sum.Invert();
      }
      else{

        // Regular inversion failed, compute pseudo-inverse

        TDecompSVD svd(Cov_Sum, 1e-10);
        LogThrowIf(!svd.Decompose(), "SVD also failed; TODO just use poisson in this case?");

        // Get U, S, V matrices from SVD
        TMatrixD U = svd.GetU();
        TVectorD S = svd.GetSig();
        TMatrixD V = svd.GetV();

        // Compute pseudo-inverse of S
        TMatrixD S_inv(S.GetNrows(), S.GetNrows());
        S_inv.Zero();

        double tolerance = 1e-10; // Tolerance for singular values considered as zero
        for (int i = 0; i < S.GetNrows(); ++i) {
          if (std::abs(S[i]) > tolerance) {
            S_inv[i][i] = 1.0 / S[i];
          }
        }

        // Pseudo-inverse of the original matrix: V * S_inv * U^T
        //TMatrixD pseudoInverse = V * S_inv * U.T();

        // Assign the pseudo-inverse to the original matrix
        Cov_Sum = V * S_inv * U.T();;

      }

      if(DoDebug){
        Cov_Sum.Print();
      }

      double chi2 = 0.;
      for(unsigned int i=0; i<_nTotalBins; i++){
        for(unsigned int j=0; j<_nTotalBins; j++){
          double this_chi2 = ( DataArr[i] - MCArr[i] ) * Cov_Sum(i, j) * ( DataArr[j] - MCArr[j] );
          chi2 += this_chi2;
        }
      }

      if(chi2<0){
        printf("[JSKIMDEBUG] Negative chi2: %e\n", chi2);
        printf("[JSKIMDEBUG] Data-MC array:\n");
        for(unsigned int i=0; i<_nTotalBins; i++){
          printf("%d\t%f\t%f\n", i, DataArr[i], MCArr[i]);
        }

        printf("[JSKIMDEBUG] Cov_Data\n");
        Cov_Data.Print();
        printf("[JSKIMDEBUG] Cov_MC\n");
        Cov_MC.Print();

        printf("[JSKIMDEUBG] chi2 calculation\n");
        double tmp_chi2 = 0.;
        for(unsigned int i=0; i<_nTotalBins; i++){
          for(unsigned int j=0; j<_nTotalBins; j++){
            double this_chi2 = ( DataArr[i] - MCArr[i] ) * Cov_Sum(i, j) * ( DataArr[j] - MCArr[j] );
            tmp_chi2 += this_chi2;
            printf("(%d, %d): Diff_i = %f, Diff_j = %f, CovInv = %f -> this_chi2 = %f; current chi2 = %f\n", i, j, DataArr[i] - MCArr[i], DataArr[j] - MCArr[j], Cov_Sum(i, j), this_chi2, tmp_chi2);
          }
        }

        Cov_Sum.Print();

        // TODO REMOVE TFile.h
        TFile *f_debug = new TFile("f_debug.root", "RECREATE");
        f_debug->cd();
        Cov_Data.Write("Cov_Data");
        Cov_MC.Write("Cov_MC");
        Cov_Sum.Write("PInvCovSum");
        f_debug->Close();

        abort();
      }

/*
      double chi2 = 0.;
      for(unsigned int i=0; i<_nTotalBins; i++){
        //double this_chi2 = ( DataArr[i] - MCArr[i] ) * 1./Cov_Data(i, i) * ( DataArr[i] - MCArr[i] );
        double this_chi2 = ( DataArr[i] - MCArr[i] ) * Cov_Sum(i, i) * ( DataArr[i] - MCArr[i] );
        std::cout << i << "\t" << Cov_Sum(i, i) << std::endl;
        chi2 += this_chi2;
      }
*/

      if(DoDebug){
        std::cout << "[JSKIMDEBUG] chi2 = " << chi2 << std::endl;
      }
      return chi2;

    }

  };

}

#endif //GUNDAM_STAT_COVARIANCE_H
