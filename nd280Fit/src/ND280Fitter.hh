#ifndef __ND280FITTER_HH__
#define __ND280FITTER_HH__

#include <algorithm>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <TFile.h>
#include <TGraph.h>
#include <TMath.h>
#include <TMatrixT.h>
#include <TMatrixTSym.h>
#include <TRandom3.h>
#include <TVectorT.h>
#include <future>

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"

#include "AnaFitParameters.hh"
#include "AnaSample.hh"
#include "ColorOutput.hh"
#include "OptParser.hh"


class ND280Fitter {

public:

    ND280Fitter();
    ~ND280Fitter();

    void Reset();

    // Pre-initialization Methods
    void SetOutputTDirectory(TDirectory* output_tdirectory_);
    void SetPrngSeed(int PRNG_seed_);
    void SetMinimizationSettings(const MinSettings& minimization_settings_);
    void SetDisableSystFit(bool disable_syst_fit_);
    void SetSaveEventTree(bool save_event_tree_);
    void SetSaveFitParams(bool save_fit_params_);
    void SetApplyStatisticalFluctuationsOnSamples(bool apply_statistical_fluctuations_on_samples_);
    void SetSaveFitParamsFrequency(int save_fit_params_frequency_);
    void SetAnaFitParametersList(std::vector<AnaFitParameters*> AnaFitParameters_list_);
    void SetAnaSamplesList(std::vector<AnaSample*> AnaSample_list_);
    void SetSelectedDataType(int selected_data_type_);

    void Initialize();

    // Post-initialization Methods
    void WritePrefitData();
    void MakeOneSigmaChecks();
    bool Fit();
    void WritePostFitData();
    void ScanParameters(unsigned int nbSteps_);

    void WriteSamplePlots(TDirectory* outputDir_);

    double EvalFit(const double* par);
    void UpdateFitParameterValues(const double* par);
    void PropagateSystematics();
    void ComputeChi2();

    void InitFitter(std::vector<AnaFitParameters*>& fitpara);
    void FixParameter(const std::string& par_name, const double& value);
    void ParameterScans(const std::vector<int>& param_list, unsigned int nsteps);

    TMatrixD* GeneratePriorCovarianceMatrix();
    TMatrixD* GeneratePosteriorCovarianceMatrix();


private:

    void InitializeThreadsParameters();
    void InitializeFitParameters();
    void InitializeDataSamples();
    void InitializeFitter();

    void PassIfInitialized(const std::string& method_name_) const;
    void PassIfNotInitialized(const std::string& method_name_) const;

    void SaveParams(const std::vector<std::vector<double>>& new_pars);
    void SaveEventHist(int fititer, bool is_final = false);
    void SaveResults(const std::vector<std::vector<double>>& parresults,
                     const std::vector<std::vector<double>>& parerrors);

    // Parallel methods
    void ReWeightEvents(int iThread_ = -1);

    // Internal
    bool _isInitialized_{};
    bool _fitHasBeenDone_{};
    bool _fitHasConverged_{};

    // Options (Fitter)
    bool _disableChi2Pulls_{};


    // Options (Debug and Monitoring)
    bool _saveFitParameters_{};
    bool _saveEventTree_{};
    bool _disableMultiThread_{};
    bool _advancedTimeMonitoring_{};

    bool _apply_statistical_fluctuations_on_samples_{};
    bool _printFitState_{};

    int _PRNG_seed_{};
    int _nb_threads_{};
    int _nb_fit_parameters_{};
    int _nbFitCalls_{};
    int _saveFitParamsFrequency_{};
    int _nbTotalEvents_{};

    double _chi2Buffer_{};
    double _chi2StatBuffer_{};
    double _chi2PullsBuffer_{};
    double _chi2RegBuffer_{};

    std::vector<std::string> _parameter_names_;
    std::vector<bool> _parameterFixedStatusList_;
    std::vector<double> _parameterPriorValues_;
    std::vector<double> _parameterSteps_;
    std::vector<double> _parameter_low_edges_;
    std::vector<double> _parameter_high_edges_;
    std::vector<AnaFitParameters*> _fitParametersGroupList_;
    std::vector<AnaSample*> _samplesList_;
    std::vector<std::vector<double>> _newParametersList_;

    std::vector<std::vector<std::vector<long long>>> _durationReweightParameters_; // [SystGroup][Thread][Iteration] = time
    std::map<int, std::vector<double>> _durationHistoryHandler_;

    DataType _selected_data_type_;
    MinSettings _minimization_settings_;

    TRandom3* _PRNG_;
    TDirectory* _outputTDirectory_{};
    ROOT::Math::Minimizer* _minimizer_;
    ROOT::Math::Functor* _functor_;

    std::chrono::high_resolution_clock::time_point _startTime_;
    std::chrono::high_resolution_clock::time_point _endTime_;

    std::vector<double> _chi2StatHistory_;
    std::vector<double> _chi2PullsHistory_;
    std::vector<double> _chi2RegHistory_;

    // Multi-thread
    std::mutex _threadMutex_;
    int _counterThread_{};
    bool _stopThreads_{};
    std::vector<bool> _triggerReweightThreads_;
    std::vector<bool> _triggerReFillMcHistogramsThreads_;
    std::vector<std::future<void>> _asyncFitThreads_;
    std::vector<std::vector<std::vector<TH1D*>>> _histThreadHandlers_;

};

#endif // __ND280FITTER_HH__
