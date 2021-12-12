//
// Created by Nadrino on 21/05/2021.
//

#ifndef GUNDAM_FITPARAMETERSET_H
#define GUNDAM_FITPARAMETERSET_H

#include "vector"
#include "string"

#include "json.hpp"
#include "TMatrixDSym.h"
#include "TVectorT.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TVectorT.h"
#include "TMatrixDSymEigen.h"

#include "Logger.h"

#include "FitParameter.h"


/*
 * FitParameterSet is a class which aims at handling a set of parameters bond together with a covariance matrix
 * User parameters:
 * - Covariance matrix (dim N)
 * - N Fit Parameters (handing dials)
 *
 * */

class FitParameterSet {

public:
  FitParameterSet();
  virtual ~FitParameterSet();

  void reset();

  // Setters
  void setConfig(const nlohmann::json &config_);
  void setSaveDir(TDirectory* saveDir_);

  // Init
  void initialize();

  // Getters
  bool isEnabled() const;

  bool isEnableThrowMcBeforeFit() const;

  bool isUseOnlyOneParameterPerEvent() const;
  const std::string &getName() const;
  std::vector<FitParameter> &getParameterList();
  const std::vector<FitParameter> &getParameterList() const;
  TMatrixDSym *getOriginalCovarianceMatrix() const;
  const nlohmann::json &getConfig() const;

  // Core
  size_t getNbParameters() const;
  FitParameter& getFitParameter( size_t iPar_ );
  double getChi2() const;

  // Throw / Shifts
  void moveFitParametersToPrior();
  void throwFitParameters(double gain_ = 1);

  // Eigen decomposition
  void setEigenParameter( int iPar_, double value_ );
  void setEigenParStepSize( int iPar_, double step_ );
  void setEigenParIsFixed( int iPar_, bool isFixed_ );

  bool isEigenParFixed( int iPar_ ) const;
  double getEigenParStepSize( int iPar_ ) const;

  bool isUseEigenDecompInFit() const;
  int getNbEnabledEigenParameters() const;
  double getEigenParameterValue(int iPar_) const;
  double getEigenSigma(int iPar_) const;
  const TMatrixD* getInvertedEigenVectors() const;
  const TMatrixD* getEigenVectors() const;
  void propagateEigenToOriginal();
  void propagateOriginalToEigen();

  // Misc
  std::string getSummary() const;

  static double toNormalizedParRange(double parRange, const FitParameter& par);
  static double toNormalizedParValue(double parValue, const FitParameter& par);
  static double toRealParValue(double normParValue, const FitParameter& par);
  static double toRealParRange(double normParRange, const FitParameter& par);

  double toNormalizedEigenParRange(double parRange, int parIndex) const;
  double toNormalizedEigenParValue(double parValue, int parIndex) const;
  double toRealEigenParValue(double normParValue, int parIndex) const;
  double toRealEigenParRange(double normParRange, int parIndex) const;

protected:
  void passIfInitialized(const std::string& methodName_) const;

  void readCovarianceMatrix();

private:
  // User parameters
  nlohmann::json _config_;

  // Internals
  bool _isInitialized_{false};
  std::vector<FitParameter> _parameterList_;
  TDirectory* _saveDir_{nullptr};

  // JSON
  std::string _name_;
  bool _isEnabled_{};
  bool _throwMcBeforeFit_{true};
  double _maxEigenFraction_{1};

  // Input file:
  std::shared_ptr<TFile> _covarianceMatrixFile_{nullptr};
  TMatrixDSym* _covarianceMatrix_{nullptr};
  TVectorD* _parameterPriorList_{nullptr};
  TVectorD* _parameterLowerBoundsList_{nullptr};
  TVectorD* _parameterUpperBoundsList_{nullptr};
  TObjArray* _parameterNamesList_{nullptr};

  double _globalParameterMinValue_{std::nan("UNSET")};
  double _globalParameterMaxValue_{std::nan("UNSET")};

  // Eigen objects
  int _nbEnabledEigen_{0};
  bool _useEigenDecompInFit_{false};
  bool _useOnlyOneParameterPerEvent_{false};
  std::shared_ptr<TMatrixDSymEigen> _eigenDecomp_{nullptr};
  std::shared_ptr<TVectorD> _eigenValues_{nullptr};
  std::shared_ptr<TVectorD> _eigenValuesInv_{nullptr};
  std::shared_ptr<TMatrixD> _eigenVectors_{nullptr};
  std::shared_ptr<TMatrixD> _eigenVectorsInv_{nullptr};
  std::shared_ptr<TMatrixD> _projectorMatrix_{nullptr};
  std::shared_ptr<TMatrixD> _inverseCovarianceMatrix_{nullptr};
  std::shared_ptr<TMatrixD> _effectiveCovarianceMatrix_{nullptr};

  std::shared_ptr<TMatrixD> _originalCorrelationMatrix_{nullptr};
  std::shared_ptr<TMatrixD> _effectiveCorrelationMatrix_{nullptr};

  std::shared_ptr<TMatrixD> _choleskyMatrix_{nullptr};

  std::shared_ptr<TVectorD> _eigenParValues_;
  std::shared_ptr<TVectorD> _originalParValues_;
  std::shared_ptr<TVectorD> _eigenParPriorValues_;
  std::shared_ptr<TVectorD> _eigenParStepSize_;
  std::vector<bool> _eigenParFixedList_;

};


#endif //GUNDAM_FITPARAMETERSET_H
