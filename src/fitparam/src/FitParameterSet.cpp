//
// Created by Nadrino on 21/05/2021.
//

#include "FitParameterSet.h"
#include "JsonUtils.h"

#include "GenericToolbox.h"
#include "GenericToolbox.Root.h"
#include "Logger.h"

LoggerInit([](){
  Logger::setUserHeaderStr("[FitParameterSet]");
} )

FitParameterSet::FitParameterSet() {
  this->reset();
}
FitParameterSet::~FitParameterSet() {
  this->reset();
}

void FitParameterSet::reset() {

  _isInitialized_ = false;

  _jsonConfig_ = nlohmann::json();
  if(_covarianceMatrixFile_ != nullptr){
    _covarianceMatrixFile_->Close(); // should delete every attached pointer
  }

  _maxEigenFraction_ = 1;
  _covarianceMatrixFile_ = nullptr;
  _originalCovarianceMatrix_ = nullptr;
  _originalCorrelationMatrix_ = nullptr;
  _parameterPriorList_ = nullptr;
  _parameterNamesList_ = nullptr;

  _parameterList_.clear();

}

void FitParameterSet::setJsonConfig(const nlohmann::json &jsonConfig) {
  _jsonConfig_ = jsonConfig;
  while(_jsonConfig_.is_string()){
    LogWarning << "Forwarding FitParameterSet config to: \"" << _jsonConfig_.get<std::string>() << "\"..." << std::endl;
    _jsonConfig_ = JsonUtils::readConfigFile(_jsonConfig_.get<std::string>());
  }
}
void FitParameterSet::setSaveDir(TDirectory* saveDir_){
  _saveDir_ = saveDir_;
}

void FitParameterSet::initialize() {

  if( _jsonConfig_.empty() ){
    LogError << "Json config not set" << std::endl;
    throw std::logic_error("json config not set");
  }

  _name_ = JsonUtils::fetchValue<std::string>(_jsonConfig_, "name", "");
  LogInfo << "Initializing parameter set: " << _name_ << std::endl;

  if( _saveDir_ != nullptr ){
    _saveDir_ = GenericToolbox::mkdirTFile(_saveDir_, _name_);
  }

  _isEnabled_ = JsonUtils::fetchValue<bool>(_jsonConfig_, "isEnabled");
  if( not _isEnabled_ ){
    LogWarning << _name_ << " parameters are disabled." << std::endl;
    return;
  }

  this->readCovarianceMatrix();

  _useOnlyOneParameterPerEvent_ = JsonUtils::fetchValue<bool>(_jsonConfig_, "useOnlyOneParameterPerEvent", false);

  if( _jsonConfig_.find("parameterLimits") != _jsonConfig_.end() ){
    auto parLimits = JsonUtils::fetchValue(_jsonConfig_, "parameterLimits", nlohmann::json());
    _globalParameterMinValue_ = JsonUtils::fetchValue(parLimits, "minValue", std::nan("UNSET"));
    _globalParameterMaxValue_ = JsonUtils::fetchValue(parLimits, "maxValue", std::nan("UNSET"));
  }

  // Optional parameters:
  std::string pathBuffer;

  // parameterPriorTVectorD
  pathBuffer = JsonUtils::fetchValue<std::string>(_jsonConfig_, "parameterPriorTVectorD", "");
  if(not pathBuffer.empty()){
    LogDebug << "Reading provided parameterPriorTVectorD..." << std::endl;
    _parameterPriorList_ = (TVectorT<double>*) _covarianceMatrixFile_->Get(pathBuffer.c_str());
    // Sanity checks
    if( _parameterPriorList_ == nullptr ){
      LogError << "Could not find \"" << pathBuffer << "\" into \"" << _covarianceMatrixFile_->GetName() << "\"" << std::endl;
      throw std::runtime_error("TObject not found.");
    }
    else if(_parameterPriorList_->GetNrows() != _originalCovarianceMatrix_->GetNrows() ){
      LogError << GET_VAR_NAME_VALUE(_parameterPriorList_->GetNrows() != _originalCovarianceMatrix_->GetNrows()) << std::endl;
      throw std::runtime_error("TObject size mismatch.");
    }
  }
  else{
    LogDebug << "No parameterPriorTVectorD provided, all parameter prior are set to 1." << std::endl;
    _parameterPriorList_ = new TVectorT<double>(_originalCovarianceMatrix_->GetNrows());
    for( int iPar = 0 ; iPar < _parameterPriorList_->GetNrows() ; iPar++ ){
      (*_parameterPriorList_)[iPar] = 1;
      if( _useEigenDecompInFit_ ) (*_originalParValues_)[iPar] = 1;
    }
  }

  pathBuffer = JsonUtils::fetchValue<std::string>(_jsonConfig_, "parameterLowerBoundsTVectorD", "");
  if( not pathBuffer.empty() ){
    _parameterLowerBoundsList_ = (TVectorT<double>*) _covarianceMatrixFile_->Get(pathBuffer.c_str());
    LogThrowIf(_parameterLowerBoundsList_ == nullptr, "Could not fetch parameterLowerBoundsTVectorD: \"" << pathBuffer)
    LogThrowIf(_parameterLowerBoundsList_->GetNrows() != _originalCovarianceMatrix_->GetNrows(),
               "parameterLowerBoundsTVectorD \"" << pathBuffer << "\" have not the right size.")
  }

  pathBuffer = JsonUtils::fetchValue<std::string>(_jsonConfig_, "parameterUpperBoundsTVectorD", "");
  if( not pathBuffer.empty() ){
    _parameterUpperBoundsList_ = (TVectorT<double>*) _covarianceMatrixFile_->Get(pathBuffer.c_str());
    LogThrowIf(_parameterUpperBoundsList_ == nullptr, "Could not fetch parameterUpperBoundsTVectorD: \"" << pathBuffer)
    LogThrowIf(_parameterUpperBoundsList_->GetNrows() != _originalCovarianceMatrix_->GetNrows(),
               "parameterUpperBoundsTVectorD \"" << pathBuffer << "\" have not the right size.")
  }

  // parameterNameTObjArray
  pathBuffer = JsonUtils::fetchValue<std::string>(_jsonConfig_, "parameterNameTObjArray", "");
  if(not pathBuffer.empty()){
    LogDebug << "Reading provided parameterNameTObjArray..." << std::endl;
    _parameterNamesList_ = (TObjArray*) _covarianceMatrixFile_->Get(pathBuffer.c_str());
    // Sanity checks
    if( _parameterNamesList_ == nullptr ){
      LogError << "Could not find \"" << pathBuffer << "\" into \"" << _covarianceMatrixFile_->GetName() << "\"" << std::endl;
      throw std::runtime_error("TObject not found.");
    }
    else if(_parameterNamesList_->GetSize() != _originalCovarianceMatrix_->GetNrows() ){
      LogError << GET_VAR_NAME_VALUE(_parameterNamesList_->GetSize() != _originalCovarianceMatrix_->GetNrows()) << std::endl;
      throw std::runtime_error("TObject size mismatch.");
    }
  }
  else{
    LogDebug << "No parameterNameTObjArray provided, parameters will be referenced with their index." << std::endl;
    _parameterNamesList_ = new TObjArray(_originalCovarianceMatrix_->GetNrows());
    for( int iPar = 0 ; iPar < _parameterPriorList_->GetNrows() ; iPar++ ){
      _parameterNamesList_->Add(new TNamed("", ""));
    }
  }

  LogDebug << "Initializing parameters..." << std::endl;
  _parameterList_.reserve(_originalCovarianceMatrix_->GetNcols()); // need to keep the memory at the same place -> FitParameter* will be used
  for(int iParameter = 0 ; iParameter < _originalCovarianceMatrix_->GetNcols() ; iParameter++ ){
    _parameterList_.emplace_back();
    _parameterList_.back().setParameterIndex(iParameter);
    _parameterList_.back().setName(_parameterNamesList_->At(iParameter)->GetName());
    _parameterList_.back().setParameterValue((*_parameterPriorList_)[iParameter]);
    _parameterList_.back().setPriorValue((*_parameterPriorList_)[iParameter]);
    _parameterList_.back().setStdDevValue(TMath::Sqrt((*_originalCovarianceMatrix_)[iParameter][iParameter]));
    _parameterList_.back().setDialsWorkingDirectory(JsonUtils::fetchValue<std::string>(_jsonConfig_, "dialSetWorkingDirectory", "./"));

    // OLD:
    if( JsonUtils::doKeyExist(_jsonConfig_, "dialSetDefinitions") ){
      _parameterList_.back().setDialSetConfig(JsonUtils::fetchValue<nlohmann::json>(_jsonConfig_, "dialSetDefinitions"));
    }

    if( JsonUtils::doKeyExist(_jsonConfig_, "parameterDefinitions") ){
      auto parsConfig = JsonUtils::fetchValue<nlohmann::json>(_jsonConfig_, "parameterDefinitions");
      JsonUtils::forwardConfig(parsConfig);
      auto parConfig = JsonUtils::fetchMatchingEntry(parsConfig, "parameterName", std::string(_parameterNamesList_->At(iParameter)->GetName()));
      if( parConfig.empty() ){
        // try with par index
        parConfig = JsonUtils::fetchMatchingEntry(parsConfig, "parameterIndex", iParameter);
      }
      _parameterList_.back().setParameterDefinitionConfig(parConfig);
    }

    _parameterList_.back().setEnableDialSetsSummary(JsonUtils::fetchValue<bool>(_jsonConfig_, "printDialSetsSummary", false));

    if( _globalParameterMinValue_ == _globalParameterMinValue_ ){
      _parameterList_.back().setMinValue(_globalParameterMinValue_);
    }
    if( _globalParameterMaxValue_ == _globalParameterMaxValue_ ){
      _parameterList_.back().setMaxValue(_globalParameterMaxValue_);
    }

    if( _parameterLowerBoundsList_ != nullptr ){
      _parameterList_.back().setMinValue((*_parameterLowerBoundsList_)[iParameter]);
    }
    if( _parameterUpperBoundsList_ != nullptr ){
      _parameterList_.back().setMaxValue((*_parameterUpperBoundsList_)[iParameter]);
    }

    _parameterList_.back().initialize();
  }


  if( _useEigenDecompInFit_ ){
    LogDebug << "Initializing eigen objects..." << std::endl;
    _originalParValues_ = std::shared_ptr<TVectorD>(new TVectorD(_parameterPriorList_->GetNrows()));
    for( int iPar = 0 ; iPar < _parameterPriorList_->GetNrows() ; iPar++ ){
      (*_originalParValues_)[iPar] = _parameterList_.at(iPar).getParameterValue();
    }
    propagateOriginalToEigen();
    _eigenParPriorValues_ = std::shared_ptr<TVectorD>( (TVectorD*) _eigenParValues_->Clone() );
    _eigenParStepSize_ = std::shared_ptr<TVectorD>(new TVectorD(_eigenParValues_->GetNrows()));
    _eigenParFixedList_.resize( _eigenParStepSize_->GetNrows(), false );
  }

  _isInitialized_ = true;
}


// Getters
bool FitParameterSet::isEnabled() const {
  return _isEnabled_;
}
const std::string &FitParameterSet::getName() const {
  return _name_;
}
std::vector<FitParameter> &FitParameterSet::getParameterList() {
  return _parameterList_;
}
const std::vector<FitParameter> &FitParameterSet::getParameterList() const{
  return _parameterList_;
}
TMatrixDSym *FitParameterSet::getOriginalCovarianceMatrix() const {
  return _originalCovarianceMatrix_;
}
const nlohmann::json &FitParameterSet::getJsonConfig() const {
  return _jsonConfig_;
}

// Core
size_t FitParameterSet::getNbParameters() const {
  return _parameterList_.size();
}
FitParameter& FitParameterSet::getFitParameter( size_t iPar_ ){
  return _parameterList_.at(iPar_);
}
double FitParameterSet::getChi2() const {

  if (not _isEnabled_) { return 0; }

  if (_inverseCovarianceMatrix_ == nullptr) {
    LogError << GET_VAR_NAME_VALUE(_inverseCovarianceMatrix_) << std::endl;
    throw std::runtime_error("inverse matrix not set");
  }

  double chi2 = 0;

  // EIGEN MATRIX PULLS ARE NOT WORKING
//  if( _useEigenDecompInFit_ ){
//    for( int iEigen = 0 ; iEigen < _eigenParValues_->GetNrows() ; iEigen++ ){
//      chi2 += std::pow( (*_eigenParValues_)[iEigen] - (*_eigenParPriorValues_)[iEigen], 2 )/(*_eigenValues_)[iEigen];
//    }
//  }
//  else{
    double iDelta, jDelta;
    for (int iPar = 0; iPar < _inverseCovarianceMatrix_->GetNrows(); iPar++) {
      if( not _parameterList_[iPar].isEnabled() ) continue;
      iDelta = (_parameterList_[iPar].getParameterValue() - _parameterList_[iPar].getPriorValue());
      for (int jPar = 0; jPar < _inverseCovarianceMatrix_->GetNrows(); jPar++) {
        if( not _parameterList_[jPar].isEnabled() ) continue;
        jDelta = iDelta;
        jDelta *= (_parameterList_[jPar].getParameterValue() - _parameterList_[jPar].getPriorValue());
        jDelta *= (*_inverseCovarianceMatrix_)(iPar, jPar);
        chi2 += jDelta;
      }
    }
//  }

  return chi2;
}


// Eigen
bool FitParameterSet::isUseEigenDecompInFit() const {
  return _useEigenDecompInFit_;
}
int FitParameterSet::getNbEnabledEigenParameters() const {
  return _nbEnabledEigen_;
}
double FitParameterSet::getEigenParameterValue(int iPar_) const{
  return (*_eigenParValues_)[iPar_];
}
double FitParameterSet::getEigenSigma(int iPar_) const{
  return TMath::Sqrt((*_eigenValues_)[iPar_]);
}
void FitParameterSet::setEigenParameter( int iPar_, double value_ ){
  LogThrowIf(_eigenParValues_ == nullptr, "_eigenParValues_ not set.");
  (*_eigenParValues_)[iPar_] = value_;
}
void FitParameterSet::setEigenParStepSize( int iPar_, double step_ ){
  LogThrowIf(_eigenParStepSize_ == nullptr, "_eigenParStepSize_ not set.");
  (*_eigenParStepSize_)[iPar_] = step_;
}
void FitParameterSet::setEigenParIsFixed( int iPar_, bool isFixed_ ){
  _eigenParFixedList_[iPar_] = isFixed_;
}

bool FitParameterSet::isEigenParFixed( int iPar_ ){
  return _eigenParFixedList_[iPar_];
}
double FitParameterSet::getEigenParStepSize( int iPar_ ){
  return (*_eigenParStepSize_)[iPar_];
}

const TMatrixD* FitParameterSet::getInvertedEigenVectors() const{
  return _invertedEigenVectors_.get();
}
const TMatrixD* FitParameterSet::getEigenVectors() const{
  return _eigenVectors_.get();
}
void FitParameterSet::propagateEigenToOriginal(){
  (*_originalParValues_) = (*_eigenParValues_);
  (*_originalParValues_) *= (*_invertedEigenVectors_);
  for( int iOrig = 0 ; iOrig < _originalParValues_->GetNrows() ; iOrig++ ){
    if( _parameterList_.at(iOrig).isFixed() ) continue;
    _parameterList_.at(iOrig).setParameterValue((*_originalParValues_)[iOrig]);
  }
}
void FitParameterSet::propagateOriginalToEigen(){
  (*_eigenParValues_) = (*_originalParValues_);
  (*_eigenParValues_) *= (*_eigenVectors_);
}

// Misc
std::string FitParameterSet::getSummary() const {
  std::stringstream ss;

  ss << "FitParameterSet: " << _name_ << " -> initialized=" << _isInitialized_ << ", enabled=" << _isEnabled_;

  if(_isInitialized_ and _isEnabled_){
    ss << ", nbParameters: " << _parameterList_.size() << "(defined)/" << _originalCovarianceMatrix_->GetNrows() << "(covariance)";
    if( not _parameterList_.empty() ){
      for( const auto& parameter : _parameterList_ ){
        ss << std::endl << GenericToolbox::indentString(parameter.getSummary(), 2);
      }
    }
  }

  return ss.str();
}


// Protected
void FitParameterSet::passIfInitialized(const std::string &methodName_) const {
  if( not _isInitialized_ ){
    LogError << "Can't do \"" << methodName_ << "\" while not initialized." << std::endl;
    throw std::logic_error("class not initialized");
  }
}
void FitParameterSet::readCovarianceMatrix(){

  _covarianceMatrixFile_ = TFile::Open(JsonUtils::fetchValue<std::string>(_jsonConfig_, "covarianceMatrixFilePath").c_str());
  if( not _covarianceMatrixFile_->IsOpen() ){
    LogError << "Could not open: _covarianceMatrixFile_: " << _covarianceMatrixFile_->GetPath() << std::endl;
    throw std::runtime_error("Could not open: _covarianceMatrixFile_");
  }

  _originalCovarianceMatrix_ = (TMatrixDSym*) _covarianceMatrixFile_->Get(
    JsonUtils::fetchValue<std::string>(_jsonConfig_, "covarianceMatrixTMatrixD").c_str()
  );
  if(_originalCovarianceMatrix_ == nullptr ){
    LogError << "Could not find: " << JsonUtils::fetchValue<std::string>(_jsonConfig_, "covarianceMatrixTMatrixD")
             << " in " << _covarianceMatrixFile_->GetPath() << std::endl;
    throw std::runtime_error("Could not find: covarianceMatrixTObjectPath");
  }

  if( _saveDir_ != nullptr ){
    _saveDir_->cd();
    _originalCovarianceMatrix_->Write("CovarianceMatrix_TMatrixDSym");
    GenericToolbox::convertToCorrelationMatrix((TMatrixD*)_originalCovarianceMatrix_)->Write("CorrelationMatrix_TMatrixD");
  }

  _originalCorrelationMatrix_ = std::shared_ptr<TMatrixDSym>(
    GenericToolbox::convertToSymmetricMatrix(GenericToolbox::convertToCorrelationMatrix((TMatrixD*) _originalCovarianceMatrix_))
    );

  _useEigenDecompInFit_ = JsonUtils::fetchValue(_jsonConfig_ , "useEigenDecompInFit", false);
  _maxEigenFraction_ = JsonUtils::fetchValue(_jsonConfig_ , "maxEigenFraction", double(1.));
  if( _maxEigenFraction_ != 1 ){
    LogInfo << "Max eigen fraction set to: " << _maxEigenFraction_*100 << "%" << std::endl;
    _useEigenDecompInFit_ = true;
  }

  LogWarning << "Computing inverse of the covariance matrix: " << _originalCovarianceMatrix_->GetNcols() << "x" << _originalCovarianceMatrix_->GetNrows() << std::endl;
  if( not _useEigenDecompInFit_ ){
    LogDebug << "Using default matrix inversion..." << std::endl;
    _inverseCovarianceMatrix_ = std::shared_ptr<TMatrixDSym>((TMatrixDSym*)(_originalCovarianceMatrix_->Clone()));
    _inverseCovarianceMatrix_->Invert();
  }
  else{
    LogInfo << "Using eigen decomposition..." << std::endl;

    LogDebug << "Computing the eigen vectors / values..." << std::endl;
    _eigenDecomp_ = std::shared_ptr<TMatrixDSymEigen>(new TMatrixDSymEigen(*_originalCovarianceMatrix_));
    LogDebug << "Eigen decomposition done." << std::endl;

    _eigenValues_ = std::shared_ptr<TVectorD>( (TVectorD*) _eigenDecomp_->GetEigenValues().Clone() );
    _eigenVectors_ = std::shared_ptr<TMatrixD>( (TMatrixD*) _eigenDecomp_->GetEigenVectors().Clone() );

    auto *eigVDec = new TMatrixDSymEigen(*GenericToolbox::convertToSymmetricMatrix(_eigenVectors_.get()));
//    eigVDec->GetEigenValues().Print();

    _invertedEigenVectors_ = std::shared_ptr<TMatrixD>( (TMatrixD*) _eigenVectors_->Clone() );
    _invertedEigenVectors_->Invert();

//    auto* idTest = (TMatrixD*) _eigenVectors_->Clone();
//    (*idTest) *= (*_invertedEigenVectors_);
////    idTest->Print();
//
//    auto* vecTest = new TVectorD(_invertedEigenVectors_->GetNrows());
//    for( int i = 0 ; i < vecTest->GetNrows() ; i++ ){
//      (*vecTest)[i] = 1;
//    }
//    vecTest->Print();
//
//    (*vecTest) *= (*_eigenVectors_); // orig -> eigen
//    vecTest->Print();
//
//    (*vecTest) *= (*_invertedEigenVectors_); // eigen -> orig
//    vecTest->Print(); // should be 1

    double eigenCumulative = 0;
    _nbEnabledEigen_ = 0;
    double eigenTotal = _eigenValues_->Sum();

    _inverseCovarianceMatrix_ = std::shared_ptr<TMatrixDSym>(new TMatrixDSym(_originalCovarianceMatrix_->GetNrows()));
    _effectiveCovarianceMatrix_ = std::shared_ptr<TMatrixDSym>(new TMatrixDSym(_originalCovarianceMatrix_->GetNrows()));
    _projectorMatrix_ = std::shared_ptr<TMatrixDSym>(new TMatrixDSym(_originalCovarianceMatrix_->GetNrows()));

    for (int iEigen = 0; iEigen < _eigenValues_->GetNrows(); iEigen++) {

      eigenCumulative += (*_eigenValues_)[iEigen];
      if( eigenCumulative / eigenTotal > _maxEigenFraction_ ){
        eigenCumulative -= (*_eigenValues_)[iEigen]; // not included
        break;
      }
      _nbEnabledEigen_++;

      for(int iDof = 0 ; iDof < _originalCovarianceMatrix_->GetNrows() ; iDof++ ){
        for(int jDof = 0 ; jDof < _originalCovarianceMatrix_->GetNcols() ; jDof++ ){
          (*_projectorMatrix_)[iDof][jDof]
            +=   (*_eigenVectors_)[iDof][iEigen]
                 * (*_eigenVectors_)[jDof][iEigen];

          (*_inverseCovarianceMatrix_)[iDof][jDof]
            += (1. / (*_eigenValues_)[iEigen] )
               * (*_eigenVectors_)[iDof][iEigen]
               * (*_eigenVectors_)[jDof][iEigen];

          (*_effectiveCovarianceMatrix_)[iDof][jDof]
            += (*_eigenValues_)[iEigen]
               * (*_eigenVectors_)[iDof][iEigen]
               * (*_eigenVectors_)[jDof][iEigen];
        }
      }

    } // iEigen

    LogWarning << "Eigen decomposition with " << _nbEnabledEigen_ << " / " << _eigenValues_->GetNrows() << " vectors" << std::endl;
    if(_nbEnabledEigen_ != _eigenValues_->GetNrows() ){
      LogInfo << "Max eigen fraction set to " << _maxEigenFraction_*100 << "%" << std::endl;
      LogInfo << "Fraction taken: " << eigenCumulative / eigenTotal*100 << "%" << std::endl;
    }

    _originalParValues_ = std::shared_ptr<TVectorD>( new TVectorD(_originalCovarianceMatrix_->GetNrows()) );
    _eigenParValues_ = std::shared_ptr<TVectorD>( new TVectorD(_originalCovarianceMatrix_->GetNrows()) );

  }

  LogInfo << "Parameter set \"" << _name_ << "\" is handling " << _originalCovarianceMatrix_->GetNcols() << " parameters." << std::endl;


}

bool FitParameterSet::isUseOnlyOneParameterPerEvent() const {
  return _useOnlyOneParameterPerEvent_;
}


