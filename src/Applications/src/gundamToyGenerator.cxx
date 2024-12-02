
#include "GundamGlobals.h"
#include "GundamApp.h"
#include "GundamUtils.h"
#include "RootUtils.h"
#include "FitterEngine.h"
#include "ConfigUtils.h"

#include "Logger.h"
#include "CmdLineParser.h"
#include "GenericToolbox.Json.h"
#include "GenericToolbox.Root.h"
#include "GenericToolbox.Utils.h"
#include "GenericToolbox.Map.h"

#include <TFile.h>
#include "TH1D.h"
#include "TH2D.h"

#include <string>
#include <vector>


LoggerInit([]{
  Logger::getUserHeader() << "[" << FILENAME << "]";
});


int main(int argc, char** argv){

  using namespace GundamUtils;

  GundamApp app{"toy generator tool"};

  // --------------------------
  // Read Command Line Args:
  // --------------------------
  CmdLineParser clParser;

  clParser.addDummyOption("Main options:");
  clParser.addOption("configFile", {"-c", "--config-file"}, "Specify path to the fitter config file");
  clParser.addOption("fitterFile", {"-f"}, "Specify the fitter output file");
  clParser.addOption("outputFile", {"-o", "--out-file"}, "Specify the ToyGenerator output file");
  clParser.addOption("useDataEntry", {"--use-data-entry"}, "Overrides \"selectedDataEntry\" in dataSet config. Second arg is to select a given dataset");
  clParser.addOption("fitSampleSetConfig", {"--fitsample-config"}, "override fitSampleSetConfig");
  clParser.addOption("plotGeneratorConfig", {"--plot-config"}, "override plotGeneratorConfig");
  clParser.addOption("nbThreads", {"-t", "--nb-threads"}, "Specify nb of parallel threads");
  clParser.addOption("nToys", {"-n"}, "Specify number of toys");
  clParser.addOption("randomSeed", {"-s", "--seed"}, "Set random seed");

  clParser.addDummyOption("Trigger options:");
  clParser.addTriggerOption("dryRun", {"-d", "--dry-run"}, "Only overrides fitter config and print it.");
  clParser.addTriggerOption("useBf", {"--use-bf"}, "Use best-fit as x-sec value instead of mean of toys.");
  clParser.addTriggerOption("usePreFit", {"--use-prefit"}, "Use prefit covariance matrices for the toy throws.");

  LogInfo << "Usage: " << std::endl;
  LogInfo << clParser.getConfigSummary() << std::endl << std::endl;

  clParser.parseCmdLine(argc, argv);

  LogThrowIf(clParser.isNoOptionTriggered(), "No option was provided.");

  LogInfo << "Provided arguments: " << std::endl;
  LogInfo << clParser.getValueSummary() << std::endl << std::endl;

  // Sanity checks
  LogThrowIf(not clParser.isOptionTriggered("configFile"), "Toy generator config file not provided.");
  LogThrowIf(not clParser.isOptionTriggered("fitterFile"), "Did not provide the output fitter file.");
  LogThrowIf(not clParser.isOptionTriggered("nToys"), "Did not provide number of toys.");

  // Global parameters
  gRandom = new TRandom3(0);     // Initialize with a UUID
  if( clParser.isOptionTriggered("randomSeed") ){
    LogAlert << "Using user-specified random seed: " << clParser.getOptionVal<ULong_t>("randomSeed") << std::endl;
    gRandom->SetSeed(clParser.getOptionVal<ULong_t>("randomSeed"));
  }
  else{
    ULong_t seed = time(nullptr);
    LogInfo << "Using \"time(nullptr)\" random seed: " << seed << std::endl;
    gRandom->SetSeed(seed);
  }

  GundamGlobals::setNumberOfThreads( clParser.getOptionVal("nbThreads", 1) );
  LogInfo << "Running the fitter with " << GundamGlobals::getNbCpuThreads() << " parallel threads." << std::endl;

  // Reading fitter file
  std::string fitterFile{clParser.getOptionVal<std::string>("fitterFile")};
  std::unique_ptr<TFile> fitterRootFile{nullptr};
  JsonType fitterConfig; // will be used to load the propagator

  if( GenericToolbox::hasExtension(fitterFile, "root") ){
    LogWarning << "Opening fitter output file: " << fitterFile << std::endl;
    fitterRootFile = std::unique_ptr<TFile>( TFile::Open( fitterFile.c_str() ) );
    LogThrowIf( fitterRootFile == nullptr, "Could not open fitter output file." );

    RootUtils::ObjectReader::throwIfNotFound = true;

    RootUtils::ObjectReader::readObject<TNamed>(fitterRootFile.get(), {{"gundam/config_TNamed"}, {"gundamFitter/unfoldedConfig_TNamed"}}, [&](TNamed* config_){
      fitterConfig = GenericToolbox::Json::readConfigJsonStr( config_->GetTitle() );
    });
  }
  else{
    LogWarning << "Reading fitter config file: " << fitterFile << std::endl;
    fitterConfig = GenericToolbox::Json::readConfigFile( fitterFile );

    clParser.getOptionPtr("usePreFit")->setIsTriggered( true );
  }

  LogAlertIf(clParser.isOptionTriggered("usePreFit")) << "Pre-fit mode enabled: will throw toys according to the prior covariance matrices..." << std::endl;

  ConfigUtils::ConfigHandler cHandler{ fitterConfig };

  // Disabling defined fit samples:
  LogInfo << "Removing defined samples..." << std::endl;
  GenericToolbox::Json::clearEntry( cHandler.getConfig(), "fitterEngineConfig/likelihoodInterfaceConfig/dataSetManagerConfig/propagatorConfig/sampleSetConfig/sampleList" );
  GenericToolbox::Json::clearEntry( cHandler.getConfig(), "fitterEngineConfig/likelihoodInterfaceConfig/dataSetManagerConfig/propagatorConfig/fitSampleSetConfig/fitSampleList" );
  GenericToolbox::Json::clearEntry( cHandler.getConfig(), "fitterEngineConfig/propagatorConfig/fitSampleSetConfig/fitSampleList" );

  // Disabling defined plots:
  LogInfo << "Removing defined plots..." << std::endl;
  GenericToolbox::Json::clearEntry( cHandler.getConfig(), "fitterEngineConfig/likelihoodInterfaceConfig/dataSetManagerConfig/propagatorConfig/plotGeneratorConfig" );
  GenericToolbox::Json::clearEntry( cHandler.getConfig(), "fitterEngineConfig/propagatorConfig/plotGeneratorConfig" );

  // Defining signal samples
  JsonType toyConfig{ ConfigUtils::readConfigFile( clParser.getOptionVal<std::string>("configFile") ) };
  cHandler.override( toyConfig );

  if( clParser.isOptionTriggered("fitSampleSetConfig") ){
    JsonType fitSampleSetConfig_new{ ConfigUtils::readConfigFile( clParser.getOptionVal<std::string>("fitSampleSetConfig") ) };
    cHandler.getConfig()["fitterEngineConfig"]["likelihoodInterfaceConfig"]["dataSetManagerConfig"]["propagatorConfig"]["fitSampleSetConfig"]["fitSampleList"] = fitSampleSetConfig_new["fitSampleList"];
  }
  if( clParser.isOptionTriggered("plotGeneratorConfig") ){
    std::vector< std::string > plotConfigKeysToCopy = {
      "varDictionnaries",
      "histogramsDefinition",
      "canvasParameters",
    };
    JsonType plotGeneratorConfig_new{ ConfigUtils::readConfigFile( clParser.getOptionVal<std::string>("plotGeneratorConfig") ) };
    for(const auto& k: plotConfigKeysToCopy){
      cHandler.getConfig()["fitterEngineConfig"]["likelihoodInterfaceConfig"]["dataSetManagerConfig"]["propagatorConfig"]["plotGeneratorConfig"][k] = plotGeneratorConfig_new[k];
    }
  }

  LogInfo << "Override done." << std::endl;

  LogInfo << "Fetching propagator config into fitter config..." << std::endl;

  // it will handle all the deprecated config options and names properly
  FitterEngine fitter{nullptr};
  fitter.configure( GenericToolbox::Json::fetchValue<JsonType>( cHandler.getConfig(), "fitterEngineConfig" ) );

  // We are only interested in our MC. Data has already been used to get the post-fit error/values
  // TODO Check
  //fitter.getLikelihoodInterface().setForceAsimovData( true );

  // Disabling eigen decomposed parameters
  fitter.getLikelihoodInterface().getModelPropagator().setEnableEigenToOrigInPropagate( false );

  // --use-data-entry
  if( clParser.isOptionTriggered("useDataEntry") ){
    auto selectedDataEntry = clParser.getOptionVal<std::string>("useDataEntry", 0);
    // Do something better in case multiple datasets are defined
    bool isFound{false};
    for( auto& dataSet : fitter.getLikelihoodInterface().getDatasetList() ){
      if( GenericToolbox::isIn( selectedDataEntry, dataSet.getDataDispenserDict() ) ){
        LogWarning << "Using data entry \"" << selectedDataEntry << "\" for dataset: " << dataSet.getName() << std::endl;
        dataSet.setSelectedDataEntry( selectedDataEntry );
        isFound = true;
      }
    }
    LogThrowIf(not isFound, "Could not find data entry \"" << selectedDataEntry << "\" among defined data sets");
  }  

  // Sample binning using parameterSetName
  for( auto& sample : fitter.getLikelihoodInterface().getModelPropagator().getSampleSet().getSampleList() ){

    if( clParser.isOptionTriggered("usePreFit") ){
      sample.setName( sample.getName() + " (pre-fit)" );
    }

    // binning already set?
    if( not sample.getBinningFilePath().empty() ){ continue; }

    LogScopeIndent;
    LogInfo << sample.getName() << ": binning not set, looking for parSetBinning..." << std::endl;
    auto associatedParSet = GenericToolbox::Json::fetchValue(
        sample.getConfig(),
        {{"parSetBinning"}, {"parameterSetName"}},
        std::string()
    );

    LogThrowIf(associatedParSet.empty(), "Could not find parSetBinning.");

    // Looking for parSet
    auto foundDialCollection = std::find_if(
        fitter.getLikelihoodInterface().getModelPropagator().getDialCollectionList().begin(),
        fitter.getLikelihoodInterface().getModelPropagator().getDialCollectionList().end(),
        [&](const DialCollection& dialCollection_){
          auto* parSetPtr{dialCollection_.getSupervisedParameterSet()};
          if( parSetPtr == nullptr ){ return false; }
          return ( parSetPtr->getName() == associatedParSet );
        });
    LogThrowIf(
        foundDialCollection == fitter.getLikelihoodInterface().getModelPropagator().getDialCollectionList().end(),
        "Could not find " << associatedParSet << " among fit dial collections: "
                          << GenericToolbox::toString(fitter.getLikelihoodInterface().getModelPropagator().getDialCollectionList(),
                                                      [](const DialCollection& dialCollection_){
                                                        return dialCollection_.getTitle();
                                                      }
                          ));

    LogThrowIf(foundDialCollection->getDialBinSet().getBinList().empty(), "Could not find binning");
    sample.setBinningFilePath( foundDialCollection->getDialBinSet().getFilePath() );

  }

  // Load everything
  fitter.getLikelihoodInterface().initialize();

  Propagator& propagator{fitter.getLikelihoodInterface().getModelPropagator()};


  if( clParser.isOptionTriggered("dryRun") ){
    std::cout << cHandler.toString() << std::endl;

    LogAlert << "Exiting as dry-run is set." << std::endl;
    return EXIT_SUCCESS;
  }


  if( not clParser.isOptionTriggered("usePreFit") and fitterRootFile != nullptr ){

    // Load post-fit parameters as "prior" so we can reset the weight to this point when throwing toys
    LogWarning << std::endl << GenericToolbox::addUpDownBars("Injecting post-fit parameters...") << std::endl;
    RootUtils::ObjectReader::readObject<TNamed>( fitterRootFile.get(), "FitterEngine/postFit/parState_TNamed", [&](TNamed* parState_){
      propagator.getParametersManager().injectParameterValues( GenericToolbox::Json::readConfigJsonStr( parState_->GetTitle() ) );
      for( auto& parSet : propagator.getParametersManager().getParameterSetsList() ){
        if( not parSet.isEnabled() ){ continue; }
        for( auto& par : parSet.getParameterList() ){
          if( not par.isEnabled() ){ continue; }
          par.setPriorValue( par.getParameterValue() );
        }
      }
    });

    // Load the post-fit covariance matrix
    LogWarning << std::endl << GenericToolbox::addUpDownBars("Injecting post-fit covariance matrix...") << std::endl;
    RootUtils::ObjectReader::readObject<TH2D>(
        fitterRootFile.get(), "FitterEngine/postFit/Hesse/hessian/postfitCovarianceOriginal_TH2D",
        [&](TH2D* hCovPostFit_){
          propagator.getParametersManager().setGlobalCovarianceMatrix(std::make_shared<TMatrixD>(hCovPostFit_->GetNbinsX(), hCovPostFit_->GetNbinsX()));
          for( int iBin = 0 ; iBin < hCovPostFit_->GetNbinsX() ; iBin++ ){
            for( int jBin = 0 ; jBin < hCovPostFit_->GetNbinsX() ; jBin++ ){
              (*propagator.getParametersManager().getGlobalCovarianceMatrix())[iBin][jBin] = hCovPostFit_->GetBinContent(1 + iBin, 1 + jBin);
            }
          }
        }
    );
  }

  // Creating output file
  std::string outFilePath{};
  if( clParser.isOptionTriggered("outputFile") ){ outFilePath = clParser.getOptionVal<std::string>("outputFile"); }
  else{
    // appendixDict["optionName"] = "Appendix"
    // this list insure all appendices will appear in the same order
    std::vector<GundamUtils::AppendixEntry> appendixDict{
        {"configFile", ""},
        {"fitterFile", "Fit"},
        {"nToys", "nToys"},
        {"randomSeed", "Seed"},
        {"usePreFit", "PreFit"},
    };

    outFilePath = "ToyGeneration_" + GundamUtils::generateFileName(clParser, appendixDict) + ".root";

    std::string outFolder{GenericToolbox::Json::fetchValue<std::string>(toyConfig, "outputFolder", "./")};
    outFilePath = GenericToolbox::joinPath(outFolder, outFilePath);
  }

  app.setCmdLinePtr( &clParser );
  app.setConfigString( ConfigUtils::ConfigHandler{toyConfig}.toString() );
  app.openOutputFile( outFilePath );
  app.writeAppInfo();

  auto* toyGenDir{ GenericToolbox::mkdirTFile(app.getOutfilePtr(), "toyGen") };
  bool useBestFitAsCentralValue{
    clParser.isOptionTriggered("useBf")
    or GenericToolbox::Json::fetchValue<bool>(toyConfig, "useBestFitAsCentralValue", false)
  };

  LogInfo << "Creating throws tree" << std::endl;
  auto* toyThrowTree = new TTree("toyThrow", "toyThrow");
  toyThrowTree->SetDirectory( GenericToolbox::mkdirTFile(toyGenDir, "throws") ); // temp saves will be done here

  auto* bestFitTree = new TTree("bestFitTree", "bestFitTree");
  bestFitTree->SetDirectory( GenericToolbox::mkdirTFile(toyGenDir, "throws") ); // temp saves will be done here

  struct ToyData{
    Sample* samplePtr{nullptr};
    JsonType config{};
    GenericToolbox::RawDataArray branchBinsData{};

    TH1D histogram{};
  };
  std::vector<ToyData> ToyDataList{};

  LogInfo << "Initializing toy generation samples..." << std::endl;
  ToyDataList.reserve(propagator.getSampleSet().getSampleList().size() );
  for( auto& sample : propagator.getSampleSet().getSampleList() ){
    ToyDataList.emplace_back();
    auto& toyDataEntry = ToyDataList.back();

    LogScopeIndent;
    LogInfo << "Defining toy entry: " << sample.getName() << std::endl;
    toyDataEntry.samplePtr = &sample;
    toyDataEntry.config = sample.getConfig();
    toyDataEntry.branchBinsData.resetCurrentByteOffset();
    std::vector<std::string> leafNameList{};
    leafNameList.reserve( sample.getHistogram().getNbBins() );
    for( int iBin = 0 ; iBin < sample.getHistogram().getNbBins(); iBin++ ){
      leafNameList.emplace_back(Form("bin_%i/D", iBin));
      toyDataEntry.branchBinsData.writeRawData( double(0) );
    }
    toyDataEntry.branchBinsData.lockArraySize();

    toyThrowTree->Branch(
        GenericToolbox::generateCleanBranchName( sample.getName() ).c_str(),
        toyDataEntry.branchBinsData.getRawDataArray().data(),
        GenericToolbox::joinVectorString(leafNameList, ":").c_str()
    );
    bestFitTree->Branch(
        GenericToolbox::generateCleanBranchName( sample.getName() ).c_str(),
        toyDataEntry.branchBinsData.getRawDataArray().data(),
        GenericToolbox::joinVectorString(leafNameList, ":").c_str()
    );

    toyDataEntry.histogram = TH1D(
        sample.getName().c_str(),
        sample.getName().c_str(),
        sample.getHistogram().getNbBins(),
        0,
        sample.getHistogram().getNbBins()
    );
  }

  int nToys{ clParser.getOptionVal<int>("nToys") };

  bool enableEventMcThrow{true};
  bool enableStatThrowInToys{true};
  auto statThrowConfig   = GenericToolbox::Json::fetchValue( cHandler.getConfig(), "statThrowConfig", JsonType() );
  enableStatThrowInToys = GenericToolbox::Json::fetchValue( statThrowConfig, "enableStatThrowInToys", enableStatThrowInToys);
  enableEventMcThrow    = GenericToolbox::Json::fetchValue( statThrowConfig, "enableEventMcThrow", enableEventMcThrow);

  auto writeBinDataFct = std::function<void()>([&]{
    for( auto& toyData : ToyDataList ){

      toyData.branchBinsData.resetCurrentByteOffset();
      for( int iBin = 0 ; iBin < toyData.samplePtr->getHistogram().getNbBins() ; iBin++ ){
        double binData{ toyData.samplePtr->getHistogram().getBinContentList()[iBin].sumWeights };

        // bin volume
        auto& bin = toyData.samplePtr->getHistogram().getBinContextList()[iBin].bin;
        double binVolume{1};

        for( auto& edges : bin.getEdgesList() ){
          if( edges.isConditionVar ){ continue; } // no volume, just a condition variable

          binVolume *= (edges.max - edges.min);
        }

        binData /= binVolume;
        toyData.branchBinsData.writeRawData( binData );
      }
    }
  });

  {
    LogWarning << "Calculating weight at best-fit" << std::endl;
    for( auto& parSet : propagator.getParametersManager().getParameterSetsList() ){ parSet.moveParametersToPrior(); }
    propagator.propagateParameters();
    writeBinDataFct();
    bestFitTree->Fill();
    GenericToolbox::writeInTFile( GenericToolbox::mkdirTFile(toyGenDir, "throws"), bestFitTree );
  }


  //////////////////////////////////////
  // THROWS LOOP
  /////////////////////////////////////
  LogWarning << std::endl << GenericToolbox::addUpDownBars( "Generating toys..." ) << std::endl;
  propagator.getParametersManager().initializeStrippedGlobalCov();

  // stats printing
  GenericToolbox::Time::AveragedTimer<1> totalTimer{};
  GenericToolbox::Time::AveragedTimer<1> throwTimer{};
  GenericToolbox::Time::AveragedTimer<1> propagateTimer{};
  GenericToolbox::Time::AveragedTimer<1> otherTimer{};
  GenericToolbox::Time::AveragedTimer<1> writeTimer{};
  GenericToolbox::TablePrinter t{};
  std::stringstream progressSs;
  std::stringstream ss; ss << LogWarning.getPrefixString() << "Generating " << nToys << " toys...";
  for( int iToy = 0 ; iToy < nToys ; iToy++ ){

    t.reset();
    t << "Total time" << GenericToolbox::TablePrinter::NextColumn;
    t << "Throw toys" << GenericToolbox::TablePrinter::NextColumn;
    t << "Propagate pars" << GenericToolbox::TablePrinter::NextColumn;
    t << "Re-normalize" << GenericToolbox::TablePrinter::NextColumn;
    t << "Write throws" << GenericToolbox::TablePrinter::NextLine;

    t << totalTimer << GenericToolbox::TablePrinter::NextColumn;
    t << throwTimer << GenericToolbox::TablePrinter::NextColumn;
    t << propagateTimer << GenericToolbox::TablePrinter::NextColumn;
    t << otherTimer << GenericToolbox::TablePrinter::NextColumn;
    t << writeTimer << GenericToolbox::TablePrinter::NextLine;

    totalTimer.stop();
    totalTimer.start();

    // loading...
    progressSs.str("");
    progressSs << t.generateTableString() << std::endl;
    progressSs << ss.str();
    GenericToolbox::displayProgressBar( iToy+1, nToys, progressSs.str() );

    // Do the throwing:
    throwTimer.start();
    propagator.getParametersManager().throwParametersFromGlobalCovariance( not GundamGlobals::isDebug() );
    throwTimer.stop();

    propagateTimer.start();
    propagator.propagateParameters();

    if( enableStatThrowInToys ){
      for( auto& toyData : ToyDataList ){
        if( enableEventMcThrow ){
          // Take into account the finite amount of event in MC
          toyData.samplePtr->getHistogram().throwEventMcError();
        }
        // JK: throwEventMcError() already treated each MC entry as Poisson event
        //     we don't want to throw again on weighted MC histogram
        // Asimov bin content -> toy data
        //toyData.samplePtr->getHistogram().throwStatError();
      }
    }
    propagateTimer.stop();

    otherTimer.start();
    // TODO: parallelize this
    writeBinDataFct();
    otherTimer.stop();

    // Write the branches
    writeTimer.start();
    toyThrowTree->Fill();
    writeTimer.stop();
  }


  LogInfo << "Writing throws..." << std::endl;
  GenericToolbox::writeInTFile( GenericToolbox::mkdirTFile(toyGenDir, "throws"), toyThrowTree );

  LogInfo << "Calculating mean & covariance matrix..." << std::endl;
  auto* meanValuesVector = GenericToolbox::generateMeanVectorOfTree(
      useBestFitAsCentralValue ? bestFitTree : toyThrowTree
  );
  auto* globalCovMatrix = GenericToolbox::generateCovarianceMatrixOfTree( toyThrowTree );

  auto* globalCovMatrixHist = GenericToolbox::convertTMatrixDtoTH2D(globalCovMatrix);
  auto* globalCorMatrixHist = GenericToolbox::convertTMatrixDtoTH2D(GenericToolbox::convertToCorrelationMatrix(globalCovMatrix));

  std::vector<TH1D> binValues{};
  binValues.reserve(propagator.getSampleSet().getSampleList().size() );
  int iBinGlobal{-1};

  for( auto& toyData : ToyDataList ){

    for( int iBin = 0 ; iBin < toyData.samplePtr->getHistogram().getNbBins() ; iBin++ ){
      iBinGlobal++;

      std::string binTitle = toyData.samplePtr->getHistogram().getBinContextList()[iBin].bin.getSummary();
      double binVolume = toyData.samplePtr->getHistogram().getBinContextList()[iBin].bin.getVolume();

      toyData.histogram.SetBinContent( 1+iBin, (*meanValuesVector)[iBinGlobal] );
      toyData.histogram.SetBinError( 1+iBin, TMath::Sqrt( (*globalCovMatrix)[iBinGlobal][iBinGlobal] ) );
      toyData.histogram.GetXaxis()->SetBinLabel( 1+iBin, binTitle.c_str() );

      globalCovMatrixHist->GetXaxis()->SetBinLabel(1+iBinGlobal, GenericToolbox::joinPath(toyData.samplePtr->getName(), binTitle).c_str());
      globalCorMatrixHist->GetXaxis()->SetBinLabel(1+iBinGlobal, GenericToolbox::joinPath(toyData.samplePtr->getName(), binTitle).c_str());
      globalCovMatrixHist->GetYaxis()->SetBinLabel(1+iBinGlobal, GenericToolbox::joinPath(toyData.samplePtr->getName(), binTitle).c_str());
      globalCorMatrixHist->GetYaxis()->SetBinLabel(1+iBinGlobal, GenericToolbox::joinPath(toyData.samplePtr->getName(), binTitle).c_str());
    }

    toyData.histogram.SetMarkerStyle(kFullDotLarge);
    toyData.histogram.SetMarkerColor(kGreen-3);
    toyData.histogram.SetMarkerSize(0.5);
    toyData.histogram.SetLineWidth(2);
    toyData.histogram.SetLineColor(kGreen-3);
    toyData.histogram.SetDrawOption("E1");
    toyData.histogram.GetXaxis()->LabelsOption("v");
    toyData.histogram.GetXaxis()->SetLabelSize(0.02);
    toyData.histogram.GetYaxis()->SetTitle( GenericToolbox::Json::fetchValue(toyData.samplePtr->getConfig(), "yAxis", "#delta#sigma").c_str() );

    GenericToolbox::writeInTFile(
        GenericToolbox::mkdirTFile(toyGenDir, "histograms"),
        &toyData.histogram, GenericToolbox::generateCleanBranchName( toyData.samplePtr->getName() )
    );

  }

  globalCovMatrixHist->GetXaxis()->SetLabelSize(0.02);
  globalCovMatrixHist->GetYaxis()->SetLabelSize(0.02);
  GenericToolbox::writeInTFile(GenericToolbox::mkdirTFile(toyGenDir, "matrices"), globalCovMatrixHist, "covarianceMatrix");

  globalCorMatrixHist->GetXaxis()->SetLabelSize(0.02);
  globalCorMatrixHist->GetYaxis()->SetLabelSize(0.02);
  globalCorMatrixHist->GetZaxis()->SetRangeUser(-1, 1);
  GenericToolbox::writeInTFile(GenericToolbox::mkdirTFile(toyGenDir, "matrices"), globalCorMatrixHist, "correlationMatrix");

  LogInfo << "Generating toyData sample plots..." << std::endl;

  for( auto& parSet : propagator.getParametersManager().getParameterSetsList() ){ parSet.moveParametersToPrior(); }
  propagator.propagateParameters();

  // manual trigger to tweak the error bars
  fitter.getLikelihoodInterface().getPlotGenerator().generateSampleHistograms( GenericToolbox::mkdirTFile(toyGenDir, "plots/histograms") );

  for( auto& histHolder : fitter.getLikelihoodInterface().getPlotGenerator().getHistHolderList(0) ){
    if( not histHolder.isData ){ continue; } // only data will print errors

    const ToyData* xsecDataPtr{nullptr};
    for( auto& xsecData : ToyDataList ){
      if( xsecData.samplePtr  == histHolder.samplePtr){
        xsecDataPtr = &xsecData;
        break;
      }
    }
    LogThrowIf(xsecDataPtr==nullptr, "corresponding data not found");

    // alright, now rescale error bars
    for( int iBin = 0 ; iBin < histHolder.histPtr->GetNbinsX() ; iBin++ ){
      // relative error should be set
      histHolder.histPtr->SetBinError(
          1+iBin,
          histHolder.histPtr->GetBinContent(1+iBin)
          * xsecDataPtr->histogram.GetBinError(1+iBin)
          / xsecDataPtr->histogram.GetBinContent(1+iBin)
      );
    }
  }

  fitter.getLikelihoodInterface().getPlotGenerator().generateCanvas(
      fitter.getLikelihoodInterface().getPlotGenerator().getHistHolderList(0),
      GenericToolbox::mkdirTFile(toyGenDir, "plots/canvas")
  );


  LogInfo << "Writing event samples in TTrees..." << std::endl;
  fitter.getLikelihoodInterface().writeEvents(
      {GenericToolbox::mkdirTFile(toyGenDir, "toyGen"),
      "events"}
  );

}
