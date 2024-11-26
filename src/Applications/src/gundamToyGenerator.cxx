
#include "GundamGlobals.h"
#include "GundamApp.h"
#include "GundamUtils.h"
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

  GundamGlobals::getParallelWorker().setNThreads( clParser.getOptionVal("nbThreads", 1) );
  LogInfo << "Running the fitter with " << GundamGlobals::getParallelWorker().getNbThreads() << " parallel threads." << std::endl;

  // Reading fitter file
  std::string fitterFile{clParser.getOptionVal<std::string>("fitterFile")};
  std::unique_ptr<TFile> fitterRootFile{nullptr};
  JsonType fitterConfig; // will be used to load the propagator

  if( GenericToolbox::hasExtension(fitterFile, "root") ){
    LogWarning << "Opening fitter output file: " << fitterFile << std::endl;
    fitterRootFile = std::unique_ptr<TFile>( TFile::Open( fitterFile.c_str() ) );
    LogThrowIf( fitterRootFile == nullptr, "Could not open fitter output file." );

    ObjectReader::throwIfNotFound = true;

    ObjectReader::readObject<TNamed>(fitterRootFile.get(), {{"gundam/config_TNamed"}, {"gundamFitter/unfoldedConfig_TNamed"}}, [&](TNamed* config_){
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
  ConfigUtils::clearEntry( cHandler.getConfig(), "fitterEngineConfig/likelihoodInterfaceConfig/dataSetManagerConfig/propagatorConfig/sampleSetConfig/sampleList" );
  ConfigUtils::clearEntry( cHandler.getConfig(), "fitterEngineConfig/likelihoodInterfaceConfig/dataSetManagerConfig/propagatorConfig/fitSampleSetConfig/fitSampleList" );
  ConfigUtils::clearEntry( cHandler.getConfig(), "fitterEngineConfig/propagatorConfig/fitSampleSetConfig/fitSampleList" );

  // Disabling defined plots:
  LogInfo << "Removing defined plots..." << std::endl;
  ConfigUtils::clearEntry( cHandler.getConfig(), "fitterEngineConfig/likelihoodInterfaceConfig/dataSetManagerConfig/propagatorConfig/plotGeneratorConfig" );
  ConfigUtils::clearEntry( cHandler.getConfig(), "fitterEngineConfig/propagatorConfig/plotGeneratorConfig" );

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
  fitter.readConfig( GenericToolbox::Json::fetchValuePath<JsonType>( cHandler.getConfig(), "fitterEngineConfig" ) );

  DataSetManager& dataSetManager{fitter.getLikelihoodInterface().getDataSetManager()};

  // TODO Check
  //dataSetManager.getPropagator().setLoadAsimovData( false );

  // Disabling eigen decomposed parameters
  dataSetManager.getPropagator().setEnableEigenToOrigInPropagate( false );

  // --use-data-entry
  if( clParser.isOptionTriggered("useDataEntry") ){
    auto selectedDataEntry = clParser.getOptionVal<std::string>("useDataEntry", 0);
    // Do something better in case multiple datasets are defined
    bool isFound{false};
    for( auto& dataSet : dataSetManager.getDataSetList() ){
      if( GenericToolbox::isIn( selectedDataEntry, dataSet.getDataDispenserDict() ) ){
        LogWarning << "Using data entry \"" << selectedDataEntry << "\" for dataset: " << dataSet.getName() << std::endl;
        dataSet.setSelectedDataEntry( selectedDataEntry );
        isFound = true;
      } 
    }   
    LogThrowIf(not isFound, "Could not find data entry \"" << selectedDataEntry << "\" among defined data sets");
  }  

  // Sample binning using parameterSetName
  for( auto& sample : dataSetManager.getPropagator().getSampleSet().getSampleList() ){

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
        dataSetManager.getPropagator().getDialCollectionList().begin(),
        dataSetManager.getPropagator().getDialCollectionList().end(),
        [&](const DialCollection& dialCollection_){
          auto* parSetPtr{dialCollection_.getSupervisedParameterSet()};
          if( parSetPtr == nullptr ){ return false; }
          return ( parSetPtr->getName() == associatedParSet );
        });
    LogThrowIf(
        foundDialCollection == dataSetManager.getPropagator().getDialCollectionList().end(),
        "Could not find " << associatedParSet << " among fit dial collections: "
                          << GenericToolbox::toString(dataSetManager.getPropagator().getDialCollectionList(),
                                                      [](const DialCollection& dialCollection_){
                                                        return dialCollection_.getTitle();
                                                      }
                          ));

    LogThrowIf(foundDialCollection->getDialBinSet().getBinList().empty(), "Could not find binning");
    sample.setBinningFilePath( foundDialCollection->getDialBinSet().getFilePath() );

  }

  // Load everything
  dataSetManager.initialize();

  Propagator& propagator{dataSetManager.getPropagator()};


  if( clParser.isOptionTriggered("dryRun") ){
    std::cout << cHandler.toString() << std::endl;

    LogAlert << "Exiting as dry-run is set." << std::endl;
    return EXIT_SUCCESS;
  }


  if( not clParser.isOptionTriggered("usePreFit") and fitterRootFile != nullptr ){

    // Load post-fit parameters as "prior" so we can reset the weight to this point when throwing toys
    LogWarning << std::endl << GenericToolbox::addUpDownBars("Injecting post-fit parameters...") << std::endl;
    ObjectReader::readObject<TNamed>( fitterRootFile.get(), "FitterEngine/postFit/parState_TNamed", [&](TNamed* parState_){
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
    ObjectReader::readObject<TH2D>(
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
    std::vector<std::pair<std::string, std::string>> appendixDict{
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
    leafNameList.reserve( sample.getMcContainer().getHistogram().nBins );
    for( int iBin = 0 ; iBin < sample.getMcContainer().getHistogram().nBins; iBin++ ){
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
        sample.getMcContainer().getHistogram().nBins,
        0,
        sample.getMcContainer().getHistogram().nBins
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
      for( int iBin = 0 ; iBin < toyData.samplePtr->getMcContainer().getHistogram().nBins ; iBin++ ){
        double binData{ toyData.samplePtr->getMcContainer().getHistogram().binList[iBin].content };

        // bin volume
        auto& bin = toyData.samplePtr->getBinning().getBinList()[iBin];
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

  std::stringstream ss; ss << LogWarning.getPrefixString() << "Generating " << nToys << " toys...";
  for( int iToy = 0 ; iToy < nToys ; iToy++ ){

    // loading...
    GenericToolbox::displayProgressBar( iToy+1, nToys, ss.str() );

    // Do the throwing:
    propagator.getParametersManager().throwParametersFromGlobalCovariance();
    propagator.propagateParameters();

    if( enableStatThrowInToys ){
      for( auto& toyData : ToyDataList ){
        if( enableEventMcThrow ){
          // Take into account the finite amount of event in MC
          toyData.samplePtr->getMcContainer().throwEventMcError();
        }
        // JK: throwEventMcError() already treated each MC entry as Poisson event
        //     we don't want to throw again on weighted MC histogram
        // Asimov bin content -> toy data
        //toyData.samplePtr->getMcContainer().throwStatError();
      }
    }

    writeBinDataFct();

    // Write the branches
    toyThrowTree->Fill();
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

    for( int iBin = 0 ; iBin < toyData.samplePtr->getMcContainer().getHistogram().nBins ; iBin++ ){
      iBinGlobal++;

      std::string binTitle = toyData.samplePtr->getBinning().getBinList()[iBin].getSummary();
      double binVolume = toyData.samplePtr->getBinning().getBinList()[iBin].getVolume();

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
  propagator.getPlotGenerator().generateSampleHistograms( GenericToolbox::mkdirTFile(toyGenDir, "plots/histograms") );

  for( auto& histHolder : propagator.getPlotGenerator().getHistHolderList(0) ){
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

  propagator.getPlotGenerator().generateCanvas(
      propagator.getPlotGenerator().getHistHolderList(0),
      GenericToolbox::mkdirTFile(toyGenDir, "plots/canvas")
  );


  LogInfo << "Writing event samples in TTrees..." << std::endl;
  dataSetManager.getTreeWriter().writeSamples(
      GenericToolbox::mkdirTFile(toyGenDir, "events"),
      dataSetManager.getPropagator()
  );

}
