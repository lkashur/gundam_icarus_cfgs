
#include "GundamGlobals.h"
#include "GundamApp.h"
#include "GundamUtils.h"
#include "Propagator.h"
#include "ConfigUtils.h"

#include "Logger.h"
#include "CmdLineParser.h"
#include "GenericToolbox.h"
#include "GenericToolbox.Json.h"
#include "GenericToolbox.Root.h"
#include "GenericToolbox.RawDataArray.h"

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

  GundamApp app{"cross-section calculator tool"};

  // --------------------------
  // Read Command Line Args:
  // --------------------------
  CmdLineParser clParser;

  clParser.addDummyOption("Main options:");
  clParser.addOption("configFile", {"-c", "--config-file"}, "Specify path to the fitter config file");
  clParser.addOption("fitterFile", {"-f"}, "Specify the fitter output file");
  clParser.addOption("outputFile", {"-o", "--out-file"}, "Specify the CalcXsec output file");
  clParser.addOption("nbThreads", {"-t", "--nb-threads"}, "Specify nb of parallel threads");
  clParser.addOption("nToys", {"-n"}, "Specify number of toys");
  clParser.addOption("randomSeed", {"-s", "--seed"}, "Set random seed");
  clParser.addOption("useDataEntry", {"--use-data-entry"}, "Overrides \"selectedDataEntry\" in dataSet config. Second arg is to select a given dataset");
  clParser.addOption("fitSampleSetConfig", {"--fitsample-config"}, "override fitSampleSetConfig");
  clParser.addOption("plotGeneratorConfig", {"--plot-config"}, "override plotGeneratorConfig");

  clParser.addDummyOption("Trigger options:");
  clParser.addTriggerOption("dryRun", {"-d", "--dry-run"}, "Only overrides fitter config and print it.");
  clParser.addTriggerOption("useBfAsXsec", {"--use-bf-as-xsec"}, "Use best-fit as x-sec value instead of mean of toys.");
  clParser.addTriggerOption("usePreFit", {"--use-prefit"}, "Use prefit covariance matrices for the toy throws.");

  LogInfo << "Usage: " << std::endl;
  LogInfo << clParser.getConfigSummary() << std::endl << std::endl;

  clParser.parseCmdLine(argc, argv);

  LogThrowIf(clParser.isNoOptionTriggered(), "No option was provided.");

  LogInfo << "Provided arguments: " << std::endl;
  LogInfo << clParser.getValueSummary() << std::endl << std::endl;


  // Sanity checks
  LogThrowIf(not clParser.isOptionTriggered("configFile"), "Xsec calculator config file not provided.");
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
  nlohmann::json fitterConfig; // will be used to load the propagator

  if( GenericToolbox::doesFilePathHasExtension(fitterFile, "root") ){
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
  ConfigUtils::applyOverrides(
      cHandler.getConfig(),
      GenericToolbox::Json::readConfigJsonStr(R"({"fitterEngineConfig":{"propagatorConfig":{"fitSampleSetConfig":{"fitSampleList":[]}}}})")
  );

  // Disabling defined plots:
  LogInfo << "Removing defined plots..." << std::endl;
  ConfigUtils::applyOverrides(
      cHandler.getConfig(),
      GenericToolbox::Json::readConfigJsonStr(R"({"fitterEngineConfig":{"propagatorConfig":{"plotGeneratorConfig":{}}}})")
  );

  // Defining signal samples
  nlohmann::json xsecConfig{ ConfigUtils::readConfigFile( clParser.getOptionVal<std::string>("configFile") ) };
  cHandler.override( xsecConfig );

  if( clParser.isOptionTriggered("fitSampleSetConfig") ){
    nlohmann::json fitSampleSetConfig_new{ ConfigUtils::readConfigFile( clParser.getOptionVal<std::string>("fitSampleSetConfig") ) };
    cHandler.getConfig()["fitterEngineConfig"]["propagatorConfig"]["fitSampleSetConfig"] = fitSampleSetConfig_new;
  }
  if( clParser.isOptionTriggered("plotGeneratorConfig") ){
    nlohmann::json plotGeneratorConfig_new{ ConfigUtils::readConfigFile( clParser.getOptionVal<std::string>("plotGeneratorConfig") ) };
    cHandler.getConfig()["fitterEngineConfig"]["propagatorConfig"]["plotGeneratorConfig"] = plotGeneratorConfig_new;
  }

  LogInfo << "Override done." << std::endl;


  LogInfo << "Fetching propagator config into fitter config..." << std::endl;
  auto configPropagator = GenericToolbox::Json::fetchValuePath<nlohmann::json>( cHandler.getConfig(), "fitterEngineConfig/propagatorConfig" );

  // Create a propagator object
  Propagator propagator;

  // Read the whole fitter config with the override parameters
  LogWarning << std::endl << GenericToolbox::addUpDownBars("Reading propagator config...") << std::endl;
  propagator.readConfig( configPropagator );

  // --use-data-entry
  if( clParser.isOptionTriggered("useDataEntry") ){
    auto selectedDataEntry = clParser.getOptionVal<std::string>("useDataEntry", 0);
    // Do something better in case multiple datasets are defined
    bool isFound{false};
    for( auto& dataSet : propagator.getDataSetList() ){
      if( GenericToolbox::doesKeyIsInMap( selectedDataEntry, dataSet.getDataDispenserDict() ) ){
        LogWarning << "Using data entry \"" << selectedDataEntry << "\" for dataset: " << dataSet.getName() << std::endl;
        dataSet.setSelectedDataEntry( selectedDataEntry );
        isFound = true;
      }
    }
    LogThrowIf(not isFound, "Could not find data entry \"" << selectedDataEntry << "\" among defined data sets");
  }

  // We are only interested in our MC. Data has already been used to get the post-fit error/values
  propagator.setLoadAsimovData( false );

  // Disabling eigen decomposed parameters
  propagator.setEnableEigenToOrigInPropagate( false );

  // Sample binning using parameterSetName
  for( auto& sample : propagator.getFitSampleSet().getFitSampleList() ){

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
        propagator.getDialCollections().begin(), propagator.getDialCollections().end(),
        [&](const DialCollection& dialCollection_){
          auto* parSetPtr{dialCollection_.getSupervisedParameterSet()};
          if( parSetPtr == nullptr ){ return false; }
          return ( parSetPtr->getName() == associatedParSet );
        });
    LogThrowIf(
        foundDialCollection == propagator.getDialCollections().end(),
        "Could not find " << associatedParSet << " among fit dial collections: "
                          << GenericToolbox::iterableToString(propagator.getDialCollections(), [](const DialCollection& dialCollection_){
                                                                return dialCollection_.getTitle();
                                                              }
                          ));

    LogThrowIf(foundDialCollection->getDialBinSet().getBinList().empty(), "Could not find binning");
    sample.setBinningFilePath( foundDialCollection->getDialBinSet().getFilePath() );

  }

  // Load everything
  propagator.initialize();

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
        {"configFile", "%s"},
        {"fitterFile", "Fit_%s"},
        {"nToys", "nToys_%s"},
        {"randomSeed", "Seed_%s"},
        {"usePreFit", "PreFit"},
    };

    outFilePath = "ToyGeneration_" + GundamUtils::generateFileName(clParser, appendixDict) + ".root";

    std::string outFolder{GenericToolbox::Json::fetchValue<std::string>(xsecConfig, "outputFolder", "./")};
    outFilePath = GenericToolbox::joinPath(outFolder, outFilePath);
  }

  app.setCmdLinePtr( &clParser );
  app.setConfigString( ConfigUtils::ConfigHandler{xsecConfig}.toString() );
  app.openOutputFile( outFilePath );
  app.writeAppInfo();

  auto* calcXsecDir{ GenericToolbox::mkdirTFile(app.getOutfilePtr(), "calcXsec") };
  bool useBestFitAsCentralValue{
    clParser.isOptionTriggered("useBfAsXsec")
    or GenericToolbox::Json::fetchValue<bool>(xsecConfig, "useBestFitAsCentralValue", false)
  };

  LogInfo << "Creating throws tree" << std::endl;
  auto* xsecThrowTree = new TTree("xsecThrow", "xsecThrow");
  xsecThrowTree->SetDirectory( GenericToolbox::mkdirTFile(calcXsecDir, "throws") ); // temp saves will be done here

  auto* xsecAtBestFitTree = new TTree("xsecAtBestFitTree", "xsecAtBestFitTree");
  xsecAtBestFitTree->SetDirectory( GenericToolbox::mkdirTFile(calcXsecDir, "throws") ); // temp saves will be done here

  struct CrossSectionData{
    Sample* samplePtr{nullptr};
    nlohmann::json config{};
    GenericToolbox::RawDataArray branchBinsData{};

    TH1D histogram{};
  };
  std::vector<CrossSectionData> crossSectionDataList{};

  LogInfo << "Initializing xsec samples..." << std::endl;
  crossSectionDataList.reserve( propagator.getFitSampleSet().getFitSampleList().size() );
  for( auto& sample : propagator.getFitSampleSet().getFitSampleList() ){
    crossSectionDataList.emplace_back();
    auto& xsecEntry = crossSectionDataList.back();

    LogScopeIndent;
    LogInfo << "Defining xsec entry: " << sample.getName() << std::endl;
    xsecEntry.samplePtr = &sample;
    xsecEntry.config = sample.getConfig();
    xsecEntry.branchBinsData.resetCurrentByteOffset();
    std::vector<std::string> leafNameList{};
    leafNameList.reserve( sample.getMcContainer().histogram->GetNbinsX() );
    for( int iBin = 0 ; iBin < sample.getMcContainer().histogram->GetNbinsX() ; iBin++ ){
      leafNameList.emplace_back(Form("bin_%i/D", iBin));
      xsecEntry.branchBinsData.writeRawData( double(0) );
    }
    xsecEntry.branchBinsData.lockArraySize();

    xsecThrowTree->Branch(
        GenericToolbox::generateCleanBranchName( sample.getName() ).c_str(),
        xsecEntry.branchBinsData.getRawDataArray().data(),
        GenericToolbox::joinVectorString(leafNameList, ":").c_str()
    );
    xsecAtBestFitTree->Branch(
        GenericToolbox::generateCleanBranchName( sample.getName() ).c_str(),
        xsecEntry.branchBinsData.getRawDataArray().data(),
        GenericToolbox::joinVectorString(leafNameList, ":").c_str()
    );

    xsecEntry.histogram = TH1D(
        sample.getName().c_str(),
        sample.getName().c_str(),
        sample.getMcContainer().histogram->GetNbinsX(),
        0,
        sample.getMcContainer().histogram->GetNbinsX()
    );
  }

  int nToys{ clParser.getOptionVal<int>("nToys") };

  bool enableEventMcThrow{true};
  bool enableStatThrowInToys{true};
  auto xsecCalcConfig   = GenericToolbox::Json::fetchValue( cHandler.getConfig(), "xsecCalcConfig", nlohmann::json() );
  enableStatThrowInToys = GenericToolbox::Json::fetchValue( xsecCalcConfig, "enableStatThrowInToys", enableStatThrowInToys);
  enableEventMcThrow    = GenericToolbox::Json::fetchValue( xsecCalcConfig, "enableEventMcThrow", enableEventMcThrow);

  auto writeBinDataFct = std::function<void()>([&]{
    for( auto& xsec : crossSectionDataList ){

      xsec.branchBinsData.resetCurrentByteOffset();
      for( int iBin = 0 ; iBin < xsec.samplePtr->getMcContainer().histogram->GetNbinsX() ; iBin++ ){
        double binData{ xsec.samplePtr->getMcContainer().histogram->GetBinContent(1+iBin) };

        // bin volume
        auto& bin = xsec.samplePtr->getMcContainer().binning.getBinList()[iBin];
        double binVolume{1};

        for( auto& edges : bin.getEdgesList() ){
          if( edges.isConditionVar ){ continue; } // no volume, just a condition variable

          binVolume *= (edges.max - edges.min);
        }

        binData /= binVolume;
        xsec.branchBinsData.writeRawData( binData );
      }
    }
  });

  {
    LogWarning << "Calculating weight at best-fit" << std::endl;
    for( auto& parSet : propagator.getParametersManager().getParameterSetsList() ){ parSet.moveFitParametersToPrior(); }
    propagator.propagateParametersOnSamples();
    writeBinDataFct();
    xsecAtBestFitTree->Fill();
    GenericToolbox::writeInTFile( GenericToolbox::mkdirTFile(calcXsecDir, "throws"), xsecAtBestFitTree );
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
    propagator.propagateParametersOnSamples();

    if( enableStatThrowInToys ){
      for( auto& xsec : crossSectionDataList ){
        if( enableEventMcThrow ){
          // Take into account the finite amount of event in MC
          xsec.samplePtr->getMcContainer().throwEventMcError();
        }
        // JK: throwEventMcError() already treated each MC entry as Poisson event
        //     we don't want to throw again on weighted MC histogram
        // Asimov bin content -> toy data
        //xsec.samplePtr->getMcContainer().throwStatError();
      }
    }

    writeBinDataFct();

    // Write the branches
    xsecThrowTree->Fill();
  }


  LogInfo << "Writing throws..." << std::endl;
  GenericToolbox::writeInTFile( GenericToolbox::mkdirTFile(calcXsecDir, "throws"), xsecThrowTree );

  LogInfo << "Calculating mean & covariance matrix..." << std::endl;
  auto* meanValuesVector = GenericToolbox::generateMeanVectorOfTree(
      useBestFitAsCentralValue ? xsecAtBestFitTree : xsecThrowTree
  );
  auto* globalCovMatrix = GenericToolbox::generateCovarianceMatrixOfTree( xsecThrowTree );

  auto* globalCovMatrixHist = GenericToolbox::convertTMatrixDtoTH2D(globalCovMatrix);
  auto* globalCorMatrixHist = GenericToolbox::convertTMatrixDtoTH2D(GenericToolbox::convertToCorrelationMatrix(globalCovMatrix));

  std::vector<TH1D> binValues{};
  binValues.reserve( propagator.getFitSampleSet().getFitSampleList().size() );
  int iBinGlobal{-1};

  for( auto& xsec : crossSectionDataList ){

    for( int iBin = 0 ; iBin < xsec.samplePtr->getMcContainer().histogram->GetNbinsX() ; iBin++ ){
      iBinGlobal++;

      std::string binTitle = xsec.samplePtr->getBinning().getBinList()[iBin].getSummary();
      double binVolume = xsec.samplePtr->getBinning().getBinList()[iBin].getVolume();

      xsec.histogram.SetBinContent( 1+iBin, (*meanValuesVector)[iBinGlobal] );
      xsec.histogram.SetBinError( 1+iBin, TMath::Sqrt( (*globalCovMatrix)[iBinGlobal][iBinGlobal] ) );
      xsec.histogram.GetXaxis()->SetBinLabel( 1+iBin, binTitle.c_str() );

      globalCovMatrixHist->GetXaxis()->SetBinLabel(1+iBinGlobal, GenericToolbox::joinPath(xsec.samplePtr->getName(), binTitle).c_str());
      globalCorMatrixHist->GetXaxis()->SetBinLabel(1+iBinGlobal, GenericToolbox::joinPath(xsec.samplePtr->getName(), binTitle).c_str());
      globalCovMatrixHist->GetYaxis()->SetBinLabel(1+iBinGlobal, GenericToolbox::joinPath(xsec.samplePtr->getName(), binTitle).c_str());
      globalCorMatrixHist->GetYaxis()->SetBinLabel(1+iBinGlobal, GenericToolbox::joinPath(xsec.samplePtr->getName(), binTitle).c_str());
    }

    xsec.histogram.SetMarkerStyle(kFullDotLarge);
    xsec.histogram.SetMarkerColor(kGreen-3);
    xsec.histogram.SetMarkerSize(0.5);
    xsec.histogram.SetLineWidth(2);
    xsec.histogram.SetLineColor(kGreen-3);
    xsec.histogram.SetDrawOption("E1");
    xsec.histogram.GetXaxis()->LabelsOption("v");
    xsec.histogram.GetXaxis()->SetLabelSize(0.02);
    xsec.histogram.GetYaxis()->SetTitle( GenericToolbox::Json::fetchValue(xsec.samplePtr->getConfig(), "yAxis", "#delta#sigma").c_str() );

    GenericToolbox::writeInTFile(
        GenericToolbox::mkdirTFile(calcXsecDir, "histograms"),
        &xsec.histogram, GenericToolbox::generateCleanBranchName( xsec.samplePtr->getName() )
    );

  }

  globalCovMatrixHist->GetXaxis()->SetLabelSize(0.02);
  globalCovMatrixHist->GetYaxis()->SetLabelSize(0.02);
  GenericToolbox::writeInTFile(GenericToolbox::mkdirTFile(calcXsecDir, "matrices"), globalCovMatrixHist, "covarianceMatrix");

  globalCorMatrixHist->GetXaxis()->SetLabelSize(0.02);
  globalCorMatrixHist->GetYaxis()->SetLabelSize(0.02);
  globalCorMatrixHist->GetZaxis()->SetRangeUser(-1, 1);
  GenericToolbox::writeInTFile(GenericToolbox::mkdirTFile(calcXsecDir, "matrices"), globalCorMatrixHist, "correlationMatrix");

  LogInfo << "Generating xsec sample plots..." << std::endl;

  for( auto& parSet : propagator.getParametersManager().getParameterSetsList() ){ parSet.moveFitParametersToPrior(); }
  propagator.propagateParametersOnSamples();

  // manual trigger to tweak the error bars
  propagator.getPlotGenerator().generateSampleHistograms( GenericToolbox::mkdirTFile(calcXsecDir, "plots/histograms") );

  for( auto& histHolder : propagator.getPlotGenerator().getHistHolderList(0) ){
    if( not histHolder.isData ){ continue; } // only data will print errors

    const CrossSectionData* xsecDataPtr{nullptr};
    for( auto& xsecData : crossSectionDataList ){
      if( xsecData.samplePtr  == histHolder.fitSamplePtr){
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
      GenericToolbox::mkdirTFile(calcXsecDir, "plots/canvas")
  );


  LogInfo << "Writing event samples in TTrees..." << std::endl;
  propagator.getTreeWriter().writeSamples( GenericToolbox::mkdirTFile(calcXsecDir, "events") );

}
