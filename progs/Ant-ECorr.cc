#include "calibration/DataManager.h"
#include "expconfig/ExpConfig.h"

#include "base/Logger.h"
#include "base/CmdLine.h"
#include "base/std_ext/string.h"
#include "base/OptionsList.h"
#include "base/WrapTFile.h"
#include "TH2.h"
#include "TH3.h"
#include "root-addons/analysis_codes/TwoPi0_MCSmearing_tools.h"
#include "TROOT.h"
#include "TRint.h"
#include "tree/TCalibrationData.h"
#include "calibration/modules/detail/TH2Storage.h"
#include "base/TH_ext.h"
#include <iostream>
#include <cstring>
#include "TDirectory.h"
#include "base/std_ext/string.h"
#include "base/TH_ext.h"
#include <array>
#include "analysis/plot/root_draw.h"
#include "base/std_ext/math.h"
#include "analysis/plot/HistogramFactory.h"
#include "base/Array2D.h"

using namespace std;
using namespace ant;
using namespace  ant::calibration::gui;

TH2* GetHist(WrapTFileInput& file, const string& histname) {
    TH2* h = nullptr;
    file.GetObject(histname, h);
    if(!h) {
        LOG(ERROR) << histname << " not found in " << file.FileNames();
        exit(EXIT_FAILURE);
    }
    return h;
}

inline string ToUpper(const string& s) {
    string out = s;
    transform(s.begin(), s.end(),out.begin(), ::toupper);
    return out;
}

inline string ToLower(const string& s) {
    string out = s;
    transform(s.begin(), s.end(),out.begin(), ::tolower);
    return out;
}

int main(int argc, char** argv) {
    SetupLogger();


    TCLAP::CmdLine cmd("Ant-mcsmearing", ' ', "0.1");
    auto cmd_verbose   = cmd.add<TCLAP::ValueArg<int>>   ("v","verbose", "Verbosity level (0..9)", false, 0,"level");
    auto cmd_batchmode = cmd.add<TCLAP::SwitchArg>       ("b","batch",   "Run in batch mode (no GUI, autosave)",false);
    auto cmd_setupname = cmd.add<TCLAP::ValueArg<string>>("s","setup",   "Setup name",       true, "", "setup");
    auto cmd_detector  = cmd.add<TCLAP::ValueArg<string>>("" ,"detector","Detector Name",    true, "", "detector");
    auto cmd_file      = cmd.add<TCLAP::ValueArg<string>>("" ,"file",    "Input file",       true, "", "file");

    const bool SaveToDatabase = false; // @todo: read from cmdline

    cmd.parse(argc, argv);

    if(cmd_verbose->isSet())
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());

    const auto det = ToUpper(cmd_detector->getValue());

    // create TRint app early in order to have valid gStyle pointer...
    int fake_argc=1;
    char* fake_argv[2];
    fake_argv[0] = argv[0];
    if(cmd_batchmode->isSet()) {
        fake_argv[fake_argc++] = strdup("-q");
    }
    auto app = new TRint("Ant-mcsmearing",&fake_argc,fake_argv,nullptr,0,true);


    const string histname = std_ext::formatter() << "MCClusterECorr/" << det << "/h_EtrueErec";
    const string statsHistName = std_ext::formatter() << "MCClusterECorr/" << det << "/h_nFills";

    WrapTFileInput infile(cmd_file->getValue());
    const auto h_ecorr  = GetHist(infile, histname);
    const auto h_nfills = GetHist(infile, statsHistName);

    const auto range = interval<double>::CenterWidth(1.0,.15);



    std_ext::IQR iqr;

    for(int x=1;x<=h_ecorr->GetNbinsX();++x) {
        for(int y=1;y<=h_ecorr->GetNbinsY();++y) {
             if(h_nfills->GetBinContent(x,y) > 500.0) {
                 const auto v = h_ecorr->GetBinContent(x,y);
                 iqr.Add(v);
             }
        }
    }

    analysis::HistogramFactory f("ECorr");

    auto h = f.makeTH1D("Factors", "ECorr Factor","", BinSettings(50, interval<double>::CenterWidth(iqr.GetMedian(), iqr.GetIQRStdDev()*3.0)));

    for(int x=1;x<=h_ecorr->GetNbinsX();++x) {
        for(int y=1;y<=h_ecorr->GetNbinsY();++y) {
             if(h_nfills->GetBinContent(x,y) > 500.0) {
                 const auto v = h_ecorr->GetBinContent(x,y);
                     h->Fill(v);
             }
        }
    }

    const auto norm = h->GetXaxis()->GetBinCenter(h->GetMaximumBin()) - 1.0;
    LOG(INFO) << norm;

    auto processed = static_cast<TH2D*>(TH_ext::Apply(h_ecorr, h_nfills, [range, norm] (const double ecorr, const double n) {
        return n > 500.0 ? range.Clip(ecorr-norm) : std_ext::NaN;
    }));

    const auto zrange = TH_ext::GetZMinMax(processed);
    processed->SetMinimum(zrange.Start());
    processed->SetMaximum(zrange.Stop());

    auto filled = TH_ext::Clone(processed, "ECorrFilled");

    Array2D_TH2D a(filled);

    a.FloodFillAverages();


    canvas("Processed") << drawoption("colz") << h << processed << filled << endc;


    shared_ptr<ExpConfig::Setup> setup = nullptr;
    shared_ptr<calibration::DataManager> manager = nullptr;

    if(SaveToDatabase) {
        const auto setup_name = cmd_setupname->getValue();
        setup = ExpConfig::Setup::Get(setup_name);
        if(setup == nullptr) {
            LOG(ERROR) << "Did not find setup instance for name " << setup_name;
            return 1;
        }

        manager = setup->GetCalibrationDataManager();
        manager->SetOverrideToDefault(true);

        const auto id = TID(0,0,{TID::Flags_t::MC});


        const string calName = std_ext::formatter() << det << "_ClusterECorr";

        TCalibrationData cdata(calName, id, id);
        calibration::detail::TH2Storage::Encode(processed, cdata);

        manager->Add(cdata,  Calibration::AddMode_t::AsDefault);
    }


    app->Run(kTRUE);
    ExpConfig::Setup::Cleanup();
    setup = nullptr;
    manager = nullptr;

    return EXIT_SUCCESS;
}
