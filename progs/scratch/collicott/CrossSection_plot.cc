#include "base/CmdLine.h"
#include "base/Logger.h"

#include "base/WrapTFile.h"

#include "tree/TEvent.h"
#include "tree/TEventData.h"

#include "TTree.h"
#include "TRint.h"
#include "TH1D.h"
#include "TH2D.h"
#include "base/std_ext/system.h"

#include "analysis/plot/HistogramFactory.h"
#include "analysis/plot/root_draw.h"
#include "analysis/physics/Physics.h"
#include "analysis/physics/scratch/collicott/gp_pi0[2g]p.h"
#include "analysis/physics/scratch/collicott/cross_section.h"
#include "analysis/physics/scratch/collicott/det_eff.h"

#include <limits>

using namespace ant;
using namespace std;
using namespace ant::analysis;
volatile bool interrupt = false;


int main( int argc, char** argv )
{
    SetupLogger();

    signal(SIGINT, [] (int) {
        cout << ">>> Interrupted" << endl;
        interrupt = true;
    });

    TCLAP::CmdLine cmd("cross_plot", ' ', "0.1");

    // ValueArg == set once
    // MultiArg == set multipl
    auto cmd_input = cmd.add<TCLAP::ValueArg<string>>("i","input","Input file",true,"","input");
    auto cmd_batchmode = cmd.add<TCLAP::MultiSwitchArg>("b","batch","Run in batch mode (no ROOT shell afterwards)",false);
    auto cmd_maxevents = cmd.add<TCLAP::MultiArg<int>>("m","maxevents","Process only max events",false,"maxevents");
    auto cmd_output = cmd.add<TCLAP::ValueArg<string>>("o","output","Output file",false,"","filename");
    //    auto cmd_setupname = cmd.add<TCLAP::ValueArg<string>>("s","setup","Override setup name", false, "Setup_2014_07_EPT_Prod", "setup");

    auto cmd_physics  = cmd.add<TCLAP::ValueArg<string>>("p","physics","Physics Class",true,"","physics");
    auto cmd_bintheta = cmd.add<TCLAP::ValueArg<int>>("","bintheta","Width of theta bins",false,5,"theta_bin_width");
    auto cmd_binphi   = cmd.add<TCLAP::ValueArg<int>>("","binphi",  "Width of phi bins",  false,10,"phi_bin_width");
    auto cmd_Tagger      = cmd.add<TCLAP::SwitchArg>("","Tagger",  "Set for tagger data");

    cmd.parse(argc, argv);

    WrapTFileInput input(cmd_input->getValue());

    // Yield tree!
    TTree* tyield = nullptr;
    if(!input.GetObject(cmd_physics->getValue()+"/CrossSection/Yield",tyield))
        LOG(ERROR) << "Cannot find yield tree in input file";
    analysis::utils::scratch_collicott_CrossSection::Yield_t yield;
    yield.LinkBranches(tyield);

    // Scaler tree!
    TTree* tscaler = nullptr;
    if(!input.GetObject(cmd_physics->getValue()+"/CrossSection/Scalers",tscaler))
        LOG(ERROR) << "Cannot find scaler tree in input file";
    analysis::utils::scratch_collicott_CrossSection::Scalers_t scalers;
    scalers.LinkBranches(tscaler);

    // DetEff trees!
    TTree* tdeteff_total = nullptr;
    if(!input.GetObject(cmd_physics->getValue()+"/DetectionEfficiency/Sig_Total",tdeteff_total))
        LOG(ERROR) << "Cannot find scaler tree in input file";
    analysis::utils::scratch_collicott_DetEff::DetEff_t deteff_total;
    deteff_total.LinkBranches(tdeteff_total);

    TTree* tdeteff_accept = nullptr;
    if(!input.GetObject(cmd_physics->getValue()+"/DetectionEfficiency/Sig_Accepted",tdeteff_accept))
        LOG(ERROR) << "Cannot find scaler tree in input file";
    analysis::utils::scratch_collicott_DetEff::DetEff_t deteff_accept;
    deteff_accept.LinkBranches(tdeteff_accept);

    unique_ptr<WrapTFileOutput> masterFile;
    if(cmd_output->isSet()) {
        masterFile = std_ext::make_unique<WrapTFileOutput>(cmd_output->getValue(),
                                                    WrapTFileOutput::mode_t::recreate,
                                                     true); // cd into masterFile upon creation
    }




    HistogramFactory HistFac("Plots");

    auto ntheta_bins = uint(180/cmd_bintheta->getValue());
    auto nphi_bins   = uint(360/cmd_binphi->getValue());
    auto nTagger = uint(47); if(cmd_Tagger->isSet()) nTagger = 352;

    LOG(INFO) << "Using theta bin width of " << (180.0/ntheta_bins) << " degrees.";
    LOG(INFO) << "Using phi   bin width of " << (360.0/nphi_bins) << " degrees.";

    BinSettings bin_theta(ntheta_bins,0,180);
    BinSettings bin_phi(nphi_bins,-180,180);
    BinSettings bin_tagger(nTagger);

    // 2D hists
    auto h_yield         = HistFac.makeTH2D("h_yield",         "tc","theta",bin_tagger,bin_theta,"h_yield", true);
    auto h_deteff_total  = HistFac.makeTH2D("h_deteff_total",  "tc","theta",bin_tagger,bin_theta,"h_deteff_total",true);
    auto h_deteff_accept = HistFac.makeTH2D("h_deteff_accept", "tc","theta",bin_tagger,bin_theta,"h_deteff_accept",true);

    // 1D hist
    auto h_flux          = HistFac.makeTH1D("h_flux",          "tc","Flux",     bin_tagger,    "h_flux",        true);
    auto h_taggeff       = HistFac.makeTH1D("h_taggeff",       "tc","TaggEff",  bin_tagger,    "h_taggeff",     true);
    auto h_livetime      = HistFac.makeTH1D("h_livetime",      "tc","Livetime", bin_tagger,    "h_livetime",    true);
    auto h_flux_corr     = HistFac.makeTH1D("h_flux_corr",     "tc","Flux Corrected for livetime and taggeff",     bin_tagger,    "h_flux_corr",   true);

//    auto t_flux          = HistFac.makeTH1D("t_flux",          "tc","Flux",     bin_tagger,    "t_flux",        true);
//    auto t_taggeff       = HistFac.makeTH1D("t_taggeff",       "tc","TaggEff",  bin_tagger,    "t_taggeff",     true);

    for(long long entry=0;entry<scalers.Tree->GetEntries();entry++) {
        scalers.Tree->GetEntry(entry);

        h_livetime->Fill(scalers.exp_livetime);

        for(auto i = 0; i < scalers.tagger_scalers().size(); i++)
        {
//            t_taggeff->SetBinContent(i,scalers.tagger_eff().at(i));
//            t_taggeff->SetBinError(i,scalers.tagger_deff().at(i));
//            t_flux->SetBinContent(i,scalers.tagger_scalers().at(i),scalers.exp_livetime);

            h_flux->Fill(i,scalers.tagger_scalers().at(i));
            h_taggeff->Fill(i,scalers.tagger_eff().at(i));
            h_flux_corr->Fill(i,scalers.tagger_scalers().at(i)*scalers.tagger_eff().at(i)*scalers.exp_livetime);
        }

//        TH1D t_flux_corr = (*t_flux)*(*t_taggeff);
//        *h_flux_corr->Add(&t_flux_corr,1);

    }

    for(long long entry=0;entry<yield.Tree->GetEntries();entry++) {
        yield.Tree->GetEntry(entry);

        if (!yield.isMC)
            h_yield->Fill(yield.tc_channel,yield. sp_theta,  yield.tc_promptrandom);
    }

    for(long long entry=0;entry<deteff_total.Tree->GetEntries();entry++) {
        deteff_total.Tree->GetEntry(entry);

        string reaction = deteff_total.reaction;
        if(reaction.substr(0, 4) != "data")
            h_deteff_total->Fill(deteff_total.tc_channel,deteff_total.sp_theta, deteff_total.tc_promptrandom);
    }


    for(long long entry=0;entry<deteff_accept.Tree->GetEntries();entry++) {
        deteff_accept.Tree->GetEntry(entry);

        string reaction = deteff_accept.reaction;
        if(reaction.substr(0, 4) != "data")
            h_deteff_accept->Fill(deteff_accept.tc_channel,deteff_accept.sp_theta, deteff_accept.tc_promptrandom);
    }

    if(!cmd_batchmode->isSet()) {
        if(!std_ext::system::isInteractive()) {
            LOG(INFO) << "No TTY attached. Not starting ROOT shell.";
        }
        else{
            argc=0; // prevent TRint to parse any cmdline
            TRint app("Ant",&argc,argv,nullptr,0,true);

            canvas("cross_plot") << drawoption("colz")
                                       << h_yield << h_deteff_total << h_deteff_accept
                                       << endc;


            canvas("cross_plot") << drawoption("colz")
                                       << h_flux << h_taggeff << h_livetime << h_flux_corr
                                       << endc;

            if(masterFile)
                LOG(INFO) << "Stopped running, but close ROOT properly to write data to disk.";

            app.Run(kTRUE); // really important to return...

            if(masterFile)
                LOG(INFO) << "Writing output file...";
            masterFile = nullptr;   // and to destroy the master WrapTFile before TRint is destroyed
        }

    }
}
