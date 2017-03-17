#include "FindCBESumThreshold.h"

using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace std;




FindCBESumThreshold::FindCBESumThreshold(const string& name, OptionsPtr opts) :
    Physics(name, opts),
    promptrandom(*ExpConfig::Setup::GetLastFound()),
    fit_model(utils::UncertaintyModels::Interpolated::makeAndLoad()),
    fitter("KinFit", 4, fit_model, true) // enable Z vertex by default
{
    fitter.SetZVertexSigma(0); // use unmeasured z vertex

    steps = HistFac.makeTH1D("Steps","","#",BinSettings(15),"steps");

    const AxisSettings axis_CBESum{"CBESum / MeV", {1600, 0, 1600}};

    h_CBESum_raw = HistFac.makeTH1D("CBESum raw ",axis_CBESum,"h_CBESum_raw");
    h_CBESum_pr  = HistFac.makeTH1D("CBESum raw prompt-random subtracted",axis_CBESum,"h_CBESum_pr");

}

void FindCBESumThreshold::ProcessEvent(const TEvent& event, manager_t&)
{
    steps->Fill("Seen",1);

    if(!triggersimu.ProcessEvent(event)) {
        steps->Fill("TriggerSimu failed", 1.0);
        return;
    }

    h_CBESum_raw->Fill(triggersimu.GetCBEnergySum());

    const TEventData& data = event.Reconstructed();

    for(const TTaggerHit& taggerhit : data.TaggerHits) {

        steps->Fill("Seen taggerhits",1.0);

        promptrandom.SetTaggerHit(taggerhit.Time);
        if(promptrandom.State() == PromptRandom::Case::Outside)
            continue;

        steps->Fill("Acc taggerhits",1.0);

        h_CBESum_pr->Fill(triggersimu.GetCBEnergySum(), promptrandom.FillWeight());

        const auto& cands = data.Candidates;
        if(cands.size() != 5)
            return;
        steps->Fill("nCands==5",1);



    }
}

void FindCBESumThreshold::ShowResult()
{
    canvas(GetName())
            << steps
            << h_CBESum_raw << h_CBESum_pr
            << endc;
}

AUTO_REGISTER_PHYSICS(FindCBESumThreshold)