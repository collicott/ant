#pragma once

#include "analysis/physics/Physics.h"

#include <string>

class TH1D;

namespace ant {
namespace analysis {
namespace physics {

class ParticleIDCheck : public Physics {
protected:

    struct branch_hists {
        branch_hists(SmartHistFactory& HistFac,const std::string& name);
        TH1D* hist;
        void Fill(const data::Event::Data& data);
    };

    branch_hists mctrue;
    branch_hists rec;

public:
    ParticleIDCheck(PhysOptPtr opts);

    void ProcessEvent(const data::Event &event) override;
    void Finish() override;
    void ShowResult() override;
};

}
}
}
