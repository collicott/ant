#include "Setup.h"

namespace ant {
namespace expconfig {
namespace setup {

/**
 * @brief Ant Setup for Rare Pion Test Beamtime in January 2015 & MC studies
 */
class Setup_2015_01_Pion : public Setup
{
    const bool MCTaggerHits;

public:

    Setup_2015_01_Pion(const std::string& name, OptionsPtr opt);

    virtual double GetElectronBeamEnergy() const override;

    virtual ExpConfig::Setup::candidatebuilder_config_t GetCandidateBuilderConfig() const override;

    bool Matches(const TID& tid) const override {
        if(!std_ext::time_between(tid.Timestamp, "2015-01-27", "2015-02-01"))
            return false;
        return true;
    }
};

}}} // namespace ant::expconfig::setup
