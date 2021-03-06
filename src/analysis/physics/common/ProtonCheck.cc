#include "physics/common/ProtonCheck.h"
#include <cmath>
#include <iostream>

#include "utils/particle_tools.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

ProtonCheck::ProtonCheck(const std::string& name,OptionsPtr opts):
    Physics(name,opts)
{
    const BinSettings e(1000);
    const BinSettings t(300,-15,15);
    const BinSettings dE(80,0,8);
    const BinSettings theta_bins(180,0,180);
    tof = HistFac.makeTH2D("Proton TOF TAPS","t [ns]","E [MeV]", t,e,"tof");
    tof_trueE = HistFac.makeTH2D("Proton TOF TAPS (MCTrue Energy)","t [ns]","E_{true} [MeV]", t,e,"tof_true");
    dEE = HistFac.makeTH2D("Proton dEE TAPS","E [MeV]","dE [MeV]", e,dE,"dEE");
    cand_mult = HistFac.makeTH1D("Candidates / Event","# Candiadates/Event","#",BinSettings(20),"mult");

    theta =  HistFac.makeTH1D("Theta","#Theta","#",BinSettings(180,0,180),"theta");

    theta_corr = HistFac.makeTH2D("Theta Corrleation","true","rec",theta_bins,theta_bins,"theta_corr");
}

void ProtonCheck::ProcessEvent(const TEvent& event, manager_t&)
{
    auto mctrue_particles = utils::ParticleTypeList::Make(event.MCTrue().ParticleTree);
    if(mctrue_particles.GetAll().size() == 1) {
        const auto& protons = mctrue_particles.Get(ParticleTypeDatabase::Proton);

        if(protons.size()==1) {
            const auto& mctrue = protons.at(0);

            if(mctrue->Theta() < std_ext::degree_to_radian(20.0)) {

                for(const auto& cand : event.Reconstructed().Candidates) {

                    if(cand.Detector & Detector_t::Any_t::TAPS_Apparatus) {
                        tof->Fill(cand.Time, cand.CaloEnergy);
                        tof_trueE->Fill(cand.Time, mctrue->Ek());
                        dEE->Fill(cand.CaloEnergy, cand.VetoEnergy);
                        constexpr auto radtodeg = std_ext::radian_to_degree(1.0);
                        theta->Fill(cand.Theta*radtodeg);
                        theta_corr->Fill(mctrue->Theta()*radtodeg, cand.Theta*radtodeg);
                    }
                }

                cand_mult->Fill(event.Reconstructed().Candidates.size());
            }

        }

    }
}


void ProtonCheck::Finish()
{

}


void ProtonCheck::ShowResult()
{
    canvas("ProtonCheck")
            << cand_mult
            << drawoption("colz")
            << tof << tof_trueE << dEE << theta_corr << endc;
}


AUTO_REGISTER_PHYSICS(ProtonCheck)
