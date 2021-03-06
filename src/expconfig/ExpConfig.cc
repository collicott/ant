#include "ExpConfig.h"

#include "setups/Setup.h"

#include "unpacker/UnpackerAcqu.h"
#include "unpacker/UnpackerA2Geant.h"

#include "tree/TID.h"
#include "base/Logger.h"

#include <type_traits>
#include <list>
#include <iostream>

using namespace std;

namespace ant { // template implementations need explicit namespace

std::string ExpConfig::Setup::manualName = ""; // default empty
std::shared_ptr<ExpConfig::Setup> ExpConfig::Setup::lastFound = nullptr; // default nothing found so far

shared_ptr<ExpConfig::Setup> ExpConfig::Setup::Get(const TID& tid)
{

    shared_ptr<Setup> config = nullptr;

    if(!manualName.empty()) {
        config = ExpConfig::Setup::Get(manualName);
        if(config == nullptr) {
            throw Exception(
                        std_ext::formatter()
                        << "Found no config matching name "
                        << manualName
                        );
        }
        return config;
    }

    // go to automatic search mode in all registered setups

    std::list< std::shared_ptr<Setup> > modules;
    for(auto setup_name : expconfig::SetupRegistry::GetNames()) {
        modules.emplace_back(expconfig::SetupRegistry::GetSetup(setup_name));
    }

    // remove the config if the config says it does not match
    modules.remove_if([&tid] (const shared_ptr<Setup>& m) {
        return !m->Matches(tid);
    });

    // check if something reasonable is left
    if(modules.empty()) {
        throw ExpConfig::Exception(std_ext::formatter()
                                   << "No setup found for TID "
                                   << tid);
    }
    if(modules.size()>1) {
        throw ExpConfig::Exception(std_ext::formatter()
                                   << "More than one setup found for TID "
                                   << tid);
    }
    Setup::lastFound = modules.back();

    return modules.back();
}



shared_ptr<ExpConfig::Setup> ExpConfig::Setup::Get(const std::string& name)
{
    lastFound = expconfig::SetupRegistry::GetSetup(name);
    return lastFound;
}

shared_ptr<ExpConfig::Setup> ExpConfig::Setup::GetLastFound()
{
    // try if we have a name
    if(lastFound==nullptr && !manualName.empty()) {
        Get(manualName);
    }
    return lastFound;
}

shared_ptr<Detector_t> ExpConfig::Setup::GetDetector(Detector_t::Type_t type)
{
    auto config = GetLastFound();
    if(config == nullptr)
        throw ExceptionNoConfig("Could not find setup to search for required detector");
    for(const auto& detector : config->GetDetectors()) {
        if(detector->Type == type)
            return detector;
    }
    throw ExceptionNoDetector("Could not find detector in given setup");
}

void ExpConfig::Setup::SetManualName(const string& setupname, bool required)
{
    manualName = setupname;
    if(required && Get(setupname) == nullptr)
        throw ExceptionNoConfig("No setup found in registry for manual name "+setupname);
}

void ExpConfig::Setup::Cleanup()
{
    lastFound = nullptr;
    manualName = "";
    expconfig::SetupRegistry::Destroy();
}

std::list<string> ExpConfig::Setup::GetNames() {
    return expconfig::SetupRegistry::GetNames();
}


} // namespace ant




