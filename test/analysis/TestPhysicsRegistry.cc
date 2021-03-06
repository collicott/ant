#include "catch.hpp"
#include "expconfig_helpers.h"

#include "analysis/physics/Physics.h"
#include "expconfig/ExpConfig.h"

#include "analysis/utils/Uncertainties.h"

#include "base/OptionsList.h"
#include "base/WrapTFile.h"
#include "base/tmpfile_t.h"

#include "TError.h"

#include <memory>

using namespace std;
using namespace ant;
using namespace ant::analysis;

void dotest();

TEST_CASE("PhysicsRegistry: Create all physics classes", "[analysis]") {
    test::EnsureSetup();
    dotest();
}

bool histogram_overwrite_detected = false;
bool duplicate_mkdir_detected = false;


void dotest() {

    // overwrite ROOT's error handler to detect some warnings
    SetErrorHandler([] (
                    int level, Bool_t abort, const char *location,
                    const char *msg) {
        // those tests are specific enough...
        if(string(location) == "TDirectory::Append")
            histogram_overwrite_detected = true;
        if(string(location) == "TDirectoryFile::Append")
            histogram_overwrite_detected = true;
        if(string(location) == "TDirectoryFile::mkdir")
            duplicate_mkdir_detected = true;
        DefaultErrorHandler(level, abort, location, msg);
    });

    // create all available physics classes
    for(auto name : PhysicsRegistry::GetList()) {

        cout << "Running physics class " << name << endl;
        // some errors only appear when some outfiles are present
        tmpfile_t tmpfile;
        auto outfile = std_ext::make_unique<WrapTFileOutput>(tmpfile.filename,
                                WrapTFileOutput::mode_t::recreate,
                                true);

        histogram_overwrite_detected = false;
        duplicate_mkdir_detected = false;
        INFO(name);
        try {
            PhysicsRegistry::Create(name);
            REQUIRE_FALSE(histogram_overwrite_detected);
            REQUIRE_FALSE(duplicate_mkdir_detected);
            auto objects = outfile->GetListOf<TNamed>();
            if(objects.size()>1) {
                for(auto o : objects)
                    cout << "Found object " << o->GetName() << endl;
            }
            REQUIRE(objects.size() == 1);
            REQUIRE(objects.front()->GetName() == name);
        }
        catch(PhysicsRegistry::Exception e) {
            FAIL(string("Physics Registry error: ")+e.what());
        }
        catch(WrapTFile::Exception) {
            // ignore silently if Physics classes can't load some files...
            continue;
        }
        catch(ExpConfig::ExceptionNoDetector) {
            // ignore silently if test setup did not provide detector
            continue;
        }
        catch(utils::UncertaintyModel::Exception) {
            // ignore silently if class cannot load uncertainty model
            continue;
        }
        catch(const std::exception& e) {
            FAIL(string("Unexpected exception: ")+e.what());
        }
        catch(...) {
            FAIL("Something weird was thrown.");
        }
        // write the file
        REQUIRE_NOTHROW(outfile = nullptr);

    }


}
