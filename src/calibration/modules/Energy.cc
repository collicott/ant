#include "Energy.h"
#include "analysis/plot/HistogramFactories.h"
#include "analysis/data/Event.h"
#include "analysis/utils/combinatorics.h"
#include "calibration/CalibrationDataManager.h"

#include "tree/TDetectorRead.h"

#include "base/Logger.h"

#include <cstdint>
#include <ctime>
#include <algorithm>
#include <sstream>
#include <vector>
#include <list>

using namespace std;
using namespace ant;
using namespace ant::calibration;

size_t Energy::CalibType::Instances = 0;

Energy::Energy(Detector_t::Type_t detectorType,
               std::shared_ptr<DataManager> calmgr,
               Calibration::Converter::ptr_t converter,
               double defaultPedestal,
               double defaultGain,
               double defaultThreshold,
               double defaultRelativeGain) :
    Calibration::Module(
        std_ext::formatter()
        << Detector_t::ToString(detectorType)
        << "_Energy"
           ),
    DetectorType(detectorType),
    calibrationManager(calmgr),
    Converter(move(converter)),
    Pedestals(defaultPedestal,"Pedestals"),
    Gains(defaultGain,"Gains"),
    Thresholds(defaultThreshold,"Thresholds"),
    RelativeGains(defaultRelativeGain,"RelativeGains")
{
    if(Converter==nullptr)
        throw std::runtime_error("Given converter should not be nullptr");
}

void Energy::ApplyTo(const readhits_t& hits, extrahits_t& extrahits)
{
    const auto& dethits = hits.get_item(DetectorType);

    // now calibrate the Energies (ignore any other kind of hits)
    for(TDetectorReadHit* dethit : dethits) {
        if(dethit->GetChannelType() != Channel_t::Type_t::Integral)
            continue;



        // Values might already be filled
        // (for example by previous calibration run, or A2Geant unpacker),
        // then we apply the threshold and the relative gain only
        std::vector<double> values(0);

        // prefer RawData if available
        if(!dethit->RawData.empty()) {
            // the Converter is smart enough to account for reference Energys!
            values = Converter->Convert(dethit->RawData);

            // for pedestal calibration, we insert extra hits here
            // containing the raw values
            extrahits.emplace_back(
                        LogicalChannel_t{
                            dethit->GetDetectorType(),
                            Channel_t::Type_t::Pedestal,
                            dethit->Channel
                        },
                        values
                        );

            // apply pedestal/gain/threshold to each of the values (might be multihit)
            for(double& value : values) {
                if(Pedestals.Values.empty())
                    value -=Pedestals.DefaultValue;
                else
                    value -= Pedestals.Values[dethit->Channel];

                if(Gains.Values.empty())
                    value *= Gains.DefaultValue;
                else
                    value *= Gains.Values[dethit->Channel];
            }

        }
        else {
            // maybe the values are already filled
            values = dethit->Values;
            dethit->Values.resize(0);
        }

        // always apply the threshold cut and the relative gains
        dethit->Values.reserve(values.size());

        for(double value : values) {
            if(RelativeGains.Values.empty())
                value *= RelativeGains.DefaultValue;
            else
                value *= RelativeGains.Values[dethit->Channel];

            const double threshold = Thresholds.Values.empty()
                                     ? Thresholds.DefaultValue
                                     : Thresholds.Values[dethit->Channel];
            if(value<threshold)
                continue;

            // only add if it passes the threshold
            dethit->Values.push_back(value);
        }

    }
}

std::vector<std::list<TID> > Energy::GetChangePoints() const
{
    vector<list<TID>> changePointLists;

    for (auto& calibType: { Pedestals.Name, Gains.Name, Thresholds.Name, RelativeGains.Name})
        changePointLists.push_back(calibrationManager->GetChangePoints(std_ext::formatter()
                                                                       << GetName()
                                                                       << "-" << calibType));
    return changePointLists;
}
void Energy::Update(size_t index, const TID& tid)
{
    for (auto calibration: {&Pedestals, &Gains, &Thresholds, &RelativeGains})
    {
        if (calibration->Index == index)
        {
            TCalibrationData cdata;
            if(calibrationManager->GetData(GUI_CalibType::ConstructName(GetName(), calibration->Name),
                                           tid, cdata)) {
                calibration->Values.clear();
                calibration->Values.reserve(cdata.Data.size());
                for (auto& val: cdata.Data)
                    calibration->Values.push_back(val.Value);
            }
            else {
                LOG(ERROR) << "Could not update calibration data for " << calibration->Name
                             << "at changepoint TID=" << tid;
            }
        }
    }
}

Energy::~Energy()
{

}

Energy::GUI_CalibType::GUI_CalibType(const string& basename, CalibType& type,
                                     const shared_ptr<DataManager>& calmgr) :
    gui::Manager_traits(basename),
    calibType(type),
    calibrationManager(calmgr)
{}

string Energy::GUI_CalibType::GetName() const
{
    // serves as the CalibrationID for the manager,
    // and as the histogram name
    return ConstructName(Manager_traits::GetName(), calibType.Name);
}

string Energy::GUI_CalibType::GetHistogramName() const
{
    return GetName();
}

void Energy::GUI_CalibType::StartRange(const interval<TID>& range)
{
    // always make sure the values are large enough
    std::vector<double>& values = calibType.Values;
    values.resize(GetNumberOfChannels(), calibType.DefaultValue);

    TCalibrationData cdata;
    if(calibrationManager->GetData(GetName(), range.Start(), cdata)) {
        for(const TKeyValue<double>& kv : cdata.Data) {
            values[kv.Key] = kv.Value;
        }
        for(const TKeyValue<vector<double>>& kv : cdata.FitParameters) {
            fitParameters.insert(make_pair(kv.Key, kv.Value));
        }
        LOG(INFO) << GetName() << ": Loaded previous values from database";
    }
    else {
        LOG(INFO) << GetName() << ": No previous values found, built new gains for all channel from default gain";
    }

    // save a copy for comparison at finish stage
    previousValues = calibType.Values;

}

void Energy::GUI_CalibType::StoreFinishRange(const interval<TID>& range)
{
    TCalibrationData cdata(
                "Unknown", /// \todo get static information about author/comment?
                "No Comment",
                time(nullptr),
                GetName(),
                range.Start(),
                range.Stop()
                );

    std::vector<double>& values = calibType.Values;

    // fill data
    cdata.Data.resize(0);
    for(unsigned ch=0;ch<values.size();ch++) {
        cdata.Data.emplace_back(ch, values[ch]);
    }

    // fill fit parameters (if any)
    cdata.FitParameters.resize(0);
    for(const auto& it_map : fitParameters) {
        const unsigned ch = it_map.first;
        const vector<double>& params = it_map.second;
        cdata.FitParameters.emplace_back(ch, params);
    }

    calibrationManager->Add(cdata);

    LOG(INFO) << "Added TCalibrationData " << cdata;
}
