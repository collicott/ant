#include "FitGausPol3.h"


#include "BaseFunctions.h"

#include "base/Logger.h"

#include "base/TF1Ext.h"

#include "TF1.h"
#include "TH1.h"

using namespace std;
using namespace ant::calibration;

void ant::calibration::gui::FitGausPol3::Sync()
{
    signal->SetParameters(&(func->GetParameters()[0]));
    bg->SetParameters(    &(func->GetParameters()[3]));
    setRange(signal,GetRange());
    setRange(bg,GetRange());
}

ant::calibration::gui::FitGausPol3::FitGausPol3()
{
    signal = functions::gaus::getTF1();
    signal->SetLineColor(kRed);

    bg = functions::pol<3>::getTF1();
    bg->SetLineColor(kBlue);

    func = functions::GausPol<3>::getTF1();
    func->SetLineColor(kGreen);

    SetRange(ant::interval<double>(100,250));
    func->SetParName(0,"A");
    func->SetParLimits(0, 0.0, 1E+12);
    func->SetParName(1,"x_{0}");
    func->SetParName(2,"#sigma");
    func->SetParLimits(2, 0.0, 1E+12);
    func->SetParName(3,"p_{0}");
    func->SetParName(4,"p_{1}");
    func->SetParName(5,"p_{2}");
    func->SetParName(6,"p_{3}");

    signal->SetNpx(1000);
    bg->SetNpx(1000);
    func->SetNpx(1000);

    AddKnob<KnobsTF1::ParameterKnob>(func->GetParName(0), func, 0, IndicatorProperties::Type_t::slider_horizontal);
    AddKnob<KnobsTF1::ParameterKnob>(func->GetParName(1), func, 1, IndicatorProperties::Type_t::slider_vertical);
    AddKnob<KnobsTF1::ReferenceParameterKnob>(func->GetParName(2), func, 2, 1, IndicatorProperties::Type_t::slider_vertical);
    AddKnob<KnobsTF1::RangeKnob>("min",func,KnobsTF1::RangeKnob::RangeEndType::lower);
    AddKnob<KnobsTF1::RangeKnob>("max",func,KnobsTF1::RangeKnob::RangeEndType::upper);
}

ant::calibration::gui::FitGausPol3::~FitGausPol3()
{
    delete signal;
    delete bg;
    delete func;
}

void ant::calibration::gui::FitGausPol3::Draw()
{
    signal->Draw("same");
    bg->Draw("same");
    func->Draw("same");
}

void ant::calibration::gui::FitGausPol3::Fit(TH1* hist)
{
    FitFunction::doFit(hist, func);
    Sync();
}

void gui::FitGausPol3::FitBackground(TH1* hist)
{
    const auto fixedPars = {0,1,2};
    FixParameters(func, fixedPars);
    FitFunction::doFit(hist, func);
    Sync();
    UnFixParameters(func, fixedPars);
}

void gui::FitGausPol3::FitSignal(TH1* hist)
{
    const auto fixedPars = {3,4,5,6};
    FixParameters(func, fixedPars);
    FitFunction::doFit(hist, func);
    Sync();
    UnFixParameters(func, fixedPars);
}

void ant::calibration::gui::FitGausPol3::SetDefaults(TH1 *hist)
{
    // defaults for taps baf2?
    // Amplitude
    func->SetParameter(0, hist->GetMaximum());
    const double max_pos = hist->GetXaxis()->GetBinCenter(hist->GetMaximumBin());

    // x0
    auto range = GetRange();
    func->SetParameter(1, range.Clip(max_pos));
    func->SetParLimits(1, range.Start(), range.Stop());

    // sigma
    func->SetParameter(2, 8);
    func->SetParLimits(2, 5, 50);

    func->SetParameter(3, 1);
    func->SetParameter(4, 1);
    func->SetParameter(5, 1);
    func->SetParameter(6, 0.1);

    Sync();
}

void ant::calibration::gui::FitGausPol3::SetRange(ant::interval<double> i)
{
    setRange(func, i);
    setRange(signal, i);
    setRange(bg, i);
    // x_0 peak position must be in range
    func->SetParLimits(1, i.Start(), i.Stop());
}

ant::interval<double> ant::calibration::gui::FitGausPol3::GetRange() const
{
    return getRange(func);
}

std::vector<double> ant::calibration::gui::FitGausPol3::Save() const
{
    SavedState_t params;
    params.reserve(2+func->GetNpar());
    saveTF1(func,params);

    return params;
}

void ant::calibration::gui::FitGausPol3::Load(const std::vector<double> &data)
{
    if(data.size() != std::size_t(2+func->GetNpar())) {
        LOG(WARNING) << "Can't load parameters";
        return;
    }

    SavedState_t::const_iterator p = data.begin();
    loadTF1(p, func);
    SetRange(getRange(func));
    Sync();
}

double ant::calibration::gui::FitGausPol3::GetPeakPosition() const
{
    return func->GetParameter(1);
}

double gui::FitGausPol3::GetPeakWidth() const
{
    return func->GetParameter(2);
}

double gui::FitGausPol3::SignalToBackground(const double x) const
{
    const auto s = func->Eval(x);
    const auto b = bg->Eval(x);

    return (s-b)/(s+b);
}

