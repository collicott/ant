#include "catch.hpp"

#include "tmpfile_t.h"

#include "tree/TCalibrationData.h"
//#include "tree/TDataRecord.h"
#include "tree/TID.h"

#include "TFile.h"
#include "TTree.h"

#include <iostream>

using namespace std;

void dotest();

TEST_CASE("TCalibrationData: Write to TTree", "[tree]") {
    dotest();
}

void dotest()
{
    ant::tmpfile_t tmpfile;


    TFile f(tmpfile.filename.c_str(),"RECREATE");

    TTree* tree = new TTree("testtree","TCalibData Test Tree");
    ant::TCalibrationData* cdata = new ant::TCalibrationData(tmpfile.filename,ant::TID(0,0u),ant::TID(0,1u,{ant::TID::Flags_t::MC}));
    tree->Branch("cdata",cdata);

    ant::TCalibrationData* cdata2 = new ant::TCalibrationData(tmpfile.filename,
                                                              ant::TID(0,1011u),
                                                              ant::TID(0,1213u,{ant::TID::Flags_t::MC})
                                                              );

    cdata2->Data.emplace_back(1,1.);
    cdata2->Data.emplace_back(2,2.1);
    tree->Branch("cdata2",cdata2);

    cdata->Author  = "Martin Wolfes";


    cdata->Data.emplace_back(1,2.1);
    tree->Fill();
    cdata->Data.emplace_back(2,3.2);
    tree->Fill();

    cout << cdata << endl;
    cout << cdata2 << endl;

    f.Write();
    f.Close();

    tree = nullptr;

    TFile f2(tmpfile.filename.c_str(),"READ");
    REQUIRE(f2.IsOpen());

    f2.GetObject("testtree",tree);

    REQUIRE(tree!=nullptr);

    ant::TCalibrationData* readcdata = nullptr;
    ant::TCalibrationData* readcdata2 = nullptr;
    tree->SetBranchAddress("cdata",&readcdata);
    tree->SetBranchAddress("cdata2",&readcdata2);

    REQUIRE(tree->GetEntries() == 2);

    tree->GetEntry(0);
    REQUIRE(readcdata->Data.size() == 1);

    tree->GetEntry(1);
    REQUIRE(readcdata->Data.size() == 2);

    tree->GetEntry(2);
    REQUIRE(readcdata2->Data.back().Value == 2.1);

}
