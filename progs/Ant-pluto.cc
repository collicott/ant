#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>

// pluto++
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#pragma GCC diagnostic ignored "-Wvla"
#pragma GCC diagnostic ignored "-Wignored-qualifiers"
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic ignored "-Wwrite-strings"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include "PParticle.h"
#include "PBeamSmearing.h"
#include "PReaction.h"
#include "PPlutoBulkDecay.h"
#include "PDataBase.h"
#pragma GCC diagnostic pop

// ROOT
#include "TClonesArray.h"
#include "TMath.h"
#include "TRandom.h"

#include "anyoption.h"



using namespace std;

// Options (default values)
static Double_t Egmax = 1.604;  // beam energy, max Tagg Energy [GeV]
static Double_t Egmin = 1.450;  // beam energy, min Tagg Energy [GeV]
static bool enableBulk = true;
static bool saveIntermediate = true;
static Int_t numEvents = 1000;
static string reactonString = "p omega";
static string prefix = "pluto_";
static string suffix = "00";
static string outfile = "";

static const Double_t me = 0.000510999;  // mass electron [GeV]

/**
 * @brief thetaCrit
 * @return
 */
inline Double_t thetaCrit() {
    return me / Egmax;
}

std::string generateFilename() {
    std::stringstream s;
    s << prefix;

    string fixedName(reactonString);
    std::replace( fixedName.begin(), fixedName.end(), ' ', '_');
    s << "-" << fixedName;

    if(enableBulk)
        s << "-bulk";

    s << "-" << suffix;

    return s.str();
}

void ReadCmdline(int argc, char** argv ) {

     AnyOption opt;

     opt.setVerbose();
     opt.autoUsagePrint(true);

     opt.addUsage(" --numEvents <n>           Number of events to simulate");
     opt.addUsage(" --saveIntermediate        Save intermedite particles");
     opt.addUsage(" --enableBulk              Enable bulk decay");

     opt.addUsage(" --reaction                Reactoin String, PLUTO notation");
     opt.addUsage(" --Emin                    Minimal photon energy [GeV]");
     opt.addUsage(" --Emax                    Maximal photon energy [GeV]");

     opt.addUsage(" --help -h");

     opt.setOption("numEvents");

     opt.setFlag("saveIntermediate");
     opt.setFlag("enableBulk");
     opt.setOption("Emin");
     opt.setOption("Emax");
     opt.setOption("prefix");
     opt.setOption("suffix");
     opt.setOption("reaction");
     opt.setOption("Output");

     opt.setFlag("help",'h');

 //    opt.processFile(CFG_FILE);

     opt.processCommandArgs( argc, argv );

     if( opt.getFlag('h') || opt.getValue("--help") ) {
          opt.printUsage();
          exit(0);
     }

     if ( opt.getValue("numEvents") )
         numEvents = atoi(opt.getValue("numEvents"));

     if ( opt.getValue("reaction") )
         reactonString = opt.getValue("reaction");

     if ( opt.getValue("prefix") )
         prefix = opt.getValue("prefix");

     if ( opt.getValue("suffix") )
         suffix = opt.getValue("suffix");

     if ( opt.getValue("Emin") )
         Egmin = atof(opt.getValue("Emin"));

     if ( opt.getValue("Emax") )
         Egmax = atof(opt.getValue("Emax"));
     if ( opt.getValue("Output") )
         outfile = opt.getValue("Output");


     saveIntermediate = ( NULL != opt.getValue("saveIntermediate") );
     enableBulk = ( NULL != opt.getValue("enableBulk") );

}

void PrintConfig() {
    cout << "\n\n Simulating " << numEvents << " events:\n\n";
    cout << "  Reactrion:  g p --> " << reactonString << "\n\n";
    cout << "  Photon Engery : " << Egmin << " to " << Egmax << " GeV\n\n";
    cout << "  Saving to " << generateFilename() << "\n\n";
    cout << "  saveIntermediate particles: ";
    if( saveIntermediate )
        cout << "yes";
    else
        cout << "no";
    cout << "\n\n";
    cout << "  enable Bulk decay: ";
    if( enableBulk )
        cout << "yes";
    else
        cout << "no";
    cout << "\n\n";

}

int main( int argc, char** argv ) {

    gRandom->SetSeed(); // Initialize ROOT's internal rng. Used for TF1s.

    ReadCmdline( argc, argv );
    PrintConfig();

    PBeamSmearing *smear = new PBeamSmearing("beam_smear", "Beam smearing");
    smear->SetReaction("g + p");

    TF1* tagger_spectrum = new TF1("bremsstrahlung","(1.0/x)", Egmin, Egmax);
    smear->SetMomentumFunction( tagger_spectrum );

    TF1* theta_smear = new TF1( "angle", "x / ( x*x + [0] )^2", 0.0, 5.0 * thetaCrit() );
    theta_smear->SetParameter( 0, thetaCrit() * thetaCrit() );

    smear->SetAngularSmearing( theta_smear );

    makeDistributionManager()->Add(smear);

    PStaticData* sdata = makeStaticData();

    // additional omega decays (from PDG Booklet 2014)
    sdata->AddDecay("w --> eta + g",        "w", "eta,g",           4.6E-4);
    sdata->AddDecay("w --> g + g + g",      "w", "g,g,g",           1.9E-4);    // upper limit
    sdata->AddDecay("w --> pi0 e+ e-",      "w", "pi0,e+,e-",       7.7E-4);
    sdata->AddDecay("w --> pi0 mu+ mu-",    "w", "pi0,mu+,mu-",     1.3E-4);
    sdata->AddDecay("w --> pi+ pi- pi0 pi0","w", "pi+,pi-,pi0,pi0", 2.0E-4);    // upper limit
    sdata->AddDecay("w --> pi+ pi- pi+ pi-","w", "pi+,pi-,pi+,pi-", 1.0E-3);    // upper limit
    sdata->AddDecay("w --> pi0 pi0 g",      "w", "pi0,pi0,g",       6.6E-5);
    sdata->AddDecay("w --> eta pi0 g",      "w", "eta,pi0,g",       3.3E-5);    // upper limit

    // Charge conjucation violating modes
    sdata->AddDecay("w --> eta pi0",        "w", "eta,pi0",         2.1E-4);    // upper limit
    sdata->AddDecay("w --> pi0 pi0",        "w", "pi0,pi0",         2.1E-4);    // upper limit
    sdata->AddDecay("w --> pi0 pi0 pi0",    "w", "pi0,pi0,pi0",     2.3E-4);    // upper limit

    if( outfile.empty() )
        outfile = generateFilename();

    // PReaction constructor requires non-const char*. so... make copies... ARGH!
    char* rs = strdup(reactonString.c_str());
    char* gf = strdup(outfile.c_str());

    cout << "filename " << gf << endl;

    PReaction* reactrion = new PReaction(Egmax, "g", "p", rs, gf, saveIntermediate, 0, 0, 0);

    if( enableBulk ) {
        PPlutoBulkDecay* p1 = new PPlutoBulkDecay();
        p1->SetRecursiveMode(1);
        p1->SetTauMax(0.001);
        reactrion->AddBulk(p1);
    }

    reactrion->Print();   //The "Print()" statement is optional

    reactrion->Loop(numEvents);

    cout << "Simulation finished." << endl;

    free(rs);
    free(gf);

    // Do not delete the reaction, otherwise: infinite loop somewhere in ROOT...
    //delete reactrion;

    return 0;
}

