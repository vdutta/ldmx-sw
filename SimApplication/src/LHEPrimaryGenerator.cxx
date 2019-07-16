#include "SimApplication/LHEPrimaryGenerator.h"

// Geant4
#include "G4Event.hh"
#include "G4IonTable.hh"

// LDMX
#include "SimApplication/UserPrimaryParticleInformation.h"

// Geant4
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

namespace ldmx {

    LHEPrimaryGenerator::LHEPrimaryGenerator( std::string lhe_file_name ) :
            reader_( lhe_file_name ) {
    }

    LHEPrimaryGenerator::~LHEPrimaryGenerator() { }

    void LHEPrimaryGenerator::GeneratePrimaryVertex(G4Event* anEvent) {

        if ( reader_.readEvent() ) {
            
            //grab package of current event info from reader
            LHEF::HEPEUP currentEventInfo = reader_.hepeup;

            //the LHE::Reader::junk string may contain a vertex definition
            std::stringstream maybeVertex;
            maybeVertex << currentEventInfo.junk;

            double vertex_x(0.0), vertex_y(0.0), vertex_z(0.0); //default vertex is the origin
            std::string flag; //parse for the #vertex flag
            while ( maybeVertex >> flag ) {
                //check for vertex flag
                if ( flag == "#vertex" ) {
                    maybeVertex >> vertex_x >> vertex_y >> vertex_z;
                }
            }

            //Open a new G4 primary vertex for this event
            G4PrimaryVertex* vertex = new G4PrimaryVertex();
            vertex->SetPosition( vertex_x , vertex_y , vertex_z );

            vertex->SetWeight( currentEventInfo.XWGTUP );

            //list of primary particles already created
            std::vector<G4PrimaryParticle*> primaryList;

            //loop through particles in event and add them to G4
            int numParticles = currentEventInfo.NUP;
            for ( int particleIndex = 0; particleIndex < numParticles; particleIndex++ ) {

                int particleStatus                      = currentEventInfo.ISTUP.at(particleIndex);
                long int particlePDG                    = currentEventInfo.IDUP.at(particleIndex);
                std::vector<double> particleMomentum    = currentEventInfo.PUP.at(particleIndex);
                double particleInvarTime                = currentEventInfo.VTIMUP.at(particleIndex);
                unsigned int primaryMotherIndex         = currentEventInfo.MOTHUP.at(particleIndex).first;
                int primaryMotherStatus                 = currentEventInfo.ISTUP.at(primaryMotherIndex);

                if ( particleStatus > 0) {
                    //status for this particle is activated

                    //create G4 primary particle for this activated particle
                    G4PrimaryParticle* primary = new G4PrimaryParticle();

                    if ( particlePDG == -623) {
                        //special importing for Tunsten Ion (W)
                        G4ParticleDefinition* tungstenIonDef = G4IonTable::GetIonTable()->GetIon(74, 184, 0.);
                        if (tungstenIonDef != NULL) {
                            primary->SetParticleDefinition(tungstenIonDef);
                        } else {
                            G4Exception("LHEPrimaryGenerator::GeneratePrimaryVertex", 
                                        "EventGenerationError", FatalException, 
                                        "Failed to find particle definition for W ion.");
                        }
                    } else {
                        primary->SetPDGcode( particlePDG );
                    }

                    primary->Set4Momentum(particleMomentum.at(0) * GeV, 
                                          particleMomentum.at(1) * GeV, 
                                          particleMomentum.at(2) * GeV, 
                                          particleMomentum.at(3) * GeV);
                    primary->SetProperTime(particleInvarTime * nanosecond);

                    UserPrimaryParticleInformation* primaryInfo = new UserPrimaryParticleInformation();
                    primaryInfo->setHepEvtStatus( particleStatus );
                    primary->SetUserInformation( primaryInfo );

                    primaryList.push_back( primary );

                    // Add primary to particle tree
                    //      add primary as daugter to its mother if the mother exists and has non-zero status
                    //      otherwise set primary as root primary for this vertex
                    if ( primaryMotherIndex > 0 and primaryMotherIndex < primaryList.size() and primaryMotherStatus > 0 ) {
                        G4PrimaryParticle* primaryMom = primaryList.at( primaryMotherIndex );
                        if (primaryMom != NULL) {
                            primaryMom->SetDaughter( primary );
                        }
                    } else {
                        vertex->SetPrimary( primary );
                    }

                } //status for current particle is activated 

            } //loop through particles (particleIndex)

            //add primary vertex to geant event with particle tree
            anEvent->AddPrimaryVertex(vertex);

        } else {
            std::cout << "[ LHEPrimaryGenerator ] : Ran out of input events so run will be aborted!" << std::endl;
            G4RunManager::GetRunManager()->AbortRun(true);
            anEvent->SetEventAborted();
        }

        return;
    }

}
