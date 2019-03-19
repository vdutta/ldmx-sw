/*
 * @file TargetBremFilter.cxx
 * @class TargetBremFilter
 * @brief Class defining a UserActionPlugin that allows a user to filter out 
 *        events that don't result in a brem within the target.
 * @author Omar Moreno, SLAC National Accelerator Laboratory
 */

#include "Biasing/TargetBremFilter.h"

//------------//
//   Geant4   //
//------------//
#include "G4RunManager.hh"

//-------------//
//   ldmx-sw   //
//-------------//
#include "Biasing/TargetBremFilterMessenger.h"
#include "SimCore/UserTrackInformation.h"

SIM_PLUGIN(ldmx, TargetBremFilter)

namespace ldmx { 

    TargetBremFilter::TargetBremFilter() {
        messenger_ = new TargetBremFilterMessenger(this);
    }

    TargetBremFilter::~TargetBremFilter() {
    }


    G4ClassificationOfNewTrack TargetBremFilter::stackingClassifyNewTrack(
            const G4Track* track, 
            const G4ClassificationOfNewTrack& currentTrackClass) {

        // Use current classification by default so values from other plugins are not overridden.
        G4ClassificationOfNewTrack classification = currentTrackClass;

        /*G4String particleName = track->GetParticleDefinition()->GetParticleName();
        G4int pdgID = track->GetParticleDefinition()->GetPDGEncoding();
        std::cout << "[ EcalProcessFilter ]: Pushing "
                    << "\tParticle " << particleName << " ( PDG ID: " << pdgID << " ) : " << "\n"
                    << "\tTrack ID: " << track->GetTrackID() << " to stack.\n" 
                    << std::endl;*/

        if (track->GetParentID() == 0) return fWaiting; 

        return classification;
    }

    void TargetBremFilter::stepping(const G4Step* step) { 

        // Get the track associated with this step.
        G4Track* track = step->GetTrack();

        // Only process the primary electron track
        if (track->GetParentID() != 0) return;

        // get the PDGID of the track.
        G4int pdgID = track->GetParticleDefinition()->GetPDGEncoding();
        
        // Make sure that the particle being processed is an electron.
        if (pdgID != 11) return;

        // Get the volume the particle is in.
        G4String volumeName = track->GetVolume()->GetName();

        // Get the particle type.
        G4String particleName = track->GetParticleDefinition()->GetParticleName();

        /* 
        // Get the kinetic energy of the particle.
        double incidentParticleEnergy = step->GetPostStepPoint()->GetTotalEnergy();
        std::cout << "*******************************\n" 
                  << "*   Step " << track->GetCurrentStepNumber() << "\n"
                  << "********************************\n"
                  << "[ TargetBremFilter ]: " << "\n" 
                  << "\tTotal energy of " << particleName      << " ( PDG ID: " << pdgID
                  << " ) : " << incidentParticleEnergy       << "\n"
                  << "\tTrack ID: " << track->GetTrackID()     << "\n" 
                  << "\tStep #: " << track->GetCurrentStepNumber() << "\n"
                  << "\tParticle currently in " << volumeName  
                  << "\tPost step process: " << step->GetPostStepPoint()->GetStepStatus() 
                  << std::endl;
        }
        */

        // If the particle isn't in the target, don't continue with the processing.
        if (volumeName.compareTo(volumeName_) != 0) return;

        // Get the secondaries produced in the step and check if which process 
        // they were created with.
        const std::vector<const G4Track*>* secondaries = step->GetSecondaryInCurrentStep(); 
        if (!secondaries->empty()) {
            
            //std::cout << "[ TargetBremFilter ]: Total secondaries: " 
            //          << secondaries->size() << std::endl;
           
            const G4Track* secondary_track = (*secondaries)[0]; 
            G4String processName = secondary_track->GetCreatorProcess()->GetProcessName();
            
            if (processName.compareTo("eBrem") == 0 
                        && secondary_track->GetKineticEnergy() > bremEnergyThreshold_) {
                    
                        
                /*std::cout << "[ TargetBremFilter ]: "                 
                          << "Tagging brem as a candidate with energy " 
                          << secondary_track->GetKineticEnergy() << std::endl;*/

                    
                    if (!secondary_track->GetUserInformation()) {
                        auto trackInfo{new UserTrackInformation()}; 
                        trackInfo->setInitialMomentum(secondary_track->GetMomentum()); 
                        const_cast<G4Track*>(secondary_track)->SetUserInformation(trackInfo); 
                    }

                static_cast<UserTrackInformation*>(secondary_track->GetUserInformation())->tagBremCandidate();
                hasBremCandidate_ = true; 
            } 
        }

        // Check if the particle is exiting the volume.  If so, check if any 
        // brem candidates were found.  If not, abort the event.
        if (step->GetPostStepPoint()->GetStepStatus() == fGeomBoundary) {
            

            /*std::cout << "[ TargetBremFilter ]: "
                      << "Particle " << particleName << "is leaving the "
                      << volumeName << " volume with momentum "
                      << track->GetMomentum().mag() << std::endl;*/

            // Check if any brem candidates were tagged.  If not, just stopped
            // processing the event.
            if (!hasBremCandidate_) {
                
                /*std::cout << "[ TargetBremFilter ]: "
                          << "Brem candidate was not found. --> Killing all tracks!"
                          << std::endl;*/
                track->SetTrackStatus(fKillTrackAndSecondaries);
                G4RunManager::GetRunManager()->AbortEvent();
                return;
            }
            // Check if the recoil electron should be killed.  If not, postpone 
            // its processing until the brem gamma has been processed.
            else if (killRecoilElectron_) { 
                track->SetTrackStatus(fStopAndKill);
            } 
            
            track->SetTrackStatus(fSuspend);  

             
        } else if (step->GetPostStepPoint()->GetKineticEnergy() == 0) { 
            
            /*std::cout << "[ TargetBremFilter ]: "
                        << "Electron never made it out of the target --> Killing all tracks!"
                        << std::endl;*/

            track->SetTrackStatus(fKillTrackAndSecondaries);
            G4RunManager::GetRunManager()->AbortEvent();
            return;
        }
    }

    void TargetBremFilter::endEvent(const G4Event*) {
        hasBremCandidate_ = false; 
    }
}

