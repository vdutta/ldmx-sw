/*
 * @file TargetBremFilter.cxx
 * @class TargetBremFilter
 * @brief Class defining a UserActionPlugin that allows a user to filter out 
 *        events that don't result in a brem within the target.
 * @author Omar Moreno, SLAC National Accelerator Laboratory
 */

#include "Biasing/TargetBremFilter.h"

#include "SimCore/UserTrackInformation.h"

SIM_PLUGIN(ldmx, TargetBremFilter)

namespace ldmx { 

    std::vector<G4Track*> TargetBremFilter::bremGammaTracks_ = {};

    TargetBremFilter::TargetBremFilter() {
        messenger_ = new TargetBremFilterMessenger(this);
    }

    TargetBremFilter::~TargetBremFilter() {
    }


    G4ClassificationOfNewTrack TargetBremFilter::stackingClassifyNewTrack(
            const G4Track* track, 
            const G4ClassificationOfNewTrack& currentTrackClass) {

        if (verbose_ > 0) { 
            std::cout << "********************************" << std::endl; 
            std::cout << "*   Track pushed to the stack  *" << std::endl;
            std::cout << "********************************" << std::endl;
        }

        // get the PDGID of the track.
        G4int pdgID = track->GetParticleDefinition()->GetPDGEncoding();

        // Get the particle type.
        G4String particleName = track->GetParticleDefinition()->GetParticleName();

        if (verbose_ > 0) { 
            std::cout << "[ TargetBremFilter ]: " << "\n" 
                      << "\tParticle " << particleName << " ( PDG ID: " << pdgID << " ) : " << "\n"
                      << "\tTrack ID: " << track->GetTrackID() << "\n"
                      << "\tEnergy: " << track->GetKineticEnergy()  
                      << std::endl;
        }


        // Use current classification by default so values from other plugins are not overridden.
        G4ClassificationOfNewTrack classification = currentTrackClass;

        if (track->GetTrackID() == 1 && pdgID == 11) {
            if (verbose_ > 0) {
                std::cout << "[ TargetBremFilter ]: Pushing track to waiting stack." << std::endl;
            }
            return fWaiting; 
        }

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
        //G4VPhysicalVolume* volume = track->GetVolume();
        G4String volumeName = track->GetVolume()->GetName();

        // Get the particle type.
        G4String particleName = track->GetParticleDefinition()->GetParticleName();

        // Get the kinetic energy of the particle.
        double incidentParticleEnergy = step->GetPostStepPoint()->GetTotalEnergy();
        if (verbose_ > 0) { 
            std::cout << "*******************************" << std::endl; 
            std::cout << "*   Step " << track->GetCurrentStepNumber() << std::endl;
            std::cout << "********************************" << std::endl;

            std::cout << "[ TargetBremFilter ]: " << "\n" 
                      << "\tTotal energy of " << particleName      << " ( PDG ID: " << pdgID
                      << " ) : " << incidentParticleEnergy       << "\n"
                      << "\tTrack ID: " << track->GetTrackID()     << "\n" 
                      << "\tStep #: " << track->GetCurrentStepNumber() << "\n"
                      << "\tParticle currently in " << volumeName  
                      << "\tPost step process: " << step->GetPostStepPoint()->GetStepStatus() 
                      << std::endl;
        }

        // If the particle isn't in the target, don't continue with the processing.
        if (volumeName.compareTo(volumeName_) != 0) return;

        // Get the secondaries produced in the step and check if which process 
        // they were created with.
        const std::vector<const G4Track*>* secondaries = step->GetSecondaryInCurrentStep(); 
        if (!secondaries->empty()) {
            
            if (verbose_ > 0) { 
                std::cout << "[ TargetBremFilter ]: Total secondaries: " 
                          << secondaries->size() << std::endl;
            }
           
            const G4Track* secondary_track = (*secondaries)[0]; 
            G4String processName = secondary_track->GetCreatorProcess()->GetProcessName();
            if (verbose_ > 0) { 
                std::cout << "[ TargetBremFilter ]: Secondary created with process: " 
                          << processName << std::endl;
            }
            
            if (processName.compareTo("eBrem") == 0 
                        && secondary_track->GetKineticEnergy() > bremEnergyThreshold_) {
                    
                    if (verbose_ > 0) { 
                        std::cout << "[ TargetBremFilter ]: " 
                                  << "Tagging brem as a candidate with energy " 
                                  << secondary_track->GetKineticEnergy() << std::endl;
                    }

                    
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
            
            if (verbose_ > 0) { 
                std::cout << "[ TargetBremFilter ]: "
                          << "Particle " << particleName << "is leaving the "
                          << volumeName << " volume with momentum "
                          << track->GetMomentum().mag() << std::endl;
            }

            /*
            if (track->GetMomentum().mag() >= recoilEnergyThreshold_) { 
                std::cout << "[ TargetBremFilter ]: "
                            << "Electron energy is above threshold --> Aborting event."
                            << std::endl;
                
                track->SetTrackStatus(fKillTrackAndSecondaries);
                G4RunManager::GetRunManager()->AbortEvent();
                return;
            
            }*/

            // Check if any brem candidates were tagged.  If not, just stopped
            // processing the event.
            if (!hasBremCandidate_) {
                
                if (verbose_ > 0) { 
                    std::cout << "[ TargetBremFilter ]: "
                              << "Brem candidate was not found. --> Killing all tracks!"
                              << std::endl;
                }
                track->SetTrackStatus(fKillTrackAndSecondaries);
                G4RunManager::GetRunManager()->AbortEvent();
                return;
            }
            // Check if the recoil electron should be killed.  If not, postpone 
            // its processing until the brem gamma has been processed.
            else if (killRecoilElectron_) { 
                if (verbose_ > 0) {
                    std::cout << "[ TargetBremFilter ]: Killing electron track."
                              << std::endl;
                }

                track->SetTrackStatus(fStopAndKill);
            } 
            
            track->SetTrackStatus(fSuspend);  

             
        } else if (step->GetPostStepPoint()->GetKineticEnergy() == 0) { 
            
            if (verbose_ > 0) { 
                std::cout << "[ TargetBremFilter ]: "
                          << "Electron never made it out of the target --> Killing all tracks!"
                          << std::endl;
            }

            track->SetTrackStatus(fKillTrackAndSecondaries);
            G4RunManager::GetRunManager()->AbortEvent();
            return;
        }
           
            /*

            

       


        } 
        else if (step->GetPostStepPoint()->GetKineticEnergy() == 0) { 
            std::cout << "[ TargetBremFilter ]: "
                        << "Electron never made it out of the target --> Killing all tracks!"
                        << std::endl;

            track->SetTrackStatus(fKillTrackAndSecondaries);
            G4RunManager::GetRunManager()->AbortEvent();
            return;
        }*/
    }

    void TargetBremFilter::endEvent(const G4Event*) {
        //bremGammaTracks_.clear();
        hasBremCandidate_ = false; 
    }
    
    void TargetBremFilter::removeBremFromList(G4Track* track) {   
        bremGammaTracks_.erase(std::remove(bremGammaTracks_.begin(), 
                    bremGammaTracks_.end(), track), bremGammaTracks_.end());
    }
}

