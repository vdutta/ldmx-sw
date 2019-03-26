/*
 * @file EcalPNProcessFilter.cxx
 * @brief User action plugin that biases Geant4 to only process events which
 *        involve a photonuclear reaction in the ECal with certain final-state
 *        kinematics.
 * @author Joshua Hiltbrand 
 *         University of Minnesota 
 */

#include "Biasing/EcalPNProcessFilter.h"

SIM_PLUGIN(ldmx, EcalPNProcessFilter)

namespace ldmx { 

    EcalPNProcessFilter::EcalPNProcessFilter() {
        messenger_ = new EcalPNProcessFilterMessenger(this);
    }

    EcalPNProcessFilter::~EcalPNProcessFilter() {
    }

    G4ClassificationOfNewTrack EcalPNProcessFilter::stackingClassifyNewTrack(
            const G4Track* track, 
            const G4ClassificationOfNewTrack& currentTrackClass) {

        if (track == currentTrack_) {
            currentTrack_ = nullptr; 
            return fWaiting; 
        }

        // Use current classification by default so values from other plugins are not overridden.
        G4ClassificationOfNewTrack classification = currentTrackClass;
        
        return classification;
    }

    void EcalPNProcessFilter::stepping(const G4Step* step) { 

        if (TargetBremFilter::getBremGammaList().empty()) return; 
        
        // Get the track associated with this step.
        G4Track* track = step->GetTrack();

        // Only process tracks whose parent is the primary particle
        if (track->GetParentID() != 1) return; 

        // get the PDGID of the track.
        G4int pdgID = track->GetParticleDefinition()->GetPDGEncoding();

        // Make sure that the particle being processed is a photon.
        if (pdgID != 22) return; 

        // Get the volume the particle is in.
        G4String volumeName = track->GetVolume()->GetName();

        // If the particle isn't in the specified volume, stop processing the event.
        std::vector<G4Track*> bremGammaList = TargetBremFilter::getBremGammaList();
        if (std::find(std::begin(volumes_), std::end(volumes_), volumeName) == std::end(volumes_)) {

            // If secondaries were produced outside of the volume of interest, 
            // and there aren't additional brems to process, abort the 
            // event.  Otherwise, suspend the track and move on to the next brem.
            if (step->GetSecondary()->size() != 0 
                    && (std::find(bremGammaList.begin(), bremGammaList.end(), track) != bremGammaList.end())) { 
                
                if (bremGammaList.size() == 1) { 
                    track->SetTrackStatus(fKillTrackAndSecondaries);
                    G4RunManager::GetRunManager()->AbortEvent();
                    currentTrack_ = nullptr;
                    return;
                } else {

                    currentTrack_ = track; 
                    track->SetTrackStatus(fSuspend);
                    TargetBremFilter::removeBremFromList(track);
                    return;
                }
            }
            return;
        }

        // The list of brems will only contain a given track/particle if it 
        // originates from the target.  If the gamma originates elsewhere, 
        // suspend it and move on to the next gamma.
        if (std::find(bremGammaList.begin(), bremGammaList.end(), track) == bremGammaList.end()) { 
            
            currentTrack_ = track; 
            track->SetTrackStatus(fSuspend);
            return;
        }
 
        // Get the particles daughters.
        const G4TrackVector* secondaries = step->GetSecondary();

        // If the particle doesn't interact, then move on to the next step.
        if (secondaries->size() == 0) { 
            
            // If the particle is exiting the bounding volume, kill it.
            if (!boundVolumes_.empty() && step->GetPostStepPoint()->GetStepStatus() == fGeomBoundary) {
                if (std::find(std::begin(boundVolumes_), std::end(boundVolumes_), volumeName) != std::end(boundVolumes_)) {

                    if (bremGammaList.size() == 1) { 

                        track->SetTrackStatus(fKillTrackAndSecondaries);
                        G4RunManager::GetRunManager()->AbortEvent();
                        currentTrack_ = nullptr;
                        return;

                    } else { 
                        currentTrack_ = track; 
                        track->SetTrackStatus(fSuspend);
                        TargetBremFilter::removeBremFromList(track);
                        return;
                    }
                }
            }

        } else {

            // If the brem gamma interacts and produces secondaries, get the 
            // process used to create them. 
            G4String processName = secondaries->at(0)->GetCreatorProcess()->GetProcessName(); 
            if (!processName.contains(BiasingMessenger::getProcess())) {

                if (bremGammaList.size() == 1) { 

                    track->SetTrackStatus(fKillTrackAndSecondaries);
                    G4RunManager::GetRunManager()->AbortEvent();
                    currentTrack_ = nullptr;
                    return;

                } else { 
                    currentTrack_ = track; 
                    track->SetTrackStatus(fSuspend);
                    TargetBremFilter::removeBremFromList(track);
                    return;
                }
            } else { // This was a photonuclear reaction!
                
                bool interesting = false;
                G4double gammaEnergy = step->GetPreStepPoint()->GetKineticEnergy();
                double totalNeutronEnergy = 0.0;
                for (auto& iSec : *secondaries) {
                    G4int pdg = abs(iSec->GetParticleDefinition()->GetPDGEncoding());
                    if (pdg != 2112 && pdg != 130 && pdg != 310 && pdg != 321) {
                        continue;
                    } else {
                        totalNeutronEnergy += iSec->GetKineticEnergy();
                        if (totalNeutronEnergy / gammaEnergy > neutronEnergyFraction_) {
                            interesting = true;
                            break;
                        }
                    }
                }

                if (!interesting) {
                
                    // Total neutron / kaon KE did not exceed threshold
                    track->SetTrackStatus(fKillTrackAndSecondaries);
                    G4RunManager::GetRunManager()->AbortEvent();
                    currentTrack_ = nullptr;
                    return;
                }
            }
            
	        TargetBremFilter::removeBremFromList(track);
            BiasingMessenger::setEventWeight(track->GetWeight());
            photonGammaID_ = track->GetTrackID(); 
        }

    }

    void EcalPNProcessFilter::postTracking(const G4Track* track) { 
       
        if (track->GetParentID() == photonGammaID_) { 
            UserTrackInformation* userInfo 
              = dynamic_cast<UserTrackInformation*>(track->GetUserInformation());
            userInfo->setSaveFlag(true); 
        }
    }

    void EcalPNProcessFilter::addVolume(std::string volume) { 
        
        std::cout << "[ EcalPNProcessFilter ]: Applying filter to volume " << volume << std::endl;
        if (volume.compare("ecal") == 0) { 
            for (G4VPhysicalVolume* physVolume : *G4PhysicalVolumeStore::GetInstance()) {
                G4String physVolumeName = physVolume->GetName();
                if ((physVolumeName.contains("W") || physVolumeName.contains("Si")) 
                        && physVolumeName.contains("phys")) {
                    volumes_.push_back(physVolumeName);
                }
            }
        } else { 
            volumes_.push_back(volume);
        } 
    }        
    
    void EcalPNProcessFilter::addBoundingVolume(std::string volume) { 
        std::cout << "[ EcalPNProcessFilter ]: Bounding particle to volume " << volume << std::endl;
        boundVolumes_.push_back(volume);
    }        
}
