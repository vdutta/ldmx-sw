/*
 * @file HighPassFilter.cxx
 * @brief User action plugin that informs Geant4 to first process particles with an 
 *        energy above a certain (settable) threshold. Also enables removing particles below 
 *        another (settable) energy threshold from event processing. 
 * @author Lene Kristian Bryngemark, Lund University (based on Omar Moreno's code)
 *         
 */

#include "Biasing/HighPassFilter.h"

SIM_PLUGIN(ldmx, HighPassFilter)

namespace ldmx { 

  HighPassFilter::HighPassFilter() {
    messenger_ = new HighPassFilterMessenger(this);
  }


  // ----------------------------------------------------------------------------------------



  HighPassFilter::~HighPassFilter() {
  }

  // ----------------------------------------------------------------------------------------



  G4ClassificationOfNewTrack HighPassFilter::stackingClassifyNewTrack(
									 const G4Track* track, 
									 const G4ClassificationOfNewTrack& currentTrackClass) {

    if (verboseRun) {
      std::cout << "********************************" << std::endl;
      std::cout << "*   Track pushed to the stack  *" << std::endl;
      std::cout << "********************************" << std::endl;

      // get the PDGID of the track.
      G4int pdgID = track->GetParticleDefinition()->GetPDGEncoding();
      // Get the particle type.
      G4String particleName = track->GetParticleDefinition()->GetParticleName();
      //      std::cout << "[ TargetBremFilter ]: " << "\n" 
      std::cout << "[ HighPassFilter ]: " << "\n" 
		<< "\tParticle " << particleName << " ( PDG ID: " << pdgID << " ) : " << "\n"
		<< "\tTrack ID: " << track->GetTrackID() << "\n" 
		<< std::endl;
    }

    if (track == currentTrack_) {
      currentTrack_ = nullptr; 
      return fWaiting; 
    }

    // Use current classification by default so values from other plugins are not overridden.
    G4ClassificationOfNewTrack classification = currentTrackClass;
        
    return classification;
  }


  // ----------------------------------------------------------------------------------------


  
  void HighPassFilter::stepping(const G4Step* step) { 
    
    if (! verboseRun) {
      std::cout << "*** verboseRun not set ********" << std::endl;
    }
    
    //akip events with no interaction in the target 
    // TODO: make more general and remove this requirement?
    if (TargetBremFilter::getBremGammaList().empty()) { 
      return;
    } 
    
    // Get the track associated with this step.
    G4Track* track = step->GetTrack();
    
    // Only process tracks whose parent is the primary particle
    if (track->GetParentID() != 1) return; 
    
    // get the PDGID of the track.
    G4int pdgID = track->GetParticleDefinition()->GetPDGEncoding();
    
    // Make sure that the particle being processed is an electron. ????? electron is 11. gamma?
    // TODO: At some point all particle types should be allowed.
    if (pdgID != 22) {
      if (verboseRun) {
	std::cout << "[ HighPassFilter ]: " << "\n"
		  << "\tSkipping particle with wrong PDG ID: " <<  pdgID 
		  << std::endl;
      }
      return; 
    }

    // Get the volume the particle is in.
    G4VPhysicalVolume* volume = track->GetVolume();
    G4String volumeName = volume->GetName();

    if (verboseRun) {
  
      std::cout << "*******************************" << std::endl; 
      std::cout << "*   Step " << track->GetCurrentStepNumber() << std::endl;
      std::cout << "********************************" << std::endl;
    }


    if (verboseRun) {
      // Get the particle type.
      G4String particleName = track->GetParticleDefinition()->GetParticleName();
      // Get the kinetic energy of the particle.
      double incidentParticleEnergy = step->GetPreStepPoint()->GetTotalEnergy();

      std::cout << "[ HighPassFilter ]:\n" 
		<< "\tTotal energy of " << particleName  << ": " << incidentParticleEnergy << " MeV"
		<< "\tPDG ID: " << pdgID 
		<< "\tTrack ID: " << track->GetTrackID() 
		<< "\tStep #: " << track->GetCurrentStepNumber()
		<< "\tParent ID: " << track->GetParentID()   << std::endl;
      //		<< "\tParticle currently in " << volumeName   << std::endl;
    }

    // If the particle isn't in the specified volume, stop processing the 
    // event.
    std::vector<G4Track*> bremGammaList = TargetBremFilter::getBremGammaList();
    if (std::find(std::begin(volumes_), std::end(volumes_), volumeName) == std::end(volumes_)) {

      if (verboseRun) {
	std::cout << "[ HighPassFilter ]: "
		  << "Brem is in " << volumeName  << std::endl;
      }
	  
      // If secondaries were produced outside of the volume of interest, 
      // and there aren't additional brems to process, abort the 
      // event.  Otherwise, suspend the track and move on to the next 
      // brem.
      if (step->GetSecondary()->size() != 0 
	  && (std::find(bremGammaList.begin(), bremGammaList.end(), track) != bremGammaList.end())) { 
	
	if (verboseRun) {
	  std::cout << "[ HighPassFilter ]: "
		    << "Reaction occured outside volume of interest " ; //--> Aborting event." 
	    //	    << std::endl;
	}
	
	if (bremGammaList.size() == 1) { 
	  track->SetTrackStatus(fKillTrackAndSecondaries);
	  G4RunManager::GetRunManager()->AbortEvent();
	  currentTrack_ = nullptr;
	  if (verboseRun) {
	    std::cout << "--> Aborting event."
                    << std::endl;
	  }
	  return;
	} else {

	  currentTrack_ = track; 
	  track->SetTrackStatus(fSuspend);
	  TargetBremFilter::removeBremFromList(track);
	  if (verboseRun) {
	    std::cout << "--> Removing brem with track ID " << track->GetTrackID() << " from list."
		      << std::endl;
          }
	  return;
	}
      }
      return;
    }
    else
      if (verboseRun) {

	std::cout << "[ HighPassFilter ]: "
		  << "No suitable volume match for track with "
		  << "\tPDG ID: " << pdgID
		  << "\tTrack ID: " << track->GetTrackID() << std::endl;
      }
    // The list of brems will only contain a given track/particle if it 
    // originates from the target.  If the gamma originates elsewhere, 
    // suspend it and move on to the next gamma.
    // TODO: At some point, this needs to be modified to include brems 
    // from downstream of the target.
    if (std::find(bremGammaList.begin(), bremGammaList.end(), track) == bremGammaList.end()) { 
           
      if (verboseRun) {
 
	std::cout << "[ HighPassFilter ]: "
		  << "Brem list doesn't contain track." << std::endl;
      }
            
      currentTrack_ = track; 
      track->SetTrackStatus(fSuspend);
      return;
    }
 
    // Get the particles daughters.
    const G4TrackVector* secondaries = step->GetSecondary();


    // If the particle doesn't interact, then move on to the next step.
    if (secondaries->size() == 0) { 
            
      if (verboseRun) {
	std::cout << "[ HighPassFilter ]: "
		  << "Brem photon did not interact --> Continue propagating track."
		  << std::endl;
      }
        
      // If the particle is exiting the bounding volume, kill it.
      if (!boundVolumes_.empty() && step->GetPostStepPoint()->GetStepStatus() == fGeomBoundary) {
	if (std::find(std::begin(boundVolumes_), std::end(boundVolumes_), volumeName) != std::end(boundVolumes_)) {

           
	  if (verboseRun) {
	    std::cout << "[ HighPassFilter ]: "
		      << "Brem photon is exiting the volume --> particle will be killed or suspended."
		      << std::endl;
	  }
                    
	  if (bremGammaList.size() == 1) { 
	    track->SetTrackStatus(fKillTrackAndSecondaries);
	    G4RunManager::GetRunManager()->AbortEvent();
	    currentTrack_ = nullptr;
	    if (verboseRun) {
	      std::cout << "[ HighPassFilter ]: " 
			<< " Brem list is empty --> Killing all tracks!"
			<< std::endl;
	    }
	    return;
	  } else { 
	    currentTrack_ = track; 
	    track->SetTrackStatus(fSuspend);
	    TargetBremFilter::removeBremFromList(track);
	    if (verboseRun) {
	      std::cout << "[ HighPassFilter ]: " 
			<< " Other tracks still need to be processed --> Suspending track!"
			<< std::endl;
	    }
	    return;
	  }
	}
      }
	    
    } else {

      // If the brem gamma interacts and produces secondaries, get the 
      // process used to create them. 
      G4String processName = secondaries->at(0)->GetCreatorProcess()->GetProcessName(); 
            

      if (verboseRun) {
	std::cout << "[ HighPassFilter ]: "
		  << "Brem photon produced " << secondaries->size() 
		  << " particles via " << processName << " process." 
		  << std::endl;
      }

      // Only record the process that is being biased
      if (!processName.contains(BiasingMessenger::getProcess())) {

	if (verboseRun) {
	  std::cout << "[ HighPassFilter ]: "
		    << "Process was not " << BiasingMessenger::getProcess() 
		    << std::endl;
	}
                
	if (bremGammaList.size() == 1) { 
	  track->SetTrackStatus(fKillTrackAndSecondaries);
	  G4RunManager::GetRunManager()->AbortEvent();
	  currentTrack_ = nullptr;
	  if (verboseRun) {
	    std::cout << "[ HighPassFilter ]: " 
		      << " Brem list is empty --> Killing all tracks!"
		      << std::endl;
	  }
	  return;
	} else { 
	  currentTrack_ = track; 
	  track->SetTrackStatus(fSuspend);
	  TargetBremFilter::removeBremFromList(track);
	  if (verboseRun) {
	    std::cout << "[ HighPassFilter ]: " 
		      << " Other tracks still need to be processed --> Suspending track!"
		      << std::endl;
	  }
	  return;
	}
      }

      //not returned yet -- means we are looking at the process of interest. so keep it            
      std::cout << "[ HighPassFilter ]: "
		<< "Brem photon produced " << secondaries->size() 
		<< " particles via " << processName << " process." 
		<< std::endl;
      std::cout << "Track ID: " << track->GetTrackID() << "\tEvent weight: " << track->GetWeight() << std::endl;

      TargetBremFilter::removeBremFromList(track);
      BiasingMessenger::setEventWeight(track->GetWeight());
      photonGammaID_ = track->GetTrackID(); 
    } 
  }


  // ----------------------------------------------------------------------------------------



  void HighPassFilter::postTracking(const G4Track* track) { 
    
    if (! verboseRun) {
      std::cout << "*** verboseRun not set ********" << std::endl;
    }
    
    //parent id is 1 for the incoming electron, 2 for the brem gamma       
    if (track->GetParentID() == photonGammaID_) { 
      UserTrackInformation* userInfo 
	= dynamic_cast<UserTrackInformation*>(track->GetUserInformation());
      userInfo->setSaveFlag(true); 
      G4ThreeVector pvec = track->GetMomentum();
      if (verboseRun) {
	// get the PDGID of the track.
	G4int pdgID = track->GetParticleDefinition()->GetPDGEncoding();
	std::cout << "[ HighPassFilter ]:\n" 
		  << "\tPDG ID: " << pdgID 
		  << "\tTrack ID: " << track->GetTrackID() 
		  << "\tStep #: " << track->GetCurrentStepNumber()
		  << "\tParent ID: " << track->GetParentID() //<< "\n"
		  << "\t p: [ " << pvec[0] << ", " << pvec[1] << ", " << pvec[2] << " ]" << std::endl;
      }
    }
  }


  // ----------------------------------------------------------------------------------------



  void HighPassFilter::setStackThreshold(std::string stackThreshold) { 
        
    std::cout << "[ HighPassFilter ]: Set filter to suspend particles below " << stackThreshold << " MeV until particles above threshold are processed" << std::endl;
    stackEnergyThreshold_=std::atof(stackThreshold.c_str());
  }        
  

  // ----------------------------------------------------------------------------------------


  void HighPassFilter::setKillThreshold(std::string killThreshold) { 
        
    std::cout << "[ HighPassFilter ]: Set filter to kill particles below " << killThreshold << " MeV " << std::endl;
    killEnergyThreshold_=std::atof(killThreshold.c_str());
    //      killThresholds_.push_back(killThreshold);
  }        

  // ----------------------------------------------------------------------------------------


  void HighPassFilter::setVerbose(std::string verbose) { 
        
    std::cout << "[ HighPassFilter ]: Setting filter verbosity to " << verbose << std::endl;
    //accept any positive number or any case version of TRUE
    if ( std::atoi(verbose.c_str()) > 0 || verbose.compare("true") == 0 || verbose.compare("TRUE") == 0 || verbose.compare("True") == 0 )
      verboseRun=true;
    else 
      verboseRun=false;
  }        
  


  void HighPassFilter::addVolume(std::string volume) {

    if (! verboseRun) {
      std::cout << "*** verboseRun not set ********" << std::endl;
    }

    std::cout << "[ HighPassFilter ]: Applying filter to volume " << volume << std::endl;
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

  void HighPassFilter::addBoundingVolume(std::string volume) {
    std::cout << "[ HighPassFilter ]: Bounding particle to volume " << volume << std::endl;
    boundVolumes_.push_back(volume);
  }




}
