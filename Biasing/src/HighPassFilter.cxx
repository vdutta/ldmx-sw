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
    delete messenger_;
  }

  // ----------------------------------------------------------------------------------------



  G4ClassificationOfNewTrack HighPassFilter::stackingClassifyNewTrack(
									 const G4Track* track, 
									 const G4ClassificationOfNewTrack& currentTrackClass) {

    //    G4int tID = track->GetTrackID();
    // Use current classification by default so values from other plugins are not overridden.                                                                                                                
    G4ClassificationOfNewTrack classification = currentTrackClass;

    //don't send incoming electron to the waiting stack...
    if ( track->GetTrackID() < 2 ) {
      //take this opportunity to reset the parent id book keeping. so this has to be done in the first if!!
      firstWaitingParent_=1000000;
      stopStacking=false;
      return classification;
    }

    if ( stopStacking )
      return classification;

    if ( track->GetParentID() > firstWaitingParent_ ) { 
      //the parent of this track has already been pushed to stack so let's assume the important stack sorting, and killing, happened in the previous step
      if (verboseRun) {
	std::cout << "The parent (ID: " << track->GetParentID() << ") of track " << track->GetTrackID() << " has already been classified before; won't touch" << std::endl;
      }
      stopStacking=true; //this is set the first time we get downstream of the stacked generation
      return classification;

    }


    // get the PDGID of the track.                                                                                                                            
    G4int pdgID = track->GetParticleDefinition()->GetPDGEncoding();
    
    // filter on kinetic energy for p and n (since they're stable particles), total otherwise
    G4double E = (pdgID == 2212 || pdgID == 2112) ? track->GetKineticEnergy() : track->GetTotalEnergy() ;


    if (verboseRun && classification != 1) {
      std::cout << "********************************" << std::endl;
      std::cout << "*   Track pushed to the stack  *" << std::endl;
      std::cout << "********************************" << std::endl;

      // Get the particle type and momentum.
      G4String particleName = track->GetParticleDefinition()->GetParticleName();
      G4ThreeVector pvec = track->GetMomentum();
      std::cout << "[ HighPassFilter ]: Decision on " //track with \n"
		<< " Particle " << particleName << " (PDG ID: " << pdgID << ")" //<< "\n"
		<< "\tstack class: " << currentTrackClass
		<< "\tTrack ID: " << track->GetTrackID()
		<< "\tStep #: " << track->GetCurrentStepNumber()
		<< "\tParent ID: " << track->GetParentID() //<< "\n"                                                                                        
		<< "\tp = (" << pvec[0] << ", " << pvec[1] << ", " << pvec[2] << "):\t" ; //<< std::endl;
    }


    //TODO: could use    processName.contains(BiasingMessenger::getProcess())
    //parent id is 1 for the incoming electron, (at least) 2 for the brem gamma
    if ( track->GetParentID() > 1 && (track->GetCreatorProcess()->GetProcessName()).contains("photonNuclear") ){

      if ( E < killEnergyThreshold_ ) {
	// remove track and daughters if energy is below threshold
	if (verboseRun) {
	  
	  if (classification == 1) { //just to make sure... perhaps some which were waiting are getting killed
	    std::cout << "********************************" << std::endl;
	    std::cout << "*   Killing track  *" << std::endl;
	    std::cout << "********************************" << std::endl;
	    
	    // Get the particle type and momentum.                                                                                                                                                                 
	    G4String particleName = track->GetParticleDefinition()->GetParticleName();
	    G4ThreeVector pvec = track->GetMomentum();
	    std::cout << "[ HighPassFilter ]: Decision on " //track with \n"                                                                                                                                       
		      << " Particle " << particleName << " (PDG ID: " << pdgID << "), " //<< "\n" 
		      << "\t stack class: " << currentTrackClass
		      << "\tTrack ID: " << track->GetTrackID()
		      << "\tStep #: " << track->GetCurrentStepNumber()
		      << "\tParent ID: " << track->GetParentID() //<< "\n"                                                                                                                                         
		      << "\t p = ( " << pvec[0] << ", " << pvec[1] << ", " << pvec[2] << " ) :\t" ; //<< std::endl;                                                                                                
	  } //if classification = 1 (fWaiting)


	  std::cout  << " --- kill (kill threshold)"  << std::endl;
	  std::cout << "[ HighPassFilter ]: " << "\n"
		    << "\tfound sth photo sth "
		    << std::endl;
	
	}

	return fKill;
      } // if below thr                                                                                                                                         

    } //if brem gamma parent                                                   


    //if we don't want to kill it, and it's already waiting, leave it at that
    if ( classification == 1 ) {
      //      if (verboseRun) std::cout  << "\t------ waiting (already set)"  << std::endl;
      return fWaiting; 
    }

    //if we don't want to kill it, could still want to put it on waiting stack
    if ( E < stackEnergyThreshold_ ) {
      if (verboseRun) std::cout  << " --- waiting (stacking threshold)"  << std::endl;

      firstWaitingParent_=track->GetParentID();
      return fWaiting;
    }


    // return current classification by default so values from other plugins are not overridden.
    // probably we rarely get this far
    if (verboseRun) {
      if ( classification == 0 ) std::cout  << " --- urgent (unchanged) "  << std::endl;
      else std::cout  << " --- " << classification << " (unchanged) "  << std::endl;

      if (track->GetParentID() > 0 ) //don't ask this for the incoming electron, which doesn't have a process that created it
	std::cout << "[ HighPassFilter ]: " << "\n"
		  << "\tWhat I found wasn't sth photo sth , it was "
		  << track->GetCreatorProcess()->GetProcessName()	
		  << std::endl;
    }

    return classification;
  }




  // ----------------------------------------------------------------------------------------



  //void HighPassFilter::preTracking(const G4Track* track) { 
    void HighPassFilter::PreTracking( G4Track* track) { 
    
    
    // get the PDGID of the track.
    G4int pdgID = track->GetParticleDefinition()->GetPDGEncoding();

    // kinetic energy for p and n, total otherwise)
    G4double E = (pdgID == 2212 || pdgID == 2112) ? track->GetKineticEnergy() : track->GetTotalEnergy() ;
    // remove track and daughters if energy is below threshold
    if ( E < stackEnergyThreshold_ ) {
      track->SetTrackStatus(fSuspend); //fWaiting); //
   
      if (verboseRun) {
	G4ThreeVector pvec = track->GetMomentum();
	std::cout << "[ HighPassFilter ]: Sent track to waiting stack: \n" 
		  << "\tPDG ID: " << pdgID 
		  << "\tTrack ID: " << track->GetTrackID() 
		  << "\tStep #: " << track->GetCurrentStepNumber()
		  << "\tParent ID: " << track->GetParentID() //<< "\n"
		  << "\t p: [ " << pvec[0] << ", " << pvec[1] << ", " << pvec[2] << " ]" << std::endl;
      } //if verbose

      return;
	
    } // if below thr
    //    else {
      if (verboseRun) {
	G4ThreeVector pvec = track->GetMomentum();
	std::cout << "[ HighPassFilter ]: Not sent track to waiting stack: \n" 
		  << "\tPDG ID: " << pdgID 
		  << "\tTrack ID: " << track->GetTrackID() 
		  << "\tStep #: " << track->GetCurrentStepNumber()
		  << "\tParent ID: " << track->GetParentID() //<< "\n"
		  << "\t p: [ " << pvec[0] << ", " << pvec[1] << ", " << pvec[2] << " ]" << std::endl;
      } //if verbose
	
      //  }// above threshold
  }



    // ----------------------------------------------------------------------------------------

     void HighPassFilter::postTracking(const G4Track* track) { 


       //       listOfPushedTrackIDs_.clear();
       // listOfPushedTrackIDs_.resize(0);


     }

    // ----------------------------------------------------------------------------------------


      void HighPassFilter::PostTracking( G4Track* track) { 
    
    
      //parent id is 1 for the incoming electron, 2 for the brem gamma       
      if ( track->GetCreatorProcess()->GetProcessName().compare("photo") == 0 ){

	//    if (track->GetParentID() == photonGammaID_) { 

	// get the PDGID of the track.
	G4int pdgID = track->GetParticleDefinition()->GetPDGEncoding();

	// kinetic energy for p and n, total otherwise)
	G4double E = (pdgID == 2212 || pdgID == 2112) ? track->GetKineticEnergy() : track->GetTotalEnergy() ;
	// remove track and daughters if energy is below threshold
	if ( E < killEnergyThreshold_ ) {
	  track->SetTrackStatus(fKillTrackAndSecondaries);
   
	  if (verboseRun) {
	    G4ThreeVector pvec = track->GetMomentum();
	    std::cout << "[ HighPassFilter ]: Killed track with \n" 
		      << "\tPDG ID: " << pdgID 
		      << "\tTrack ID: " << track->GetTrackID() 
		      << "\tStep #: " << track->GetCurrentStepNumber()
		      << "\tParent ID: " << track->GetParentID() //<< "\n"
		      << "\t p: [ " << pvec[0] << ", " << pvec[1] << ", " << pvec[2] << " ]" << std::endl;
	  } //if verbose

	} // if below thr
	else { 
	  if (verboseRun) {
	    G4ThreeVector pvec = track->GetMomentum();
	    std::cout << "[ HighPassFilter ]: Didn't kill track with \n" 
		      << "\tPDG ID: " << pdgID 
		      << "\tTrack ID: " << track->GetTrackID() 
		      << "\tStep #: " << track->GetCurrentStepNumber()
		      << "\tParent ID: " << track->GetParentID() //<< "\n"
		      << "\t p: [ " << pvec[0] << ", " << pvec[1] << ", " << pvec[2] << " ]" << std::endl;
	  } //if verbose
	} //if above threshold

      } //if brem gamma parent

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
