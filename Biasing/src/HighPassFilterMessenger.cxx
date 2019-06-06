/**
 * @file HighPassFilterMessenger.cxx
 * @brief Messenger for setting parameters on HighPassFilter.
 * @author Lene Kristian Bryngemark, Lund University
 */

#include "Biasing/HighPassFilterMessenger.h"
#include <iostream>

namespace ldmx { 
    
    HighPassFilterMessenger::HighPassFilterMessenger(HighPassFilter* filter) :
        UserActionPluginMessenger(filter), filter_(filter) {

            killThresholdCmd_ 
                = new G4UIcmdWithAString{std::string(getPath() + "killThreshold").c_str(), this};
            killThresholdCmd_->AvailableForStates(G4ApplicationState::G4State_PreInit,
                    G4ApplicationState::G4State_Idle);
            killThresholdCmd_->SetGuidance("All particles below this energy will be discarded from processing. Defaults to 0."); 


            stackThresholdCmd_ = new G4UIcmdWithAString{std::string(getPath() + "stackThreshold").c_str(), this};
            stackThresholdCmd_->AvailableForStates(G4ApplicationState::G4State_PreInit,
                                           G4ApplicationState::G4State_Idle);
            stackThresholdCmd_->SetGuidance("Particles below this energy threshold are sent to the stack until particles above threshold are processed. Defaults to 0.");     

            verboseCmd_ = new G4UIcmdWithAString{std::string(getPath() + "verbosity").c_str(), this};
            verboseCmd_->AvailableForStates(G4ApplicationState::G4State_PreInit,
                                           G4ApplicationState::G4State_Idle);
            verboseCmd_->SetGuidance("Set verbose output level by passing positive value or \"true\". Defaults to false.");     //TODO: default to false

	    
	    boundCmd_
	      = new G4UIcmdWithAString{std::string(getPath() + "bound_volume").c_str(), this};
	    boundCmd_->AvailableForStates(G4ApplicationState::G4State_PreInit,
					  G4ApplicationState::G4State_Idle);
	    boundCmd_->SetGuidance("Bound a particle to the given volume.");
	    
	    
	    volumeCmd_ = new G4UIcmdWithAString{std::string(getPath() + "volume").c_str(), this};
	    volumeCmd_->AvailableForStates(G4ApplicationState::G4State_PreInit,
					   G4ApplicationState::G4State_Idle);
	    volumeCmd_->SetGuidance("Volume to apply the filter to. Note that multiple volumes may be added.");
	    


    }

    HighPassFilterMessenger::~HighPassFilterMessenger() {
        delete killThresholdCmd_; 
        delete stackThresholdCmd_; 
        delete verboseCmd_; 
	//	delete volumeCmd_;
	//	delete boundCmd_;
 }

    void HighPassFilterMessenger::SetNewValue(G4UIcommand* command, G4String newValue) {

        // Handles passed commands.
        UserActionPluginMessenger::SetNewValue(command, newValue);

        if (command == stackThresholdCmd_) filter_->setStackThreshold(newValue);
	else if (command == killThresholdCmd_) filter_->setKillThreshold(newValue);
        else if (command == verboseCmd_) filter_->setVerbose( newValue );
	else if (command == volumeCmd_) filter_->addVolume(newValue);
        else if (command == boundCmd_) filter_->addBoundingVolume(newValue);
 
  }
}
