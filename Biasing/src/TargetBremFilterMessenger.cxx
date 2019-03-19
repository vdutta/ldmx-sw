/**
 * @file TargetBremFilterMessenger.cxx
 * @brief Messenger for setting parameters on TargetBremFilter.
 * @author Omar Moreno, SLAC National Accelerator Laboratory
 */

#include "Biasing/TargetBremFilterMessenger.h"

//------------//
//   Geant4   //
//------------//
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"

//----------------//
//   C++ StdLib   //
//----------------//
#include <iostream>

//-------------//
//   ldmx-sw   //
//-------------//
#include "Biasing/TargetBremFilter.h"

namespace ldmx { 
    
    TargetBremFilterMessenger::TargetBremFilterMessenger(TargetBremFilter* filter) :
        UserActionPluginMessenger(filter), filter_(filter) {

            bremEnergyThresholdCmd_ 
                = std::make_unique<G4UIcmdWithAString>(std::string(getPath() + "brem_threshold").c_str(), this);
            bremEnergyThresholdCmd_->AvailableForStates(G4ApplicationState::G4State_PreInit,
                                              G4ApplicationState::G4State_Idle);
            bremEnergyThresholdCmd_->SetGuidance("Minium energy that the brem electron should have."); 

            killRecoilCmd_ 
                = std::make_unique<G4UIcmdWithoutParameter>(std::string(getPath() + "kill_recoil").c_str(), this);
            killRecoilCmd_->AvailableForStates(G4ApplicationState::G4State_PreInit,
                    G4ApplicationState::G4State_Idle);
            killRecoilCmd_->SetGuidance("Enable killing of the electron track that produces the brem."); 

            volumeCmd_ = std::make_unique<G4UIcmdWithAString>(std::string(getPath() + "volume").c_str(), this);
            volumeCmd_->AvailableForStates(G4ApplicationState::G4State_PreInit,
                                           G4ApplicationState::G4State_Idle);
            volumeCmd_->SetGuidance("Volume to apply the filter to.");     
    }

    TargetBremFilterMessenger::~TargetBremFilterMessenger() {
    }

    void TargetBremFilterMessenger::SetNewValue(G4UIcommand* command, G4String newValue) {
        
        // Handles verbose command.
        UserActionPluginMessenger::SetNewValue(command, newValue);

        if (command == killRecoilCmd_.get()) { 
            filter_->setKillRecoilElectron(true); 
        } else if (command == volumeCmd_.get()) {
            filter_->setVolume(newValue);
        } else if (command == bremEnergyThresholdCmd_.get()) {
            filter_->setBremEnergyThreshold(G4UIcommand::ConvertToDouble(newValue));
        }
    }
}
