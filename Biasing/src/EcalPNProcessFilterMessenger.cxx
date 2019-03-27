/**
 * @file EcalPNProcessFilterMessenger.cxx
 * @brief Messenger for setting parameters on EcalPNProcessFilter.
 * @author Joshua Hiltbrand
 *         University of Minnesota 
 */

#include "Biasing/EcalPNProcessFilterMessenger.h"
#include <iostream>

namespace ldmx { 
    
    EcalPNProcessFilterMessenger::EcalPNProcessFilterMessenger(EcalPNProcessFilter* filter) :
        UserActionPluginMessenger(filter), filter_(filter) {

            boundCmd_ 
                = new G4UIcmdWithAString{std::string("/ldmx/plugins/EcalPNProcessFilter/bound_volume").c_str(), this};
            boundCmd_->AvailableForStates(G4ApplicationState::G4State_PreInit,
                    G4ApplicationState::G4State_Idle);
            boundCmd_->SetGuidance("Bound a particle to the given volume."); 


            volumeCmd_ = new G4UIcmdWithAString{std::string("/ldmx/plugins/EcalPNProcessFilter/volume").c_str(), this};
            volumeCmd_->AvailableForStates(G4ApplicationState::G4State_PreInit,
                                           G4ApplicationState::G4State_Idle);
            volumeCmd_->SetGuidance("Volume to apply the filter to. Note that multiple volumes may be added.");     
            energyFractionThresholdCmd_ = new G4UIcmdWithAString{std::string("/ldmx/plugins/EcalPNProcessFilter/energyFractionThreshold").c_str(), this};
            energyFractionThresholdCmd_->AvailableForStates(G4ApplicationState::G4State_PreInit, G4ApplicationState::G4State_Idle);
            energyFractionThresholdCmd_->SetGuidance("Fraction of incident KE that neutron and kaon daughters must carry.");
    }

    EcalPNProcessFilterMessenger::~EcalPNProcessFilterMessenger() {
        delete boundCmd_;
        delete volumeCmd_; 
        delete energyFractionThresholdCmd_;
    }

    void EcalPNProcessFilterMessenger::SetNewValue(G4UIcommand* command, G4String newValue) {

        // Handles verbose command.
        UserActionPluginMessenger::SetNewValue(command, newValue);

        if (command == volumeCmd_) filter_->addVolume(newValue);
        else if (command == boundCmd_) filter_->addBoundingVolume(newValue); 
        else if (command == energyFractionThresholdCmd_) filter_->setEnergyFractionThreshold(G4UIcommand::ConvertToDouble(newValue));
    }
}
