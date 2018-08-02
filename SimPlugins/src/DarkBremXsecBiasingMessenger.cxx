#include "SimPlugins/DarkBremXsecBiasingMessenger.h"

// LDMX
#include "SimPlugins/DarkBremXsecBiasingPlugin.h"

// STL
#include <string>

namespace ldmx {

    DarkBremXsecBiasingMessenger::DarkBremXsecBiasingMessenger(DarkBremXsecBiasingPlugin* biasingPlugin) :
            UserActionPluginMessenger(biasingPlugin), biasingPlugin_(biasingPlugin) {

        xsecFactorCmd_ = new G4UIcommand(std::string(getPath() + "xsecFactor").c_str(), this);
        G4UIparameter* modulus = new G4UIparameter("xsecFactor", 'd', false);
        xsecFactorCmd_->SetParameter(modulus);
        xsecFactorCmd_->SetGuidance("Set the cross section biasing factor for the Dark Brem process.");
        xsecFactorCmd_->AvailableForStates(G4ApplicationState::G4State_PreInit, G4ApplicationState::G4State_Idle);
    }

    void DarkBremXsecBiasingMessenger::SetNewValue(G4UIcommand *command, G4String newValue) {

        // Handles verbose command.
        UserActionPluginMessenger::SetNewValue(command, newValue);

        if (command == xsecFactorCmd_) {
            biasingPlugin_->setXsecBiasingFactor(std::stod(newValue));
        }
    }

} // namespace sim
