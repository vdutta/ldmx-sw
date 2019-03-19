/**
 * @file TargetBremFilterMessenger.h
 * @brief Messenger for setting parameters on TargetBremFilter.
 * @author Omar Moreno, SLAC National Accelerator Laboratory
 */

#ifndef BIASING_TARGETBREMFILTERMESSENGER_H
#define BIASING_TARGETBREMFILTERMESSENGER_H

//-------------//
//   ldmx-sw   //
//-------------//
#include "SimPlugins/UserActionPluginMessenger.h"

// Forward declarations
class G4UIcmdWithoutParameter; 
class G4UIcmdWithAString; 

namespace ldmx { 
   
    // Forward declare to avoid circular dependency in headers
    class TargetBremFilter;

    class TargetBremFilterMessenger : public UserActionPluginMessenger {
        
        public: 

            /** 
             * Constructor
             *
             * @param Filter associated with this messenger.
             */
            TargetBremFilterMessenger(TargetBremFilter* filter); 

            /** Destructor */
            ~TargetBremFilterMessenger();

            /**
             */
            void SetNewValue(G4UIcommand * command, G4String newValue);

        private:

            /** The filter associated with this messenger. */
            TargetBremFilter* filter_{nullptr};
            
            /**
             * Command allowing a user to specify the minimum energy that the
             * brem gamma must have. 
             */
            std::unique_ptr<G4UIcmdWithAString> bremEnergyThresholdCmd_;

            /** Command dictating whether the electron track gets killed. */
            std::unique_ptr<G4UIcmdWithoutParameter> killRecoilCmd_;

            /** 
             * Command allowing a user to specify what volume the filter 
             * should be applied to.
             */
            std::unique_ptr<G4UIcmdWithAString> volumeCmd_;

    }; // TargetBremFilterMessenger

} // ldmx

#endif // BIASING_TARGETBREMFILTERMESSENGER_H

