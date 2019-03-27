/**
 * @file EcalPNProcessFilterMessenger.h
 * @brief Messenger for setting parameters on EcalPNProcessFilter.
 * @author Joshua Hiltbrand
 *         University of Minnesota
 */

#ifndef BIASING_ECALPNPROCESSFILTERMESSENGER_H
#define BIASING_ECALPNPROCESSFILTERMESSENGER_H

//------------//
//   Geant4   //
//------------//
#include "G4UIcmdWithAString.hh"

//-------------//
//   ldmx-sw   //
//-------------//
#include "Biasing/EcalPNProcessFilter.h"
#include "SimPlugins/UserActionPluginMessenger.h"

namespace ldmx { 
   
    // Forward declare to avoid circular dependency in headers
    class EcalPNProcessFilter;

    class EcalPNProcessFilterMessenger : UserActionPluginMessenger {
        
        public: 

            /** 
             * Constructor
             *
             * @param Filter associated with this messenger.
             */
            EcalPNProcessFilterMessenger(EcalPNProcessFilter* filter); 

            /** Destructor */
            ~EcalPNProcessFilterMessenger();

            /**
             */
            void SetNewValue(G4UIcommand* command, G4String newValue);

        private:

            /** The filter associated with this messenger. */
            EcalPNProcessFilter* filter_{nullptr};
            
            /** 
             * Command allowing a user to specify what volume the filter 
             * should be applied to.
             */
            G4UIcmdWithAString* volumeCmd_{nullptr};

            /** 
             * Command allowing a user to specify whether a particle should 
             * be bound to the specified volume.  If so, once the particle
             * exits the volume it will be killed. 
             */
            G4UIcmdWithAString* boundCmd_{nullptr};

            /**
             * Command allowing a user to specifiy what fraction of the
             * incident pn gamma KE must be carried away by neutrons
             * and kaons.
             */
            G4UIcmdWithAString* energyFractionThresholdCmd_{nullptr};

    }; // EcalPNProcessFilterMessenger
}

#endif // BIASING_ECALPNPROCESSFILTERMESSENGER_H
