/**
 * @file HighPassFilterMessenger.h
 * @brief Messenger for setting parameters on HighPassFilter.
 * @author Lene Kristian Bryngemark, Lund University
 */

#ifndef BIASING_HIGHPASSFILTERMESSENGER_H
#define BIASING_HIGHPASSFILTERMESSENGER_H

//------------//
//   Geant4   //
//------------//
#include "G4UIcmdWithAString.hh"

//-------------//
//   ldmx-sw   //
//-------------//
#include "Biasing/HighPassFilter.h"
#include "SimPlugins/UserActionPluginMessenger.h"

namespace ldmx { 
   
    // Forward declare to avoid circular dependency in headers
    class HighPassFilter;

    class HighPassFilterMessenger : UserActionPluginMessenger {
        
        public: 

            /** 
             * Constructor
             *
             * @param Filter associated with this messenger.
             */
            HighPassFilterMessenger(HighPassFilter* filter); 

            /** Destructor */
            ~HighPassFilterMessenger();

            /**
             */
            void SetNewValue(G4UIcommand * command, G4String newValue);

        private:

            /** The filter associated with this messenger. */
            HighPassFilter* filter_{nullptr};
            
            /** 
             * Command allowing a user to specify a lower energy threshold
             * below which the particle will be killed.
             */
            G4UIcmdWithAString* killThresholdCmd_{nullptr};

            /** 
             * Command allowing a user to specify a lower energy threshold
             * below which the particle will be suspended (sent to stack) until 
             * particles above threshold have been processed.
             */
            G4UIcmdWithAString* stackThresholdCmd_{nullptr};

            /** 
             * Command allowing a user to specify verbose output level.
             */
            G4UIcmdWithAString* verboseCmd_{nullptr};

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




    }; // HighPassFilterMessenger
}

#endif // BIASING_HIGHPASSFILTERMESSENGER_H

