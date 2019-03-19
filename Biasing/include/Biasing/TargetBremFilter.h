/**
 * @file TargetBremFilter.h
 * @class TargetBremFilter
 * @brief Class defining a UserActionPlugin that allows a user to filter out 
 *        events that don't result in a brem within the target.
 * @author Omar Moreno, SLAC National Accelerator Laboratory
 */

#ifndef _BIASING_TARGET_BREM_FILTER_H_
#define _BIASING_TARGET_BREM_FILTER_H_

//----------------//
//   C++ StdLib   //
//----------------//
#include <algorithm>

//-------------//
//   ldmx-sw   //
//-------------//
#include "SimPlugins/UserActionPlugin.h"

namespace ldmx {

    class TargetBremFilterMessenger; 

    class TargetBremFilter : public UserActionPlugin {

        public:

            /** Class constructor. */
            TargetBremFilter();

            /** Class destructor. */
            ~TargetBremFilter();

            /**
             * Get the name of the plugin.
             * @return The name of the plugin.
             */
            inline std::string getName() { return "TargetBremFilter"; }

            /**
             * Get whether this plugin implements the event action.
             * @return True if the plugin implements the event action.
             */
            inline bool hasEventAction() { return true; }

            /**
             * Get whether this plugin implements the stepping action.
             * @return True to indicate this plugin implements the stepping action.
             */
            inline bool hasSteppingAction() { return true; }

            /**
             * Get whether this plugin implements the stacking aciton.
             * @return True to indicate this plugin implements the stacking action.
             */
            inline bool hasStackingAction() { return true; }

            /**
             * Implement the stepping action which performs the target volume biasing.
             * @param step The Geant4 step.
             */
            void stepping(const G4Step* step);

            /**
             * End of event action.
             * @param event The Geant4 event.
             */
            void endEvent(const G4Event* event);

            /**
             * Classify a new track which postpones track processing.
             * Track processing resumes normally if a target PN interaction occurred.
             * @param aTrack The Geant4 track.
             * @param currentTrackClass The current track classification.
             */
            G4ClassificationOfNewTrack stackingClassifyNewTrack(const G4Track* aTrack, 
                    const G4ClassificationOfNewTrack& currentTrackClass);

            /** 
             * Enable/disable killing of the recoil electron track.  If the 
             * recoil track is killed, only the brem gamma is propagated.
             */
            void setKillRecoilElectron(bool killRecoilElectron) { 
                killRecoilElectron_ = killRecoilElectron; 
            }

            /** 
             * @param volume Set the volume that the filter will be applied to. 
             */
            void setVolume(std::string volumeName) { volumeName_ = volumeName; }; 

            /**
             * Set the minimum energy that the brem gamma must have.
             */
            void setBremEnergyThreshold(double bremEnergyThreshold) { 
                bremEnergyThreshold_ = bremEnergyThreshold; 
            }

        private:
            
            /** Messenger used to pass arguments to this class. */
            TargetBremFilterMessenger* messenger_{nullptr};

            /** The volume that the filter will be applied to. */
            G4String volumeName_{"target_PV"};

            /** Brem gamma energy treshold. */
            double bremEnergyThreshold_{0}; 

            /** Flag indicating if the recoil electron track should be killed. */
            bool killRecoilElectron_{false};

            /** Flag denoting that an event has a brem candidate. */
            bool hasBremCandidate_{false};

            /** Verbosity level. */
            bool verbose_{0};  

    }; // TargetBremFilter
}

#endif // BIASING_TARGETBREMFILTER_H__
