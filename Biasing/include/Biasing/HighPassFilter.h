/**
 * @file HighPassFilter.h
 * @brief User action plugin that informs Geant4 to first process particles with an 
 *        energy above a certain (settable) threshold. Also enables removing particles below 
 *        another (settable) energy threshold from event processing.
 * @author Lene Kristian Bryngemark, Lund University (based on Omar Moreno's code)
 */

#ifndef BIASING_HIGHPASSFILTER_H_
#define BIASING_HIGHPASSFILTER_H_

//----------------//
//   C++ StdLib   //
//----------------//
#include <algorithm>

//------------//
//   Geant4   //
//------------//
#include "G4VPhysicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RunManager.hh"

//----------//
//   LDMX   //
//----------//
#include "SimPlugins/UserActionPlugin.h"
#include "Biasing/BiasingMessenger.h"
#include "Biasing/HighPassFilterMessenger.h"
#include "Biasing/TargetBremFilter.h"
#include "SimCore/UserTrackInformation.h"

class UserTrackingAction;

namespace ldmx {

    
    // Forward declaration to avoid circular dependecies
    class HighPassFilterMessenger; 

    class HighPassFilter : public UserActionPlugin {

        public:

            /** Default Ctor */
            HighPassFilter();

            /** Destructor */
            ~HighPassFilter();

            /** @return A std::string descriptor of the class. */
            virtual std::string getName() {
                return "HighPassFilter";
            }

            /**
             * Get whether this plugin implements the stepping action.
             * @return True to indicate this plugin implements the stepping action.
             */
            bool hasSteppingAction() {
                return true;
            }

            /**
             * Get whether this plugin implements the tracking action.
             * @return True if the plugin implements the tracking action.
             */
            bool hasTrackingAction() { 
                return true;
            }

            /**
             * Get whether this plugin implements the stacking aciton.
             * @return True to indicate this plugin implements the stacking action.
             */
            bool hasStackingAction() { 
                return true;
            }

            void stepping(const G4Step* step);

            /**
             * Post-tracking action.
             */
            void postTracking(const G4Track*); 

            /**
             * Classify a new track which postpones track processing.
             * Track processing resumes normally if an energy threshold is reached.
             * @param aTrack The Geant4 track.
             * @param currentTrackClass The current track classification.
             */
            G4ClassificationOfNewTrack stackingClassifyNewTrack(const G4Track* aTrack, const G4ClassificationOfNewTrack& currentTrackClass);

            /** 
             * Add a stackThreshold for the filter.
             *
             * @param StackThreshold name
             */
            void setStackThreshold(std::string stackThreshold);

            /** 
             * Add a killThreshold for the filter.
             *
             * @param KillThreshold name
             */
            void setKillThreshold(std::string killThreshold);

            /** 
             * Set verbose level for the filter.
             *
             * @param Verbose name
             */
            void setVerbose(std::string verbose);

            /** 
             * Add a volume to apply the filter to.
             *
             * @param Volume name
             */
            void addVolume(std::string volume);

            /** 
             * Add a volume to bound the particle of interest to.
             *
             * @param Volume name
             */
            void addBoundingVolume(std::string volume);

        private:

            /** Messenger used to pass arguments to this class. */
            HighPassFilterMessenger* messenger_{nullptr};

            /** Pointer to the current track being processed. */
            G4Track* currentTrack_{nullptr};

            /** List of volumes to apply filter to. */
            std::vector<std::string> volumes_; 

            /** List of volumes to bound the particle to. */
            std::vector<std::string> boundVolumes_;

            /** Filter action energy thresholds */
            double stackEnergyThreshold_{0.}; // MeV
            double killEnergyThreshold_{0.}; // MeV

            /** PN gamma parent ID. */
            double photonGammaID_{-1}; 

	    bool verboseRun{true}; //Lene: print loads and loads of output from every step of the filter. 
	    //should totally be false, normally, and also not hardwired like this.


    }; // HighPassFilter 
}

#endif // BIASING_HIGHPASSFILTER_H__
