/**
 * @file EcalPNProcessFilter.h
 * @brief User action plugin that biases Geant4 to only process events which
 *        involve a photonuclear reaction in the ECal with certain final-state
 *        kinematics.
 * @author Joshua Hiltbrand 
 *         University of Minnesota 
 */

#ifndef BIASING_ECALPNPROCESSFILTER_H_
#define BIASING_ECALPNPROCESSFILTER_H_

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
#include "Biasing/EcalPNProcessFilterMessenger.h"
#include "Biasing/TargetBremFilter.h"
#include "SimCore/UserTrackInformation.h"

class UserTrackingAction;

namespace ldmx {
    
    // Forward declaration to avoid circular dependecies
    class EcalPNProcessFilterMessenger; 

    class EcalPNProcessFilter : public UserActionPlugin {

        public:

            /** Default Ctor */
            EcalPNProcessFilter();

            /** Destructor */
            ~EcalPNProcessFilter();

            /** @return A std::string descriptor of the class. */
            virtual std::string getName() {
                return "EcalPNProcessFilter";
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
             * Pre-tracking action.
             */
            void postTracking(const G4Track*); 

            /**
             * Classify a new track which postpones track processing.
             * Track processing resumes normally if a target process interaction occurred.
             * @param aTrack The Geant4 track.
             * @param currentTrackClass The current track classification.
             */
            G4ClassificationOfNewTrack stackingClassifyNewTrack(const G4Track* aTrack, const G4ClassificationOfNewTrack& currentTrackClass);

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

            /**
             * Add the fraction KE threshold for daughter neutrons and kaons.
             *
             * @param KE fraction threshold
             */
            void setEnergyFractionThreshold(double energyFractionThreshold);

        private:

            /** Messenger used to pass arguments to this class. */
            EcalPNProcessFilterMessenger* messenger_{nullptr};

            /** Pointer to the current track being processed. */
            G4Track* currentTrack_{nullptr};

            /** List of volumes to apply filter to. */
            std::vector<std::string> volumes_; 

            /** List of volumes to bound the particle to. */
            std::vector<std::string> boundVolumes_;

            /** KE fraction threshold for PN neutron and kaon daughters */
            double energyFractionThreshold_{0.5}; // MeV

            /** PNProcess gamma parent ID. */
            double photonGammaID_{-1}; 

    }; // EcalPNProcessFilter 
}

#endif // BIASING_ECALPNPROCESSFILTER_H
