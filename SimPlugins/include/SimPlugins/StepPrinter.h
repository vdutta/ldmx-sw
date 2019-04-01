/**
 * @file StepPrinter.h
 * @class StepPrinter
 * @brief Sim plugin that prints the details of a step taken by a particle.
 * @author Omar Moreno, SLAC National Accelerator Laboratory
 */

#ifndef __SIM_PLUGINS_STEP_PRINTER_H__
#define __SIM_PLUGINS_STEP_PRINTER_H__

//------------//
//   Geant4   //
//------------//
#include "SimPlugins/UserActionPlugin.h"

class G4Step; 

namespace ldmx {

    class StepPrinter : public UserActionPlugin {

        public:

            /** @return The name of the plugin. */
            inline std::string getName() { return "StepPrinter"; }
            
            /**
             * @return True to indicate this plugin implements the stepping
             *         action.
             */
            inline bool hasSteppingAction() { return true; }
           
            /**
             * @return True to indicate this plugin implements the stacking 
             *         action.
             */
            inline bool hasStackingAction() { return true; }
            /**
             * Method used to classify a new track. 
             * @param aTrack The Geant4 track.
             * @param currentTrackClass The current track classification.
             */
            G4ClassificationOfNewTrack stackingClassifyNewTrack(const G4Track* aTrack, 
                    const G4ClassificationOfNewTrack& currentTrackClass);

            /**
             * Implement the stepping action that prints the step. 
             * @param step The Geant4 step.
             */
            void stepping(const G4Step* step);
        
        private:

    }; // StepPrinter

} // ldmx

#endif // __SIM_PLUGINS_STEP_PRINTER_H__
