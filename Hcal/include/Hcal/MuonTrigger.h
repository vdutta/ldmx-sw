/**
 * @file MuonTrigger.h
 * @brief Header file for MuonTrigger producer
 * @author Tom Eichlersmith, University of Minnesota
 */

#ifndef HCAL_MUONTRIGGER_H
#define HCAL_MUONTRIGGER_H

//LDMX Framework
#include "Event/Event.h"
#include "Framework/EventProcessor.h" //Needed to declare processor
#include "Framework/ParameterSet.h" // Needed to import parameters from configuration file
#include "Event/TriggerResult.h" //Store trigger result in event bus

namespace ldmx {
    
    /**
     * @class MuonTrigger
     * @brief Producer to trigger on  muons in order to separate them
     *  from interactions originating from the beam.
     *
     * Currently, this trigger is cutting on the number of consecutive layers and strips hit in any section.
     * The user has control over the minimums for each section separately, so that the different sizes of the
     * sections are represented.
     * The time it takes for this trigger to run scales roughly as three times the number of hits.
     * One loop through the hits is to sort them by section, layer, strip information.
     * Another checks for consecutive strips (split according to section).
     * And finally the last loop checks for consecutive layers (split according to section).
     * This is an upper bound estimate because if layers and strips are not counted more than once when
     * sorting the hits, so any overlap decreases the number of loops in the analysis.
     */
    class MuonTrigger : public ldmx::Producer {
        public:

            MuonTrigger(const std::string& name, ldmx::Process& process) : ldmx::Producer(name, process) {}

            virtual void configure(const ldmx::ParameterSet& ps);

            virtual void produce(ldmx::Event& event);

            virtual void onFileOpen() { }

            virtual void onFileClose() { }

            virtual void onProcessStart() { } 

            virtual void onProcessEnd() { }

        private:
            
            /**
             * Counts the number of consecutive numbers in the input STL set.
             */
            int consecCount( const std::set &numbers ) const;

            /** Collection name of HcalHits to look in */
            std::string hcalHitCollName_;

            /** Pass name that created the HcalHits collection */
            std::string hcalHitPassName_;

            /** Name of this trigger object */
            std::string triggerObjectName_;

            /** Array of threshholds for consecutive layers to hit (one for each section) */
            int minConsecLayersHit_[5];

            /** Array of threshholds for consecutive strips to hit (one for each section) */
            int minConsecStripsHit_[5];

            /** Result to be put in event bus */
            TriggerResult result_;
    };
}

#endif /* HCAL_MUONTRIGGER_H */
