/**
 * @file HcalMipTriggerProducer.h
 * @brief Header file for HcalMipTriggerProducer
 * @author Tom Eichlersmith, University of Minnesota
 */

#ifndef HCAL_HCALMIPTRIGGERPRODUCER_H
#define HCAL_HCALMIPTRIGGERPRODUCER_H

//STL
#include <set> //listing layers that have already been counted
#include <map> //logging hits
#include <vector> //creating track

//LDMX Framework
#include "Event/Event.h" //getting event bus
#include "Framework/EventProcessor.h" //Needed to declare processor
#include "Framework/ParameterSet.h" // Needed to import parameters from configuration file
#include "DetDescr/HcalID.h" //For HcalSection enum and HcalID creation
#include "Event/HcalHit.h" //getting hits

namespace ldmx {
    
    /**
     * @class HcalMipTriggerProducer
     * @brief Constructs Trigger Result depending on Section, Layer, Strip, and 
     *  Amplitude information of hits in the Hcal.
     */
    class HcalMipTriggerProducer : public ldmx::Producer {
        public:

            HcalMipTriggerProducer(const std::string& name, ldmx::Process& process) 
                : ldmx::Producer(name, process) {}

            virtual void configure(const ldmx::ParameterSet& ps);

            virtual void produce(ldmx::Event& event);

            virtual void onFileOpen() { }

            virtual void onFileClose() { }

            virtual void onProcessStart() { } 

            virtual void onProcessEnd() { }

        private:
            
            /**
             * Determine if an HcalHit could be a MIP passing through.
             * Checks on mininimum PE and maximum Energy.
             *
             * This function could be changed to check a two bit energy
             * flag instead of checking against actual PE and Energy
             * measurements.
             *
             * @param hit HcalHit* that needs to be checked
             * @bool true if plausible MIP
             */
            bool isPlausibleMip( ldmx::HcalHit* hit ) const;

            /**
             * Find end points that haven't been tested before.
             * 
             * @param orientation HcalOrientation of orientation to search
             * @return true if new end points are found
             */
            bool findEndPoints( int orientation );

            /**
             * Remove a given raw ID value from the hit log.
             *
             * @param orientation HcalOrienation
             * @param rawID raw ID value to be removed
             */
            void removeHcalID( int orientation , DetectorID::RawValue rawID );

            /** struct to help organize hitLog */
            struct HitLogNode {
                /** ID of hit in log */
                HcalID id;

                /** flag if hit could be a good end point */
                bool isGood;
            };

            /** Hits sorted by their orientations and stored as their raw IDs */
            std::map< DetectorID::RawValue , HitLogNode > hitLog_[6];

            /** Current Starting Point */
            HcalID* startPt_;

            /** Current Finish Point */
            HcalID* finishPt_;
            
            /** array of layers in each orientation */
            int nLayersPerOrientation_[6]; 
            
            /** Name of HcalHit collection */
            std::string hitCollName_;

            /** Name of HcalHit pass */
            std::string hitPassName_;
            
            /** Maximum Difference Between a hit and the center line of track */
            float trackRadius_;

            /** Minimum Fraction of Layers Hit to be considered a MIP */
            float minFracLayersHit_;

            /** Maximum Energy of a HcalHit to be considered a MIP */
            float maxEnergy_;

            /** Minimum PE of a HcalHit to be considered a real hit (not noise) */
            float minPE_;

            /** Name of this trigger object */
            std::string triggerObjectName_;
    };
}

#endif /* HCAL_HCALMIPTRIGGERPRODUCER_H */
