/**
 * @file HcalMipTriggerProducer.h
 * @brief Header file for HcalMipTriggerProducer
 * @author Tom Eichlersmith, University of Minnesota
 */

#ifndef HCAL_HCALMIPTRIGGERPRODUCER_H
#define HCAL_HCALMIPTRIGGERPRODUCER_H

//LDMX Framework
#include "Event/Event.h"
#include "Framework/EventProcessor.h" //Needed to declare processor
#include "Framework/ParameterSet.h" // Needed to import parameters from configuration file
#include "DetDescr/HcalID.h" //For HcalSection enum and HcalID creation

namespace ldmx {
    
    /**
     * @class HcalMipTriggerProducer
     * @brief Constructs Trigger Result depending on Section, Layer, Strip, and Amplitude information of hits in the Hcal.
     */
    class HcalMipTriggerProducer : public ldmx::Producer {
        public:

            HcalMipTriggerProducer(const std::string& name, ldmx::Process& process) : ldmx::Producer(name, process) {}

            virtual void configure(const ldmx::ParameterSet& ps);

            virtual void produce(const ldmx::Event& event);

            virtual void onFileOpen();

            virtual void onFileClose();

            virtual void onProcessStart(); 

            virtual void onProcessEnd();

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
            bool isPlausibleMip( HcalHit* hit ) const;

            /**
             * Find end points that haven't been tested before.
             * 
             * @param orientation HcalOrientation of orientation to search
             * @return true if new end points are found
             */
            bool findEndPoints( HcalOrientation orientation );

            /**
             * Remove a given raw ID value from the hit log.
             *
             * @param orientation HcalOrienation
             * @param rawID raw ID value to be removed
             */
            void removeHcalID( HcalOrienation orientation , DetectorID::RawValue rawID );

            /** enum for layers orientated in the same way */
            enum HcalOrientation {
                BACK_EVEN = 0,
                BACK_ODD = 1,
                TOP = 2,
                BOTTOM = 3,
                LEFT = 4,
                RIGHT = 5
            };
            
            /** list of orientation enum */
            static const std::vector< HcalOrientation > HcalOrientationList_({ BACK_EVEN , BACK_ODD , TOP , BOTTOM , LEFT , RIGHT })
            
            /** Hits sorted by their orientations and stored as their raw IDs */
            std::map< HcalOrientation , std::map< DetectorID::RawValue , HcalID > > hitLog_;

            /** Bad End Points for each HcalOrientation */
            std::map< HcalOrientation , std::set< DetectorID::RawValue > > badEndPts_;

            /** Current Starting Point */
            HcalID* startPt_;

            /** Current Finish Point */
            HcalID* finishPt_;
            
            /** Name of HcalHit collection and pass */
            std::string hitCollName_;
            std::string hitPassName_;

            /** Minimum Number of Layers Hit to be considered a MIP */
            int minLayersHit_;

            /** Maximum Energy of a HcalHit to be considered a MIP */
            float maxEnergy_;

            /** Minimum PE of a HcalHit to be considered a real hit (not noise) */
            float minPE_;
    };
}

#endif /* HCAL_HCALMIPTRIGGERPRODUCER_H */
