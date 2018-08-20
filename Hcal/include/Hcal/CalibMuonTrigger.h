/**
 * @file CalibMuonTrigger.h
 * @brief Header file for CalibMuonTrigger
 * @author Tom Eichlersmith, University of Minnesota
 */

#ifndef HCAL_CALIBMUONTRIGGER_H
#define HCAL_CALIBMUONTRIGGER_H

//STL
#include <cmath> //abs value
#include <set> //listing layers that have already been counted
#include <map> //logging hits
#include <vector> //creating track

//LDMX Framework
#include "Event/Event.h" //getting event bus
#include "Framework/EventProcessor.h" //Needed to declare processor
#include "Framework/ParameterSet.h" // Needed to import parameters from configuration file
#include "DetDescr/HcalID.h" //For HcalSection enum and HcalID creation
#include "Event/HcalHit.h" //getting hits
#include "Event/TriggerResult.h" //storing result in event bus

namespace ldmx {
    
    /**
     * @class CalibMuonTrigger
     * @brief Constructs Trigger Result depending on Section, Layer, Strip, PE, and 
     *  Energy information of hits in the Hcal.
     *
     * In the interest of speed (and at the expense of accuracy), this producer looks for
     * tracks in each orientation of the Hcal (i.e. layers in the Hcal that are orientated
     * the same way). This means if there is one track through the back Hcal, this producer
     * may count that as two tracks: one through the horizontal layers and one through the
     * vertical layers. Additionally, the sections in the side Hcal are isolated similarly.
     * No attempt at combining these tracks is made, if a track is found anywhere, the event
     * passes this trigger.
     */
    class CalibMuonTrigger : public ldmx::Producer {
        public:

            CalibMuonTrigger(const std::string& name, ldmx::Process& process) 
                : ldmx::Producer(name, process) {}

            virtual void configure(const ldmx::ParameterSet& ps);

            virtual void produce(ldmx::Event& event);

            virtual void onFileOpen() { }

            virtual void onFileClose() { }

            virtual void onProcessStart() { } 

            /**
             * Prints performance trackers
             */
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
            bool isPlausibleMip( ldmx::HcalHit* hit ) const;

            /**
             * Find end points that haven't been tested before.
             * 
             * @param orientation orientation to search
             * @return true if new end points are found
             */
            bool findEndPoints( int orientation );

            /** struct to help organize hitLog */
            struct HitLogNode {
                /** Layer of HcalHit */
                int layer;

                /** Strip of HcalHit */
                int strip;

                /** flag if hit has been used as an end point */
                bool isUsed;
            };
                   
            // INPUT PARAMETERS

            /** Name of HcalHit collection */
            std::string hitCollName_;

            /** Name of HcalHit pass */
            std::string hitPassName_;

            /** Origin of Muons (Target or Cosmic) */
            std::string muonOrigin_;
            
            /** Maximum Difference Between a hit and the center line of track [in N layers/strips] */
            double trackRadius_;

            /** Minimum Fraction of Hits included in track to accept track */
            double minFracHit_;

            /** Absolute minimum number of hits to attempt to construct a track */
            int absoluteMinHits_;

            /** Maximum Energy of a HcalHit to be considered a MIP */
            double maxEnergy_;

            /** Minimum PE of a HcalHit to be considered a real hit (not noise) */
            double minPE_;

            /** Name of this trigger object */
            std::string triggerObjectName_;

            /** Trigger object to add to event */
            TriggerResult result_;

            // HELPER MEMBER VARIABLES
         
            /** Hits sorted by their orientations and stored as their raw IDs */
            std::map< unsigned int , HitLogNode > hitLog_[6];

            /** Current Starting Point */
            std::map< unsigned int , HitLogNode >::iterator startPt_;

            /** Current Finish Point */
            std::map< unsigned int , HitLogNode >::iterator finishPt_;

            /** Set that contains list of orientations that should use the layer
             * as the independent variable */
            std::set< int > layerUsers_;
    
            // PERFORMANCE TRACKERS

            /** Number of Events Passed */
            int numPass_;
            
            /** Number events that have a certain number of tracks */
            std::map< int , int > numTracksPerEvent_;
    };
}

#endif /* HCAL_CALIBMUONTRIGGER_H */
