/**
 * @file HcalMipTriggerResult.h
 * @brief Header file for HcalMipTriggerResult class.
 * @author Tom Eichlersmith, University of Minnesota
 */

#ifndef EVENT_HCALMIPTRIGGERRESULT_H
#define EVENT_HCALMIPTRIGGERRESULT_H

//STL
#include <vector> //storage of vectors

//LDMX Framework
#include "DetDescr/HcalID.h" //RawValue typedef
#include "Event/TriggerResult.h" //Base class

namespace ldmx {
    
    /**
     * @class HcalMipTriggerResult
     * @brief Storage class for result from HcalMipTriggerProducer.
     */
    class HcalMipTriggerResult : public TriggerResult {
        public:
            /**
             * Default Constructor
             */
            HcalMipTriggerResult() : TriggerResult() { }

            /**
             * Destructor
             */
            ~HcalMipTriggerResult() {
                Clear();
            }
            
            /**
             * Clear the HcalMipTriggerResult
             */
            void Clear( Option_t *opt = "" );

            /**
             * Copy the HcalMipTriggerResult
             */
            void Copy( Option_t *opt = "" );
            /**
             * Add a track to the vector of tracks.
             */
            void addTrack( const std::vector< DetectorID::RawValue > &track );

            /**
             * Get number of tracks in this result.
             */
            int getNumTracks() const { return trackVec_.size(); }

            /**
             * Get the whole vector of tracks.
             */
            std::vector< std::vector< DetectorID::RawValue > > &getTrackVec() const { return trackVec_; }
            
            /**
             * Set the fraction of layers hit threshold.
             */
            void setFracLayersHit( const float fracLayersHit ) { fracLayersHit_ = fracLayersHit; return; }

            /**
             * Get the fraction of layers hit threshold.
             */
            float getFracLayersHit() const { return fracLayersHit_; }

            /**
             * Set the radius of the track cylinder.
             */
            void setTrackRadius( const int trackRadius ) { trackRadius_ = trackRadius; return; }

            /**
             * Get the radius of the track cylinder.
             */
            int getTrackRadius() const { return trackRadius_; }

        private:

            /** The list of tracks found (may be empty) */
            std::vector< std::vector< DetectorID::RawValue > > trackVec_;

    };
}

#endif /* EVENT_HCALMIPTRIGGERRESULT_H */
