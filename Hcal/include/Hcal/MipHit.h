/**
 * @file MipHit.h
 * @author Tom Eichlersmith, University of Minnesota
 * @brief Header file for class MipHit
 */

#ifndef HCAL_MIPHIT_H
#define HCAL_MIPHIT_H

//Standard Libraries
#include <vector> //for storage fector of HcalHit pointers

//ROOT Libs

//LDMX Framework
#include "Event/HcalTrack.h" //for typedef of hcalhits
#include "DetDescr/HcalID.h" //for the HcalSection enum

namespace ldmx {
    
    /**
     * @class MipHit
     * @brief Stores pointers to HcalHits that are considered a single MIP hit (usually due to proximity).
     */
    class MipHit {
        public:
            /**
             * Default Constructor
             */
            MipHit();

            /**
             * Preferred Constructor
             * Uses input hit list to initialize mip.
             */
            MipHit( std::vector< HitPtr > hcalHits );

            /**
             * Add an HcalHit to the MipHit
             */
            void addHit( HitPtr hit );
            
            /**
             * Set Up this MipHit. Once this MipHit has the HcalHits added to it,
             *  running this function will calculate the other member variables from the
             *  HcalHits stored in this class.
             * 
             * @return true if successful
             */
            bool setUp();

            /**
             * Get the section.
             */
            int getSection() const { return section_; }

            /**
             * Get the layer.
             */
            int getLayer() const { return layer_; }

            /**
             * Get the lower strip.
             */
            int getLowStrip() const { return lowstrip_; }

            /**
             * Get the upper strip.
             */
            int getUpStrip() const { return upstrip_; }
            
            /**
             * Get the total energy of the MipHit
             */
            float getEnergy() const { return totalEnergy_; }

            /**
             * Get the number of HcalHits in this MipHit
             */
            int getNumHits() const { return hcalHits_.size(); }

            /**
             * Get MipHit box center
             */
            const std::vector<float> &getBoxCenter() const { return boxCenter_; }

            /**
             * Get vector of HcalHits in this MipHit
             */
            const std::vector< HitPtr > &getHcalHits() const { return hcalHits_; }

        private:
            
            /**
             * Determine the Section, Layer, Strip information and check if it is consistent.
             * All hits are in the same section and strips don't differ by more than one.
             */
            bool setSLS();

            /**
             * Determine the six planes of the MipHit box in detector coordinates.
             *
             * @return true if successful
             */
            bool setBoxPlanes();

            /**
             * Calculate the center of the MipHit box in detector coordinates.
             */
            void setBoxCenter();

            /**
             * Calculate the total energy of the MipHit.
             */
            void setTotalEnergy();

            /** Section of the Hits */
            int section_;

            /** Layer of the Hits */
            int layer_;

            /** low and up strip of Hits (may be equal) */
            int lowstrip_, upstrip_;

            /** The total energy of the MipHit */
            float totalEnergy_;

            /** The center of the MipHit box in detector coordinates */
            std::vector< float > boxCenter_;

            /** The planes of the MipHit box in detector coordinates */
            float xMin_, xMax_ , yMin_, yMax_, zMin_, zMax_;
 
            /** Storage vector of pointers to HcalHits */
            std::vector< HitPtr > hcalHits_;
   };

    /** Pointer to MipHit instance */
    typedef MipHit * MipHitPtr;

}

#endif /* HCAL_MIPHIT_H */
