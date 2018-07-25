/**
 * @file MIPHit.h
 * @author Tom Eichlersmith, University of Minnesota
 * @brief Header file for class MIPHit
 */

#ifndef HCAL_MIPHIT_H
#define HCAL_MIPHIT_H

//Standard Libraries
#include <vector> //for storage fector of HcalHit pointers

//ROOT Libs

//LDMX Framework
#include "Event/HcalTrack.h" //for typedef of hcalhits
#include "DetDescr/HcalID.h" //for the HcalSection enum
#include "EventDisplay/DetectorGeometry.h" //for the HCAL Geometry constants

namespace ldmx {
    
    /**
     * @class MIPHit
     * @brief Stores pointers to HcalHits that are considered a single MIP hit (usually due to proximity).
     */
    class MIPHit {
        public:
            /**
             * Default Constructor
             */
            MIPHit();

            /**
             * Add an HcalHit to the MIPHit
             */
            void addHit( HitPtr hit );
            
            /**
             * Evaluate this MIPHit. Once this MIPHit has the HcalHits added to it,
             *  running this function will calculate the other member variables from the
             *  HcalHits stored in this class.
             * 
             * @return true if successful
             */
            bool Eval();

            /**
             * Get the section.
             */
            HcalSection getSection() const { return section_; }

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
             * Get the total energy of the MIPHit
             */
            float getEnergy() const { return totalEnergy_; }

            /**
             * Get the number of HcalHits in this MIPHit
             */
            int getNumHits() const { return hcalHits_.size(); }

            /**
             * Get MIPHit box center
             */
            const std::vector<float> &getBoxCenter() const { return boxCenter_; }

        private:
            
            /**
             * Determine the Section, Layer, Strip information and check if it is consistent.
             * All hits are in the same section and strips don't differ by more than one.
             */
            bool setSLS();

            /**
             * Determine the six planes of the MIPHit box in detector coordinates.
             *
             * @return true if successful
             */
            bool setBoxPlanes();

            /**
             * Calculate the center of the MIPHit box in detector coordinates.
             */
            void setBoxCenter();

            /**
             * Calculate the total energy of the MIPHit.
             */
            void setTotalEnergy();

            /** Storage vector of pointers to HcalHits */
            std::vector< HitPtr > hcalHits_;

            /** Section of the Hits */
            HcalSection section_;

            /** Layer of the Hits */
            int layer_;

            /** low and up strip of Hits (may be equal) */
            int lowstrip_, upstrip_;

            /** The total energy of the MIPHit */
            float totalEnergy_;

            /** The center of the MIPHit box in detector coordinates */
            std::vector< float > boxCenter_;

            /** The planes of the MIPHit box in detector coordinates */
            float xMin_, xMax_ , yMin_, yMax_, zMin_, zMax_;
    };
}

#endif /* HCAL_MIPHIT_H */
