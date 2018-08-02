/**
 * @file HcalDetectorGeometry.h
 * @author Tom Eichlersmith, University of Minnesota
 * @brief Header file for class HcalDetectorGeometry
 */

#ifndef HCAL_HCALDETECTORGEOMETRY_H
#define HCAL_HCALDETECTORGEOMETRY_H

//STL
#include <map> //storage maps
#include <cmath> //sqrt

//LDMX Framework
#include "DetDescr/HcalID.h" //HcalSection enum
#include "Event/HcalHit.h" //hit pointer

namespace ldmx {
    
    /**
     * @class HcalDetectorGeometry
     * @brief Class to translated between detector location (section, layer, strip) and real space.
     */
    class HcalDetectorGeometry {
        public:
            /**
             * Constructor
             * This is where all the detector constants are set.
             */
            HcalDetectorGeometry()

            /**
             * Calculate real space coordinates from detector location.
             *
             * @param hit HcalHit to find real space hit for
             * @param point vector that will contain real space point
             * @param errs vector that will contain errors in each coordinate
             */
            void transformDet2Real( const HcalHit* hit ,
                std::vector< float > &point , std::vector< float > &errs ) const;
            
            /**
             * Calculate real space coordinates of a cluster of hits.
             *
             * Determines cluster's coordinates by a weighted mean of the individuals.
             * 
             * @param hitVec vector of HcalHits to find a "center" for
             * @param point vector that will contain real space point
             * @param errs vector that will contain errors in each coordinate
             */
            void transformDet2Real( const std::vector< const HcalHit* > &hitVec ,
                std::vector< float > &point , std::vector< float > &errs ) const;
        
        private:
            /** Number of layers in each section */
            std::map< HcalID::HcalSection , int > nLayers_;

            /** Number of strips per layer in each section */
            std::map< HcalID::HcalSection , int > nStrips_;

            /** Length of Scintillator Strip [mm] */
            std::map< HcalID::HcalSection , float > lengthScint_;

            /** The plane of the zero'th layer of each section [mm] */
            std::map< HcalID::HcalSection , float > zeroLayer_;
            
            /** The plane of the zero'th strip of each section [mm] */
            std::map< HcalID::HcalSection , float > zeroStrip_;
 
            /** an example layer number of a vertical layer */
            int parityVertical_;

            /** Uncertainty in timing position along a bar/strip [mm] */
            float uncertaintyTimingPos_;

            /** Thickness of Scintillator Strip [mm] */
            float thicknessScint_;

            /** Width of Scintillator Strip [mm] */
            float widthScint_;
 
            /** Thickness of a whole layer  [mm] */
            float thicknessLayer_;
           
    };
}
