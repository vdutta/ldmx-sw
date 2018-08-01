/**
 * @file HcalDetectorGeometry.h
 * @author Tom Eichlersmith, University of Minnesota
 * @brief Header file for class HcalDetectorGeometry
 */

#ifndef HCAL_HCALDETECTORGEOMETRY_H
#define HCAL_HCALDETECTORGEOMETRY_H

//STL
#include <map> //storage maps

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
             * @param point vector containing real space point
             * @param errs vector containing errors in each coordinate
             */
            void transformDet2Real( const HcalHit* hit , std::vector< float > &point ,
                std::vector< float > &errs ) const;

        private:
            /** Number of layers in each section */
            std::map< HcalID::HcalSection , int > nLayers_;

            /** Number of strips per layer in each section */
            std::map< HcalID::HcalSection , int > nStrips_;

            /** Length of Scintillator Strip [mm] */
            std::map< HcalID::HcalSection , float > lengthScint_;

            /** Front of each section, the plane of the zero'th layer of each section [mm] */
            std::map< HcalID::HcalSection , float > frontSection_;
            
            /** Floor of each section, the plane of the zero'th strip of each section [mm] */
            std::map< HcalID::HcalSection , float > floorSection_;
 
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
