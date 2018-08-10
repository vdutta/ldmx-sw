/**
 * @file HcalDetectorGeometry.h
 * @author Tom Eichlersmith, University of Minnesota
 * @brief Header file for class HcalDetectorGeometry
 */

#ifndef TOOLS_HCALDETECTORGEOMETRY_H
#define TOOLS_HCALDETECTORGEOMETRY_H

//STL
#include <map> //storage maps
#include <cmath> //sqrt
#include <iostream> //cerr

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
            HcalDetectorGeometry();

            /**
             * Calculate real space coordinates from detector location.
             *
             * @param hit HcalHit to find real space hit for
             * @param point vector that will contain real space point
             * @param errs vector that will contain errors in each coordinate
             */
            void transformDet2Real( HcalHit* hit ,
                std::vector< double > &point , std::vector< double > &errs ) const;
            
            /**
             * Calculate real space coordinates of a cluster of hits.
             *
             * Determines cluster's coordinates by a weighted mean of the individuals.
             * 
             * @param hitVec vector of HcalHits to find a "center" for
             * @param point vector that will contain real space point
             * @param errs vector that will contain errors in each coordinate
             */
            void transformDet2Real( const std::vector< HcalHit* > &hitVec ,
                std::vector< double > &point , std::vector< double > &errs ) const;
        
        private:
            /** Number of layers in each section */
            std::map< HcalSection , int > nLayers_;

            /** Number of strips per layer in each section */
            std::map< HcalSection , int > nStrips_;

            /** Length of Scintillator Strip [mm] */
            std::map< HcalSection , double > lengthScint_;

            /** The plane of the zero'th layer of each section [mm] */
            std::map< HcalSection , double > zeroLayer_;
            
            /** The plane of the zero'th strip of each section [mm] */
            std::map< HcalSection , double > zeroStrip_;
 
            /** an example layer number of a vertical layer */
            int parityVertical_;

            /** Uncertainty in timing position along a bar/strip [mm] */
            double uncertaintyTimingPos_;

            /** Thickness of Scintillator Strip [mm] */
            double thicknessScint_;

            /** Width of Scintillator Strip [mm] */
            double widthScint_;
 
            /** Thickness of a whole layer  [mm] */
            double thicknessLayer_;
           
    };
}

#endif /* TOOLS_HCALDETECTORGEOMETRY_H */
