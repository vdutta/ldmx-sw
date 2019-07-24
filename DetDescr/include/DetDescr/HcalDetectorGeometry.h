/**
 * @file HcalDetectorGeometry.h
 * @author Tom Eichlersmith, University of Minnesota
 * @brief Header file for class HcalDetectorGeometry
 */

#ifndef TOOLS_HCALDETECTORGEOMETRY_H
#define TOOLS_HCALDETECTORGEOMETRY_H

//STL
#include <map> //storage maps
#include <vector> //BoundingBox
#include <utility> //BoundingBox
#include <cmath> //sqrt
#include <iostream> //cerr

//LDMX Framework
#include "DetDescr/HcalID.h" //HcalSection enum
#include "Event/HcalHit.h" //hit pointer

namespace ldmx {
    
    /**
     * @type BoundingBox
     * @brief Stores the minimum and maximum of each coordinate for a box.
     *
     * This has all of the information needed to define an axis-aligned rectangular prism.
     */
    typedef std::vector< std::pair< double, double > > BoundingBox;

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
             * @return BoundingBox in real space
             */
            BoundingBox transformDet2Real( HcalHit* hit ) const;
            
            /**
             * Calculate real space coordinates of a cluster of hits.
             *
             * Determines cluster's coordinates by a weighted mean of the individuals.
             * 
             * @param hitVec vector of HcalHits to find a "center" for
             * @return BoundingBox in real space
             */
            BoundingBox transformDet2Real( const std::vector< HcalHit* > &hitVec ) const;

            /**
             * Get the depth of an input section.
             *
             * Depth is the direction perpendicular to the plane of the layers.
             *
             * @param section HcalSection that you want the depth of
             * @return double the depth of section
             */
            double getDepth( HcalSection section ) const { return ( nLayers_.at( section ) * thicknessLayer_ ); }

            /**
             * Get the width of an input section.
             *
             * Width is defined as the direction in the plane of the layers,
             * but perpendicular to the direction of the scintillator strips/bars.
             *
             * @param section HcalSection that you want the width of
             * @return double the widht of section
             */
            double getWidth( HcalSection section ) const { return ( nStrips_.at( section ) * widthScint_ ); }

            /**
             * Get the length of an input section.
             *
             * Length is defined as the direction along the strips/bars.
             *
             * @param section HcalSection that you want the length of
             * @return double the length of section
             */
            double getLength( HcalSection section ) const { return lengthScint_.at(section); }

            /**
             * Get the plane of zero'th layer for the input section.
             *
             * @return double plane of zero'th layer
             */
            double getZeroLayer( HcalSection section ) const { return zeroLayer_.at( section ); }

            /**
             * Get bounding box for the input section.
             *
             * @param section HcalSection
             * @return BoundingBox that bounds section
             */
            BoundingBox getBoundingBox( HcalSection section ) const;
        
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

    /**
     * Namespace Wide Instance
     */
    const HcalDetectorGeometry HCAL_DETECTOR_GEOMETRY;
}

#endif /* TOOLS_HCALDETECTORGEOMETRY_H */
