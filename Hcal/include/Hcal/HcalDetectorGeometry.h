/**
 * @file HcalDetectorGeometry.h
 * @author Tom Eichlersmith, University of Minnesota
 * @brief Header file for class HcalDetectorGeometry
 */

#ifndef HCAL_HCALDETECTORGEOMETRY_H
#define HCAL_HCALDETECTORGEOMETRY_H

//INCLUDES

namespace ldmx {
    
    /**
     * @class HcalDetectorGeometry
     * @brief Class to translated between detector location (section, layer, strip) and real space.
     */
    class HcalDetectorGeometry {
        public:
            /**
             * Constructor
             */
            HcalDetectorGeometry() { }

            /**
             * Calculate real space coordinates from detector location.
             */
            void transformDet2Real( /* Args */ );

        private:
            /** Number of layers in each section */
            static const int nLayersVerticalBackHcal_( 40 );
            static const int nLayersHorizontalBackHcal_( 41 );
            static const int nLayersTopHcal_( 17 );
            static const int nLayersBottomHcal_( 17 );
            static const int nLayersRightHcal_( 17 );
            static const int nLayersLeftHcal_( 17 );

            /** Number of strips per layer in each section */
            static const int nStripsBackHcal_( 31 );
            static const int nStripsTopHcal_( 31 );
            static const int nStripsBottomHcal_( 31 );
            static const int nStripsRightHcal_( 31 );
            static const int nStripsLeftHcal_( 31 );

            /** Flag to say which parity the vertical layers are */
            static const bool isVerticalOdd_( false );

            /** Uncertainty in timing position along a bar/strip [mm] */
            static const float uncertaintyTimingPos_( 200.0 );

            /** Thickness of Scintillator Strip [mm] */
            static const float thicknessScint_( 20.0 );

            /** Width of Scintillator Strip [mm] */
            static const float widthScint_( 100.0 );

            /** Length of Scintillator Strip [mm] */
            static const float lengthBackHcal_( 3100. );
            static const float lengthTopHcal_( (3100.+525.)/2 );
            static const float lengthBottomHcal_( (3100.+525.)/2 );
            static const float lengthLeftHcal_( (3100.+525.)/2 );
            static const float lengthRightHcal_( (3100.+525.)/2 );
 
            /** Thickness of a whole layer (scint+air+absorber) [mm] */
            static const float thicknessLayer_( 50. + 20. + 2*2. );

            /** Front of each section, the plane of the zero'th layer of each section [mm] */
            static const float frontBackHcal_( 200. + 290.0 ); //front in z
            static const float frontTopHcal_( 525./2 ); //bottom in y
            static const float frontBottomHcal_( -525./2 ); //top in y
            static const float frontLeftHcal_( 525./2 ); //rightside in x
            static const float frontRightHcal_( -525./2 ); //leftside in x
            
            /** Floor of each section, the plane of the zero'th strip of each section [mm] */
            static const float floorBackHcal_( -3100.0/2 ); //bottom in y
            static const float floorTopHcal_( 200. ); //front in z
            static const float floorBottomHcal_( 200. ); //front in z
            static const float floorLeftHcal_( 200. ); //front in z
            static const float floorRightHcal_( 200. ); //front in z
            
    };
}
