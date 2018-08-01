/**
 * @file HcalDetectorGeometry.cxx
 * @brief Implementation file for class HcalDetectorGeometry
 */

#include "Hcal/HcalDetectorGeometry.h"

namespace ldmx {
    
    HcalDetectorGeometry::HcalDetectorGeometry() {
        
        nLayers_[ HcalSection::BACK   ] = 81
        nLayers_[ HcalSection::TOP    ] = 17
        nLayers_[ HcalSection::BOTTOM ] = 17
        nLayers_[ HcalSection::LEFT   ] = 17
        nLayers_[ HcalSection::RIGHT  ] = 17
        
        nStrips_[ HcalSection::BACK   ] = 31
        nStrips_[ HcalSection::TOP    ] = 31
        nStrips_[ HcalSection::BOTTOM ] = 31
        nStrips_[ HcalSection::LEFT   ] = 31
        nStrips_[ HcalSection::RIGHT  ] = 31
         
        lengthScint_[ HcalSection::BACK   ] = 3100.
        lengthScint_[ HcalSection::TOP    ] = (3100.+525.)/2.
        lengthScint_[ HcalSection::BOTTOM ] = (3100.+525.)/2.
        lengthScint_[ HcalSection::LEFT   ] = (3100.+525.)/2.
        lengthScint_[ HcalSection::RIGHT  ] = (3100.+525.)/2.
         
        zeroLayer_[ HcalSection::BACK   ] = 200. + 290.
        zeroLayer_[ HcalSection::TOP    ] = 525./2.
        zeroLayer_[ HcalSection::BOTTOM ] = -525./2.
        zeroLayer_[ HcalSection::LEFT   ] = 525./2.
        zeroLayer_[ HcalSection::RIGHT  ] = -525./2.
         
        zeroStrip_[ HcalSection::BACK   ] = -3100./2.
        zeroStrip_[ HcalSection::TOP    ] = 200.
        zeroStrip_[ HcalSection::BOTTOM ] = 200.
        zeroStrip_[ HcalSection::LEFT   ] = 200.
        zeroStrip_[ HcalSection::RIGHT  ] = 200.

        parityVertical_ = 0;

        uncertaintyTimingPos_ = 200.0;

        thicknessScint_ = 20.0;

        widthScint_ = 100.0;

        thicknessLayer_ = 50. + 20. + 2*2.;
    }
    void HcalDetectorGeometry::transformDet2Real( const HcalHit* hit , std::vector< float > &point ,
        std::vector< float &errs ) const {
        
        point.clear();
        errs.clear();
        for ( int i = 0; i < 3; i++ ) {
            point.push_back( 0.0 );
            errs.push_back( 0.0 );
        }

        float x,y,z; //point coordinates
        float ex,ey,ez; //errors in coordinates
        
        HcalSection section = (HcalSection)( hit->getSection() );
        int layer = hit->getLayer();
        int strip = hit->getStrip();

        //calculate center of layer,strip with respect to detector section
        float layercenter = (layer + 0.5)*thicknessLayer_;
        float stripcenter = (strip + 0.5)*thicknessScint_;

        //calculate error in layer,strip position
        float elayer = 0.5*thicknessLayer_;
        float estrip = 0.5*thicknessScint_;

        if ( section == HcalSection::BACK ) {
            
            z = zeroLayer_[ section ] + layercenter;
            ez = elayer
            
            if ( ( (layer ^ parityVertical_) & 1) == 0 ) { //checks for same parity
                //Vertical Layers
                
                x = zeroStrip_[ section ] + stripcenter;
                ex = estrip;

                y = hit->getY();
                ey = uncertaintyTimingPos_;

            } else {
                //Horizontal Layers
                
                x = hit->getX();
                ex = uncertaintyTimingPos_;

                y = zeroStrip_[ section ] + stripcenter;
                ey = estrip;

            } //calculate depending on layer

        } else {
            
            z = zeroStrip_[ section ] + stripcenter;
            ez = estrip;

            if ( section == HcalSection::TOP or section == HcalSection::BOTTOM ) {
                
                x = hit->getX();
                ex = uncertaintyTimingPos_;
                
                if ( section == HcalSection::TOP ) {
                    y = zeroLayer_[ section ] + layercenter;
                } else {
                    y = zeroLayer_[ section ] - layercenter;
                } //top or bottom hcal
                ey = elayer;
                
            } else if ( section == HcalSection::LEFT or section == HcalSection::RIGHT ) {
                
                y = hit->getY();
                ey = uncertaintyTimingPos_;

                if ( section == HcalSection::LEFT ) {
                    x = zeroLayer_[ section ] + layercenter;
                } else {
                    x = zeroLayer_[ section ] - layercenter;
                } //left or right hcal
                ex = elayer;
    
            } else {
                std::cerr << "[ HcalDetectorGeometry::transformDet2Real ] : Unknow Hcal Section!" << std::endl;
                return;
            } //side hcal
        
        } //calculate depending on section

        point[0] = x;
        point[1] = y;
        point[3] = z;

        errs[0] = ex;
        errs[1] = ey;
        errs[2] = ez;

        return;
    }
}
