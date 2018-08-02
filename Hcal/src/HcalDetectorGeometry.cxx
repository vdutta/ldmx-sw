/**
 * @file HcalDetectorGeometry.cxx
 * @brief Implementation file for class HcalDetectorGeometry
 */

#include "Hcal/HcalDetectorGeometry.h"

namespace ldmx {
    
    HcalDetectorGeometry::HcalDetectorGeometry() {
        
        nLayers_[ HcalSection::BACK   ] = 81;
        nLayers_[ HcalSection::TOP    ] = 17;
        nLayers_[ HcalSection::BOTTOM ] = 17;
        nLayers_[ HcalSection::LEFT   ] = 17;
        nLayers_[ HcalSection::RIGHT  ] = 17;
        
        nStrips_[ HcalSection::BACK   ] = 31;
        nStrips_[ HcalSection::TOP    ] = 31;
        nStrips_[ HcalSection::BOTTOM ] = 31;
        nStrips_[ HcalSection::LEFT   ] = 31;
        nStrips_[ HcalSection::RIGHT  ] = 31;
         
        lengthScint_[ HcalSection::BACK   ] = 3100.;
        lengthScint_[ HcalSection::TOP    ] = (3100.+525.)/2.;
        lengthScint_[ HcalSection::BOTTOM ] = (3100.+525.)/2.;
        lengthScint_[ HcalSection::LEFT   ] = (3100.+525.)/2.;
        lengthScint_[ HcalSection::RIGHT  ] = (3100.+525.)/2.;
         
        zeroLayer_[ HcalSection::BACK   ] = 200. + 290.;
        zeroLayer_[ HcalSection::TOP    ] = 525./2.;
        zeroLayer_[ HcalSection::BOTTOM ] = -525./2.;
        zeroLayer_[ HcalSection::LEFT   ] = 525./2.;
        zeroLayer_[ HcalSection::RIGHT  ] = -525./2.;
         
        zeroStrip_[ HcalSection::BACK   ] = -3100./2.;
        zeroStrip_[ HcalSection::TOP    ] = 200.;
        zeroStrip_[ HcalSection::BOTTOM ] = 200.;
        zeroStrip_[ HcalSection::LEFT   ] = 200.;
        zeroStrip_[ HcalSection::RIGHT  ] = 200.;

        parityVertical_ = 0;

        uncertaintyTimingPos_ = 200.0;

        thicknessScint_ = 20.0;

        widthScint_ = 100.0;

        thicknessLayer_ = 50. + 20. + 2*2.;
    }

    void HcalDetectorGeometry::transformDet2Real( const HcalHit* hit ,
        std::vector< float > &point , std::vector< float > &errs ) const {
        
        point.clear();
        errs.clear();
        for ( int i = 0; i < 3; i++ ) {
            point.push_back( 0.0 );
            errs.push_back( 0.0 );
        }    
        
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
            
            point[2] = zeroLayer_.at( section ) + layercenter;
            errs[2] = elayer;
            
            if ( ( (layer ^ parityVertical_) & 1) == 0 ) { //checks for same parity
                //Vertical Layers
                
                point[0] = zeroStrip_.at( section ) + stripcenter;
                errs[0] = estrip;

                point[1] = hit->getY();
                errs[1] = uncertaintyTimingPos_;

            } else {
                //Horizontal Layers
                
                point[0] = hit->getX();
                errs[0] = uncertaintyTimingPos_;

                point[1] = zeroStrip_.at( section ) + stripcenter;
                errs[1] = estrip;

            } //calculate depending on layer

        } else {
            
            point[2] = zeroStrip_.at( section ) + stripcenter;
            errs[2] = estrip;

            if ( section == HcalSection::TOP or section == HcalSection::BOTTOM ) {
                
                point[0] = hit->getX();
                errs[0] = uncertaintyTimingPos_;
                
                if ( section == HcalSection::TOP ) {
                    point[1] = zeroLayer_.at( section ) + layercenter;
                } else {
                    point[1] = zeroLayer_.at( section ) - layercenter;
                } //top or bottom hcal
                errs[1] = elayer;
                
            } else if ( section == HcalSection::LEFT or section == HcalSection::RIGHT ) {
                
                point[1] = hit->getY();
                errs[1] = uncertaintyTimingPos_;

                if ( section == HcalSection::LEFT ) {
                    point[0] = zeroLayer_.at( section ) + layercenter;
                } else {
                    point[0] = zeroLayer_.at( section ) - layercenter;
                } //left or right hcal
                errs[0] = elayer;
    
            } else {
                std::cerr << "[ HcalDetectorGeometry::transformDet2Real ] : Unknow Hcal Section!" << std::endl;
                return;
            } //side hcal
        
        } //calculate depending on section

        return;
    }
    
    void HcalDetectorGeometry::transformDet2Real( const std::vector<const HcalHit*>  &hitVec ,
        std::vector< float > &point , std::vector< float > &errs ) const {
        
        point.clear();
        errs.clear();
        for ( int i = 0; i < 3; i++ ) {
            point.push_back( 0.0 );
            errs.push_back( 0.0 );
        }
        
        std::vector<float> pointSum( 3 , 0.0 ); //sums of weighted coordinates
        std::vector<float> weightSum( 3 , 0.0 ); //sums of weights for each coordinate
        
        //calculate real space point for each hit
        for ( const HcalHit* hit : hitVec ) {
            
            std::vector<float> cpt, cer;

            transformDet2Real( hit , cpt , cer );
            
            //Add weighted values to sums
            float weight;
            for ( unsigned int iC = 0; iC < 3; iC++ ) {
                weight = 1.0 / ( cer[iC]*cer[iC] );
                weightSum[ iC ] += weight;
                pointSum[ iC ] += weight*cpt[iC];
            }
        } //go through hitVec
        
        //calculate final weighted mean
        for ( unsigned int iC = 0; iC < 3; iC++ ) {
            point[iC] = pointSum[iC]/weightSum[iC];
            errs[iC] = sqrt( 1.0 / weightSum[iC] );
        }

        return;
    }
}
