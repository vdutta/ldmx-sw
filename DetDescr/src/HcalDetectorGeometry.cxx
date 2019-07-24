/**
 * @file HcalDetectorGeometry.cxx
 * @brief Implementation file for class HcalDetectorGeometry
 */

#include "Tools/HcalDetectorGeometry.h"

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

        parityVertical_ = 1;

        uncertaintyTimingPos_ = 200.0;

        thicknessScint_ = 20.0; 

        widthScint_ = 100.0;

        thicknessLayer_ = 50. + thicknessScint_ + 2*2.; //absorber + scint +2*air
    }

    HitBox HcalDetectorGeometry::transformDet2Real( HcalHit* hit ) const {
        
        //pairs that will go into HitBox
        std::pair<double,double> X,Y,Z;

        HcalSection section = (HcalSection)( hit->getSection() );
        int layer = hit->getLayer();
        int strip = hit->getStrip();

        //calculate center of layer,strip with respect to detector section
        double layercenter = layer*thicknessLayer_ + 0.5*thicknessScint_;
        double stripcenter = (strip + 0.5)*widthScint_;

        //calculate error in layer,strip position
        double elayer = 0.5*thicknessScint_;
        double estrip = 0.5*widthScint_;
        
        double x,y,z;
        if ( section == HcalSection::BACK ) {
            
            z = zeroLayer_.at( section ) + layercenter;
            Z.first  = z-elayer;
            Z.second = z+elayer;
            
            //only horizontal layers implemented currently
            if ( false ) { //( (layer ^ parityVertical_) & 1) == 0 ) { //checks for same parity
                //Vertical Layers
                
                x = zeroStrip_.at( section ) + stripcenter;
                X.first  = x - estrip;
                X.second = x + estrip;
                
                y = hit->getY();
                Y.first  = y - uncertaintyTimingPos_;
                Y.second = y + uncertaintyTimingPos_;

            } else {
                //Horizontal Layers
                
                x = hit->getX();
                X.first  = x - uncertaintyTimingPos_;
                X.second = x + uncertaintyTimingPos_;

                y = zeroStrip_.at( section ) + stripcenter;
                Y.first  = y - estrip;
                Y.second = y + estrip;

            } //calculate depending on layer

        } else {
            
            z = zeroStrip_.at( section ) + stripcenter;
            Z.first  = z - estrip;
            Z.second = z + estrip;

            if ( section == HcalSection::TOP or section == HcalSection::BOTTOM ) {
                
                x = hit->getX();
                X.first  = x - uncertaintyTimingPos_;
                X.second = x + uncertaintyTimingPos_;
                
                if ( section == HcalSection::TOP ) {
                    y = zeroLayer_.at( section ) + layercenter;
                } else {
                    y = zeroLayer_.at( section ) - layercenter;
                } //top or bottom hcal

                Y.first  = y - elayer;
                Y.second = y + elayer;
                
            } else if ( section == HcalSection::LEFT or section == HcalSection::RIGHT ) {
                
                y = hit->getY();
                Y.first  = y - uncertaintyTimingPos_;
                Y.second = y + uncertaintyTimingPos_;

                if ( section == HcalSection::LEFT ) {
                    x = zeroLayer_.at( section ) + layercenter;
                } else {
                    x = zeroLayer_.at( section ) - layercenter;
                } //left or right hcal

                X.first  = x - elayer;
                X.second = x + elayer;
    
            } else {
                std::cerr << "[ HcalDetectorGeometry::transformDet2Real ] : Unknown Hcal Section!" << std::endl;
                std::cerr << "    Returning a valid HitBox but with values that are all zero." << std::endl;
                std::pair<double,double> zeroes(0.0,0.0);
                HitBox hbox( 3 , zeroes );
                return hbox;
            } //side hcal
        
        } //calculate depending on section

        HitBox hbox;
        hbox.push_back( X );
        hbox.push_back( Y );
        hbox.push_back( Z );
        return hbox;
    }
    
    HitBox HcalDetectorGeometry::transformDet2Real( const std::vector<HcalHit*>  &hitVec ) const {
        
        std::vector<double> pointSum ( 3 , 0.0 ); //sums of weighted coordinates
        std::vector<double> weightSum( 3 , 0.0 ); //sums of weights for each coordinate
        
        //calculate real space point for each hit
        for ( HcalHit* hit : hitVec ) {
            
            HitBox box = transformDet2Real( hit );
            
            //Add weighted values to sums
            double weight;
            for ( unsigned int iC = 0; iC < 3; iC++ ) {
                
                double cer = abs(box[iC].second - box[iC].first)/2.0;

                weight = 1.0 / ( cer*cer );
                weightSum[ iC ] += weight;
                pointSum[ iC ] += weight*( ( box[iC].second + box[iC].first )/2.0 );
            }
        } //go through hitVec
        
        //Construct final HitBox
        HitBox hbox;
        for ( int iC = 0; iC < 3; iC++ ) {
            double c = pointSum[ iC ] / weightSum[ iC ];
            double ec = 1.0 / sqrt( weightSum[ iC ] );
            hbox.emplace_back( c - ec , c + ec );
        }

        return hbox;
    }
}
