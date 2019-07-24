/**
 * @file HcalDetectorGeometry.cxx
 * @brief Implementation file for class HcalDetectorGeometry
 */

#include "DetDescr/HcalDetectorGeometry.h"

namespace ldmx {
    
    HcalDetectorGeometry::HcalDetectorGeometry() {
        
        parityVertical_ = 1;
        thicknessScint_ = 15.0; 

        widthScint_ = 100.0;

        thicknessLayer_ = 25. + thicknessScint_ + 2*2.; //absorber + scint +2*air

        nLayers_[ HcalSection::BACK   ] = 100;
        nLayers_[ HcalSection::TOP    ] = 28;
        nLayers_[ HcalSection::BOTTOM ] = 28;
        nLayers_[ HcalSection::LEFT   ] = 28;
        nLayers_[ HcalSection::RIGHT  ] = 28;
        
        nStrips_[ HcalSection::BACK   ] = 30;
        nStrips_[ HcalSection::TOP    ] = 3;
        nStrips_[ HcalSection::BOTTOM ] = 3;
        nStrips_[ HcalSection::LEFT   ] = 3;
        nStrips_[ HcalSection::RIGHT  ] = 3;
         
        double ecal_z  = 290.;
        double ecal_xy = 525.;
        double ecal_front = 200.;

        lengthScint_[ HcalSection::BACK   ] = 3000.;
        lengthScint_[ HcalSection::TOP    ] = (3000.+ecal_xy)/2.;
        lengthScint_[ HcalSection::BOTTOM ] = (3000.+ecal_xy)/2.;
        lengthScint_[ HcalSection::LEFT   ] = (3000.+ecal_xy)/2.;
        lengthScint_[ HcalSection::RIGHT  ] = (3000.+ecal_xy)/2.;
         
        zeroLayer_[ HcalSection::BACK   ] = ecal_front + nStrips_[ HcalSection::TOP ] * widthScint_;
        zeroLayer_[ HcalSection::TOP    ] = ecal_xy/2.;
        zeroLayer_[ HcalSection::BOTTOM ] = ecal_xy/2.;
        zeroLayer_[ HcalSection::LEFT   ] = ecal_xy/2.;
        zeroLayer_[ HcalSection::RIGHT  ] = ecal_xy/2.;
         
        zeroStrip_[ HcalSection::BACK   ] = 3000./2.; 
        zeroStrip_[ HcalSection::TOP    ] = 200.;
        zeroStrip_[ HcalSection::BOTTOM ] = 200.;
        zeroStrip_[ HcalSection::LEFT   ] = 200.;
        zeroStrip_[ HcalSection::RIGHT  ] = 200.;
    }

    BoundingBox HcalDetectorGeometry::getBoundingBox( HcalHit* hit ) const {
        
        //pairs that will go into BoundingBox
        std::pair<double,double> X(0,0), Y(0,0), Z(0,0);

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
                std::cerr << "[ HcalDetectorGeometry::getBoundingBox ] : Unknown Hcal Section!" << std::endl;
                std::cerr << "    Returning a valid BoundingBox but with values that are all zero." << std::endl;
            } //side hcal
        
        } //calculate depending on section

        BoundingBox hbox;
        hbox.push_back( X );
        hbox.push_back( Y );
        hbox.push_back( Z );
        return hbox;
    }
    
    BoundingBox HcalDetectorGeometry::getBoundingBox( const std::vector<HcalHit*>  &hitVec ) const {
        
        std::vector<double> pointSum ( 3 , 0.0 ); //sums of weighted coordinates
        std::vector<double> weightSum( 3 , 0.0 ); //sums of weights for each coordinate
        
        //calculate real space point for each hit
        for ( HcalHit* hit : hitVec ) {
            
            BoundingBox box = getBoundingBox( hit );
            
            //Add weighted values to sums
            double weight;
            for ( unsigned int iC = 0; iC < 3; iC++ ) {
                
                double cer = abs(box[iC].second - box[iC].first)/2.0;

                weight = 1.0 / ( cer*cer );
                weightSum[ iC ] += weight;
                pointSum[ iC ] += weight*( ( box[iC].second + box[iC].first )/2.0 );
            }
        } //go through hitVec
        
        //Construct final BoundingBox
        BoundingBox hbox;
        for ( int iC = 0; iC < 3; iC++ ) {
            double c = pointSum[ iC ] / weightSum[ iC ];
            double ec = 1.0 / sqrt( weightSum[ iC ] );
            hbox.emplace_back( c - ec , c + ec );
        }

        return hbox;
    }

    BoundingBox HcalDetectorGeometry::getBoundingBox( HcalSection section ) const {

        std::pair< double, double > X(0,0), Y(0,0), Z(0,0);

        double total_strip_width = nStrips_.at( section ) * widthScint_;
        double total_thickness = nLayers_.at( section ) * thicknessLayer_;
        if ( section == HcalSection::BACK ) {
           
            X.first  = -zeroStrip_.at( HcalSection::BACK );
            X.second = X.first + total_strip_width;

            Y.first  = -lengthScint_.at( HcalSection::BACK )/2.0;
            Y.second =  lengthScint_.at( HcalSection::BACK )/2.0;

            Z.first  = zeroLayer_.at( HcalSection::BACK );
            Z.second = Z.first + total_thickness;

        } else {

            Z.first  = zeroStrip_.at( section );
            Z.second = Z.first + total_strip_width;

            if ( section == HcalSection::LEFT ) {
                
                X.first  = zeroLayer_.at( HcalSection::LEFT );
                X.second = X.first + total_thickness;

                Y.second = zeroLayer_.at( HcalSection::TOP );
                Y.first  = Y.second - lengthScint_.at( HcalSection::LEFT );

            } else if ( section == HcalSection::RIGHT ) {

                X.second = -zeroLayer_.at( HcalSection::RIGHT );
                X.first  = X.second - total_thickness;

                Y.first  = -zeroLayer_.at( HcalSection::BOTTOM );
                Y.second = Y.first + lengthScint_.at( HcalSection::RIGHT );

            } else if ( section == HcalSection::TOP ) {

                Y.first  = zeroLayer_.at( HcalSection::TOP );
                Y.second = Y.first + total_thickness;

                X.first  = -zeroLayer_.at( HcalSection::RIGHT );
                X.second = X.first + lengthScint_.at( HcalSection::TOP );

            } else if ( section == HcalSection::BOTTOM ) {

                Y.second = -zeroLayer_.at( HcalSection::BOTTOM );
                Y.first  = Y.second - total_thickness;

                X.second = zeroLayer_.at( HcalSection::LEFT );
                X.first  = X.second - lengthScint_.at( HcalSection::BOTTOM );

            } else {
                std::cerr << "[ Warning ] : Unrecognized HcalSection in HcalDetectorGeometry::getBoundingBox." << std::endl;
                std::cerr << "    Will return an incorrect geometry description!" << std::endl;
            }
        }

        BoundingBox boundingbox;
        boundingbox.push_back( X );
        boundingbox.push_back( Y );
        boundingbox.push_back( Z );

        return boundingbox;
    }
}
