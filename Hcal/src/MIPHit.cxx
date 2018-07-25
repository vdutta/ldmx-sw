/**
 * @file MIPHit.cxx
 * @brief Implementation file for class MIPHit
 * @author Tom Eichlersmith, University of Minnesota
 */

#include "Hcal/MIPHit.h"

namespace ldmx {

    MIPHit::MIPHIt() : hcalHits_(), section_(0), layer_(0), lowstrip_(0), upstrip_(0)_, totalEnergy_(0.0),
        boxCenter_( 3 , -1.0 ), xMin_(0.0), xMax_(0.0), yMin_(0.0), yMax_(0.0), zMin_(0.0), zMax_(0.0) 
        {  }

    void MIPHit::addHit( HitPtr hit ) {
        hcalHits_.push_back( hit );
        return;
    }

    bool MIPHit::Eval() {
        
        setTotalEnergy();

        if ( setBoxPlanes() ) {
            setBoxCenter();
        } else {
            return false;
        }
        
        return true;
    }

    bool MIPHit::setSLS() {
        
        layer_ = static_cast<int>(hcalHits_[0]->getLayer());
        section_ = hcalHits_[0]->getSection();
        lowstrip_ = hcalHits_[0]->getStrip();
        upstrip_ = hcalHits_[0]->getStrip();
        for ( auto it = hcalHits_.begin(); it != hcalHits_.end(); ++it ) {
            if ( static_cast<int>((*it)->getLayer()) != layer ) {
                std::cout << "[ MIPHit::setBoxCorners ] : Different layers were grouped into the same MIPHit." << std::endl;
                return false;
            } else if ( (*it)->getSection() != section ) {
                std::cout << "[ MIPHit::setBoxCorners ] : Different sections were grouped into the same MIPHit." << std::endl;
                return false;
            }
            
            int cstrip = (*it)->getStrip();
            if ( cstrip < lowstrip_ ) {
                lowstrip_ = cstrip;
            } else if ( cstrip > upstrip_ ) {
                upstrip_ = cstrip;
            } //calculate low,up strip
        } //check that all hits are in the same section,layer (it) and set low,up strip
        
        return true;
    }
    
    bool MIPHit::setBoxPlanes() {
        
        bool success = false;
        if ( !hcalHits_.empty() ) {
            
            //Determine the Section,Layer,Strip information
            if ( !setSLS() ) {
                return false;
            }

            //all hits are in the same section,layer
            // calculate planes depending on section,layer,strip information
            float frontlayer = (layer_-1)*hcal_layer_thick; //front of layer w.r.t. hcal section
            float backlayer = layer_*hcal_layer_thick; //back of layer w.r.t. hcal section
            float leftstrip = lowstrip_*bar_width; //left side of strip w.r.t. layer
            float rightstrip = (upstrip_+1)*bar_width; //right side of strip w.r.t. layer
            if ( section_ == ldmx::HcalSection::BACK ) {
                
                zMin_ = hcal_front_z + frontlayer;
                zMax_ = hcal_front_z + backlayer;

                if ( layer_ % 2 == 0 ) { //even (horizontal) layer
                    xMin_ = -hcal_x_width/2;
                    xMax_ = hcal_x_width/2;
                    yMin_ = -hcal_y_width/2 + leftstrip;
                    yMax_ = -hcal_y_width/2 + rightstrip;
                } else { //odd (vertical) layer
                    xMin_ = -hcal_x_width/2 + leftstrip;
                    xMax_ = -hcal_x_width/2 + rightstrip;
                    yMin_ = -hcal_y_width/2;
                    yMax_ = hcal_y_width/2;
                } //sort by layer

            } else { //side hcal
                
                zMin_ = ecal_front_z + leftstrip;
                zMax_ = ecal_front_z + rightstrip;

                if ( section_ == ldmx::HcalSection::TOP ) {
                    xMin_ = -hcal_ecal_xy/2;
                    xMax_ = hcal_x_width/2;
                    yMin_ = hcal_ecal_xy/2 + frontlayer;
                    yMax_ = hcal_ecal_xy/2 + backlayer;
                } else if ( section_ == ldmx::HcalSection::BOTTOM ) {
                    xMin_ = -hcal_x_width/2;
                    xMax_ = hcal_ecal_xy/2;
                    yMin_ = -hcal_ecal_xy/2 - backlayer;
                    yMax_ = -hcal_ecal_xy/2 - frontlayer;
                } else if ( section_ == ldmx::HcalSection::LEFT ) {
                    xMin_ = hcal_ecal_xy/2 + frontlayer;
                    xMax_ = hcal_ecal_xy/2 + backlayer;
                    yMin_ = -hcal_y_width/2;
                    yMax_ = hcal_ecal_xy/2;
                } else if ( section_ == ldmx::HcalSection::RIGHT ) {
                    xMin_ = -hcal_ecal_xy/2 - backlayer;
                    xMax_ = -hcal_ecal_xy/2 - frontlayer;
                    yMin_ = -hcal_ecal_xy/2;
                    yMax_ = hcal_y_width/2;
                } else {
                    std::cout << "[ MIPHit::setBoxPlanes ] : Unknown HcalSection!" << std::endl;
                    return false;
                }
            
            } //if-else tree to sort by section
            
            success = true;

        } //vector of hits is non empty
        
        return success;
    }

    void MIPHit::setBoxCenter() {
        boxCenter_[0] = (xMin_ + xMax_)/2;
        boxCenter_[1] = (yMin_ + yMax_)/2;
        boxCenter_[2] = (zMin_ + zMax_)/2;
        return;
    }

    void MIPHit::setTotalEnergy() {
        totalEnergy_ = 0.0;
        for ( unsigned int i = 0; i < hcalHits_.size(); i++ ) {
            totalEnergy_ += hcalHits_[i]->getEnergy();
        } //iterate through list of hits (i)
        return;
    }
}
