/**
 * @file MipHit.cxx
 * @brief Implementation file for class MipHit
 * @author Tom Eichlersmith, University of Minnesota
 */

#include "Hcal/MipHit.h"

namespace ldmx {

    MipHit::MipHit() : section_(0), layer_(0), lowstrip_(0), upstrip_(0), totalEnergy_(0.0),
        boxCenter_( 3 , -1.0 ), xMin_(0.0), xMax_(0.0), yMin_(0.0), yMax_(0.0), zMin_(0.0), zMax_(0.0) { }

    MipHit::MipHit( std::vector< HitPtr > hcalHits ) : MipHit() {
        for ( HcalHit* hit : hcalHits ) {
            addHit( hit );
        }
    }

    void MipHit::addHit( HitPtr hit ) {
        hcalHits_.push_back( hit );
        return;
    }

    bool MipHit::setUp() {
        
        setTotalEnergy();

        if ( !setBoxPlanes() ) {
            return false;
        }

        setBoxCenter();
        
        return true;
    }

    bool MipHit::setSLS() {
        
        layer_ = static_cast<int>(hcalHits_[0]->getLayer());
        section_ = hcalHits_[0]->getSection();
        lowstrip_ = hcalHits_[0]->getStrip();
        upstrip_ = hcalHits_[0]->getStrip();
        for ( auto it = hcalHits_.begin(); it != hcalHits_.end(); ++it ) {
            if ( static_cast<int>((*it)->getLayer()) != layer_ ) {
                std::cout << "[ MipHit::setBoxCorners ] : Different layers were grouped into the same MipHit." << std::endl;
                return false;
            } else if ( (*it)->getSection() != section_ ) {
                std::cout << "[ MipHit::setBoxCorners ] : Different sections were grouped into the same MipHit." << std::endl;
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
    
    bool MipHit::setBoxPlanes() {
        
        bool success = false;
        
        return success;
    }

    void MipHit::setBoxCenter() {
        boxCenter_[0] = (xMin_ + xMax_)/2;
        boxCenter_[1] = (yMin_ + yMax_)/2;
        boxCenter_[2] = (zMin_ + zMax_)/2;
        return;
    }

    void MipHit::setTotalEnergy() {
        totalEnergy_ = 0.0;
        for ( unsigned int i = 0; i < hcalHits_.size(); i++ ) {
            totalEnergy_ += hcalHits_[i]->getEnergy();
        } //iterate through list of hits (i)
        return;
    }
}
