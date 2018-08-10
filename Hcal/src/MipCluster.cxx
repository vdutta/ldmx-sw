/**
 * @file MipCluster.cxx
 * @brief Implementation file for class MipCluster
 * @author Tom Eichlersmith, University of Minnesota
 */

#include "Hcal/MipCluster.h"

namespace ldmx {

    MipCluster::MipCluster() { }

    void MipCluster::addHit( HcalHit* hit ) {
        hcalHits_.push_back( hit );
        return;
    }

    void MipCluster::mergeCluster( const MipCluster &cluster ) {
        
        for ( unsigned int iH = 0; iH < cluster.getNumHits(); iH++ ) {
            this->addHit( cluster.getHcalHit( iH ) );
        }

        return;
    }

    void MipCluster::set() {
        
        //HCAL_DETECTOR_GEOMETRY
        // is meant to contain all relationship between detector coordinates (section, layer, strip)
        // and real space coordinates (x,y,z). This allows for there to be only one place where updates
        // to the detector geometry need to be input. Investigate the HcalDetectorGeoemetry class to see
        // about changing how the real space point is calculated from HcalHits.
        box_ = HCAL_DETECTOR_GEOMETRY.transformDet2Real( hcalHits_ );

        setTotalEnergy();
        
        //cluster hasn't been checked as a seed yet
        wasBadSeed( false );

        return;
    }

    HitBox MipCluster::getBox() const {
        return box_;
    }

    void MipCluster::setTotalEnergy() {
        totalEnergy_ = 0.0;
        for ( unsigned int i = 0; i < hcalHits_.size(); i++ ) {
            totalEnergy_ += hcalHits_[i]->getEnergy();
        } //iterate through list of hits (i)
        return;
    }

}
