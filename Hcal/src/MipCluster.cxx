/**
 * @file MipCluster.cxx
 * @brief Implementation file for class MipCluster
 * @author Tom Eichlersmith, University of Minnesota
 */

#include "Hcal/MipCluster.h"

namespace ldmx {

    MipCluster::MipHit() { }

    void MipCluster::addHit( ldmx::HcalHit* hit ) {
        hcalHits_.push_back( hit );
        return;
    }

    void MipCluster::setTotalEnergy() {
        totalEnergy_ = 0.0;
        for ( unsigned int i = 0; i < hcalHits_.size(); i++ ) {
            totalEnergy_ += hcalHits_[i]->getEnergy();
        } //iterate through list of hits (i)
        return;
    }
}
