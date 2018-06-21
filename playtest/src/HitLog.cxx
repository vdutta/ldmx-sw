/**
 * @file HitLog.cxx
 * @brief Implementation file for class HitLog
 * @author Tom Eichlersmith, University of Minnesota
 */

#include "playtest/HitLog.h"

namespace ldmx {
    
    HitLog::HitLog() : minPE_(0.0), conedepth_(2), coneangle_(3), origin_(17.5),
                       lowside_(12.5), upside_(22.5), nstrips_(34), nlayers_(81), layermod_(1000) {
        
    }

    HitLog::HitLog( const float minPE , const int conedepth , const int coneangle , const float origin ,
            const float lowside , const float upside , const int nstrips , const int nlayers ) :
            HitLog(), minPE_(minPE), conedepth_(conedepth), coneangle_(coneangle), origin_(origin),
            lowside_(lowside), upside_(upside), nstrips_(nstrips), nlayers_(nlayers)  {

    }

    void HitLog::AddHit( HitPtr hit ) {
        
        int key = KeyGen( hit );
        log_[ key ] = hit;

        return;
    }

    int HitLog::KeyGen( HitPtr hit ) const {
        return static_cast<int>( hit->getLayer()*layermod_ + hit->getStrip() );
    }

    void HitLog::SetSearchCone( const int seedlayer , const int seedstrip ) {
        
    }

    std::pair< HitPtr , HitPtr > HitLog::search( const HitLog log , const int lowkey , const int upkey ) const {
        
        std::pair< HitPtr , HitPtr > ret( nullptr , nullptr );
        
        if ( lowkey > upkey ) { //Mis-use correction
            std::cout << "Input hit keys to search in wrong order: " << lowkey << " " << upkey << std::endl;
            std::cout << "Returning an empty search" << std::endl;
            return ret;
        } else { //inputs are correct form

            auto lowbound = log.lower_bound( lowkey ); //points to first key that is not before lowkey (equivalent or after) (map::end if all are before lowkey)
            auto upbound = log.upper_bound( upkey ); //points to first key after upkey (map::end if nothing after upkey)
        
            if ( lowbound != upbound ) { //check to see if range has any thickness
                //there is at least one hit in the range
                //lowbound points to the first one
                //see if there is a hit on either side of lowerbound
                auto beforeside = std::prev( lowbound ); 
                auto afterside = std::next( lowbound );
                
                int beforekeydif = lowbound->first - beforeside->first;
                int afterkeydif = afterside->first - lowbound->first;

                if ( beforekeydif != 1 or afterkeydif != 1 ) {
                    //lowbound has at most one neighbor
                    ret.first = lowbound->second;
                    
                    if ( beforekeydif == 1 ) {
                        //beforeside is the neighbor for lowbound
                        ret.second = beforeside->second;
                    } else if ( afterkeydif == 1 ) {
                        //afterside is the neighbor for lowerbound
                        ret.second = afterside->second;
                    } //else: lowbound is truly isolated

                } //check if lowbound could be isolated

            } //check to see if range has any thickness

        } //make sure inputs are correct
            
        return ret;
    }

    bool HitLog::stripbounds( const int seedlayer , const int seedstrip , const int layer , int &lowstrip , int &upstrip ) const {
        
        float slope = (seedstrip - origin_)/seedlayer; //slope of line between seed and origin

        lowstrip = static_cast<int>(floor( layer*slope + lowside_ ));
        upstrip = static_cast<int>(ceil( layer*slope + upside_ ));
        
        //Change if out of bounds
        if ( lowstrip < 1 )
            lowstrip = 1;
        if ( upstrip > nStrips_ )
            upstrip = nStrips_;

        if ( lowstrip > nStrips_ or upstrip < 1 ) {
            return false; //projected track is outside of Hcal
        }

        return true;
    }

}
