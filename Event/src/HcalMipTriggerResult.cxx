/**
 * @file HcalMipTriggerResult.cxx
 * @brief Implementation file for HcalMipTriggerResult class.
 * @author Tom Eichlersmith, University of Minnesota
 */

#include "Event/HcalMipTriggerResult.h"

namespace ldmx {
    
    HcalMipTriggerResult::HcalMipTriggerResult() : TriggerResult() { } 

    HcalMipTriggerResult::~HcalMipTriggerResult() {
        Clear();
    }

    void HcalMipTriggerResult::Print(Option_t *option) const {
        
        TriggerResult::Print( option );

        return;
    }

    void HcalMipTriggerResult::Clear(Option_t *option) {
        
        TriggerResult::Clear( option );

        trackVec_.clear();

        return;
    }

    void HcalMipTriggerResult::Copy(TObject& o) const {
        
        TriggerResult::Copy( o );

        HcalMipTriggerResult& hmtr = (HcalMipTriggerResult&)( o );
        hmtr.trackVec_ = trackVec_;

        return;
    }

    void HcalMipTriggerResult::set( bool pass ) {
        TriggerResult::set( "HcalMipTrigger" , pass , 4 );
        return;
    }
    
    void HcalMipTriggerResult::addTrack( const std::vector< unsigned int > &track ) {
        trackVec_.push_back( track );
        return;
    }

    int HcalMipTriggerResult::getNumTracks() const {
        return trackVec_.size();
    }

    std::vector< std::vector< unsigned int > > HcalMipTriggerResult::getTrackVec() const {
        return trackVec_;
    }
    
    void HcalMipTriggerResult::setFracLayersHit( const float fracLayersHit ) {
        TriggerResult::setAlgoVar( 0 , static_cast<double>(fracLayersHit) );
        return;
    }

    float HcalMipTriggerResult::getFracLayersHit() const {
        return static_cast<float>( TriggerResult::getAlgoVar0() );
    }

    void HcalMipTriggerResult::setTrackRadius( const float trackRadius ) {
        TriggerResult::setAlgoVar( 1 , static_cast<double>( trackRadius ) );
        return;
    }

    float HcalMipTriggerResult::getTrackRadius() const {
        return static_cast<float>( TriggerResult::getAlgoVar1() );
    }

    void HcalMipTriggerResult::setMaxEnergy( const float maxEnergy ) {
        TriggerResult::setAlgoVar( 2 , static_cast<double>( maxEnergy ) );
        return;
    }

    float HcalMipTriggerResult::getMaxEnergy() const {
        return static_cast<float>( TriggerResult::getAlgoVar2() );
    }

    void HcalMipTriggerResult::setMinPE( const float maxEnergy ) {
        TriggerResult::setAlgoVar( 3 , static_cast<double>( minPE ) );
        return;
    }

    float HcalMipTriggerResult::getMinPE() const {
        return static_cast<float>( TriggerResult::getAlgoVar3() );
    }
}
