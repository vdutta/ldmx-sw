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
}
