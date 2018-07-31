/**
 * @file HcalMipTriggerResult.cxx
 * @brief Implementation file for HcalMipTriggerResult class.
 * @author Tom Eichlersmith, University of Minnesota
 */

namespace ldmx {
   
    void HcalMipTriggerResult::addTrack( const std::vector< DetectorID::RawValue > &track ) {
        
        trackVec_.push_back( track );
        
        return;
    }
}
