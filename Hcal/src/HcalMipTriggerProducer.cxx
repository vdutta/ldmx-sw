/**
 * @file HcalMipTriggerProducer.cxx
 * @brief Implementaiton of HcalMipTriggerProducer class
 * @author Tom Eichlersmith, University of Minnesota
 */

#include "Hcal/HcalMipTriggerProducer.h"

namespace ldmx {

    void HcalMipTriggerProducer::configure(const ldmx::ParameterSet& ps) {
        
        hitCollName_ = ps.getString( "HcalHitCollectionName" );
        hitPassName_ = ps.getString( "HcalHitPassName" );
        
        nBackLayers = ps.getInteger( "NumLayersBackHcal" );
        nLayersPerOrientation_[ HcalOrientation::BACK_EVEN ] = nBackLayers/2;
        nLayersPerOrientation_[ HcalOrientation::BACK_ODD ] = nBackLayers/2;
        if ( nBackLayers % 2 == 1 ) { //odd number back layers
            nLayersPerOrientation_[ HcalOrientation::BACK_ODD ] ++;
        }
        nLayersPerOrientation_[ HcalOrientation::TOP ] = ps.getInteger( "NumLayersTopHcal" );
        nLayersPerOrientation_[ HcalOrientation::BOTTOM ] = ps.getInteger( "NumLayersBottomHcal" );
        nLayersPerOrientation_[ HcalOrientation::LEFT ] = ps.getInteger( "NumLayersLeftHcal" );
        nLayersPerOrientation_[ HcalOrientation::RIGHT ] = ps.getInteger( "NumLayersRightHcal" );

        maxStripDif_ = ps.getDouble( "MaximumStripDifference" );

        minFracLayersHit_ = ps.getDouble( "MinimumFractionLayersHit" );

        maxEnergy_ = ps.getDouble( "MaximumEnergy" );

        minPE_ = ps.getDouble( "MinimumPE" );

        triggerObjectName_ = ps.getString( "HcalMipTriggerObjectName" );

        return;
    }

    void HcalMipTriggerProducer::produce(const ldmx::Event& event) {
        
        //initialize event containers
        hitLog_.clear();
        
        //obtain list of raw hits
        const TClonesArray *rawhits = event.getCollection( hitCollName_ , hitPassName_ );

        //add only plausible mip hits to hitLog_
        for ( size_t iH = 0; iH < rawhits->GetEntriesFast(); iH++ ) {
            HcalHit* chit = (HcalHit*)(rawhits->At(iH));
            if ( isPlausibleMip( chit ) ) { //could be a MIP
                
                HcalSection csect = chit->getSection();
                int clayer = chit->getLayer();
                HcalOrientation corient = csect;
                
                //The underlying ints for the HcalSection and HcalOrientation enums are shifted by one
                // except for EVEN layers in the BACK Hcal.
                if ( csect != HcalSection::BACK or (csect = HcalSection::BACK and clayer % 2 == 1 ) ) {
                    corient++;
                }
                
                //Each hit in hitLog_[ orientation ] has the same section number, so the standard sort
                // will sort them correctly according to layer,strip information
                HcalID cID;
                cID.setRawValue( chit->getID() );
                cID.unpack();
                hitLog_[ corient ][ chit->getID() ] = HitLogNode( cID , true );
            } //could be a MIP
        } //iterate through rawhits (iH)

        //for each HcalOrientation
        for ( HcalOrientation corient : HcalOrientationList_ ) {
            
            //while good end points are being found
            int trackcnt = 0;
            while ( findEndPoints( corient ) ) {
                
                int laycnt = 0;
                if ( startPt_->getLayerID() != finishPt_->getLayerID() ) {
                    
                    //calculate slope
                    int startLayer = startPt_->getLayerID();
                    int startStrip = startPt_->getStrip();
                    int dstrip = ( finishPt_->getStrip() - startStrip );
                    int dlayer = ( finishPt_->getLayerID() - startLayer );
                    float slope = static_cast<float>(dstrip)/static_cast<float>(dlayer);
                    
                    //count hits in this orientation that are in the cylinder
                    std::set<int> countedLayers; //layers that already have been counted
                    for ( auto it : hitLog_[ corient ] ) {
                        
                        HcalID* chit = &( (it->second).id );
    
                        int cstrip = chit->getStrip();
                        int clayer = chit->getLayerID();
    
                        float trackstrip = slope*( clayer - startLayer ) + startStrip;
    
                        float stripdif = std::abs( trackstrip - cstrip );
    
                        if ( stripdif < maxStripDif_ ) {
                            //hit is in track cylinder
                            if ( countedLayers.find( clayer ) != countedLayers.end() )
                                laycnt++;
    
                            track.push_back( chit->getRawValue() );
                        } //check if hit is in track cylinder
    
                    } //iterate through hitLog of current orientation (chit)
                
                } //startPt and finishPt are not on the same layer
                
                float layfrac = static_cast<float>( laycnt ) / static_cast<float>( nLayersPerOrientation_.at( corient ) );

                if ( layfrac > minFracLayersHit_ ) { 
                    //good track found
                    for ( DectorID::RawValue rawID : track ) {
                        removeHcalID( corient , rawID );
                    } //iterate through track (rawID)
                    
                    trackcnt++;

                } else {
                    //mark the end points as not good end points
                    hitLog_[ corient ][ startPt_->getRawValue() ].isGood = false;
                    hitLog_[ corient ][ finishPt_->getRawValue() ].isGood = false;
                } //what to do if track is good

            } //while good end points are still being found

        } //for each orientation

        return;
    }

    bool HcalMipTriggerProducer::isPlausibleMip( HcalHit* hit ) const {
        return ( hit->getPE() > minPE_ and hit->getEnergy() < maxEnergy_ );
    }

    bool HcalMipTriggerProducer::findEndPoints( HcalOrientation orientation ) {
        
        startPt_ = nullptr;
        finishPt_ = nullptr;

        if ( hitLog_[orientation].size() > 1 ) {
            
            for ( auto itH : hitLog_[ orientation ] ) {

                HcalID *chit = &( (itH->second).id );
                int clayer = chit->getLayerID();

                if ( startPt_ ) {
                    
                    if ( clayer < startPt_->getLayerID() and (itH->second).isGood ) {
                        startPt_ = chit;
                    } //if chit could be a good start point

                } else {
                    startPt_ = chit;
                } //if startPt_ is nullptr

                if ( finishPt_ ) {
                    
                    if ( clayer > finishPt_->getLayerID() and (itH->second).isGood ) {
                        finishPt_ = chit;
                    } //if chit could be a good finish point

                } else {
                    finishPt_ = chit;
                } //if finishPt_ is nullptr
            
            } //iterate through hits in this orientation hitlog (chit)
        
        } //two ore more hits in this part of the log

        return ( startPt_ and finishPt_ );

    }

    void HcalMipTriggerProducer::removeHcalID( HcalOrientation orientation , DetectorID::RawValue rawID ) {
        
        hitLog_[ orientation ].erase( rawID );

        return;
    }

}

DECLARE_PRODUCER_NS(ldmx, HcalMipTriggerProducer);