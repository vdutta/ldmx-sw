/**
 * @file HcalMuonTriggerProducer.cxx
 * @brief Implementaiton of HcalMuonTriggerProducer class
 * @author Tom Eichlersmith, University of Minnesota
 */

#include "Hcal/HcalMuonTriggerProducer.h"

namespace ldmx {

    void HcalMuonTriggerProducer::configure(const ldmx::ParameterSet& ps) {
        
        hitCollName_ = ps.getString( "HcalHitCollectionName" );
        hitPassName_ = ps.getString( "HcalHitPassName" );
        
        muonOrigin_ = ps.getString( "HcalMuonOrigin" );
        
        int minL = 0, maxL = 2;
        if ( muonOrigin_ == "Cosmic" ) {
            //Cosmic Muons
            // Use layers in side hcal and strips in Back
            minL = 2;
            maxL = 6;
        } else if ( muonOrigin_ != "Target" ) {
            std::cerr << "WARNING [ HcalMuonTrigger ] : Unknown Muon Origin. Defaulting to Target" << std::endl;
        }
        for ( int l = minL; l < maxL; l++ )
            layerUsers_.insert( l );

        trackRadius_ = ps.getDouble( "TrackRadius" );

        minFracHit_ = ps.getDouble( "MinFractionHit" );
        
        absoluteMinHits_ = ps.getInteger( "AbsoluteMinNumberHits" );

        maxEnergy_ = ps.getDouble( "MaximumEnergy" );

        minPE_ = ps.getDouble( "MinimumPE" );

        triggerObjectName_ = ps.getString( "HcalMuonTriggerObjectName" );

        numPass_ = 0;

        return;
    }

    void HcalMuonTriggerProducer::produce(ldmx::Event& event) {
        
        //initialize event containers
        for ( int iO = 0; iO < 6; iO++ ) {
            hitLog_[ iO ].clear();
        }
        
        //obtain list of raw hits
        const TClonesArray *rawhits = event.getCollection( hitCollName_ , hitPassName_ );

        //add only plausible mip hits to hitLog
        for ( size_t iH = 0; iH < rawhits->GetEntriesFast(); iH++ ) {
            ldmx::HcalHit* chit = (ldmx::HcalHit*)(rawhits->At(iH));
            if ( isPlausibleMip( chit ) ) { //could be a MIP
                
                int corient = chit->getSection();
                int clayer = chit->getLayer();
                int cstrip = chit->getStrip();
                
                //The underlying ints for the HcalSection and HcalOrientation enums are shifted by one
                // except for EVEN layers in the BACK Hcal.
                if ( corient != HcalSection::BACK or clayer % 2 == 1 ) {
                    corient++;
                }
                
                HitLogNode cNode;
                cNode.layer = clayer;
                cNode.strip = cstrip;
                cNode.isUsed = false;
                hitLog_[ corient ][ chit->getID() ] = cNode;
            } //could be a MIP
        } //iterate through rawhits (iH)
        
        //Find (and count) tracks
        int trackcnt = 0;
        for ( int corient = 0; corient < 6; corient++ ) {
            
            //lowest number of hits to accept track
            int hitCntFloor = absoluteMinHits_;
            if ( minFracHit_*hitLog_[ corient ].size() > hitCntFloor )
                hitCntFloor = static_cast<int>( minFracHit_*hitLog_[ corient ].size() );

            //while good end points are being found
            int longtracklen = 0; //num hits in longest track
            while ( findEndPoints( corient ) ) {
                
                int hitcnt = 0; //will not double count hits in same layer (if shallow) or same strip (if steep)
                int startLayer = (startPt_->second).layer;
                int startStrip = (startPt_->second).strip;
                int dstrip = ( (finishPt_->second).strip - startStrip );
                int dlayer = ( (finishPt_->second).layer - startLayer );
                
                bool shouldUseLayers = true;
                if ( layerUsers_.find( corient ) == layerUsers_.end() )
                    shouldUseLayers = false;

                if ( shouldUseLayers and abs(dlayer) > 0 ) {
                    //BACK HCAL, layer is independent variable for target muons
                    
                    //calculate slope
                    float slope = static_cast<float>(dstrip)/static_cast<float>(dlayer);
                    
                    //count hits in this orientation that are in the cylinder
                    std::set<int> countedLayers; //layers that already have been counted
                    for ( auto node : hitLog_[ corient ] ) {
                        
                        int cstrip = (node.second).strip;
                        int clayer = (node.second).layer;
    
                        float trackstrip = slope*( clayer - startLayer ) + startStrip;
    
                        float stripdif = std::abs( trackstrip - cstrip );
    
                        if ( stripdif < trackRadius_ ) {
                            //hit is in track cylinder
                            hitcnt++;    
                        } //check if hit is in track cylinder
    
                    } //iterate through hitLog of current orientation (node)
                
                } else if ( abs(dstrip) > 0 ) {
                    //SIDE HCAL, strip is independent variable for target muons

                    //calculate slope
                    float slope = static_cast<float>(dlayer)/static_cast<float>(dstrip);

                    //count hits in this orientation that are in the track cylinder
                    std::set<int> countedStrips; //strips that have already been counted
                    for ( auto node : hitLog_[ corient ] ) {
                        
                        int cstrip = (node.second).strip;
                        int clayer = (node.second).layer;

                        float tracklayer = slope*( cstrip - startStrip ) + startLayer;

                        float layerdif = std::abs( tracklayer - clayer );

                        if ( layerdif < trackRadius_ ) {
                            //hit is in track cylinder
                            hitcnt++;
                         } //check if hit is in the track cylinder

                    } //iterate through hitLog of current orientation (node)
                
                } //steep/shallow slope
                //skips end points if -1 <= dlayer <= 1 and -1 <= dstrip <= 1

                if ( hitcnt > longtracklen ) {
                    //good track found
                    longtracklen = hitcnt;
                } 

                //mark the end points as used end points
                (startPt_->second).isUsed = true;
                (finishPt_->second).isUsed = true;
        
            } //while good end points are still being found
            
            //check if longest "track" found can be considered a track
            if ( longtracklen > hitCntFloor )
                trackcnt++;

        } //for each orientation (corient)
        
        bool pass = false;
        if ( trackcnt > 0 ) { 
            pass = true;
            numPass_++;
        }

        result_.set( triggerObjectName_ , pass , 5 );
        result_.setAlgoVar( 0 , minPE_ );
        result_.setAlgoVar( 1 , maxEnergy_ );
        result_.setAlgoVar( 2 , minFracHit_ );
        result_.setAlgoVar( 3 , trackRadius_ );
        result_.setAlgoVar( 4 , static_cast<double>( trackcnt ) );

        event.addToCollection( "Trigger" , result_ );
        
        numTracksPerEvent_[ trackcnt ] ++;

        return;
    }
    
    void HcalMuonTriggerProducer::onProcessEnd() {
        
        printf( "\n" );
        printf( " ============================================\n" );
        printf( " | HcalMuonTriggerProducer | %14s |\n" , muonOrigin_.c_str() );
        printf( " |==========================================|\n" );
        printf( " |          Num Passed : %-19d|\n" , numPass_ );
        printf( " |==========================================|\n" );
        printf( " |            N Tracks : N Events           |\n" );
        for ( auto track_event : numTracksPerEvent_ ) {
            printf( " |%20d : %-19d|\n" , track_event.first , track_event.second );
        }
        printf( " ============================================\n" );
        
        return;
    }

    bool HcalMuonTriggerProducer::isPlausibleMip( ldmx::HcalHit* hit ) const {
        return ( hit->getPE() > minPE_ and hit->getEnergy() < maxEnergy_ );
    }

    bool HcalMuonTriggerProducer::findEndPoints( int orientation ) {
        
        startPt_ = hitLog_[ orientation ].end();
        finishPt_ = hitLog_[ orientation ].end();

        if ( hitLog_[ orientation ].size() >= absoluteMinHits_ ) {
            
            for ( std::map< unsigned int , HitLogNode >::iterator it = hitLog_[ orientation ].begin();
                it != hitLog_[ orientation ].end(); ++it ) {

                if ( !(it->second).isUsed ) {
                
                    int clayer = (it->second).layer;
                    int cstrip = (it->second).strip;
    
                    if ( startPt_ == hitLog_[ orientation ].end() //haven't chosen a starting point yet
                         or clayer < (startPt_->second).layer //lower layer
                         or (clayer == (startPt_->second).layer and cstrip < (startPt_->second).strip) //same layer but lower strip
                       ) {
                        startPt_ = it;
                    } //if it->second could be a good start point
    
                    if ( finishPt_ == hitLog_[ orientation ].end() //haven't chosen a finishing point yet
                         or clayer > (finishPt_->second).layer //higher layer
                         or (clayer == (finishPt_->second).layer and cstrip > (finishPt_->second).strip) //same layer but higher strip
                       ) {
                        finishPt_ = it;
                    } //if it->second could be a good finish point
                
                } //is current hit good

            } //iterate through hits in this orientation hitlog (it)
        
        } //two ore more hits in this part of the log

        return ( startPt_ != hitLog_[ orientation ].end() and finishPt_ != hitLog_[ orientation ].end() );

    }

}

DECLARE_PRODUCER_NS( ldmx , HcalMuonTriggerProducer );
