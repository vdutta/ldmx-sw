/**
 * @file HcalTrackProducer.cxx
 * @brief Implementation file for class HcalTrackProducer
 * @author Tom Eichlersmith, University of Minnesota
 */

#include "Hcal/HcalTrackProducer.h"

namespace ldmx {
    
    HcalTrackProducer::~HcalTrackProducer() {
        //clean up class used for processing
        for ( auto itS = mipLog_.begin(); itS != mipLog_.end(); ++itS ) {
            for ( auto itH = (itS->second).begin(); itH != (itS->second).end(); ++itH ) {
                delete itH->second;
            } //iterate through hits in current section (itH)
        } //iterate through sections in mipLog (itS)
    }
    
    void HcalTrackProducer::configure( const ParameterSet& ps ) {
        
        hitcollname_ = ps.getString( "HitCollectionName" );
        hitpassname_ = ps.getString( "HitPassName" );
        
        nLayersBack_ = ps.getInteger( "NumBackHcalLayers" );
        nStripsBack_ = ps.getInteger( "NumBackHcalStrips" );
        nLayersTopBot_ = ps.getInteger( "NumTopBottomHcalLayers" );
        nStripsTopBot_ = ps.getInteger( "NumTopBottomHcalStrips" );
        nLayersLeftRight_ = ps.getInteger( "NumLeftRightHcalLayers" );
        nStripsLeftRight_ = ps.getInteger( "NumLeftRightHcalStrips" );
            
        //setting moduli to next highest order of 10
        layerMod_ = (int)(pow( 10 , std::ceil( log10( nStripsBack_ ) ) ));
        sectionMod_ = (int)(pow( 10 , std::ceil( log10( nLayersBack_ ) ) ));

        minPE_ = static_cast<float>( ps.getDouble( "MinimumPE" ) );
        maxEnergy_ = static_cast<float>( ps.getDouble( "MaximumEnergy" ) );

        trackWidth_ = ps.getInteger( "TrackWidth" );
        
        minTrackLayHits_ = ps.getInteger( "MinTrackLayerHits" );
        
        maxTrackCnt_ = ps.getInteger( "MaxTrackCount" );

        hcalTracksName_ = ps.getString( "HcalTrackCollectionName" );
        hcalTracks_ = new TClonesArray( "ldmx::HcalTrack" , 1000 );
        
        return; 
    }

    void HcalTrackProducer::produce( Event& event ) {
   
        //initialize event containters
        nonoiseLog_.clear();
        mipLog_.clear();
        trackLog_.clear();
    
        //obtain list of hits
        const TClonesArray *rawhits = event.getCollection( hitcollname_ , hitpassname_ );
     
        //add rawhits to rawLog_ so that they get sorted
        for ( size_t i = 0; i < rawhits->GetEntriesFast(); i++ ) {
            HitPtr chit = (HitPtr)(rawhits->At(i));
            int ckey = KeyGen( chit );
            if ( chit->getPE() > minPE_ ) { //chit is not noise
                nonoiseLog_[ ckey ] = chit;
            } //chit is not noise
        } //iterate through rawhits (i)
        
        //filter rawLog_ into mip hits
        std::vector< HitPtr > cgroup;
        for ( auto itH = rawLog_.begin(); itH != rawLog_.end(); ++itH ) {
            
            if ( itH != rawLog_.begin() ) {
                //itH has a previous iterator     
                auto previtH = std::prev( itH );

                int keydif = itH->first - previtH->first;

                if ( keydif > 1 ) {
                    //different group of hits
                    if ( isMIP( cgroup ) ) {
                        //cgroup is a mip
                        //add new MipHit to mipLog
                        MipHitPtr cmip = new MipHit( cgroup );
                        if ( cmip->SetUp() ) {
                            HcalSection csect = cmip->getSection();
                            int ckey = KeyGen( cmip );
                            mipLog_[ ckey ] = cmip;
                        } else {
                            std::cout << "Setting up MIP Hit went wrong" << std::endl;
                            delete cmip;
                        }
                    } //is cgroup a mip

                    cgroup.clear(); //wipe current group

                } //are we in a different group

            } //does itH have a previous iterator

            cgroup.push_back( itH->second );
        } //iterate through rawLog_, grouping hits into mip hits (itH)
      
        //search for tracks in each section
        //merge tracks from each section into full detector tracks
        //add full detector tracks to collection
        
        //add collection to event bus
        event.add( hcaltracksname_ , hcaltracks_ );
        
        /* Skim Rule depending on number of tracks
        if ( trackcnt > 1 ) {
            setStorageHint( hint_mustKeep );
        } else {
            setStorageHint( hint_mustDrop );
        }
        */        
        return;
    }
    
    bool HcalTrackProducer::TrackSearch( HcalSection section ) {
        
        //Determine two end points
        //Find all mips within traversing cylinder
        //Check if track is acceptable
        // YES: add track to trackLog and remove hits from mipLog
        // NO : add end points to badendpoints set
        //Repeat until out of possible end point pairs

        return false;
    }

    int HcalTrackProducer::KeyGen( const int section , const int layer , const int strip ) const {
        return( section*sectionMod_*layerMod_ + layer*layerMod_ + strip );
    }
    
    int HcalTrackProducer::KeyGen( HitPtr hit ) const {
        int section = static_cast<int>( hit->getSection() );
        int layer = static_cast<int>( hit->getLayer() );
        int strip = static_cast<int>( hit->getStrip() );
        return KeyGen( section , layer , strip );
    }

    int HcalTrackProducer::KeyGen( MipHitPtr mip ) const {
        int section = static_cast<int>( mip->getSection() );
        int layer = mip->getLayer();
        int lowstrip = mip->getLowStrip();
        return KeyGen( section , layer , lowstrip );
    }

    bool HcalTrackProducer::findEndPts( HcalSection section , std::pair< MipHitPtr > &endpts ) {
        
        bool success = false;

        if ( !mipLog_.empty() ) { //mip log for this section is nonempty
            //Get first mip
            //Iterate through log finding mip furthest from start that isn't in badEndPts
            //Call this mip end
            //Iterate through log finding mip furthest from end that isn't in badEndPts
            //Call this mip start
            //Return pair
        } //miplog for this section is non emtpy

        return success;
    }
    
    bool HcalTrackProducer::connectTrack( HcalSection section , const std::pair< MipHitPtr > &endpts , std::vector< MipHitPtr > &possibletrack ) {
        
        //draw line between end points
        //iterate through log collecting mips that are within trackwidth of line

        return isAcceptableTrack( possibletrack );
    }

    bool HcalTrackProducer::isAcceptableTrack( const std::vector< MipHitPtr > &track ) const {
        //For now, keeping all tracks
        return ( true );
    }

    bool HcalTrackProducer::isMIP( const std::vector< HitPtr > &group ) const {
        
        unsigned int groupsize = group.size();
        //check if group is too small or too big in number
        if ( groupsize < 1 or groupsize > 2 ) {
            return false;
        }
        
        //calculate total energy of group
        float groupE = 0.0;
        for ( int i = 0; i < groupsize; i++ ) {
            groupE += group[i]->getEnergy();
        }

        return ( groupE < maxEnergy_ );
    }

}

DECLARE_PRODUCER_NS( ldmx , HcalTrackProducer );
