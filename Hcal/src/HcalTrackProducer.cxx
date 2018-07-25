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
                            mipLog_[ csect ][ ckey ] = cmip;
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
      
        //search for tracks
        HcalTrack *track = new HcalTrack();
        int trackcnt = 0;
        while( TrackSearch( track ) and trackcnt < maxtrackcnt_ ) {
            //add track to collection
            HcalTrack *toadd = (HcalTrack *)(hcaltracks_->ConstructedAt(trackcnt));
            *toadd = *track;
            
            //Remove track from log
            RemoveTrack( track );

            track->Clear(); //re-initialize track
            seedlayer_ = *layercheck_.begin(); //change seedlayer
            trackcnt++;
        } //keep searching for tracks until can't find anymore
        
        //add collection to event bus
        event.add( hcaltracksname_ , hcaltracks_ );
        
        //memory clean up
        delete track;
        
        /* Skim Rule depending on number of tracks
        if ( trackcnt > 1 ) {
            setStorageHint( hint_mustKeep );
        } else {
            setStorageHint( hint_mustDrop );
        }
        */        
        return;
    }
    
    bool HcalTrackProducer::TrackSearch( HcalTrack *track ) {
        
        seedstrip_ = 0;
        while ( FindSeed() ) { //seed found
            
            track->setSeed( seedlayer_ , seedstrip_ );
            
            SetSearchCone();
            
            //Checks if track is started successfully and then tries to
            // extend the track.
            if ( BeginPartialTrack( track ) and ExtendTrack( track ) ) { //track successfully constructed
                return true;
            }
            
            //bad seed
            badseeds_.insert( KeyGen( 0 , seedlayer_ , seedstrip_ ) );

        } //while possible good seed found

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
    
    bool HcalTrackProducer::FindSeed() {
        
        if ( layercheck_.empty() ) { //no more layers to check
            return false;
        }
       
        //check if seedlayer_ has been searched before
        std::set<int>::iterator seedlayer_it = layercheck_.find( seedlayer_ );

        if ( seedlayer_it != layercheck_.end() ) { //seedlayer_ hasn't been searched before
           
            //set keys to cover entire layer
            int lowkey = seedlayer_*layermod_;
            int upkey = (seedlayer_+1)*layermod_ - 1;
            
            HcalTrack *trackseed = new HcalTrack(); //temp HcalTrack for SearchByKey
            int seedkey = lowkey - 1; //preserve lowkey when entering loop
            bool badseedfound, mipfound;
           
            do { //SearchByKey found a mip and it is listed as a bad seed
                lowkey = seedkey+1;//move range up to start after bad seed
                mipfound = SearchByKey( lowkey , upkey , trackseed );//search for another mip 
                if ( mipfound ) { //mip found
                    seedstrip_ = static_cast<int>( trackseed->getHit(0)->getStrip() ); //reset seedstrip
                    seedkey = KeyGen( 0 , seedlayer_ , seedstrip_ ); //reset seedkey
                    badseedfound = ( badseeds_.find( seedkey ) != badseeds_.end() ); //see if new seed is bad
                } else {
                    badseedfound = false; //no mip found so exit loop
                }
                trackseed->Clear();
            } while ( badseedfound ); //SearchByKey found a mip and it is listed as a bad seed
            
            if ( mipfound ) { //loop exited with a mipfound, then it is not a bad seed
                return true;
            }
            
            //entire layer searched, no good seeds
            layercheck_.erase( seedlayer_it );

        } //seedlayer_ hasn't been searched before

        //Function would exit by now if it hadn't found a seed
        seedlayer_ = *layercheck_.begin();
        return ( FindSeed() );
    }
    bool HcalTrackProducer::BeginPartialTrack( HcalTrack *track ) {
        
        //make sure track is empty
        if ( track->getNHits() > 0 ) {
            track->Clear();
        }

        while ( !cone_.empty() ) { //loop through cone to find mips
            
            SearchByKey( cone_.front().first , cone_.front().second , track );
            cone_.pop_front();

        } //loop through cone to find mips

        //check if enough hits in cone
        return ( track->getNHits() >= minconehits_ );
    }
            
    bool HcalTrackProducer::ExtendTrack( HcalTrack *track ) {
        
        //Iterate through list of layers to search for hits
        bool addednewhit = true;
        while ( !layerlist_.empty() ) { //loop through elements of layerlist_
            
            int layer = layerlist_.front();
            layerlist_.pop_front();
            
            //Linearly extrapolate points to layer
            // HcalTrack::evalFit calculates fit and then evaluates
            float centerstrip = track->evalFit( layer );
            
            //Define lowstrip and upstrip
            int lowstrip = static_cast<int>(std::floor( centerstrip - trackwidth_/2.0 ));
            int upstrip = static_cast<int>(std::ceil( centerstrip + trackwidth_/2.0 ));
            CorrectStrip( lowstrip );
            CorrectStrip( upstrip );

            //SearchByKey
            int lowkey = KeyGen( 0 , layer , lowstrip );
            int upkey = KeyGen( 0 , layer , upstrip );

            addednewhit = SearchByKey( lowkey , upkey , track , centerstrip );
            
        } //loop through elements of layerlist

        return isAcceptableTrack( track );
    }
    
    bool HcalTrackProducer::isAcceptableTrack( const HcalTrack *track ) const {
        //For now, only going by number of layer hits in track
        return ( track->getNLayHits() > mintracklayhits_ );
    }

    bool HcalTrackProducer::SearchByKey( const std::map< int , MipHitPtr > log , const int lowkey , const int upkey , HcalTrack *track , const float prefstrip ) {
        
        bool success = false;

        if ( lowkey > upkey ) { //Mis-use correction
            std::cout << "Input hit keys to HcalTrackProducer::SearchByKey in wrong order: " << lowkey << "->" << upkey << std::endl;
            std::cout << "Returning an empty search" << std::endl;
        } else { //inputs are correct form
            
            auto lowbound = log.lower_bound( lowkey ); //points to first key that is not before lowkey (equivalent or after) (map::end if all are before lowkey)
            auto upbound = log.upper_bound( upkey ); //points to first key after upkey (map::end if nothing after upkey)

            //count number of mips and determine which one is closest to prefstrip
            int nmips = 0;
            MipHitPtr bestmip;
            float beststripdif = 10000.0;
            for( auto it = lowbound; it != upbound; ++it ) {
                nmips++;
                MipHitPtr cmip = it->second;
                float cstrip = static_cast<float>(cmip->getLowStrip() + cmip->getUpStrip())/2.0;

                float cstripdif = std::abs( prefstrip - cstrip );

                if ( cstripdif < beststripdif ) {
                    beststripdif = cstripdif;
                    bestmip = cmip;
                } //is cmip closer than bestmip to prefstrip
            } //iterate through [lowbound, upbound) range of log (it)
        
            //cases based on how many mips were found
            if ( nmips > 0 ) { //there are some mips
                
                track->incLayHit();
                track->addGroup( *bestmip );
                success = true;

            } //some mips were found

        } //make sure inputs are correct
        
        return success;
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
