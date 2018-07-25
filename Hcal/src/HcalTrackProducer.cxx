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
                } //are we in a different group
            } //does itH have a previous iterator

            cgroup.push_back( itH->second );
        } //iterate through rawLog_, grouping hits into mip hits (itH)
      
        //search for tracks
        HcalTrack *track = new HcalTrack();
        seedlayer_ = firstseedlayer_;
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

    void HcalTrackProducer::RemoveTrack( const HcalTrack *track ) {
        
        bool alreadywarned = false;
        for ( int i = 0; i < track->getNHits(); i++ ) {
            
            std::map< int , HitPtr >::iterator toremove = log_.find( KeyGen( track->getHit(i) ) );
            if ( toremove == log_.end() ) {
                if ( !alreadywarned ) {
                    std::cout << "[ HcalTrackProducer::RemoveTrack ]: ";
                    std::cout << "Unable to locate hit to be removed from log." << std::endl;
                    std::cout << "                                    ";
                    std::cout << "This bodes ill for how this producer was defined." << std::endl;
                    alreadywarned = true;
                }
            } else {
                log_.erase( toremove );
            }
        
        } //iterate through track and remove each hit (i)

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
    
    void HcalTrackProducer::CorrectStrip( int &strip ) const {
        
        if ( strip < 0 ) {
            strip = 0;
        } else if ( strip > nstrips_ ) {
            strip = nstrips_;
        }
        
        return;
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

    void HcalTrackProducer::SetSearchCone() {
        
        //reset lists
        cone_.clear();
        layerlist_.clear();
        
        //calculate slope of cone
        float slope = static_cast<float>( coneangle_ )/( conedepth_*2.0 );

        //Set Cone and Layer Lists 
        for ( int l = 1; l < nlayers_+1; l++ ) {
        //for ( std::set<int>::iterator it = layercheck_.begin(); it != layercheck_.end(); ++it ) {
            
        //    int l = *it;
            if ( l < seedlayer_ - conedepth_ or l > seedlayer_ + conedepth_ ) { //layer outside of cone
                layerlist_.push_back( l );
            } else { //layer inside cone
                
                int centerstrip, lowstrip, upstrip;
                if ( true ) { //!(( seedlayer ^ l ) & 1 ) ) { 
                    //current layer and seedlayer have same parity (same orientation)
                    centerstrip = seedstrip_; 
                
                    //calculate current halfwidth of cone
                    float halfwidth = std::abs(slope*(l - seedlayer_))/2.0;
                
                    lowstrip = static_cast<int>(std::floor( centerstrip - halfwidth ));
                    upstrip = static_cast<int>(std::ceil( centerstrip + halfwidth ));
                
                } else { 
                    //current layer has different orientation
                    //no way to determine centerstrip because we have no track direction yet
                    //search entire layer
                    lowstrip = 0;
                    upstrip = nstrips_;
                }
 
                CorrectStrip( lowstrip );
                CorrectStrip( upstrip );

                //add keys to cone
                cone_.push_back( std::pair<int,int>( KeyGen( 0 , l , lowstrip ), KeyGen( 0 , l , upstrip ) ) );

            } //layer in or out of cone

        } //iterate through layercheck_ (it)

        //sort layerlist according to proximity to seed RUNTIME ERROR
//        std::sort( layerlist_.begin() , layerlist_.end() , 
//            [&]( const int &a , const int &b ) -> bool {
//                return ( std::abs(a - this->seedlayer_) < std::abs(b - this->seedlayer_) );
//            });

        return;
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

    bool HcalTrackProducer::SearchByKey( const int lowkey , const int upkey , HcalTrack *track , const float prefstrip ) {
        
        bool success = false;

        if ( lowkey > upkey ) { //Mis-use correction
            std::cout << "Input hit keys to HcalTrackProducer::SearchByKey in wrong order: " << lowkey << "->" << upkey << std::endl;
            std::cout << "Returning an empty search" << std::endl;
        } else { //inputs are correct form
            
            auto lowbound = log_.lower_bound( lowkey ); //points to first key that is not before lowkey (equivalent or after) (map::end if all are before lowkey)
            auto upbound = log_.upper_bound( upkey ); //points to first key after upkey (map::end if nothing after upkey)
            std::vector< std::vector<HitPtr> > mipvec; //list of mips in this key range

            //group hits based on key separation and then check if a group is a mip
            std::vector<HitPtr> curr_group;
            for( auto it = lowbound; it != upbound; ++it ) {
                
               if ( it == log_.begin() ) {
                    //beginning of log, isolated on left
                    if ( curr_group.size() > 0 ) {
                        std::cout << "[ HcalTrackProducer::SearchByKey ] Iterator reached begining of log after adding to a group.\n";
                        std::cout << "                                   Something was jumbled." << std::endl;
                        curr_group.clear();
                    }

                } else {
                    //inside log, find previous hit
                    auto previt = std::prev( it );
 
                    int keydif = it->first - previt->first;
                    
                    if ( keydif != 1 ) {
                        //different group
                        if ( isMIP( curr_group ) ) {
                            mipvec.push_back( curr_group );
                        }
                        curr_group.clear();
                    } //checking key difference 
                } //if it is log_.begin()
                
                curr_group.push_back( it->second );
            
            } //iterate through [lowbound, upbound) range of log (it)
            
            //check last group (may include upbound)
            // improve? need to scan past upbound to finish constructing last group?
            if ( isMIP( curr_group ) ) {
                mipvec.push_back( curr_group );
            }
        
            //mipvec has group(s) considered mips
            //cases based on how many mips were found
            int nmips = mipvec.size();
            if ( nmips > 0 ) { //there are some mips
                
                auto bestmip = mipvec.begin();
                if ( prefstrip > 0  and nmips > 1 ) {
                    
                    float beststripdif = static_cast<float>( nstrips_+1 ), currstripdif;
                    for ( auto it = mipvec.begin(); it != mipvec.end(); ++it ) {
                        //calculate currstrip (energy weighted average of strips in group)
                        float currstrip = 0.0;
                        float currenergy = 0.0;
                        for ( int i = 0; i < it->size(); i++ ) {
                            float en = it->at(i)->getEnergy();
                            currstrip += it->at(i)->getStrip()*en;
                            currenergy += en;
                        } //iterate throug currmip (i)
                        currstrip = currstrip/currenergy;

                        currstripdif = std::abs( prefstrip - currstrip );
                        
                        if ( currstripdif < beststripdif ) {
                            beststripdif = currstripdif;
                            bestmip = it;
                        } //check if currmip is better
                 
                    } //iterate through all mips found (it)
                
                } //check if there is an preferred key and a decision needs to be made
                
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
