/**
 * @file HcalTrackProducer.cxx
 * @brief Implementation file for class HcalTrackProducer
 * @author Tom Eichlersmith, University of Minnesota
 */

#include "Hcal/HcalTrackProducer.h"

namespace ldmx {
    
    void HcalTrackProducer::configure( const ParameterSet& ps ) {
        
        hitcollname_ = ps.getString( "HitCollectionName" , "hcalDigis" );
        hitpassname_ = ps.getString( "HitPassName" , "recon" );

        nlayers_ = ps.getInteger( "NumHcalLayers" , 81 );
        nstrips_ = ps.getInteger( "NumHcalStrips" , 34 );

        //setting moduli to next highest order of 10
        layermod_ = (int)(pow( 10 , std::ceil( log10( nstrips_ ) ) ));
        sectionmod_ = (int)(pow( 10 , std::ceil( log10( nlayers_ ) ) ));

        minPE_ = static_cast<float>( ps.getDouble( "MinimumPE" , 0.0 ) );
        maxEnergy_ = static_cast<float>( ps.getDouble( "MaximumEnergy" , 4000.0 ) );

        firstseedlayer_ = ps.getInteger( "FirstSeedLayer" , 1 );
        
        conedepth_ = ps.getInteger( "SearchConeDepth" , 3 );
        coneangle_ = ps.getInteger( "SearchConeAngle" , 3 );
        minconehits_ = ps.getInteger( "MinConeHits" , 3 );
        
        trackwidth_ = ps.getInteger( "TrackWidth" , 3 );
        
        mintracklayhits_ = ps.getInteger( "MinTrackLayerHits" , 20 );
        
        maxtrackcnt_ = ps.getInteger( "MaxTrackCount", 100 );

        hcaltracksname_ = ps.getString( "HcalTrackCollectionName" , "HcalTracks" );
        hcaltracks_ = new TClonesArray( "ldmx::HcalTrack" , 1000 );
        
        return; 
    }

    void HcalTrackProducer::produce( Event& event ) {

        //initialize event containters
        log_.clear();
        
        layercheck_.clear();
        for ( int i = 1; i < nlayers_+1; i++ ) {
            layercheck_.insert( i );
        }

        badseeds_.clear();

        //obtain list of hits
        const TClonesArray *rawhits = event.getCollection( hitcollname_ , hitpassname_ );

        //pre-process hits and add to log
        for ( size_t i = 0; i < rawhits->GetEntriesFast(); i++ ) {
            HitPtr curr_hit = (HitPtr)(rawhits->At(i));
            if ( curr_hit->getPE() > minPE_ ) { //curr_hit is not noise
                AddHit( curr_hit );
            } //curr_hit is not noise
        } //iterate through rawhits (i)

        //search for tracks
        HcalTrack *track = new HcalTrack();
        int seedlayer = firstseedlayer_;
        int trackcnt = 0;
        while( TrackSearch( seedlayer , track ) and trackcnt < maxtrackcnt_ ) {
            //add track to collection
            HcalTrack *toadd = (HcalTrack *)(hcaltracks_->ConstructedAt(trackcnt));
            *toadd = *track; 

            //Remove track from log
            RemoveTrack( track );

            track->Clear(); //re-initialize track
            seedlayer = *layercheck_.begin(); //change seedlayer
            trackcnt++;
        } //keep searching for tracks until can't find anymore

        //add collection to event bus
        event.add( hcaltracksname_ , hcaltracks_ );
        
        return;
    }

    void HcalTrackProducer::AddHit( HitPtr hit ) {
        
        int key = KeyGen( hit );
        
        log_[ key ] = hit;

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
    
    bool HcalTrackProducer::TrackSearch( int seedlayer , HcalTrack *track ) {
        
        seedlayer_ = seedlayer;
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
        return( section*sectionmod_*layermod_ + layer*layermod_ + strip );
    }
    
    int HcalTrackProducer::KeyGen( HitPtr hit ) const {
        int section = static_cast<int>( hit->getSection() );
        int layer = static_cast<int>( hit->getLayer() );
        int strip = static_cast<int>( hit->getStrip() );
        return KeyGen( section , layer , strip );
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
        for ( std::set<int>::iterator it = layercheck_.begin(); it != layercheck_.end(); ++it ) {
            
            int l = *it;
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

        //sort layerlist according to proximity to seed
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
        
        float leftslope, rightslope;
        std::pair< HitPtr , HitPtr > leftmost( track->getHit(0), track->getHit(1) ), rightmost( track->getHit(0) , track->getHit(1) );
 
        while ( !layerlist_.empty() ) { //loop through elements of layerlist_
            
            int layer = layerlist_.front();
            layerlist_.pop_front();
            
            //calculate extremes for given layer oritentation
            //Find leftmost, secondleftmost, rightmost, secondrightmost (left and right sides could be equal)
            for ( int i = 0; i < track->getNHits(); i++ ) {
            
                HitPtr curr_hit = track->getHit( i );
                float curr_strip = curr_hit->getStrip();
                float curr_layer = curr_hit->getLayer();
                if ( true ) { //!(( layer ^ curr_layer ) & 1 ) ) { 
                //current layer and layer have same parity (same orientation)
                 
                    bool isdifleftlayer = ( std::abs(curr_layer - leftmost.first->getLayer()) >= 1.0 );
                    bool isdifrightlayer = ( std::abs(curr_layer - rightmost.first->getLayer()) >= 1.0 );

                    //Check if curr_strip is first or second most left
                    if ( curr_strip < leftmost.first->getStrip() ) {
                        if ( isdifleftlayer ) {
                            leftmost.second = leftmost.first;
                        } //check to make sure second most left is not on same layer
                        leftmost.first = curr_hit;
                    } else if ( curr_strip < leftmost.second->getStrip() and isdifleftlayer ) {
                        leftmost.second = curr_hit;
                    } //check to update second most (and make sure it isn't on the same level)
            
                    //Check if curr_strip is first or second most right
                    if ( curr_strip > rightmost.first->getStrip() ) {
                        if ( isdifrightlayer ) {
                            rightmost.second = rightmost.first;
                        }
                        rightmost.first = curr_hit;
                    } else if ( curr_strip > rightmost.second->getStrip() and isdifrightlayer ) {
                        rightmost.second = curr_hit;
                    }
 
                } //skip curr_hit if not same orientation as layer
           
            } //iterate through partial track (it)
        
            //Extend from leftmost to layer
            if ( leftmost.first->getLayer() - leftmost.second->getLayer() > 0.01 ) {
                leftslope = (leftmost.first->getStrip() - leftmost.second->getStrip())/(leftmost.first->getLayer() - leftmost.second->getLayer());
            } else {
                leftslope = 0.0;
            }

            //Extend from rightmost to layer
            if ( rightmost.first->getLayer() - rightmost.second->getLayer() > 0.01 ) {
                rightslope = (rightmost.first->getStrip() - rightmost.second->getStrip())/(rightmost.first->getLayer() - rightmost.second->getLayer());
            } else {
                rightslope = 0.0;
            }

            float leftedge = (layer - leftmost.first->getLayer())*leftslope + leftmost.first->getStrip();
            float rightedge = (layer - rightmost.first->getLayer())*rightslope + rightmost.first->getStrip();

            //Arithmetic mean of right edge and left edge
            float centertrack = (leftedge+rightedge)/2;
        
            //Define lowstrip and upstrip
            int lowstrip = static_cast<int>(std::floor( centertrack - trackwidth_/2.0 ));
            int upstrip = static_cast<int>(std::ceil( centertrack + trackwidth_/2.0 ));
            CorrectStrip( lowstrip );
            CorrectStrip( upstrip );

            //SearchByKey
            int lowkey = KeyGen( 0 , layer , lowstrip );
            int upkey = KeyGen( 0 , layer , upstrip );

            SearchByKey( lowkey , upkey , track );

        } //loop through elements of layerlist

        return isAcceptableTrack( track );
    }
    
    bool HcalTrackProducer::isAcceptableTrack( const HcalTrack *track ) const {
        //For now, only going by number of layer hits in track
        return ( track->getNLayHits() > mintracklayhits_ );
    }

    bool HcalTrackProducer::SearchByKey( const int lowkey , const int upkey , HcalTrack *track ) {
        
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
            switch( mipvec.size() ) {
                case 0:
                    //no mips in this key range
                    break;
                case 1:
                    //one mip in this key range
                    track->incLayHit(); //think about moving this to a more general spot (avoid double counting a layer)
                    track->addGroup( mipvec[0] );
                    success = true;
                    break;
                default:
                    //more than one mip found
                    std::cout << "[ HcalTrackProducer::SearchByKey ] More than one MIP found in key range " << lowkey << "->" << upkey << "\n";
                    std::cout << "                                   Adding first one found." << std::endl;
                    track->incLayHit();
                    track->addGroup( mipvec[0] );
                    success = true;
                    break;
            } //cases based on how many mips were found

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
