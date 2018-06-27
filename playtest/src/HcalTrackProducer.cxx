/**
 * @file HcalTrackProducer.cxx
 * @brief Implementation file for class HcalTrackProducer
 * @author Tom Eichlersmith, University of Minnesota
 */

#include "playtest/HcalTrackProducer.h"

namespace ldmx {
    
    void HcalTrackProducer::configure( const ParameterSet& ps ) {
        
        hitcollname_ = ps.getString( "HitCollectionName" );

        nlayers_ = ps.getInteger( "NumHcalLayers" , 81 );
        nstrips_ = ps.getInteger( "NumHcalStrips" , 34 );

        layermod_ = ps.getInteger( "LayerModulus" , 1000 );

        minPE_ = static_cast<float>( ps.getDouble( "MinimumPE" , 0.0 ) );

        conedepth_ = ps.getInteger( "SearchConeDepth" , 3 );
        coneangle_ = ps.getInteger( "SearchConeAngle" , 3 );
        minconehits_ = ps.getInteger( "MinConeHits" , 3 );
        
        trackwidth_ = ps.getInteger( "TrackWidth" , 3 );
        
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

        while ( !cone_.empty() ) {
            cone_.pop();
        }

        while ( !layerlist_.empty() ) {
            layerlist_.pop();
        }

        badseeds_.clear();

        //obtain list of hits
        const TClonesArray *rawhits = event.getCollection( hitcollname_ );

        //pre-process hits and add to log
        for ( size_t i = 0; i < rawhits->GetEntriesFast(); i++ ) {
            HitPtr curr_hit = (HitPtr)(rawhits->At(i));
            if ( curr_hit->getPE() > minPE_ ) { //curr_hit is not noise
                AddHit( curr_hit );
            } //curr_hit is not noise
        } //iterate through rawhits (i)

        //search for tracks
        HcalTrack track;
        int seedlayer = 0;
        int trackcnt = 0;
        while( TrackSearch( seedlayer , track ) and trackcnt < 5 ) {
            //add track to collection
            HcalTrack *toadd = (HcalTrack *)(hcaltracks_->At(trackcnt));
            std::cout << "Point to new location in TClonesArray" << std::endl;
            toadd = new HcalTrack( track ); 
            std::cout << "Copy over track information" << std::endl;

            track.Clear(); //re-initialize track
            seedlayer = *layercheck_.begin(); //change seedlayer
            std::cout << seedlayer << std::endl;
            trackcnt++;
        } //keep searching for tracks until can't find anymore

        //add collection to event bus
        event.add( "HcalTracks" , hcaltracks_ );

        return;
    }

    void HcalTrackProducer::onFileOpen() {
         
    }

    void HcalTrackProducer::onFileClose() {
    
    }

    void HcalTrackProducer::onProcessStart() {

    }

    void HcalTrackProducer::onProcessEnd() {

    }

    void HcalTrackProducer::AddHit( HitPtr hit ) {
        
        int key = KeyGen( hit );
        log_[ key ] = hit;

        return;
    }
    
    bool HcalTrackProducer::TrackSearch( int seedlayer , HcalTrack &track ) {
        std::cout << "HcalTrackProducer::TrackSearch" << std::endl;
        int seedstrip = 0;
        while ( FindSeed( seedlayer , seedstrip ) ) { //seed found
            
            SetSearchCone( seedlayer , seedstrip );
            
            //Checks if track is started successfully and then tries to
            // extend the track.
            if ( BeginPartialTrack( track ) and ExtendTrack( track ) ) { //track successfully constructed
                
                return true;
            
            } else { //bad seed
                
                badseeds_.insert( seedlayer*layermod_ + seedstrip );
            
            } //track has good or bad seed

        } //while possible good seed found

        return false;
    }

    int HcalTrackProducer::KeyGen( HitPtr hit ) const {
        return static_cast<int>( hit->getLayer()*layermod_ + hit->getStrip() );
    }
    
    void HcalTrackProducer::CorrectStrip( int &strip ) const {
        
        if ( strip < 0 ) {
            strip = 0;
        } else if ( strip > nstrips_ ) {
            strip = nstrips_;
        }
        
        return;
    }
    
    bool HcalTrackProducer::FindSeed( int &seedlayer , int &seedstrip ) {
        std::cout << "HcalTrackProducer::FindSeed" << std::endl; 
        if ( layercheck_.empty() ) { //no more layers to check
            return false;
        }
       
        //check if seedlayer has been searched before
        std::set<int>::iterator seedlayer_it = layercheck_.find( seedlayer );

        if ( seedlayer_it != layercheck_.end() ) { //seedlayer hasn't been searched before
           
            //set keys to cover entire layer
            int lowkey = seedlayer*layermod_;
            int upkey = (seedlayer+1)*layermod_ - 1;
            HcalTrack trackseed;
            int seedkey = seedlayer*layermod_ + seedstrip;
            
            bool mipfound = SearchByKey( lowkey , upkey , trackseed );
            while ( badseeds_.find( seedkey ) != badseeds_.end() and mipfound ) { //SearchByKey found a mip and it is listed as a bad seed
                seedstrip = static_cast<int>( trackseed.getHit(0)->getStrip() );
                seedkey = seedlayer*layermod_ + seedstrip;
                lowkey = seedkey+1;
                trackseed.Clear();
                mipfound = SearchByKey( lowkey , upkey , trackseed );
            } //SearchByKey found a mip and it is listed as a bad seed

            if ( mipfound ) { //loop exited with a mipfound, then it is not a bad seed
                seedstrip = static_cast<int>( trackseed.getHit(0)->getStrip() );
                return true;
            } else { //entire layer searched, no good seeds
                layercheck_.erase( seedlayer_it );
            }

        } //seedlayer hasn't been searched before

        //Function would exit by now if it hadn't found a seed
        seedlayer = *layercheck_.begin();
        return ( FindSeed( seedlayer , seedstrip ) );
    }

    void HcalTrackProducer::SetSearchCone( const int seedlayer , const int seedstrip ) {
        
        //reset lists
        while ( !cone_.empty() ) {
            cone_.pop();
        }
        
        while ( !layerlist_.empty() ) {
            layerlist_.pop();
        } 
        
        //calculate slope of cone
        float slope = static_cast<float>( coneangle_ )/( conedepth_*2.0 );

        for ( int l = seedlayer - conedepth_; l < seedlayer + conedepth_ + 1; l++ ) {
            
            std::set<int>::iterator l_it = layercheck_.find( l );
            
            if ( l_it != layercheck_.end() ) { //current layer hasn't been check entirely
                //calculate strip numbers for layer
                float left = (l - seedlayer)*slope + seedstrip;
                float right = (l - seedlayer)*(-1)*slope + seedstrip;
                
                if ( left > right ) {
                    float tmp = left;
                    left = right;
                    right = tmp;
                } //they are in wrong order

                int lowstrip = static_cast<int>(std::floor( left ));
                int upstrip = static_cast<int>(std::ceil( right ));
                CorrectStrip( lowstrip );
                CorrectStrip( upstrip );

                //add keys to cone
                cone_.push( std::pair<int,int>( l*layermod_+lowstrip , l*layermod_+upstrip ) );

            } //current layer hasn't been exhaustively checked
        } //iterate through layers in cone (l)

        //Set Layer List (layers outside cone that haven't been searched exhaustively
        for ( std::set<int>::iterator it = layercheck_.begin(); it != layercheck_.end(); ++it ) {
            
            if ( *it < seedlayer - conedepth_ or *it > seedlayer + conedepth_ ) { //layer outside of cone
                layerlist_.push( *it );
            }

        } //iterate through layercheck_ (it)

        return;
    }
    
    bool HcalTrackProducer::BeginPartialTrack( HcalTrack &track ) {
        std::cout << "HcalTrackProducer::BeginPartialTrack" << std::endl;
        while ( !cone_.empty() ) { //loop through cone to find mips
            
            SearchByKey( cone_.front().first , cone_.front().second , track );
            cone_.pop();

        } //loop through cone to find mips

        //check if enough hits in cone
        return ( track.getNHits() >= minconehits_ );
    }
            
    bool HcalTrackProducer::ExtendTrack( HcalTrack &track ) {
        std::cout << "HcalTrackProducer::ExtendTrack" << std::endl;
        //check to see if track has been changed
        bool addednewhit = true;
        float leftslope, rightslope;
        std::pair< HitPtr , HitPtr > leftmost( track.getHit(0), track.getHit(1) ), rightmost( track.getHit(0) , track.getHit(1) );
 
        while ( !layerlist_.empty() ) { //loop through elements of layerlist_
            
            int layer = layerlist_.front();
            layerlist_.pop();
            
            if ( addednewhit ) { //track has been changed, so edit slope

                //Find leftmost, secondleftmost, rightmost, secondrightmost (left and right sides could be equal)
                for ( int i = 0; i < track.getNHits(); i++ ) {
                
                    HitPtr curr_hit = track.getHit( i );
                    float curr_strip = curr_hit->getStrip();
                    float curr_layer = curr_hit->getLayer();
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
                } //iterate through partial track (it)
            
                //Extend from leftmost to layer
                leftslope = (leftmost.first->getStrip() - leftmost.second->getStrip())/(leftmost.first->getLayer() - leftmost.second->getLayer());

                //Extend from rightmost to layer
                rightslope = (rightmost.first->getStrip() - rightmost.second->getStrip())/(rightmost.first->getLayer() - rightmost.second->getLayer());
            
            } //track has been changed, so edit slope

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
            int lowkey = layer*layermod_ + lowstrip;
            int upkey = layer*layermod_ + upstrip;

            addednewhit = SearchByKey( lowkey , upkey , track );

        } //loop through elements of layerlist_

        return isAcceptableTrack( track );
    }
    
    bool HcalTrackProducer::isAcceptableTrack( const HcalTrack track ) const {
        //For now, accepting all tracks
        return true;
    }

    bool HcalTrackProducer::SearchByKey( const int lowkey , const int upkey , HcalTrack &track ) const {
       
        bool success = false;

        if ( lowkey > upkey ) { //Mis-use correction
            std::cout << "Input hit keys to HcalTrackProducer::SearchByKey in wrong order: " << lowkey << " " << upkey << std::endl;
            std::cout << "Returning an empty search" << std::endl;
        } else { //inputs are correct form

            auto lowbound = log_.lower_bound( lowkey ); //points to first key that is not before lowkey (equivalent or after) (map::end if all are before lowkey)
            auto upbound = log_.upper_bound( upkey ); //points to first key after upkey (map::end if nothing after upkey)
        
            if ( lowbound != upbound ) { //check to see if range has any thickness
                //there is at least one hit in the range
                //lowbound points to the first one
                //see if there is a hit on either side of lowerbound
                
                //prev and next will return iterators outside of container
                //must check if lowbound is on an edge
                auto beforeside = std::prev( lowbound ); 
                auto afterside = std::next( lowbound );
                
                int beforekeydif = lowbound->first - beforeside->first;
                int afterkeydif = afterside->first - lowbound->first;
                
                if ( lowbound == log_.begin() ) { //there cannot be a beforeside
                    beforekeydif = 1000;
                }

                if ( lowbound == log_.end() ) { //there cannot be an afterside
                    afterkeydif = 1000;
                }

                if ( beforekeydif != 1 or afterkeydif != 1 ) {
                    //lowbound has at most one neighbor
                    track.addHit( lowbound->second );

                    if ( beforekeydif == 1 ) {
                        //beforeside is the neighbor for lowbound
                        track.addHit( beforeside->second );
                    } else if ( afterkeydif == 1 ) {
                        //afterside is the neighbor for lowerbound
                        track.addHit( afterside->second );
                    } //else: lowbound is truly isolated

                    success = true;

                } //check if lowbound could be isolated

            } //check to see if range has any thickness

        } //make sure inputs are correct
        
        return success;
    }

}

DECLARE_PRODUCER_NS( ldmx , HcalTrackProducer );
