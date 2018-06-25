/**
 * @file HitLog.cxx
 * @brief Implementation file for class HitLog
 * @author Tom Eichlersmith, University of Minnesota
 */

#include "playtest/HitLog.h"

namespace ldmx {
    
    HitLog::HitLog() : nlayers_(81), nstrips_(34), layermod_(1000), minPE_(0.0), conedepth_(2), coneangle_(3),
                       minconehits_(3), trackwidth_(3), origin_(17.5), lowside_(12.5), upside_(22.5) {
        for ( int i = 1; i < nlayers_+1; i++) {
            layercheck_.insert( i );
        }
    }

    HitLog::HitLog( const int nlayers , const int nstrips , const float minPE , const int conedepth , const int coneangle ,
                    const int minconehits , const int trackwidth , const float origin , const float lowside , const float upside ) :
            nlayers_(nlayers), nstrips_(nstrips), layermod_(1000), minPE_(minPE), conedepth_(conedepth), coneangle_(coneangle),
            minconehits_(minconehits), trackwidth_(trackwidth), origin_(origin), lowside_(lowside), upside_(upside) {
        for ( int i = 1; i < nlayers_+1; i++) {
            layercheck_.insert( i );
        }
    }


    void HitLog::AddHit( HitPtr hit ) {
        
        int key = KeyGen( hit );
        log_[ key ] = hit;

        return;
    }
    
    bool HitLog::TrackSearch( int seedlayer , std::vector< HitPtr > &track ) {
        
        int seedstrip = 0;
        while ( FindSeed( seedlayer , seedstrip ) ) { //seed found
            
            SetSearchCone( seedlayer , seedstrip );
            
            if ( BeginPartialTrack( track ) ) { //track successfully started
                
                ExtendTrack( track );
                return ( isTrackAcceptable( track ) );
            
            } else { //bad seed
                
                badseeds_.insert( seedlayer*layermod_ + seedstrip );
            
            } //track has good or bad seed

        } //while possible good seed found

        return false;
    }

    int HitLog::KeyGen( HitPtr hit ) const {
        return static_cast<int>( hit->getLayer()*layermod_ + hit->getStrip() );
    }
    
    void CorrectStrip( int &strip ) const {
        
        if ( strip < 0 ) {
            strip = 0;
        } else if ( strip > nstrips_ ) {
            strip = nstrips_;
        }
        
        return;
    }
    
    bool HitLog::FindSeed( int &seedlayer , int &seedstrip ) const {
        
        //check if seedlayer has been searched before
        std::set<int>::iterator seedlayer_it = layercheck_.find( seedlayer );

        if ( seedlayer_it != layercheck_.end() ) { //seedlayer hasn't been searched before
           
            //set keys to cover entire layer
            int lowkey = seedlayer*layermod_;
            int upkey = (seedlayer+1)*layermod_ - 1;
            std::vector< HitPtr > trackseed;
            int seedkey = seedlayer*layermod_ + seedstrip;
            
            bool mipfound = SearchByKey( lowkey , upkey , trackseed );
            while ( badseeds_.find( seedkey ) != badseeds_.end() and mipfound ) { //SearchByKey found a mip and it is listed as a bad seed
                seedstrip = static_cast<int>( trackseed[0]->getStrip() );
                seedkey = seedlayer*layermod_ + seedstrip;
                lowkey = seedkey+1;
                trackseed.clear();
                mipfound = SearchByKey( lowkey , upkey , trackseed );
            } //SearchByKey found a mip and it is listed as a bad seed

            if ( mipfound ) { //loop exited with a mipfound, then it is not a bad seed
                seedstrip = static_cast<int>( trackseed[0]->getStrip() );
                return true;
            } else { //entire layer searched, no good seeds
                layercheck_.erase( seedlayer_it );
            }

        } //seedlayer hasn't been searched before

        if ( layercheck_.empty() ) { //no more layers to check
            return false;
        }

        //Function would exit by now if it hadn't found a seed
        seedlayer = static_cast<int>( (*layercheck_.begin())->getLayer() );
        return ( FindSeed( seedlayer , seedstrip ) );
    }

    void HitLog::SetSearchCone( const int seedlayer , const int seedstrip ) {
        
        //reset lists
        cone_.clear();
        layerlist_.clear();
        
        //calculate slope of cone
        float slope = static_cast<float>( conewidth_ )/( conedepth_*2.0 );

        for ( int l = seedlayer - conedepth_; l < seedlayer + conedepth_ + 1; l++ ) {
            
            std::set<int>::iterator l_it = layercheck_.find( l );
            
            if ( l_it != layercheck_.end() ) { //current layer hasn't been check entirely
                //make keys
                float left = (l - seedlayer)*slope + seedstrip;
                float right = (l - seedlayer)*(-1)*slope + seedstrip;
                
                if ( left > right ) {
                    float tmp = left;
                    left = right;
                    right = tmp;
                } //they are in wrong order

                int lowstrip = static_cast<int>(std::floor( left ));
                int upstrip = static_cast<int>(std::ceil( right ));
                CorrectStrip( lowkey );
                CorrectStrip( upkey );

                //add them to cone
                cone_.push( std::pair< l*layermod_+lowstrip , l*layermod_+upstrip > );

            }
        } //iterate through layers in cone (l)

        //Set Layer List (layers outside cone that haven't been searched exhaustively
        for ( std::set<int>::iterator it = layercheck_.begin(); it != layercheck_.end(); ++it ) {
            
            if ( *it < seedlayer - conedepth_ or *it > seedlayer + conedepth_ ) { //layer outside of cone
                layerlist_.push( *it );
            }

        } //iterate through layercheck_ (it)

        return;
    }
    
    bool HitLog::BeginPartialTrack( std::vector< HitPtr > &track ) const {
        
        while ( !cone_.empty() ) { //loop through cone to find mips
            
            SearchByKey( cone_.front()->first , cone.front()->second , track );
            cone_.pop();

        } //loop through cone to find mips

        //check if enough hits in cone
        return ( track.size() >= minconehits_ );
    }
            
    bool HitLog::ExtendTrack( std::vector< HitPtr > &track ) const {
        
        //check to see if track has been changed
        bool addednewhit = true;
        float leftslope, rightslope;
        std::pair< HitPtr , HitPtr > leftmost( track[0], track[1] ), rightmost( track[0] , track[1] );
 
        while ( !layerlist_.empty() ) { //loop through elements of layerlist_
            
            int layer = layerlist_.front();
            layerlist_.pop();
            
            if ( addednewhit ) { //track has been changed, so edit slope

                //Find leftmost, secondleftmost, rightmost, secondrightmost (left and right sides could be equal)
                for ( std::vector< HitPtr >::iterator it = track.begin(); it != track.end(); ++it ) {
                
                    float curr_strip = (*it)->getStrip();
            
                    //Check if curr_strip is first or second most left
                    if ( curr_strip < leftmost.first->getStrip() ) {
                        leftmost.second = leftmost.first;
                        leftmost.first = (*it);
                    } else if ( curr_strip < leftmost.second->getStrip() ) {
                        leftmost.second = (*it);
                    }   
            
                    //Check if curr_strip is first or second most right
                    if ( curr_strip > rightmost.first->getStrip() ) {
                        rightmost.second = rightmost.first;
                        rightmost.first = (*it);
                    } else if ( curr_strip > rightmost.second->getStrip() ) {
                        rightmost.second = (*it);
                    }
                } //iterate through partial track (it)
            
                //Extend from leftmost to layer
                leftslope = (leftmost.first->getStrip() - leftmost.second->getStrip())/(leftmost.first->getLayer() - leftmost.second->getLayer());

                //Extend from rightmost to layer
                rightslope = (rightmost.first->getStrip() - rightmost.second->getStrip())/(rightmost.first->getLayer() - rightmost.second->getLayer());
            
            } //track has been changed, so edit slope

            float leftedge = (layer - leftmost.first->getLayer())*slope + leftmost.first->getStrip();
            float rightedge = (layer - rightmost.first->getLayer())*slope + rightmost.first->getStrip();

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

        return true;
    }
    
    bool HitLog::isAcceptableTrack( const std::vector< HitPtr > track ) const {
        //For now, accepting all tracks
        return true;
    }

    bool HitLog::SearchByKey( const int lowkey , const int upkey , std::vector< HitPtr > &track ) const {
       
        bool success = false;

        if ( lowkey > upkey ) { //Mis-use correction
            std::cout << "Input hit keys to HitLog::SearchByKey in wrong order: " << lowkey << " " << upkey << std::endl;
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
                    track.push_back( lowbound->second );

                    if ( beforekeydif == 1 ) {
                        //beforeside is the neighbor for lowbound
                        track.push_back( beforeside->second );
                    } else if ( afterkeydif == 1 ) {
                        //afterside is the neighbor for lowerbound
                        track.push_back( afterside->second );
                    } //else: lowbound is truly isolated

                    success = true;

                } //check if lowbound could be isolated

            } //check to see if range has any thickness

        } //make sure inputs are correct
            
        return success;
    }

}
