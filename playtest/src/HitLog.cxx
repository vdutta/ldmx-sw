/**
 * @file HitLog.cxx
 * @brief Implementation file for class HitLog
 * @author Tom Eichlersmith, University of Minnesota
 */

#include "playtest/HitLog.h"

namespace ldmx {
    
    HitLog::HitLog() : nlayers_(81), nstrips_(34), layermod_(1000), minPE_(0.0), conedepth_(2), coneangle_(3),
                       trackwidth_(3), layercheck_( nlayers_ , false ), origin_(17.5), lowside_(12.5), upside_(22.5) { }

    HitLog::HitLog( const int nlayers , const int nstrips , const float minPE , const int conedepth , const int coneangle ,
                    const int trackwidth , const float origin , const float lowside , const float upside ) :
            nlayers_(nlayers), nstrips_(nstrips), layermod_(1000), minPE_(minPE), conedepth_(conedepth), coneangle_(coneangle),
            trackwidth_(trackwidth), layercheck_( nlayers_ , false ), origin_(origin), lowside_(lowside), upside_(upside) { }

    void HitLog::AddHit( HitPtr hit ) {
        
        int key = KeyGen( hit );
        log_[ key ] = hit;

        return;
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

    }

    void HitLog::SetSearchCone( const int seedlayer , const int seedstrip ) {
        
    }
    
    bool HitLog::BeginPartialTrack( std::vector< HitPtr > &track ) const {
    
    }
            
    bool HitLog::SearchLayer( const int layer , std::vector< HitPtr > &track ) const {
        //Find leftmost, secondleftmost, rightmost, secondrightmost (left and right sides could be equal)
        std::pair< HitPtr , HitPtr > leftmost( track[0], track[1] ), rightmost( track[0] , track[1] );
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
        float leftslope = (leftmost.first->getStrip() - leftmost.second->getStrip())/(leftmost.first->getLayer() - leftmost.second->getLayer());
        float leftedge = (layer - leftmost.first->getLayer())*slope + leftmost.first->getStrip();

        //Extend from rightmost to layer
        float rightslope = (rightmost.first->getStrip() - rightmost.second->getStrip())/(rightmost.first->getLayer() - rightmost.second->getLayer());
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
        
        return ( SearchByKey( lowkey , upkey , track ) );
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
