/**
 * @file HcalTrack.cxx
 * @brief Class implementation of a track through the Hcal
 * @author Tom Eichlersmith, University of Minnesota
 */

#include "Event/HcalTrack.h"

ClassImp( ldmx::HcalTrack );

namespace ldmx {

    HcalTrack::~HcalTrack() {
        Clear();
        hits_->Delete();
    }
    
    HcalTrack& HcalTrack::operator= ( const HcalTrack &track ) {
        
        if ( this != &track ) { //self-assignment guard
            
            this->hits_ = new TRefArray( *track.hits_ );
            this->setSeed( track.getSeedLayer() , track.getSeedStrip() );
            this->nlayhits_ = track.getNLayHits();

        }

        return *this;
    }
             
    void HcalTrack::Clear( Option_t *opt ) {

        TObject::Clear();

        hits_->Clear(); 
        nlayhits_ = 0;

        seedlayer_ = 0;
        seedstrip_ = 0;

        return;
    }

    void HcalTrack::addHit( HitPtr hit ) {
        hits_->Add( hit );
        return;
    }

    void HcalTrack::incLayHit() {
        nlayhits_++;
        return;
    }

    void HcalTrack::setSeed(int seedlayer , int seedstrip) {
        this->seedlayer_ = seedlayer;
        this->seedstrip_ = seedstrip;
        return;
    }

    void HcalTrack::addGroup( const std::vector<HitPtr> group ) {
        for( auto it = group.begin(); it != group.end(); ++it ) {
            addHit( *it );
        }
        return;
    }

    int HcalTrack::getNHits() const {
        return hits_->GetEntriesFast();
    }

    int HcalTrack::getNLayHits() const {
        return nlayhits_;
    }
    
    int HcalTrack::getSeedLayer() const {
        return seedlayer_;
    } 
    
    int HcalTrack::getSeedStrip() const {
        return seedstrip_;
    }

    HitPtr HcalTrack::getHit( int i ) const {
        return ( (HitPtr)(hits_->At(i)) );
    }

}
