/**
 * @file HcalMipTrack.cxx
 * @brief Class implementation of a track through the Hcal
 * @author Tom Eichlersmith, University of Minnesota
 */

#include "Event/HcalMipTrack.h"

ClassImp( ldmx::HcalMipTrack );

namespace ldmx {
    
    HcalMipTrack::HcalMipTrack() 
        : TObject(), hits_(new TRefArray())
          { } 

    HcalMipTrack::~HcalMipTrack() {
        Clear();
        hits_->Delete();
        
    }
    
    HcalMipTrack& HcalMipTrack::operator= ( const HcalMipTrack &track ) {
        
        if ( this != &track ) { //self-assignment guard
            
            this->hits_ = new TRefArray( *track.hits_ );
        }

        return *this;
    }
             
    void HcalMipTrack::Clear( Option_t *opt ) {

        TObject::Clear();

        hits_->Clear(); 
        
        return;
    }

    void HcalMipTrack::addHit( HitPtr hit ) {
        //add hit to track
        hits_->Add( hit );

        return;
    }
    
    int HcalMipTrack::getNHits() const {
        return hits_->GetEntriesFast();
    }

    HcalHit* HcalMipTrack::getHit( int i ) const {
        return ( dynamic_cast<HcalHit*>(hits_->At(i)) );
    }
    
    bool HcalMipTrack::isEmpty() const {
        return ( hits_->IsEmpty() );
    }

    bool HcalMipTrack::isBroken() const {
        for ( int iH = 0; iH < getNHits(); iH++ ) {
            if ( getHit( iH ) == nullptr )
                return true;
        } //iterate through hits checking for nullptr (iH)

        return false;
    }
}
