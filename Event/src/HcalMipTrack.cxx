/**
 * @file HcalMipTrack.cxx
 * @brief Class implementation of a track through the Hcal
 * @author Tom Eichlersmith, University of Minnesota
 */

#include "Event/HcalMipTrack.h"

ClassImp( ldmx::HcalMipTrack );

namespace ldmx {
    
    HcalMipTrack::HcalMipTrack() 
        : TObject(), hcalHits_(new TRefArray()), zxGr_(), zyGr_()
          { } 

    HcalMipTrack::~HcalMipTrack() {
        Clear();
        hcalHits_->Delete();
    }
    
    HcalMipTrack& HcalMipTrack::operator= ( const HcalMipTrack &track ) {
        
        if ( this != &track ) { //self-assignment guard
            
            this->hcalHits_ = new TRefArray( *track.hcalHits_ );
            this->zxGr_ = track.zxGr_;
            this->zyGr_ = track.zyGr_;
        }

        return *this;
    }
             
    void HcalMipTrack::Clear( Option_t *opt ) {

        TObject::Clear();

        hcalHits_->Clear();

        zxGr_.Set( 0 );
        zyGr_.Set( 0 );
        
        return;
    }
 
    void HcalMipTrack::addHit( HcalHit* hit ) {
        hcalHits_->Add( hit );
        return;
    }

    void HcalMipTrack::addPoint( const float x , const float y , const float z,
                   const float ex, const float ey, const float ez ) {
        
        zxGr_.SetPoint( zxGr_.GetN() , z , x );
        zxGr_.SetPointError( zxGr_.GetN()-1 , ez , ex );

        zyGr_.SetPoint( zyGr_.GetN() , z , y );
        zyGr_.SetPointError( zyGr_.GetN()-1 , ez , ey );

        return;       
    }
       
    int HcalMipTrack::getNHits() const {
        return hcalHits_->GetEntriesFast();
    }

    HcalHit* HcalMipTrack::getHit( int i ) const {
        return ( dynamic_cast<HcalHit*>(hcalHits_->At(i)) );
    }

    void HcalMipTrack::evalFit( const float z , float &x , float &y ) {
        
        zxGr_.Fit( "pol1" , "Q" );
        zyGr_.Fit( "pol1" , "Q" );
        
        TF1 *fitresult;
        fitresult = zxGr_.GetFunction( "pol1" );
        x = fitresult->Eval( z );
        fitresult = zyGr_.GetFunction( "pol1" );
        y = fitresult->Eval( z );
        
        return;
    }
    
    bool HcalMipTrack::isEmpty() const {
        return ( hcalHits_->IsEmpty() );
    }

    bool HcalMipTrack::isBroken() const {
        for ( int iH = 0; iH < getNHits(); iH++ ) {
            if ( getHit( iH ) == nullptr )
                return true;
        } //iterate through hits checking for nullptr (iH)

        return false;
    }
}
