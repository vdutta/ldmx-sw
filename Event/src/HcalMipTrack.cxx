/**
 * @file HcalMipTrack.cxx
 * @brief Class implementation of a track through the Hcal
 * @author Tom Eichlersmith, University of Minnesota
 */

#include "Event/HcalMipTrack.h"

ClassImp( ldmx::HcalMipTrack );

namespace ldmx {
    
    HcalMipTrack::HcalMipTrack() 
        : TObject(), hits_(new TRefArray()), zxGr_(), zyGr_()
          { } 

    HcalMipTrack::~HcalMipTrack() {
        Clear();
        hits_->Delete();
    }
    
    HcalMipTrack& HcalMipTrack::operator= ( const HcalMipTrack &track ) {
        
        if ( this != &track ) { //self-assignment guard
            
            this->hits_ = new TRefArray( *track.hits_ );
            this->zxGr_ = track.zxGr_;
            this->zyGr_ = track.zyGr_;
        }

        return *this;
    }
             
    void HcalMipTrack::Clear( Option_t *opt ) {

        TObject::Clear();

        hits_->Clear();

        zxGr_.Set( 0 );
        zyGr_.Set( 0 );
        
        return;
    }

    void HcalMipTrack::addCluster( const MipCluster &cluster ) {
       
        for ( unsigned int iH = 0; iH < cluster.getNumHits(); iH++ ) {
            hcalHits_->Add( cluster.getHcalHit( iH ) );
        } //add all hits in cluster (iH)
       
        //put the real space point in the graph for fitting
        float x,y,z,ex,ey,ez;
        
        cluster.getPoint( x , y , z , ex , ey , ez );

        zxGr_.SetPoint( zxGr_->GetN() , z , x );
        zxGr_.SetPointError( zxGr_->GetN()-1 , ez , ex );

        zyGr_.SetPoint( zyGr_->GetN() , z , y );
        zyGr_.SetPointError( zyGr_->GetN()-1 , ez , ey );

        return;
    }
    
    int HcalMipTrack::getNHits() const {
        return hits_->GetEntriesFast();
    }

    HcalHit* HcalMipTrack::getHit( int i ) const {
        return ( dynamic_cast<const HcalHit*>(hits_->At(i)) );
    }

    void HcalMipTrack::evalFit( const float z , float &x , float &y ) {
        
        zxGr_.Fit( "pol1" , "Q" );
        zyGr_.Fit( "pol1" , "Q" );

        x = (zxGr_.GetFunction( "pol1" ))->Eval( z );
        y = (zyGr_.GetFunction( "pol1" ))->Eval( z );
        
        return;
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
