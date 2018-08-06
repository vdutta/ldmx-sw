/**
 * @file HcalMipTrack.cxx
 * @brief Class implementation of a track through the Hcal
 * @author Tom Eichlersmith, University of Minnesota
 */

#include "Event/HcalMipTrack.h"

ClassImp( ldmx::HcalMipTrack );

namespace ldmx {
    
    HcalMipTrack::HcalMipTrack() 
        : TObject(), hcalHits_(new TRefArray()),
          zxGr_(new TGraphErrors()), zyGr_(new TGraphErrors())
          { } 

    HcalMipTrack::~HcalMipTrack() {
        Clear();
        hcalHits_->Delete();
        
        if ( zxGr_ )
            delete zxGr_;
        
        if ( zyGr_ )
            delete zyGr_;
    }
    
    HcalMipTrack& HcalMipTrack::operator= ( const HcalMipTrack &track ) {
        
        if ( this != &track ) { //self-assignment guard
            
            //delete existing objects
            if ( this->hcalHits_ )
                delete hcalHits_;
            if ( this->zxGr_ )
                delete this->zxGr_;
            if ( this->zyGr_ )
                delete this->zyGr_;
            
            this->hcalHits_ = new TRefArray( *track.hcalHits_ );
            this->zxGr_ = new TGraphErrors( *track.zxGr_ );
            this->zyGr_ = new TGraphErrors( *track.zyGr_ );
        }

        return *this;
    }
             
    void HcalMipTrack::Clear( Option_t *opt ) {

        TObject::Clear(opt);

        hcalHits_->Clear(opt);

        if ( zxGr_ )
            delete zxGr_;
        
        if ( zyGr_ )
            delete zyGr_;

        zxGr_ = new TGraphErrors();
        zyGr_ = new TGraphErrors();

        return;
    }
 
    void HcalMipTrack::addHit( HcalHit* hit ) {
        hcalHits_->Add( hit );
        return;
    }

    void HcalMipTrack::addPoint( const std::vector<double> &point , const std::vector<double> &errors ) {
        
        zxGr_->SetPoint( zxGr_->GetN() , point[2] , point[0] );
        zxGr_->SetPointError( zxGr_->GetN()-1 , errors[2] , errors[0] );

        zyGr_->SetPoint( zyGr_->GetN() , point[2] , point[1] );
        zyGr_->SetPointError( zyGr_->GetN()-1 , errors[2] , errors[1] );

        return;       
    }
       
    int HcalMipTrack::getNHits() const {
        return hcalHits_->GetEntriesFast();
    }

    HcalHit* HcalMipTrack::getHit( int i ) const {
        return ( dynamic_cast<HcalHit*>(hcalHits_->At(i)) );
    }

    void HcalMipTrack::evalFit( const float z , float &x , float &y ) {
        
        zxGr_->Fit( "pol1" , "Q" );
        zyGr_->Fit( "pol1" , "Q" );
        
        TF1 *fitresult;
        fitresult = zxGr_->GetFunction( "pol1" );
        x = fitresult->Eval( z );
        fitresult = zyGr_->GetFunction( "pol1" );
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