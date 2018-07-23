/**
 * @file HcalTrack.cxx
 * @brief Class implementation of a track through the Hcal
 * @author Tom Eichlersmith, University of Minnesota
 */

#include "Event/HcalTrack.h"

ClassImp( ldmx::HcalTrack );

namespace ldmx {
    
    HcalTrack::HcalTrack() 
        : TObject(), hits_(new TRefArray()), nlayhits_(0),
          seedlayer_(0), seedstrip_(0),
          oddfitres_(nullptr), oddgr_(new TGraph()), evenfitres_(nullptr), evengr_(new TGraph()) { }
            

    HcalTrack::~HcalTrack() {
        Clear();
        hits_->Delete();
        
        delete oddgr_;
        delete evengr_;
    }
    
    HcalTrack& HcalTrack::operator= ( const HcalTrack &track ) {
        
        if ( this != &track ) { //self-assignment guard
            
            this->hits_ = new TRefArray( *track.hits_ );
            this->setSeed( track.getSeedLayer() , track.getSeedStrip() );
            this->nlayhits_ = track.getNLayHits();
            this->oddfitres_ = nullptr; //re-set whenever evalFit is called on an odd layer 
            this->oddgr_ = new TGraph( *track.oddgr_ );
            this->evenfitres_ = nullptr; //re-set whenever evalFit is called on an even layer
            this->evengr_ = new TGraph( *track.evengr_ );

        }

        return *this;
    }
             
    void HcalTrack::Clear( Option_t *opt ) {

        TObject::Clear();

        hits_->Clear(); 
        nlayhits_ = 0;

        seedlayer_ = 0;
        seedstrip_ = 0;

        delete oddgr_;
        delete evengr_;
        
        oddgr_ = new TGraph();
        evengr_ = new TGraph();
        
        return;
    }

    void HcalTrack::addHit( HitPtr hit ) {
        //add hit to track
        hits_->Add( hit );

        //add hit information to graphs
        float layer = hit->getLayer();
        if ( true ) { //static_cast<int>(layer) % 2 == 1 ) { //odd layer
            oddgr_->SetPoint( oddgr_->GetN() , layer , hit->getStrip() );
        } else { //even layer
            evengr_->SetPoint( evengr_->GetN() , layer , hit->getStrip() );
        }

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

    float HcalTrack::evalFit( const int layer ) {
        float strip;
        if ( true ) { //layer % 2 == 1 ) { //odd layer
            oddgr_->Fit( "pol1" , "Q" );
            oddfitres_ = oddgr_->GetFunction( "pol1" );
            strip = oddfitres_->Eval( layer );
        } else { //even layer
            evengr_->Fit( "pol1" , "Q" );
            evenfitres_ = evengr_->GetFunction( "pol1" );
            strip = evenfitres_->Eval( layer );
        }

        return strip;
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
