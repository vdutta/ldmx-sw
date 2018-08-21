/**
 * @file MuonTrigger.cxx
 * @brief Implementation file for MuonTrigger producer
 * @author Tom Eichlersmith, University of Minnesota
 */

#include "Hcal/MuonTrigger.h"

namespace ldmx {

    void MuonTrigger::configure(const ldmx::ParameterSet& ps) {
        
        hcalHitCollName_ = ps.getString( "HcalHitCollectionName" );
        hcalHitPassName_ = ps.getString( "HcalHitPassName" );

        triggerObjectName_ = ps.getString( "TriggerObjectName" );

        minConsecLayersHit_[0] = ps.getInteger( "MinConsecutiveLayersHitBackHcal" );
        minConsecStripsHit_[0] = ps.getInteger( "MinConsecutiveStripsHitBackHcal" );

        minConsecLayersHit_[1] = ps.getInteger( "MinConsecutiveLayersHitTopHcal" );
        minConsecStripsHit_[1] = ps.getInteger( "MinConsecutiveStripsHitTopHcal" );

        minConsecLayersHit_[2] = ps.getInteger( "MinConsecutiveLayersHitBottomHcal" );
        minConsecStripsHit_[2] = ps.getInteger( "MinConsecutiveStripsHitBottomHcal" );

        minConsecLayersHit_[3] = ps.getInteger( "MinConsecutiveLayersHitLeftHcal" );
        minConsecStripsHit_[3] = ps.getInteger( "MinConsecutiveStripsHitLeftHcal" );

        minConsecLayersHit_[4] = ps.getInteger( "MinConsecutiveLayersHitRightHcal" );
        minConsecStripsHit_[4] = ps.getInteger( "MinConsecutiveStripsHitRightHcal" );

        return;
    }

    void MuonTrigger::produce(ldmx::Event& event) {
        
        //get hcal hits
        const TClonesArray *hcalHits = event.getCollection( hcalHitCollName_ , hcalHitPassName_ );

        //sort layers hit and strips hit into sections
        std::set<int> layersHit[5], stripsHit[5];
        int nHits = hcalHits->GetEntriesFast();
        HcalHit *first(nullptr), *zero(nullptr);
        double dist2first(0);
        if ( nHits > 0 ) {
            zero = (HcalHit *)(hcalHits->At(0));
        }

        for ( int iH = 0; iH < nHits; iH++ ) {
            HcalHit *chit = (HcalHit *)(hcalHits->At(iH));
            int section = chit->getSection();
            if ( section >= 0 and section < 5 ) {
                layersHit[ section ].insert( chit->getLayer() );
                stripsHit[ section ].insert( chit->getStrip() );
                
                if ( !chit->getNoise() and zero ) {
                    
                    double dx = zero->getX() - chit->getX();
                    double dy = zero->getY() - chit->getY();
                    double dz = zero->getZ() - chit->getZ();
                    double dist = sqrt( dx*dx + dy*dy + dz*dz );

                    if ( dist > dist2first ) {
                        dist2first = dist;
                        first = chit;
                    }

                } //chit is not noise

            } else {
                std::cerr << "WARNING [ MuonTrigger::produce ] : Unknown HcalSection!" << std::endl;
            }
        } //add each hit information to set

        //find last
        HcalHit *last(nullptr);
        if ( first ) {
            double dist2last(0);
            for ( int iH = 0; iH < nHits; iH++ ) {
                HcalHit *chit = (HcalHit *)(hcalHits->At(iH));
    
                if ( !chit->getNoise() ) {
    
                    double dx = first->getX() - chit->getX();
                    double dy = first->getY() - chit->getY();
                    double dz = first->getZ() - chit->getZ();
                    
                    double dist = sqrt( dx*dx + dy*dy + dz*dz );

                    if ( dist > dist2last ) {
                        dist2last = dist;
                        last = chit;
                    } //improve on distance
                } //ignore noise hits
            } //loop through all hits
        } //first exists

        bool pass = false;
        int consecLayersHit[5], consecStripsHit[5];
        for ( int s = 0; s < 5; s++ ) {
            
            consecLayersHit[s] = consecCount( layersHit[s] );
            consecStripsHit[s] = consecCount( stripsHit[s] );

            if ( consecLayersHit[s] > minConsecLayersHit_[s] and
                 consecStripsHit[s] > minConsecStripsHit_[s]
               ) {
                pass = true;
            } //if any section passes
            
        } //check each section
        
        double pathUnc(-1);
        if ( first != last ) {
            double dx = abs(last->getX() - first->getX());
            double dy = abs(last->getY() - first->getY());
            double dz = abs(last->getZ() - first->getZ());
            double dx2 = dx*dx;
            double dy2 = dy*dy;
            double dz2 = dz*dz;
            double scint_width = 100.0;
            double scint_thick = 6.0;
            double w2 = scint_width*scint_width;
            double t2 = scint_thick*scint_thick;
                
            if ( first->getSection() == 0 or last->getSection() == 0 ) {
                //back - layers along z
                double path = sqrt( t2*( 2 + dx2/dz2 + dy2/dz2 ) );
                pathUnc = sqrt( ( t2*(dx2+dy2)*(t2*(dx2+dy2)+w2*dz2) )/( sqrt(3)*dz2*dz2*(dx2+dy2+2*dz2) ) );
                pathUnc /= path;
            } else if ( first->getSection() < 3 and last->getSection() < 3 ) {
                //Top and/or bottom - layers along y
                double path = sqrt( t2*( 2 + dx2/dy2 + dz2/dy2 ) );
                pathUnc = sqrt( ( t2*(dx2+dz2)*(t2*(dx2+dz2)+w2*dy2) )/( sqrt(3)*dy2*dy2*(dx2+dz2+2*dy2) ) );
                pathUnc /= path;
            } else if ( first->getSection() > 2 and last->getSection() > 2 ) {
                //Left and/or right - layers along x
                double path = sqrt( t2*( 2 + dz2/dx2 + dy2/dx2 ) );
                pathUnc = sqrt( ( t2*(dz2+dy2)*(t2*(dz2+dy2)+w2*dx2) )/( sqrt(3)*dx2*dx2*(dz2+dy2+2*dx2) ) );
                pathUnc /= path;
            }
        }

        //build result object
        result_.set( triggerObjectName_ , pass , 21 );
        for ( int s = 0; s < 5; s++ ) {
            result_.setAlgoVar( 4*s   , minConsecLayersHit_[s] );
            result_.setAlgoVar( 4*s+1 , minConsecStripsHit_[s] );
            result_.setAlgoVar( 4*s+2 , consecLayersHit[s] );
            result_.setAlgoVar( 4*s+3 , consecStripsHit[s] );
        }
        result_.setAlgoVar( 4*5 , pathUnc );

        event.addToCollection( "Trigger" , result_ );

        return;
    }
    
    int MuonTrigger::consecCount( const std::set<int> &numbers ) const {
        
        int maxConsec(-1);
        
        if ( !numbers.empty() ) {
            
            int consec(0);
            int prev = *numbers.cbegin();
            for ( int curr : numbers ) {
                   
                if ( curr - prev > 1 ) {
                    if ( consec > maxConsec ) //improve upon max
                        maxConsec = consec;
    
                    consec = 0; //reset counter
                } //check if not consecutive anymore
     
                consec++; //count this current number

                prev = curr; //update prev
            } //go through ordered numbers
    
            if ( consec > maxConsec )
                maxConsec = consec;
    
        } //numbers is non-empty

        return maxConsec;

    }   
}

DECLARE_PRODUCER_NS(ldmx, MuonTrigger);
