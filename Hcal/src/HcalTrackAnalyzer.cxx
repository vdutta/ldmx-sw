/**
 * @file HcalTrackAnalyzer.cxx
 * @brief Implementation file for HcalTrackAnalyzer
 * @author Tom Eichlersmith, University of Minnesota
 */

#include "Hcal/HcalTrackAnalyzer.h"

#include "TROOT.h"

namespace ldmx {
    
    void HcalTrackAnalyzer::configure(const ldmx::ParameterSet& ps) {
        
        hcalMipTracksCollName_ = ps.getString( "HcalMipTracksCollectionName" );
        hcalMipTracksPassName_ = ps.getString( "HcalMipTracksPassName" );

        return;
    }

    void HcalTrackAnalyzer::analyze(const ldmx::Event& event) {
       
        const TClonesArray* tracks = event.getCollection( hcalMipTracksCollName_ , hcalMipTracksPassName_ );
        
        int ntracks = tracks->GetEntriesFast();
        hTracksPerEvent_->Fill( ntracks );
        for( int iT = 0; iT < ntracks; iT++ ) {
            HcalMipTrack *track = (HcalMipTrack *)( tracks->At( iT ));
            hClustersPerTrack_->Fill( track->getNClusters() );
            
            if ( track->isEmpty() or track->isBroken() ) {
                std::cerr << "Lost the HcalHits" << std::endl;
            } else {
                for ( int iH = 0; iH < track->getNHits(); iH++ ) {
                    HcalHit* chit = track->getHit( iH );
                    hStripsInTracks_->Fill( chit->getStrip() );
                } //hits in track (iH)
            }
        } //tracks in event (iT)
        
        //Drop non-interesting events
        if ( ntracks != 1 ) {
            setStorageHint( hint_mustKeep );
        } else {
            setStorageHint( hint_mustDrop );
        }
        
        return;
    }

    void HcalTrackAnalyzer::onProcessStart() {
 
        getHistoDirectory();
        
        hTracksPerEvent_ = new TH1F( "hTracksPerEvent_" , "Tracks Per Event" ,
            5 , -0.5 , 4.5 );
        
        hClustersPerTrack_ = new TH1F( "hClustersPerTrack_" , "MIP Clusters Per Track" ,
            100 , 0.0 , 200.0 );

        hStripsInTracks_ = new TH1F( "hStripsInTracks_" , "Strips In Each Track" ,
            50 , 0.0 , 50.0 );

        return;
    }
}

DECLARE_ANALYZER_NS(ldmx, HcalTrackAnalyzer);
