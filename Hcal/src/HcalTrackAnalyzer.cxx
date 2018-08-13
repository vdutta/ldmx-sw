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
        
        for ( int actual = 0; actual < 4; actual++ ) {
            for ( int pred = 0; pred < 4; pred++ ) {
                numTracks_[actual][pred] = 0;
            }
        }

        return;
    }

    void HcalTrackAnalyzer::analyze(const ldmx::Event& event) {
       
        const TClonesArray *simparticles = event.getCollection( "SimParticles" , "sim" );
        const TClonesArray* tracks = event.getCollection( hcalMipTracksCollName_ , hcalMipTracksPassName_ );
        
        int nmuons = 0;
        for ( int iP = 0; iP < simparticles->GetEntriesFast(); iP++ ) {
            SimParticle *simp = (SimParticle *)(simparticles->At(iP));
            if ( simp->getPdgID() == 13 or simp->getPdgID() == -13 ) {
                nmuons++;
            }
        }

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
        
        numTracks_[ nmuons ][ ntracks ] ++;

//        //Drop non-interesting events
//        if ( ntracks != 1 ) {
//            setStorageHint( hint_mustKeep );
//        } else {
//            setStorageHint( hint_mustDrop );
//        }
        
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

    void HcalTrackAnalyzer::onProcessEnd() {
        
        double reconaccuracy = 0.0;
        unsigned int numEvents = 0;
        for ( int i = 0; i < 4; i++ ) {
            reconaccuracy += (double)( numTracks_[i][i] );
            for ( int j = 0; j < 4; j++ ) {
                numEvents += numTracks_[i][j];
            }
        }
        reconaccuracy /= (double)(numEvents);
        
        printf( "\n" );
        printf( " ======================================================\n" );
        printf( " |      Mip Track Reconstruction Confusion Table      |\n" );
        printf( " | Predicted ||            Actual N Tracks            |\n" );
        printf( " | N Tracks  ||    0    |    1    |    2    |    3    |\n" );
        for ( int pred = 0; pred < 4; pred++ ) {
            
            printf( " |%10d ||" , pred );
            for ( int actual = 0; actual < 4; actual++ ) {
                printf( " %7d |" , numTracks_[ actual ][ pred ] );
            } //actual number of tracks (actual)
            printf( "\n" );
        } //predicted number of tracks (pred)
        printf( " |====================================================|\n" );
        printf( " | Accuracy  || %-37f |\n" , reconaccuracy );
        printf( " ======================================================\n" );
        
        return;
    }
}

DECLARE_ANALYZER_NS(ldmx, HcalTrackAnalyzer);
