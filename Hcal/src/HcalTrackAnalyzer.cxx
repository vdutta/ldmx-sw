/**
 * @file HcalTrackAnalyzer.cxx
 * @brief Implementation file for HcalTrackAnalyzer
 * @author Tom Eichlersmith, University of Minnesota
 */

#include "Hcal/HcalTrackAnalyzer.h"

namespace ldmx {
    
    void HcalTrackAnalyzer::configure(const ldmx::ParameterSet& ps) {
        trackcollname_ = ps.getString( "HcalTrackCollectionName" );
       
        return;
    }

    void HcalTrackAnalyzer::analyze(const ldmx::Event& event) {
        std::cout << "In Analyzer" << std::endl;
        const TClonesArray* tracks = event.getCollection(trackcollname_);
        std::cout << "Got tracks" << std::endl;
        int ntracks = tracks->GetEntriesFast();
        h_tracksperevent_->Fill( ntracks );
        for( int iT = 0; iT < ntracks; iT++ ) {
            HcalTrack *curr_track = (HcalTrack *)( tracks->At(iT) );
            if (iT < 3 ) {
                h_layhitspertrack_[iT]->Fill( curr_track->getNLayHits() );
            } else {
                std::cout << "[ HcalTrackAnalyzer::analyze ]: More than 3 tracks!" << std::endl;
            }
            //Log hits in track
            int nhits = curr_track->getNHits();
            std::cout << "NHits: " << nhits << std::endl;
            for ( int iH = 0; iH < nhits; iH++ ) {
                HitPtr curr_hit = curr_track->getHit( iH );
                float curr_Layer = curr_hit->getLayer();
            } //hits in track (iH)

        } //tracks in event (iT)
        std::cout << "Made track histograms" << std::endl;
        //const TClonesArray* hits = event.getCollection("hcalDigis");
        std::cout << "Touched via TClonesArray" << std::endl;
        for( int iT = 0; iT < ntracks; iT++ ) {
            HcalTrack *curr_track = (HcalTrack *)( tracks->At(iT) );
            //Log hits in track
            int nhits = curr_track->getNHits();
            std::cout << "NHits: " << nhits << std::endl;
            for ( int iH = 0; iH < nhits; iH++ ) {
                HitPtr curr_hit = curr_track->getHit( iH );
                
                if ( curr_hit == nullptr ) {
                    std::cout << "nullptr after touch by TClonesArray" << std::endl;
                    break;
                }
                std::cout << iH << "\r";
            } //hits in track (iHi)
            std::cout << "Done with Track" << std::endl;
        }
        std::cout << "Setting Storage Hint" << std::endl;
        //Drop non-interesting events
        if ( ntracks != 1 ) {
            setStorageHint( hint_mustKeep );
        } else {
            setStorageHint( hint_mustDrop );
        }
        std::cout << "Exiting Analysis" << std::endl;
        return;
    }

    void HcalTrackAnalyzer::onProcessStart() {
 
        getHistoDirectory();
        
        h_tracksperevent_ = new TH1F( "h_tracksperevent_" , "Tracks Per Event" ,
            11 , -0.5 , 10.5 );
        
        for ( int i = 0; i < 3; i++ ) {
            h_layhitspertrack_[i] = new TH1F( ("h_layhitspertrack_"+std::to_string(i)).c_str() ,
                ("Layer Hits Per Track "+std::to_string(i)).c_str() ,
                201 , -0.5 , 200.5 );
            h_layhitspertrack_[i]->SetLineColor( i+1 );
        }

        return;
    }

}

DECLARE_ANALYZER_NS(ldmx, HcalTrackAnalyzer);
