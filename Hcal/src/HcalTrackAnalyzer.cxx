/**
 * @file HcalTrackAnalyzer.cxx
 * @brief Implementation file for HcalTrackAnalyzer
 * @author Tom Eichlersmith, University of Minnesota
 */

#include "Hcal/HcalTrackAnalyzer.h"

namespace ldmx {
    
    void HcalTrackAnalyzer::configure(const ldmx::ParameterSet& ps) {
        trackcollname_ = ps.getString( "HcalTrackCollectionName" , "HcalTracks" );
       
        return;
    }

    void HcalTrackAnalyzer::analyze(const ldmx::Event& event) {
        
        const TClonesArray* tracks = event.getCollection(trackcollname_);
        
        int ntracks = tracks->GetEntriesFast();
        h_tracksperevent_->Fill( ntracks );
        for( int i = 0; i < ntracks; i++ ) {
            HcalTrack *curr_track = (HcalTrack *)( tracks->At(i) );
            if (i < 9 ) {
                h_layhitspertrack_[i]->Fill( curr_track->getNLayHits() );
            } else {
                std::cout << "[ HcalTrackAnalyzer::analyze ]: More than 9 tracks!" << std::endl;
            }
        }
        
        //Drop non-interesting events
        if ( ntracks > 1 )
            setStorageHint( hint_mustKeep );
        else
            setStorageHint( hint_mustDrop );
        
        return;
    }

    void HcalTrackAnalyzer::onProcessStart() {
 
        getHistoDirectory();
        
        h_tracksperevent_ = new TH1F( "h_tracksperevent_" , "Tracks Per Event" ,
            11 , -0.5 , 10.5 );
        
        for ( int i = 0; i < 9; i++ ) {
            h_layhitspertrack_[i] = new TH1F( ("h_layhitspertrack_"+std::to_string(i)).c_str() ,
                ("Layer Hits Per Track "+std::to_string(i)).c_str() ,
                201 , -0.5 , 200.5 );
            h_layhitspertrack_[i]->SetLineColor( i+1 );
        }

        return;
    }

}

DECLARE_ANALYZER_NS(ldmx, HcalTrackAnalyzer);
