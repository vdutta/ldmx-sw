/**
 * @file HcalTrackAnalyzer.cxx
 * @brief Implementation file for HcalTrackAnalyzer
 * @author Tom Eichlersmith, University of Minnesota
 */

#include "Hcal/HcalTrackAnalyzer.h"

namespace ldmx {
    
    void HcalTrackAnalyzer::configure(const ldmx::ParameterSet& ps) {
        trackcollname_ = ps.getString("TrackCollectionName");
       
        return;
    }

    void HcalTrackAnalyzer::analyze(const ldmx::Event& event) {
        
        const TClonesArray* tracks = event.getCollection(trackcollname_);
        
        int ntracks = tracks->GetEntriesFast();
        h_tracksperevent_->Fill( ntracks );
        for( int i = 0; i < ntracks; i++ ) {
            HcalTrack *curr_track = (HcalTrack *)( tracks->At(i) );
            h_hitspertrack_->Fill( curr_track->getNHits() );
        }
        
        return;
    }

    void HcalTrackAnalyzer::onProcessStart() {
 
        getHistoDirectory();
        
        h_tracksperevent_ = new TH1F( "h_tracksperevent_" , "Tracks Per Event" ,
            11 , -0.5 , 10.5 );
        h_hitspertrack_ = new TH1F( "h_hitspertrack_" , "Hits Per Track" ,
            151 , -0.5 , 150.5 );

        return;
    }

    void HcalTrackAnalyzer::onProcessEnd() {
    }

}

DECLARE_ANALYZER_NS(ldmx, HcalTrackAnalyzer);
