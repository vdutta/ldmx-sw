/**
 * @file HcalTrackAnalyzer.cxx
 * @brief Implementation file for HcalTrackAnalyzer
 * @author Tom Eichlersmith, University of Minnesota
 */

#include "playtest/HcalTrackAnalyzer.h"

namespace ldmx {
    
    void HcalTrackAnalyzer::configure(const ldmx::ParameterSet& ps) {
        trackcollname_ = ps.getString("TrackCollectionName");
       
        return;
    }

    void HcalTrackAnalyzer::analyze(const ldmx::Event& event) {
        
        const TClonesArray* tracks = event.getCollection(trackcollname_);

        h_tracksperevent_->Fill( tracks->GetEntriesFast() );

        return;
    }

    void HcalTrackAnalyzer::onProcessStart() {
 
        getHistoDirectory();
        
        h_tracksperevent_ = new TH1F( "h_tracksperevent_" , "Tracks Per Event" ,
            10 , -0.5 , 10.5 );

        return;
    }

    void HcalTrackAnalyzer::onProcessEnd() {
    }

}

DECLARE_ANALYZER_NS(ldmx, HcalTrackAnalyzer);
