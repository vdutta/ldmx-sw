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
        const TClonesArray* hits = event.getCollection("hcalDigis");

        int nhits = hits->GetEntriesFast();
        for( int i = 0; i < nhits; i++) {
            HitPtr curr_hit = (HitPtr)( hits->At(i) );
            h_peperhit_[0]->Fill( curr_hit->getPE() );
            h_energyperhit_[0]->Fill( curr_hit->getEnergy() );
            h_layerofhit_[0]->Fill( curr_hit->getLayer() );
            h_stripofhit_[0]->Fill( curr_hit->getStrip() );
        } //hits in event (i)
        
        int ntracks = tracks->GetEntriesFast();
        h_tracksperevent_->Fill( ntracks );
        for( int i = 0; i < ntracks; i++ ) {
            HcalTrack *curr_track = (HcalTrack *)( tracks->At(i) );
            if (i < 3 ) {
                h_layhitspertrack_[i]->Fill( curr_track->getNLayHits() );
            } else {
                std::cout << "[ HcalTrackAnalyzer::analyze ]: More than 3 tracks!" << std::endl;
            }
            for( int j = 0; j < curr_track->getNHits(); j++) {
                HitPtr curr_hit = curr_track->getHit( j );
                h_peperhit_[1]->Fill( curr_hit->getPE() );
                h_energyperhit_[1]->Fill( curr_hit->getEnergy() );
                h_layerofhit_[1]->Fill( curr_hit->getLayer() );
                h_stripofhit_[1]->Fill( curr_hit->getStrip() );
            } //hits in track (j)
        } //tracks in event (i)
        
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
        
        for ( int i = 0; i < 3; i++ ) {
            h_layhitspertrack_[i] = new TH1F( ("h_layhitspertrack_"+std::to_string(i)).c_str() ,
                ("Layer Hits Per Track "+std::to_string(i)).c_str() ,
                201 , -0.5 , 200.5 );
            h_layhitspertrack_[i]->SetLineColor( i+1 );
        }

    
        h_peperhit_[0] = new TH1F( "h_peperhit_allhits" , "PE Per Hit for all hits" ,
                        2000 , 0 , 2000 );
        h_energyperhit_[0] = new TH1F( "h_energyperhit_allhits" , "Energy Per Hit for all hits" ,
                        200 , 0 , 200 );
        h_layerofhit_[0] = new TH1F( "h_layerofhit_allhits" , "Layer of Hit for all hits" ,
                        100 , 0 , 100 );
        h_stripofhit_[0] = new TH1F( "h_stripofhit_allhits" , "Strip of Hit for all hits" ,
                        100 , 0 , 100 );


        h_peperhit_[1] = new TH1F( "h_peperhit_trackhits" , "PE Per Hit for track hits" ,
                        2000 , 0 , 2000 );
        h_energyperhit_[1] = new TH1F( "h_energyperhit_trackhits" , "Energy Per Hit for track hits" ,
                        200 , 0 , 200 );
        h_layerofhit_[1] = new TH1F( "h_layerofhit_trackhits" , "Layer of Hit for track hits" ,
                        100 , 0 , 100 );
        h_stripofhit_[1] = new TH1F( "h_stripofhit_trackhits" , "Strip of Hit for track hits" ,
                        100 , 0 , 100 );

        return;
    }

}

DECLARE_ANALYZER_NS(ldmx, HcalTrackAnalyzer);
