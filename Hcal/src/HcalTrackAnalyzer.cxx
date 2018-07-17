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
        
        const TClonesArray* tracks = event.getCollection(trackcollname_);
//        std::set<int> trackhitkeys;
        
        int ntracks = tracks->GetEntriesFast();
        h_tracksperevent_->Fill( ntracks );
        for( int iT = 0; iT < ntracks; iT++ ) {
            HcalTrack *curr_track = (HcalTrack *)( tracks->At(iT) );
            if (iT < 3 ) {
                h_layhitspertrack_[iT]->Fill( curr_track->getNLayHits() );
            } else {
                std::cout << "[ HcalTrackAnalyzer::analyze ]: More than 3 tracks!" << std::endl;
            }
/*            //Log hits in track
            int nhits = curr_track->getNHits();
            for ( int iH = 0; iH < nhits; iH++ ) {
                HitPtr curr_hit = curr_track->getHit( iH );
                trackhitkeys.insert( static_cast<int>( curr_hit->getLayer()*100 + curr_hit->getStrip() ) );
            } //hits in track (iH)
*/
        } //tracks in event (iT)
/*
        const TClonesArray* hits = event.getCollection("hcalDigis");
        int nhits = hits->GetEntriesFast();
        for ( int iH = 0; iH < nhits; iH++ ) {
            
            const HitPtr curr_hit = (const HitPtr)(hits->At(iH));
            float curr_PE = curr_hit->getPE();
            float curr_Energy = curr_hit->getEnergy();
            float curr_Strip = curr_hit->getStrip();
           
            if ( curr_PE > 5.5 ) {
                h_pe_nonnoise_->Fill( curr_PE );
                h_energy_nonnoise_->Fill( curr_Energy );
                h_strip_nonnoise_->Fill( curr_Strip );
                
                int curr_Key = curr_hit->getLayer()*100+curr_Strip;
                if ( trackhitkeys.count( curr_Key ) < 1 ) {
                    h_pe_notrack_->Fill( curr_PE );
                    h_energy_notrack_->Fill( curr_Energy );
                    h_strip_notrack_->Fill( curr_Strip );
                } //check if not in track

            } //check in nonnoise
            
        } //hits in event (iH)
*/        
        //Drop non-interesting events
        if ( ntracks != 1 )
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

/* Histograms to study track behavior
        h_pe_nonnoise_ = new TH1F( "h_pe_nonnoise_" , "PE Per Non-Noise Hit" , 600 , 0.0 , 600.0 ); 
        h_pe_notrack_ = new TH1F( "h_pe_notrack_" , "PE Per No Track Hit" , 600 , 0.0 , 600.0 ); 

        h_energy_nonnoise_ = new TH1F( "h_energy_nonnoise_" , "Energy Per Non-Noise Hit" , 100 , 0.0 , 100.0 );
        h_energy_notrack_ = new TH1F( "h_energy_notrack_" , "Energy Per No Track Hit" , 100 , 0.0 , 100.0 );

        h_strip_nonnoise_ = new TH1F( "h_strip_nonnoise_" , "Strip Per Non-Noise Hit" , 40 , 0.0 , 40.0 );
        h_strip_notrack_ = new TH1F( "h_strip_notrack_" , "Strip Per No Track Hit" , 40 , 0.0 , 40.0 );
*/
        return;
    }

}

DECLARE_ANALYZER_NS(ldmx, HcalTrackAnalyzer);
