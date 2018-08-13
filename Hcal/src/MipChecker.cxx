/**
 * @file MipChecker.cxx
 * @brief Implementation file for MipChecker 
 * @author Tom Eichlersmith, University of Minnesota
 */

#include "Hcal/MipChecker.h"

namespace ldmx {

    void MipChecker::configure(const ldmx::ParameterSet& ps) {
        
        numFalsePass_ = 0;
        numTruePass_ = 0;
        numFalseFail_ = 0;
        numTrueFail_ = 0;

        for ( int actual = 0; actual < 4; actual++ ) {
            for ( int pred = 0; pred < 4; pred++ ) {
                numTracks_[actual][pred] = 0;
            }
        }

        numEvents_ = 0;

        return;
    }

    void MipChecker::analyze(const ldmx::Event& event) {
        
        //get necessary collections
        const TClonesArray *simparticles = event.getCollection( "SimParticles" , "sim" );
        const TClonesArray *miptracks = event.getCollection( "hcalMipTracks" , "recon" );
        const TClonesArray *triggers = event.getCollection( "Trigger" , "recon" );
       
        //get trigger object
        const TriggerResult *hcalMipTrigger;
        for ( int iT = 0; iT < triggers->GetEntriesFast(); iT++ ) {
            hcalMipTrigger = (const TriggerResult *)(triggers->At(iT));
            if ( hcalMipTrigger->getName() == "hcalMipTrigger" )
                break;
        }

        //count real number of muons in hcal
        int nmuons = 0;
        for ( size_t iH = 0; iH < simparticles->GetEntriesFast(); iH++ ) {
            SimParticle* simparticle = (SimParticle *)(simparticles->At(iH));
            if ( simparticle->getPdgID() == 13 or simparticle->getPdgID() == -13 ) {
                //muon or anti-muon
                nmuons++;
            } //check if muon
        } //go through sim particles
        
        //check accuracy of mip track recon
        if ( nmuons < 4 and miptracks->GetEntriesFast() < 4 )
            numTracks_[ nmuons ][ miptracks->GetEntriesFast() ] ++;

        //check accuracy of mip trigger
        bool triggerpass = hcalMipTrigger->passed();
        bool realpass = ( nmuons > 0 );
        if ( triggerpass and realpass)
            numTruePass_++;
        else if ( triggerpass and !realpass )
            numFalsePass_++;
        else if ( !triggerpass and realpass )
            numFalseFail_++;
        else
            numTrueFail_++;
        
        numEvents_++;

        return;
    }
    
    void MipChecker::onProcessEnd() {
        
        double triggeraccuracy = (numTruePass_ + numTrueFail_)/(double)(numEvents_);
        printf( "\n" );
        printf( " ===============================\n" );
        printf( " | Mip Trigger Confusion Table |\n" );
        printf( " | Mip     ||    Sim Particle  |\n" );
        printf( " | Trigger ||   Pass | Fail    |\n" );
        printf( " |    Pass ||%7d | %-7d |\n" , numTruePass_ , numFalsePass_ );
        printf( " |    Fail ||%7d | %-7d |\n" , numFalseFail_ , numTrueFail_ );
        printf( " |=============================|\n" );
        printf( " | Accuracy | %-16f |\n" , triggeraccuracy );
        printf( " ===============================\n" );
        
        double reconaccuracy = 0.0;
        for ( int i = 0; i < 4; i++ ) {
            reconaccuracy += (double)( numTracks_[i][i] );
        }
        reconaccuracy /= (double)(numEvents_);
        
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

DECLARE_ANALYZER_NS(ldmx, MipChecker);
