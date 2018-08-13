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

        numWrong_ = 0;
        numRight_ = 0;

        numEvents_ = 0;

        return;
    }

    void MipChecker::analyze(const ldmx::Event& event) {
        
        //get necessary collections
        const TClonesArray *simparticles = event.getCollection( "SimParticles" , "sim" );
        const TClonesArray *miptracks = event.getCollection( "hcalMipTracks" , "recon" );
        const TClonesArray *triggers = event.getCollection( "Trigger" , "recon" );
        
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
        if ( miptracks->GetEntriesFast() == nmuons ) {
            numRight_++;
        } else {
            numWrong_++;
        }
        
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
        
        //Percentage
        double events = numEvents_/100.0;

        double truepass  = (double)(numTruePass_)/events;
        double falsepass = (double)(numFalsePass_)/events;
        double falsefail = (double)(numFalseFail_)/events;
        double truefail  = (double)(numTrueFail_)/events;

        printf( "\n==================================\n" );
        printf( "|                |  Mip Trigger    |\n" );
        printf( "|________________| Pass   | Fail   |\n" );
        printf( "|Sim        Pass | %5.2f%% | %5.2f%% |\n" , truepass , falsefail );
        printf( "|Particles  Fail | %5.2f%% | %5.2f%% |\n" , falsepass , truefail );
        printf( "==================================\n" );
        
        double right = (double)(numRight_)/events;
        double wrong = (double)(numWrong_)/events;

        printf( "\n=============================\n" );
        printf( "|      Mip Track Recon      |\n" );
        printf( "| Right Count | Wrong Count |\n" );
        printf( "| %10.6f%% | %10.6f%% |\n" , right , wrong );
        printf( "=============================\n" );
        
        return;
    }

}

DECLARE_ANALYZER_NS(ldmx, MipChecker);
