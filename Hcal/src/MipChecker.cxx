/**
 * @file MipChecker.cxx
 * @brief Implementation file for MipChecker 
 * @author Tom Eichlersmith, University of Minnesota
 */

#include "Hcal/MipChecker.h"

namespace ldmx {

    void MipChecker::configure(const ldmx::ParameterSet& ps) {
        
        numFalsePos_ = 0;
        numTruePos_ = 0;
        numFalseNeg_ = 0;
        numTrueNeg_ = 0;

        numWrong_ = 0;
        numRight_ = 0;

        numEvents_ = 0;

        return;
    }

    void MipChecker::analyze(const ldmx::Event& event) {
        
        //get necessary collections
        const TClonesArray *simhits = event.getCollection( "HcalSimHits" , "sim" );
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
        std::vector<int> muon_pdgs( 2 , 13 );
        muon_pdgs[1] *= -1;

        for ( size_t iH = 0; iH < simhits->GetEntriesFast(); iH++ ) {
            SimCalorimeterHit* simhit = (SimCalorimeterHit *)(simhits->At(iH));
            nmuons += simhit->getNumberOfContribs( muon_pdgs );
        }
        
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
            numTruePos_++;
        else if ( triggerpass and !realpass )
            numFalsePos_++;
        else if ( !triggerpass and realpass )
            numFalseNeg_++;
        else
            numTrueNeg_++;

        return;
    }
    
    void MipChecker::onProcessEnd() {
        
        printf( "\n==================================\n" );
        printf( "|                | Mip Trigger   |\n" );
        printf( "|________________| Pass  | Fail  |\n" );
        printf( "|Sim        Pass | %5.4f | %5.4f |\n" );
        printf( "|Particles  Fail | %5.4f | %5.4f |\n" );
        printf( "==================================\n" );
        
        printf( "\n===========================\n" );
        printf( "|     Mip Track Recon     |\n" );
        printf( "|Right Count | Wrong Count|\n" );
        printf( "|%11.6f | %11.6f|\n" );
        printf( "===========================\n" );
        
        return;
    }

}

DECLARE_ANALYZER_NS(ldmx, MipChecker);
