/**
 * @file MipTriggerAnalyzer.cxx
 * @brief Implementation file for MipTriggerAnalyzer class
 * @author Tom Eichlersmith, University of Minnesota
 */

#include "Hcal/MipTriggerAnalyzer.h"

namespace ldmx {

    void MipTriggerAnalyzer::configure(const ldmx::ParameterSet& ps) {
        
        hcalMipTriggerObjectName_ = ps.getString( "HcalMipTriggerObjectName" );

        hcalMipTriggerPassName_ = ps.getString( "HcalMipTriggerPassName" );

        return;
    }

    void MipTriggerAnalyzer::analyze(const ldmx::Event& event) {
        
        //get list of triggers
        const TClonesArray *triggers = event.getCollection( "Trigger" , hcalMipTriggerPassName_ );
        
        //get hcal mip trigger
        int ntriggers = triggers->GetEntriesFast();
        const TriggerResult *hcalMipTrigger;
        for ( int i = 0; i < ntriggers; i++ ) {
            hcalMipTrigger = (const TriggerResult *)(triggers->At(i));
            if ( hcalMipTrigger->getName() == hcalMipTriggerObjectName_ )
                break;
        }
        
        if ( hcalMipTrigger ) {
            hTracksPerEvent_->Fill( hcalMipTrigger->getAlgoVar4() );
        } else {
            std::cerr << hcalMipTriggerObjectName_ << " was not found in Trigger Collection in pass ";
            std::cerr << hcalMipTriggerPassName_ << std::endl;
        }

        return;
    }

    void MipTriggerAnalyzer::onProcessStart() {
        
        getHistoDirectory();

        hTracksPerEvent_ = new TH1F( "hTracksPerEvent_" , "Tracks Found Per Event" ,
            11 , -0.5 , 10.5 );

        return;
    }

}

DECLARE_ANALYZER_NS(ldmx, MipTriggerAnalyzer);
