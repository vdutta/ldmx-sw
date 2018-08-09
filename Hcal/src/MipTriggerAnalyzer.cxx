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
        const TriggerResult *hcalMipTrigger = 
            (const TriggerResult *)(triggers->FindObject( hcalMipTriggerObjectName_.c_str() ));
        
        hTracksPerEvent_->Fill( hcalMipTrigger->getAlgoVar4() );
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
